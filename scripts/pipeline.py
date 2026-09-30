"""In-fold pipeline, models, metrics and a cached runner (PLAN §4-§7).

Every data-dependent step (expression filter, limma selection, region scaling, graph
construction, covariate residualisation, module clustering) is fitted on the training
samples of the fold only.
"""
import hashlib
import json
import os
import subprocess
import sys
import time

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (confusion_matrix, f1_score, matthews_corrcoef, precision_recall_fscore_support,
                             roc_auc_score)
from sklearn.svm import LinearSVC

from common import CACHE, CLASSES, PROC, RESULTS, ROOT, SPLITS, provenance, write_json

K = 500                       # genes kept by limma (fixed for all models)
EXPR_TSV = os.path.join(PROC, "dev_genes.tsv")
SAMPLES_TSV = os.path.join(PROC, "samples.tsv")
LIMMA_DIR = os.path.join(CACHE, "limma")
COVARS = ["ph", "rin", "pmi", "age", "sex_male"]

_EXPR = None
_META = None


def data():
    global _EXPR, _META
    if _EXPR is None:
        _EXPR = pd.read_csv(EXPR_TSV, sep="\t", index_col=0).T          # samples x genes
        m = pd.read_csv(SAMPLES_TSV, sep="\t", index_col=0)
        m["sex_male"] = (m.sex_expr == "M").astype(float)
        test = set(json.load(open(os.path.join(SPLITS, "test_lock.json")))["test_gsm"])
        assert not test & set(_EXPR.index), "test arrays present in development expression matrix"
        _META = m.loc[_EXPR.index]
    return _EXPR, _META


# ---------------------------------------------------------------- limma (R, cached)
def limma_key(train, mode, labels=None):
    s = mode + "|" + ",".join(sorted(train)) + ("|" + ",".join(labels) if labels is not None else "")
    return hashlib.sha1(s.encode()).hexdigest()[:16]


def perm_labels(train, i):
    """Donor-level label permutation among training donors (C7), deterministic in (train set, i)."""
    _, M = data()
    seed = int(hashlib.sha1((",".join(sorted(train)) + f"|{i}").encode()).hexdigest()[:8], 16)
    rng = np.random.default_rng(seed)
    d = M.loc[train].groupby("donor").diagnosis.first()
    perm = dict(zip(d.index, rng.permutation(d.values)))
    return [perm[M.loc[g, "donor"]] for g in train]


def limma_spec(train, mode):
    """(train, limma_mode, labels_or_None) actually used for selection."""
    if mode.startswith("perm:"):
        return train, "plain", perm_labels(train, int(mode.split(":")[1]))
    if mode == "global":
        return global_train(), "plain", None
    return train, ("cov" if mode == "cov" else "plain"), None


def global_train():
    """C8: all clean development samples (mimics a top table computed before splitting)."""
    cv = json.load(open(os.path.join(SPLITS, "cv_multiclass.json")))
    f0 = cv["folds"][0]
    return sorted(f0["train"] + f0["val"])


def ensure_limma(jobs):
    """jobs: list of (train_gsm_list, mode, labels). Runs missing ones in a single R session."""
    os.makedirs(LIMMA_DIR, exist_ok=True)
    todo = {}
    for train, mode, labels in jobs:
        k = limma_key(train, mode, labels)
        if not os.path.exists(os.path.join(LIMMA_DIR, k + ".tsv")) and k not in todo:
            order = sorted(range(len(train)), key=lambda i: train[i])
            todo[k] = dict(key=k, mode=mode, gsm=",".join(train[i] for i in order),
                           labels=",".join(labels[i] for i in order) if labels is not None else "")
    if not todo:
        return
    jf = os.path.join(LIMMA_DIR, f"jobs_{os.getpid()}.tsv")
    pd.DataFrame(list(todo.values())).to_csv(jf, sep="\t", index=False)
    r = subprocess.run(["Rscript", os.path.join(ROOT, "scripts", "limma_batch.R"), EXPR_TSV, SAMPLES_TSV, jf,
                        LIMMA_DIR], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(r.stderr[-3000:])
    os.remove(jf)


def limma_table(train, mode, labels=None):
    return pd.read_csv(os.path.join(LIMMA_DIR, limma_key(train, mode, labels) + ".tsv"), sep="\t", comment="#")


# ---------------------------------------------------------------- fold preparation
class Fold:
    """Region-standardised top-K genes for one train/val split."""

    def __init__(self, train, val, mode="plain", k=K):
        X, M = data()
        self.train, self.val = list(train), list(val)
        spec = limma_spec(self.train, mode)
        lt = limma_table(*spec)
        lt = lt.sort_values("P", kind="mergesort")
        self.limma = lt
        self.genes = list(lt.gene[:k])
        self.pvals = lt.P.values[:k]
        xtr, xva = X.loc[self.train, self.genes].values, X.loc[self.val, self.genes].values
        rtr, rva = M.loc[self.train, "region"].values, M.loc[self.val, "region"].values
        ztr, zva = np.empty_like(xtr), np.empty_like(xva)
        for r in np.unique(rtr):
            mu, sd = xtr[rtr == r].mean(0), xtr[rtr == r].std(0, ddof=1) + 1e-8
            ztr[rtr == r] = (xtr[rtr == r] - mu) / sd
            zva[rva == r] = (xva[rva == r] - mu) / sd
        assert set(np.unique(rva)) <= set(np.unique(rtr))
        self.Xtr, self.Xva = ztr, zva
        ytr_lab = spec[2] if mode.startswith("perm:") else list(M.loc[self.train, "diagnosis"])
        self.ytr = np.array([CLASSES.index(c) for c in ytr_lab])
        self.yva = np.array([CLASSES.index(c) for c in M.loc[self.val, "diagnosis"]])
        self.donor_va = M.loc[self.val, "donor"].values
        self.Ctr = M.loc[self.train, COVARS].values.astype(float)
        self.Cva = M.loc[self.val, COVARS].values.astype(float)


# ---------------------------------------------------------------- models
def softmax(s):
    s = s - s.max(1, keepdims=True)
    e = np.exp(s)
    return e / e.sum(1, keepdims=True)


def _lr():
    return LogisticRegression(C=1.0, class_weight="balanced", max_iter=5000)


def m_logreg(f, seed):
    clf = _lr().fit(f.Xtr, f.ytr)
    return clf.predict_proba(f.Xva)


def m_logreg_donor_vote(f, seed):
    p = m_logreg(f, seed)
    out = p.copy()
    for d in np.unique(f.donor_va):
        idx = f.donor_va == d
        out[idx] = p[idx].mean(0)
    return out


def m_covariates(f, seed):
    mu, sd = f.Ctr.mean(0), f.Ctr.std(0) + 1e-8
    clf = _lr().fit((f.Ctr - mu) / sd, f.ytr)
    return clf.predict_proba((f.Cva - mu) / sd)


def m_residualised(f, seed):
    mu, sd = f.Ctr.mean(0), f.Ctr.std(0) + 1e-8
    Ctr = np.c_[np.ones(len(f.Ctr)), (f.Ctr - mu) / sd]
    Cva = np.c_[np.ones(len(f.Cva)), (f.Cva - mu) / sd]
    beta, *_ = np.linalg.lstsq(Ctr, f.Xtr, rcond=None)
    rtr, rva = f.Xtr - Ctr @ beta, f.Xva - Cva @ beta
    s = rtr.std(0) + 1e-8
    clf = _lr().fit(rtr / s, f.ytr)
    return clf.predict_proba(rva / s)


def m_modules(f, seed, n_modules=30):
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform
    r = np.corrcoef(f.Xtr.T)
    d = 1 - np.abs(r)
    np.fill_diagonal(d, 0)
    lab = fcluster(linkage(squareform(d, checks=False), "average"), n_modules, "maxclust")
    ftr, fva = [], []
    for m in np.unique(lab):
        idx = lab == m
        A = f.Xtr[:, idx]
        if idx.sum() == 1:
            ftr.append(A[:, 0]); fva.append(f.Xva[:, idx][:, 0]); continue
        mu = A.mean(0)
        _, _, vt = np.linalg.svd(A - mu, full_matrices=False)
        v = vt[0] * np.sign(vt[0].sum() or 1)
        ftr.append((A - mu) @ v); fva.append((f.Xva[:, idx] - mu) @ v)
    ftr, fva = np.array(ftr).T, np.array(fva).T
    m, s = ftr.mean(0), ftr.std(0) + 1e-8
    clf = _lr().fit((ftr - m) / s, f.ytr)
    return clf.predict_proba((fva - m) / s)


def m_rf(f, seed):
    clf = RandomForestClassifier(n_estimators=500, max_features="sqrt", class_weight="balanced",
                                 random_state=seed, n_jobs=1).fit(f.Xtr, f.ytr)
    return clf.predict_proba(f.Xva)


def m_svm(f, seed):
    clf = LinearSVC(C=1.0, class_weight="balanced", max_iter=20000).fit(f.Xtr, f.ytr)
    s = clf.decision_function(f.Xva)
    if s.ndim == 1:
        s = np.c_[-s, s]
    return softmax(s)


# ---------------------------------------------------------------- metrics
def metrics(y, prob, labels):
    pred = prob.argmax(1)
    L = list(range(len(labels)))
    p, r, _, _ = precision_recall_fscore_support(y, pred, labels=L, zero_division=0)
    out = dict(macro_f1=f1_score(y, pred, labels=L, average="macro", zero_division=0),
               mcc=matthews_corrcoef(y, pred),
               precision={labels[i]: p[i] for i in L}, recall={labels[i]: r[i] for i in L},
               confusion=confusion_matrix(y, pred, labels=L).tolist())
    try:
        if len(labels) == 2:
            out["auc"] = roc_auc_score(y, prob[:, 1])
        else:
            out["auc"] = roc_auc_score(y, prob, multi_class="ovr", average="macro", labels=L)
    except ValueError:
        out["auc"] = None
    return out


# ---------------------------------------------------------------- runner
def run(config, model_fn, cv_name="multiclass", seeds=(0,), mode="plain", train_mod=None, workers=1,
        folds_filter=None, extra=None):
    """Evaluate one configuration over all folds; one JSON per (rep, fold, seed), skipped if present.

    train_mod: optional function(train_gsm_list) -> new train list (C6 sensitivity).
    """
    cv = json.load(open(os.path.join(SPLITS, f"cv_{cv_name}.json")))
    folds = cv["folds"] if folds_filter is None else [f for f in cv["folds"] if folds_filter(f)]
    X, M = data()
    labels = CLASSES if cv_name != "binary_scz_ctrl" else ["Control", "SCZ"]
    outdir = os.path.join(RESULTS, "raw", config)
    os.makedirs(outdir, exist_ok=True)
    jobs = []
    for f in folds:
        tr = train_mod(f["train"]) if train_mod else f["train"]
        jobs.append(limma_spec(tr, mode))
    ensure_limma(jobs)
    cells = [(f, s) for f in folds for s in seeds
             if not os.path.exists(os.path.join(outdir, f"r{f['rep']}_f{f['fold']}_s{s}.json"))]
    if workers > 1 and cells:
        import multiprocessing as mp
        with mp.get_context("fork").Pool(workers) as pool:
            pool.starmap(_cell, [(config, model_fn, f, s, mode, train_mod, labels, outdir, extra) for f, s in cells])
    else:
        for f, s in cells:
            _cell(config, model_fn, f, s, mode, train_mod, labels, outdir, extra)
    return outdir


def _cell(config, model_fn, f, seed, mode, train_mod, labels, outdir, extra):
    t0 = time.time()
    tr = train_mod(f["train"]) if train_mod else f["train"]
    fold = Fold(tr, f["val"], mode=mode)
    if labels != CLASSES:           # binary: remap indices to the 2-label space
        fold.ytr = np.array([labels.index(CLASSES[i]) for i in fold.ytr])
        fold.yva = np.array([labels.index(CLASSES[i]) for i in fold.yva])
    fold.n_classes = len(labels)
    prob = model_fn(fold, seed)
    X, M = data()
    rec = dict(config=config, rep=f["rep"], fold=f["fold"], seed=seed, labels=labels,
               metrics=metrics(fold.yva, prob, labels),
               val_gsm=f["val"], y_true=fold.yva, prob=np.round(prob, 6),
               fold_sizes=dict(f["sizes"], train_samples_used=len(tr),
                               train_donors_used=int(M.loc[tr, "donor"].nunique())),
               n_genes=len(fold.genes), seconds=round(time.time() - t0, 2), provenance=provenance())
    if extra:
        rec["extra"] = extra(fold)
    write_json(os.path.join(outdir, f"r{f['rep']}_f{f['fold']}_s{seed}.json"), rec)
