"""Single final evaluation on the locked test set (rule 9). Refuses to run twice.

Usage: python 08_final_test.py <best_model_config> <best_baseline_config>

Protocol: the frozen configuration is refit on ALL clean development samples (in-sample limma,
region scaling, graph, etc. exactly as inside a CV fold, with the full dev set as "train"),
then predicts the 49 test samples once. fRMA is per-array, so the test arrays are processed
independently. Stochastic models: each of the 5 seeds is evaluated; the headline is the mean
over seeds (matching how CV was seed-averaged), and every seed's result is stored.
CI: donor-level cluster bootstrap (10,000 resamples) of the headline predictions.
"""
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd

if len(sys.argv) != 3:
    sys.exit(__doc__)
import pipeline as P
from analysis import summarise
from common import CLASSES, PROC, RAW, RESULTS, ROOT, SPLITS, provenance, write_json

OUT = os.path.join(RESULTS, "final_test.json")
if os.path.exists(OUT):
    sys.exit(f"{OUT} exists: the locked test set has already been used once.")
REG = __import__("06_run").REG
configs = sys.argv[1:3]

lock = json.load(open(os.path.join(SPLITS, "test_lock.json")))
meta = pd.read_csv(os.path.join(PROC, "samples.tsv"), sep="\t", index_col=0)
test = meta.loc[lock["test_gsm"]]

# 1. fRMA on the test arrays (per-array; frozen vectors)
tg = os.path.join(PROC, "test_genes.tsv")
if not os.path.exists(tg):
    d = os.path.join(RAW, "cel_test")
    os.makedirs(d, exist_ok=True)
    lst = os.path.join(PROC, "test_cel_list.txt")
    open(lst, "w").write("\n".join(test.cel) + "\n")
    subprocess.run(["tar", "-xf", os.path.join(RAW, "GSE53987_RAW.tar"), "-C", d, "-T", lst], check=True)
    subprocess.run(["Rscript", os.path.join(ROOT, "scripts", "03_frma.R"), d, os.path.join(PROC, "frma_test_probes.tsv")], check=True)
    subprocess.run([sys.executable, os.path.join(ROOT, "scripts", "05_genes.py"), os.path.join(PROC, "frma_test_probes.tsv"),
                    tg, "qa_frma_test.json"], check=True)
Xte_all = pd.read_csv(tg, sep="\t", index_col=0).T.loc[test.index]

# 2. training set = all clean development samples
train = P.global_train()
P.ensure_limma([(train, "plain", None), (train, "cov", None)])
X, M = P.data()


class FinalFold(P.Fold):
    """Same preprocessing as a CV fold, with the locked test samples as the held-out set."""

    def __init__(self, mode):
        spec = P.limma_spec(train, mode if mode != "global" else "plain")
        lt = P.limma_table(*spec).sort_values("P", kind="mergesort")
        self.limma, self.genes, self.pvals = lt, list(lt.gene[:P.K]), lt.P.values[:P.K]
        self.train, self.val = train, list(test.index)
        xtr = X.loc[train, self.genes].values
        xva = Xte_all.reindex(columns=self.genes).values
        assert not np.isnan(xva).any(), "gene missing in test matrix"
        rtr, rva = M.loc[train, "region"].values, test.region.values
        self.Xtr, self.Xva = np.empty_like(xtr), np.empty_like(xva)
        for r in np.unique(rtr):
            mu, sd = xtr[rtr == r].mean(0), xtr[rtr == r].std(0, ddof=1) + 1e-8
            self.Xtr[rtr == r] = (xtr[rtr == r] - mu) / sd
            self.Xva[rva == r] = (xva[rva == r] - mu) / sd
        self.ytr = np.array([CLASSES.index(c) for c in M.loc[train, "diagnosis"]])
        self.yva = np.array([CLASSES.index(c) for c in test.diagnosis])
        self.donor_va = test.donor.values
        m2 = test.copy()
        m2["sex_male"] = (m2.sex_expr == "M").astype(float)
        self.Ctr = M.loc[train, P.COVARS].values.astype(float)
        self.Cva = m2[P.COVARS].values.astype(float)
        self.n_classes = 4


def cluster_bootstrap(y, pred, donors, B=10000, seed=20260930):
    from sklearn.metrics import f1_score
    rng = np.random.default_rng(seed)
    ud = np.unique(donors)
    idx_by = {d: np.where(donors == d)[0] for d in ud}
    vals = []
    for _ in range(B):
        pick = rng.choice(ud, size=len(ud), replace=True)
        ii = np.concatenate([idx_by[d] for d in pick])
        vals.append(f1_score(y[ii], pred[ii], labels=range(4), average="macro", zero_division=0))
    return [float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))]


results = {}
for c in configs:
    spec = REG[c]
    f = FinalFold(spec.get("mode", "plain"))
    seeds = spec.get("seeds", (0,))
    per_seed, probs = [], []
    for s in seeds:
        p = spec["fn"](f, s)
        probs.append(p)
        per_seed.append(P.metrics(f.yva, p, CLASSES))
    head = dict(macro_f1=float(np.mean([m["macro_f1"] for m in per_seed])),
                mcc=float(np.mean([m["mcc"] for m in per_seed])),
                auc=float(np.mean([m["auc"] for m in per_seed])))
    cis = [cluster_bootstrap(f.yva, p.argmax(1), f.donor_va) for p in probs]
    cv = summarise(c)
    results[c] = dict(test=dict(headline_mean_over_seeds=head, macro_f1_ci95_donor_bootstrap_per_seed=cis,
                                per_seed=per_seed, n_test_samples=len(f.yva), n_test_donors=int(len(set(f.donor_va))),
                                chance_macro_f1=0.25, predictions=dict(gsm=f.val, y_true=f.yva, prob=[np.round(p, 6) for p in probs])),
                      cv=dict(macro_f1_mean=cv["macro_f1_mean"], macro_f1_fold_sd=cv["macro_f1_fold_sd"],
                              macro_f1_seed_sd=cv["macro_f1_seed_sd"], mcc_mean=cv["mcc_mean"], auc_mean=cv["auc_mean"],
                              n_folds=cv["n_folds"], n_seeds=cv["n_seeds"]))
write_json(OUT, dict(configs=configs, results=results, n_train_samples=len(train), provenance=provenance(),
                     test_lock_sha256=open(os.path.join(SPLITS, "test_lock.json.sha256")).read().split()[0]))
print(json.dumps({c: dict(test=r["test"]["headline_mean_over_seeds"], ci=r["test"]["macro_f1_ci95_donor_bootstrap_per_seed"],
                          cv=r["cv"]) for c, r in results.items()}, indent=1))
