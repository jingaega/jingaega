"""Aggregation and statistics (PLAN §5, rule 8)."""
import glob
import json
import os

import numpy as np
from scipy import stats

from common import RESULTS


def load(config):
    recs = [json.load(open(p)) for p in sorted(glob.glob(os.path.join(RESULTS, "raw", config, "*.json")))]
    if not recs:
        raise FileNotFoundError(config)
    return recs


def per_fold(recs, metric="macro_f1"):
    """dict (rep, fold) -> list of seed values."""
    d = {}
    for r in recs:
        d.setdefault((r["rep"], r["fold"]), []).append(r["metrics"][metric])
    return d


def summarise(config, chance=None):
    recs = load(config)
    labels = recs[0]["labels"]
    chance = chance if chance is not None else 1.0 / len(labels)
    pf = per_fold(recs)
    keys = sorted(pf)
    fold_means = np.array([np.mean(pf[k]) for k in keys])
    seed_sd = np.array([np.std(pf[k], ddof=1) if len(pf[k]) > 1 else 0.0 for k in keys])
    out = dict(config=config, labels=labels, chance_macro_f1=chance, n_folds=len(keys),
               n_seeds=len(pf[keys[0]]),
               macro_f1_mean=fold_means.mean(), macro_f1_fold_sd=fold_means.std(ddof=1),
               macro_f1_seed_sd=seed_sd.mean(), fold_macro_f1=fold_means.tolist())
    for m in ("mcc", "auc"):
        v = per_fold(recs, m)
        vals = [np.mean([x for x in v[k] if x is not None]) for k in keys if any(x is not None for x in v[k])]
        out[f"{m}_mean"], out[f"{m}_fold_sd"] = float(np.mean(vals)), float(np.std(vals, ddof=1))
    for part in ("precision", "recall"):
        out[part] = {c: float(np.mean([r["metrics"][part][c] for r in recs])) for c in labels}
    # confusion: sum over folds within (rep, seed), then average over (rep, seed)
    cm = {}
    for r in recs:
        cm.setdefault((r["rep"], r["seed"]), np.zeros((len(labels),) * 2))
        cm[(r["rep"], r["seed"])] += np.array(r["metrics"]["confusion"])
    out["confusion_rows_true_cols_pred"] = np.mean(list(cm.values()), 0).round(2).tolist()
    tr = np.mean([r["fold_sizes"]["train_samples_used"] for r in recs])
    va = np.mean([r["fold_sizes"]["val_samples"] for r in recs])
    out["mean_train_samples"], out["mean_val_samples"] = float(tr), float(va)
    out["vs_chance"] = nb_test(fold_means - chance, va / tr)
    out["provenance"] = recs[0]["provenance"]
    out["fold_sizes_example"] = recs[0]["fold_sizes"]
    return out


def nb_test(d, ratio):
    """Nadeau-Bengio corrected resampled t-test on J fold-level differences.

    Assumes the correlation between fold estimates is rho = n_test/(n_train+n_test), i.e.
    variance inflation (1/J + n_test/n_train). Uncorrected paired t reported beside it.
    """
    d = np.asarray(d, float)
    J, m, v = len(d), d.mean(), d.var(ddof=1)
    if v == 0:
        return dict(mean_diff=m, J=J, t_corrected=None, p_corrected=None, t_uncorrected=None, p_uncorrected=None)
    tc = m / np.sqrt((1 / J + ratio) * v)
    tu = m / np.sqrt(v / J)
    return dict(mean_diff=float(m), J=J, n_test_over_n_train=float(ratio),
                t_corrected=float(tc), p_corrected=float(2 * stats.t.sf(abs(tc), J - 1)),
                t_uncorrected=float(tu), p_uncorrected=float(2 * stats.t.sf(abs(tu), J - 1)))


def compare(config_a, config_b, metric="macro_f1"):
    """a - b on paired folds (seed-averaged within fold)."""
    ra, rb = load(config_a), load(config_b)
    pa, pb = per_fold(ra, metric), per_fold(rb, metric)
    keys = sorted(set(pa) & set(pb))
    assert len(keys) == len(pa) == len(pb), "folds are not paired"
    d = np.array([np.mean(pa[k]) - np.mean(pb[k]) for k in keys])
    tr = np.mean([r["fold_sizes"]["train_samples_used"] for r in ra])
    va = np.mean([r["fold_sizes"]["val_samples"] for r in ra])
    return dict(a=config_a, b=config_b, metric=metric, **nb_test(d, va / tr))


def holm(pvals):
    p = np.asarray(pvals, float)
    order = np.argsort(p)
    adj = np.empty_like(p)
    running = 0
    for i, idx in enumerate(order):
        running = max(running, min(1, (len(p) - i) * p[idx]))
        adj[idx] = running
    return adj.tolist()
