"""Run one named configuration over the development CV (cached, resumable).

Usage: python 06_run.py <CONFIG> [--workers N]
Every configuration name is registered below; nothing runs that is not listed here.
"""
import argparse
import json
import os

import numpy as np

import gnn
import pipeline as P
from common import SPLITS

SEEDS5 = (0, 1, 2, 3, 4)
C6 = json.load(open(os.path.join(SPLITS, "c6_sensitivity.json")))
_, M = P.data()
DONOR = M.donor.to_dict()


def add_mixed(train):
    return sorted(set(train) | set(C6["mixed_sex_gsm"]))


def add_mixed_drop(i):
    drop = set(C6["random_drop_draws"][i])
    return lambda train: sorted(g for g in add_mixed(train) if DONOR[g] not in drop)


def pvalb_extra(f):
    lt = f.limma.reset_index(drop=True)
    row = lt.index[lt.gene == "PVALB"]
    return dict(pvalb_rank=int(row[0]) + 1 if len(row) else None,
                pvalb_logFC_SCZ_minus_Control=float(lt.loc[row[0], "logFC"]) if len(row) and "logFC" in lt else None,
                pvalb_P=float(lt.loc[row[0], "P"]) if len(row) else None, n_tested=len(lt))


REG = {
    # ---- required controls (not counted)
    "C1_logreg": dict(fn=P.m_logreg),
    "C4_logreg_sample_level": dict(fn=P.m_logreg, cv="multiclass_sample_level"),
    "C5_binary_logreg": dict(fn=P.m_logreg, cv="binary_scz_ctrl", extra=pvalb_extra),
    "C6b_logreg_with_mixed": dict(fn=P.m_logreg, train_mod=add_mixed),
    "C8_logreg_global_selection": dict(fn=P.m_logreg, mode="global"),
    "C3_gcn_nograph": dict(fn=gnn.make_gcn(gnn.g_none), seeds=SEEDS5, extra=gnn.edge_count),
    "C2_V1_gcn_uniform": dict(fn=gnn.make_gcn(gnn.g_coexpr, gnn.randomise_uniform), seeds=SEEDS5, extra=gnn.edge_count),
    "C2b_V1_gcn_degree": dict(fn=gnn.make_gcn(gnn.g_coexpr, gnn.randomise_degree), seeds=SEEDS5, extra=gnn.edge_count),
    "C2_V2_gcn_uniform": dict(fn=gnn.make_gcn(gnn.g_string, gnn.randomise_uniform), seeds=SEEDS5, extra=gnn.edge_count),
    "C2b_V2_gcn_degree": dict(fn=gnn.make_gcn(gnn.g_string, gnn.randomise_degree), seeds=SEEDS5, extra=gnn.edge_count),
    # ---- exploratory variants (counted)
    "V1_gcn_coexpr": dict(fn=gnn.make_gcn(gnn.g_coexpr), seeds=SEEDS5, extra=gnn.edge_count),
    "V2_gcn_string": dict(fn=gnn.make_gcn(gnn.g_string), seeds=SEEDS5, extra=gnn.edge_count),
    "V3_logreg_donor_vote": dict(fn=P.m_logreg_donor_vote),
    "V4_covariates_only": dict(fn=P.m_covariates),
    "V5_logreg_residualised": dict(fn=P.m_residualised, mode="cov"),
    "V6_modules_logreg": dict(fn=P.m_modules),
    "V7_rf": dict(fn=P.m_rf, seeds=SEEDS5),
    "V8_svm": dict(fn=P.m_svm),
    "V9_mlp": dict(fn=gnn.m_mlp, seeds=SEEDS5),
}
for i in range(20):
    REG[f"C6c_logreg_drop{i:02d}"] = dict(fn=P.m_logreg, train_mod=add_mixed_drop(i))
    # C7: 20 donor-level label permutations, repetition 0 only (5 folds each)
    REG[f"C7_logreg_perm{i:02d}"] = dict(fn=P.m_logreg, mode=f"perm:{i}", folds_filter=lambda f: f["rep"] == 0)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("config")
    ap.add_argument("--workers", type=int, default=1)
    a = ap.parse_args()
    c = REG[a.config]
    out = P.run(a.config, c["fn"], cv_name=c.get("cv", "multiclass"), seeds=c.get("seeds", (0,)),
                mode=c.get("mode", "plain"), train_mod=c.get("train_mod"), workers=a.workers,
                folds_filter=c.get("folds_filter"), extra=c.get("extra"))
    print("done", a.config, out)
