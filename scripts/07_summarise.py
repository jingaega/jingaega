"""Build results/summary.json, results/variants_table.md and results/confusion/*.json from raw cells."""
import glob
import json
import os

import numpy as np

from analysis import compare, holm, load, nb_test, per_fold, summarise
from common import RESULTS, provenance, write_json

VARIANTS = ["V1_gcn_coexpr", "V2_gcn_string", "V3_logreg_donor_vote", "V4_covariates_only",
            "V5_logreg_residualised", "V6_modules_logreg", "V7_rf", "V8_svm", "V9_mlp"]
RESERVE = sorted(os.path.basename(p) for p in glob.glob(os.path.join(RESULTS, "raw", "R*_*")))
CONTROLS = ["C1_logreg", "C3_gcn_nograph", "C2_V1_gcn_uniform", "C2b_V1_gcn_degree", "C2_V2_gcn_uniform",
            "C2b_V2_gcn_degree", "C4_logreg_sample_level", "C5_binary_logreg", "C6b_logreg_with_mixed",
            "C8_logreg_global_selection"]


def have(c):
    return os.path.isdir(os.path.join(RESULTS, "raw", c)) and glob.glob(os.path.join(RESULTS, "raw", c, "*.json"))


def complete(c, n_expected):
    return have(c) and len(per_fold(load(c))) == n_expected


summ, comps = {}, {}
for c in CONTROLS + VARIANTS + RESERVE:
    if have(c):
        summ[c] = summarise(c)
        write_json(os.path.join(RESULTS, "confusion", c + ".json"),
                   dict(config=c, labels=summ[c]["labels"], rows_true_cols_pred=summ[c]["confusion_rows_true_cols_pred"],
                        note="summed over the 5 folds of each (repetition, seed), then averaged"))

# variants vs baseline C1 (+ Holm over the family actually run)
fam = [v for v in VARIANTS + RESERVE if v in summ and complete(v, 25) and complete("C1_logreg", 25)]
for v in fam:
    comps[f"{v} - C1_logreg"] = compare(v, "C1_logreg")
adj = holm([comps[f"{v} - C1_logreg"]["p_corrected"] or 1 for v in fam])
for v, a in zip(fam, adj):
    comps[f"{v} - C1_logreg"]["p_corrected_holm"] = a

# graph controls
for g, ctrls in (("V1_gcn_coexpr", ["C2_V1_gcn_uniform", "C2b_V1_gcn_degree", "C3_gcn_nograph"]),
                 ("V2_gcn_string", ["C2_V2_gcn_uniform", "C2b_V2_gcn_degree", "C3_gcn_nograph"]),
                 ("V9_mlp", ["C3_gcn_nograph"])):
    for c in ctrls:
        if g in summ and c in summ:
            comps[f"{g} - {c}"] = compare(g, c)
for a, b in (("C4_logreg_sample_level", "C1_logreg"), ("C8_logreg_global_selection", "C1_logreg"),
             ("C6b_logreg_with_mixed", "C1_logreg")):
    if a in summ and b in summ:
        comps[f"{a} - {b}"] = compare(a, b)

# C6: mixed-sex sensitivity vs random donor drops (all share C1's validation folds)
c6 = {}
drops = sorted(os.path.basename(p) for p in glob.glob(os.path.join(RESULTS, "raw", "C6c_logreg_drop*")))
drops = [d for d in drops if complete(d, 25)]
if "C6b_logreg_with_mixed" in summ and drops:
    b = summ["C6b_logreg_with_mixed"]["macro_f1_mean"]
    a = summ["C1_logreg"]["macro_f1_mean"]
    dd = [summarise(d)["macro_f1_mean"] - b for d in drops]
    c6 = dict(primary_minus_with_mixed=a - b, random_drop_minus_with_mixed=dd,
              random_drop_mean=float(np.mean(dd)), random_drop_sd=float(np.std(dd, ddof=1)) if len(dd) > 1 else None,
              empirical_quantile_of_primary=float(np.mean(np.array(dd) <= a - b)), n_draws=len(dd))

# C7: permutation null (repetition 0 only)
c7 = {}
perms = sorted(os.path.basename(p) for p in glob.glob(os.path.join(RESULTS, "raw", "C7_logreg_perm*")))
perms = [q for q in perms if complete(q, 5)]
if perms and "C1_logreg" in summ:
    null = [float(np.mean([np.mean(v) for v in per_fold(load(p)).values()])) for p in perms]
    c1r0 = float(np.mean([np.mean(v) for k, v in per_fold(load("C1_logreg")).items() if k[0] == 0]))
    c7 = dict(null_macro_f1_rep0=null, null_mean=float(np.mean(null)), null_sd=float(np.std(null, ddof=1)),
              c1_macro_f1_rep0=c1r0, perm_p=float((1 + sum(x >= c1r0 for x in null)) / (1 + len(null))),
              n_permutations=len(null))

n_used = sum(1 for v in VARIANTS + RESERVE if v in summ)
write_json(os.path.join(RESULTS, "summary.json"), dict(
    provenance=provenance(), variants_run=n_used, variant_budget=12, summaries=summ, comparisons=comps,
    c6_mixed_sex_sensitivity=c6, c7_permutation_null=c7))

# markdown table
L = ["| # | Config | Kind | Macro-F1 mean (chance) | fold SD | seed SD | MCC | AUC | Δ vs C1 | NB p (corr.) | p (uncorr.) | Holm | Folds×seeds |",
     "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
count = 0
for c in CONTROLS + VARIANTS + RESERVE:
    if c not in summ:
        continue
    s = summ[c]
    kind = "variant" if c in VARIANTS + RESERVE else "control"
    if kind == "variant":
        count += 1
    k = f"{c} - C1_logreg"
    d = comps.get(k)
    fmt = lambda x: "–" if x is None else f"{x:.3f}"
    L.append(f"| {count if kind == 'variant' else ''} | {c} | {kind} | {s['macro_f1_mean']:.3f} ({s['chance_macro_f1']:.2f}) | "
             f"{s['macro_f1_fold_sd']:.3f} | {s['macro_f1_seed_sd']:.3f} | {s['mcc_mean']:.3f} | {fmt(s['auc_mean'])} | "
             f"{fmt(d['mean_diff']) if d else '–'} | {fmt(d['p_corrected']) if d else '–'} | "
             f"{fmt(d['p_uncorrected']) if d else '–'} | {fmt(d.get('p_corrected_holm')) if d else '–'} | "
             f"{s['n_folds']}×{s['n_seeds']} |")
L.append(f"\nRunning exploratory-variant count: **{count} / 12**. Chance macro-F1 = 0.25 (4-class), 0.50 (binary).")
L.append("NB p = Nadeau–Bengio corrected resampled t-test on 25 paired fold-level differences (seed-averaged within fold); "
         "assumes between-fold correlation ρ = n_test/(n_train+n_test).")
open(os.path.join(RESULTS, "variants_table.md"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
