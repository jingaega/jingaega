"""Two-track final comparison table (test set + repeated CV) -> results/final_comparison.{md,json}."""
import json
import os

from analysis import summarise
from common import RESULTS, provenance, write_json

ft = json.load(open(os.path.join(RESULTS, "final_test.json")))
role = {ft["configs"][0]: "Best model", ft["configs"][1]: "Best simple baseline"}
rows, L = [], []
L.append("# Final comparison — both evaluations side by side\n")
L.append("Chance macro-F1 = **0.25** (4 balanced classes). Test set = 49 samples from 18 donors; "
         "**a test set this small makes the confidence interval wide**: the donor-bootstrap 95% CI spans "
         "about 0.3 macro-F1, and removing any single test donor moves the estimate by up to ±0.04. CV = donor-grouped stratified 5 folds × 5 repeats on "
         "143 development samples / 52 donors.\n")
L.append("| Role | Config | **Test** macro-F1 [95% CI, donor bootstrap] | Test MCC | Test AUC | "
         "**CV** macro-F1 mean ± SD (25 folds) | CV seed SD | CV MCC | CV AUC |")
L.append("|---|---|---|---|---|---|---|---|---|")
for c in ft["configs"]:
    t = ft["results"][c]["test"]
    m = t["per_seed"][0]
    ci = t["macro_f1_ci95_donor_bootstrap_per_seed"][0]
    s = summarise(c)
    rows.append(dict(role=role[c], config=c, test_macro_f1=m["macro_f1"], test_ci95=ci, test_mcc=m["mcc"],
                     test_auc=m["auc"], test_precision=m["precision"], test_recall=m["recall"],
                     test_confusion=m["confusion"], cv_macro_f1_mean=s["macro_f1_mean"],
                     cv_macro_f1_fold_sd=s["macro_f1_fold_sd"], cv_macro_f1_seed_sd=s["macro_f1_seed_sd"],
                     cv_mcc=s["mcc_mean"], cv_auc=s["auc_mean"], cv_precision=s["precision"], cv_recall=s["recall"],
                     cv_confusion_mean=s["confusion_rows_true_cols_pred"], cv_vs_chance=s["vs_chance"]))
    L.append(f"| {role[c]} | `{c}` | **{m['macro_f1']:.3f}** [{ci[0]:.3f}, {ci[1]:.3f}] | {m['mcc']:.3f} | "
             f"{m['auc']:.3f} | **{s['macro_f1_mean']:.3f} ± {s['macro_f1_fold_sd']:.3f}** | "
             f"{s['macro_f1_seed_sd']:.3f} | {s['mcc_mean']:.3f} | {s['auc_mean']:.3f} |")
L.append("\nBoth models are deterministic, so seed SD = 0 by construction.\n")
for r in rows:
    L.append(f"## {r['role']}: `{r['config']}`\n")
    L.append("| Class | Test precision | Test recall | CV precision | CV recall |")
    L.append("|---|---|---|---|---|")
    for k in r["test_precision"]:
        L.append(f"| {k} | {r['test_precision'][k]:.2f} | {r['test_recall'][k]:.2f} | "
                 f"{r['cv_precision'][k]:.2f} | {r['cv_recall'][k]:.2f} |")
    L.append("\nTest confusion matrix (rows = true, cols = predicted; BD, Control, MDD, SCZ):\n")
    L.append("```")
    for lab, row in zip(["BD", "Control", "MDD", "SCZ"], r["test_confusion"]):
        L.append(f"{lab:8s} " + " ".join(f"{x:4d}" for x in row))
    L.append("```\nCV confusion (summed over the 5 folds of a repetition, averaged over 5 repetitions):\n")
    L.append("```")
    for lab, row in zip(["BD", "Control", "MDD", "SCZ"], r["cv_confusion_mean"]):
        L.append(f"{lab:8s} " + " ".join(f"{x:5.1f}" for x in row))
    L.append("```")
    vc = r["cv_vs_chance"]
    L.append(f"\nCV vs chance: Nadeau–Bengio corrected p = {vc['p_corrected']:.3f} "
             f"(uncorrected {vc['p_uncorrected']:.3g}).\n")
write_json(os.path.join(RESULTS, "final_comparison.json"), dict(rows=rows, provenance=provenance()))
open(os.path.join(RESULTS, "final_comparison.md"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
