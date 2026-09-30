"""Step 0 QA (PLAN §3): verify which donors are sex-discordant by expression.

Uses only sex-chromosome genes and metadata sex; diagnosis is never used.
Runs on the submitters' series matrix before any split.
"""
import os

import numpy as np
import pandas as pd

from common import MIXED_SEX, PROC, RESULTS, SEX_GENES_Y, load_metadata, load_series_values, provenance, write_json

meta = load_metadata()
expr = np.log2(load_series_values())  # series values are linear-scale RMA
annot = pd.read_csv(os.path.join(PROC, "probe_annot.tsv"), sep="\t", index_col=0)


def gene_mean(symbol):
    probes = annot.index[annot["SYMBOL"] == symbol]
    return expr.loc[expr.index.intersection(probes)].mean(axis=0)


y = pd.concat([gene_mean(g) for g in SEX_GENES_Y], axis=1).mean(axis=1)
xist = gene_mean("XIST")
score = (y - xist).loc[meta.index]
# Threshold: midpoint between the two modes (2-means on the 1-D score)
lo, hi = score.min(), score.max()
for _ in range(100):
    thr = (lo + hi) / 2
    lo_new, hi_new = score[score < thr].mean(), score[score >= thr].mean()
    if np.isclose(lo_new, lo) and np.isclose(hi_new, hi):
        break
    lo, hi = lo_new, hi_new
thr = (lo + hi) / 2
meta["expr_sex"] = np.where(score >= thr, "M", "F")
meta["sex_score"] = score
meta["discordant"] = meta["expr_sex"] != meta["sex"]

donor_sexes = meta.groupby("donor")["expr_sex"].nunique()
mixed = sorted(donor_sexes.index[donor_sexes > 1])
per_sample = meta.loc[meta.donor.isin(mixed) | meta.discordant,
                      ["donor", "region", "sex", "expr_sex", "sex_score"]]
out = dict(
    provenance=provenance(),
    method="score = mean log2 expr of " + ",".join(SEX_GENES_Y) + " minus XIST; 2-means threshold",
    threshold=float(thr), male_mode=float(hi), female_mode=float(lo),
    min_gap_to_threshold=float((score - thr).abs().min()),
    n_discordant_samples=int(meta.discordant.sum()),
    donors_with_both_sexes=mixed,
    expected=sorted(MIXED_SEX),
    matches_expected=mixed == sorted(MIXED_SEX),
    donors_with_discordant_metadata=sorted(meta.loc[meta.discordant, "donor"].unique()),
    samples=per_sample.reset_index().to_dict(orient="records"),
)
write_json(os.path.join(RESULTS, "qa_sex.json"), out)
print({k: v for k, v in out.items() if k not in ("samples", "provenance")})
print(per_sample.sort_values("donor").to_string())
