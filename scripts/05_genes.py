"""Probe sets -> genes (fixed annotation, per-sample; PLAN §1) and fRMA QA.

Usage: python 05_genes.py <frma_probe_tsv> <out_gene_tsv> <qa_json_name>
"""
import os
import sys

import numpy as np
import pandas as pd

from common import PROC, RESULTS, load_series_values, provenance, write_json

src, dst, qa_name = sys.argv[1], sys.argv[2], sys.argv[3]
probes = pd.read_csv(src, sep="\t", index_col=0)
annot = pd.read_csv(os.path.join(PROC, "probe_annot.tsv"), sep="\t", index_col=0)
annot = annot.loc[annot.index.intersection(probes.index)]
ok = annot.SYMBOL.notna() & ~annot.SYMBOL.astype(str).str.contains(r"\|")
chrom = annot.CHR.astype(str)
sexchr_specific = chrom.eq("Y") | annot.SYMBOL.eq("XIST")
keep = annot.index[ok & ~sexchr_specific]
genes = probes.loc[keep].groupby(annot.loc[keep, "SYMBOL"]).mean()
genes.to_csv(dst, sep="\t", float_format="%.5f")

# QA: agreement with the submitters' matrix (same arrays), identity check only
series = np.log2(load_series_values())[probes.columns]
common_p = probes.index.intersection(series.index)
r = [np.corrcoef(probes.loc[common_p, g], series.loc[common_p, g])[0, 1] for g in probes.columns]
# is each array most correlated with its own counterpart? (identity check)
ps = probes.loc[common_p].sub(probes.loc[common_p].mean(1), axis=0)
ss = series.loc[common_p].sub(series.loc[common_p].mean(1), axis=0)
cross = np.corrcoef(ps.T.values, ss.T.values)[: len(r), len(r):]
self_best = int((cross.argmax(1) == np.arange(len(r))).sum())
out = dict(provenance=provenance(), n_arrays=probes.shape[1], n_probes=probes.shape[0],
           n_probes_unique_symbol=int(ok.sum()), n_removed_sex_specific=int((ok & sexchr_specific).sum()),
           n_genes=genes.shape[0], per_array_r_with_submitters_min=float(np.min(r)),
           per_array_r_with_submitters_median=float(np.median(r)),
           arrays_best_matched_to_own_submitted_profile=self_best)
write_json(os.path.join(RESULTS, qa_name), out)
print(out)
