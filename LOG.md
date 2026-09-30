# LOG — decisions, pre-registrations, results

Running exploratory-variant count: **0 / 12**

## 2026-09-30 — Setup (no modelling)

- **Network.** NCBI GEO, CRAN, Bioconductor and download.pytorch.org were initially blocked by
  the environment's network policy. I did not substitute an unofficial data copy. The policy
  was opened by the user, and all data then came directly from ftp.ncbi.nlm.nih.gov.
- **Data downloaded.** GSE53987_RAW.tar (md5 67576fd879ea2ef079aec41f723b7f09, 205 CEL.gz),
  GSE53987_series_matrix.txt.gz (md5 26b164ffb27bb115709a06b6f1e9d635).
  Submitters' processing (series matrix): "RMA … Brain regions were normalized separately."
- **Metadata parse.** 205 samples, 75 donors (BD 19, Ctrl 19, SCZ 19, MDD 18); samples
  BD 52 / Ctrl 55 / MDD 50 / SCZ 48. Donor sample counts: 3 (58 donors), 2 (14), 1 (3).
  Metadata sex is consistent within every donor, so the mixed-sex problem is only visible in
  expression.
- **Software installed.** R 4.3.3, limma 3.58.1, affy 1.80.0 (Ubuntu archive);
  hgu133plus2cdf 2.18.0 (Bioconductor source, md5 verified);
  hgu133plus2.db 3.13.0, frma 1.54.0, hgu133plus2frmavecs 1.5.0, preprocessCore 1.64.0 (Bioconductor 3.18, verified to load); scikit-learn 1.9.1, torch 2.14.0+cpu (PyPI / pytorch CPU index).
- **Decision: fRMA instead of whole-dataset RMA.** Reason: whole-dataset RMA fits its quantile
  reference and probe effects on the test arrays, which conflicts with rule 3. Argued in
  PLAN.md §1 and §10. The fallback is RMA run separately on the dev and test arrays.
- **Decision: exclude the 5 mixed-sex donors from the primary analysis and the test pool.**
  Reasons are in PLAN.md §3 (cross-fold leakage risk, label noise, the foreign array can't be
  identified).
- **Papers read.** Relevant passages are summarised in PLAN.md §8. The warfarin DDI paper is
  off-topic for expression classification.

PLAN.md written; **waiting for approval**. Nothing has been split, normalised or trained.

## 2026-09-30 — Supplied script and literature search (no modelling)

- User supplied `gen_prep-v7.py`. Assessed in PLAN.md §12. **Not adopted**: gene selection
  from a global GEO2R top table (selection bias), per-sample labels with no donor/region,
  probe-level features, submitters' all-array RMA. Its log2 step is correct for this dataset
  (series-matrix values are linear, max 30,750). It is kept as control **C8** to measure the
  selection bias.
- Literature search (web). Read: Lanz 2015 (source study; 19 tetrads matched on age, sex,
  PMI; PVALB decrease in SCZ confirmed), Brouard 2024 (GNN benchmark; GNNs rarely beat
  simple models; configuration-model random graph), Ambroise & McLachlan 2002, Varma &
  Simon 2006, Tomita 2004, McCall 2010 (abstract only). WebFetch cannot reach NCBI; full
  texts were fetched with curl from PMC / Europe PMC.
- **Plan changes:** added C2b (degree-preserving random graph, from Brouard 2024), C8
  (selection-bias demo), and the PVALB positive control; runtime estimate 14.6 → 17.7 h.
  The variant budget is unchanged (0 / 12 used).

## 2026-09-30 — Plan approved. Step 0: sex QA (`scripts/01_sex_qa.py` → `results/qa_sex.json`)

- Score = mean log2 of RPS4Y1, DDX3Y, KDM5D, USP9Y, EIF1AY minus XIST, computed on the
  submitters' matrix. Two clean modes (male ≈ +3.9, female ≈ −5.2; every sample ≥ 3.9 away
  from the threshold).
- Donors with arrays of both sexes = **BD_18, MDD_11, MDD_2, SCZ_4, SCZ_7**. This matches
  the protocol exactly.
- **All 8 discordant arrays of those donors are hippocampus arrays.** No PFC or striatum array
  is discordant, which suggests the upstream problem sits in hippocampus sample handling.
- **New finding: Control_7.** All 3 arrays (hip, pfc, str) are male by expression, but the
  metadata says F. The arrays agree with each other, so this is not a mixed-sex donor. Either
  the metadata sex is wrong, or the whole donor is someone other than recorded.
- Per PLAN §3 ("if the list differs … stop and report before splitting"), **stopped before
  the split** and asked the user how to handle Control_7.
