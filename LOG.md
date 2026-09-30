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

## 2026-09-30 — Split, preprocessing, baseline

- **Control_7 decision (user):** keep it; use expression-derived sex as the sex covariate for
  every donor (`samples.tsv: sex_expr`). Only V4/V5 use sex.
- **Test set locked** (`splits/test_lock.json`, sha256 c822ca51…9e9f, scikit-learn 1.9.1,
  random_state 20260930): 18 donors / 49 samples. BD 5/14, Control 5/14, MDD 4/12, SCZ 4/9
  (donors/samples). Development: 52 clean donors / 143 samples (BD 13/36, Ctrl 14/41,
  MDD 12/33, SCZ 13/33) + 13 mixed-sex samples, which are quarantined for C6. Control_7 is in
  the development set.
- **CV folds** saved (`splits/cv_*.json`): 5×5 StratifiedGroupKFold; validation folds hold
  27–32 samples. C6 design fixed: validation folds identical to the primary; only the
  training set changes.
- **fRMA** on the 156 development arrays only (test arrays have not been extracted from the
  tar). QA: per-array r with the submitters' RMA values 0.961–0.98, and every array's
  best-matching submitted profile is its own. 20,794 genes after probe→gene collapse; 67
  sex-specific probe sets removed.
- **limma per fold:** consensus within-donor correlation ≈ 0.24; about 38 s per fit
  (`duplicateCorrelation`). Parallelised across 4 R processes; the method is unchanged.
- **Runtime re-estimate:** GCN 56 s per cell single-threaded → about 0.5 h per 125-cell
  configuration with 4 workers. 7 GCN configurations ≈ 3.5 h; about 700 limma fits for the
  controls ≈ 1.9 h. Total ≈ 6 h, well under 20 h, so no cuts are needed.

### C1 baseline (control, not counted) — `results/raw/C1_logreg`
Macro-F1 **0.311** (chance 0.25), fold SD 0.088, seed SD 0 (deterministic); MCC 0.096;
AUC 0.588. Versus chance: NB corrected t = 1.28, **p = 0.21** (uncorrected p = 0.002).
**It does not beat chance under the corrected test.** Per class recall/precision:
BD 0.24/0.26, Control 0.29/0.24, MDD 0.21/0.24, **SCZ 0.56/0.60**. Whatever signal there is,
is in SCZ.

## 2026-09-30 — Queue incident
- The first launch of every GCN configuration failed immediately (the multiprocessing
  closure couldn't be pickled). No cells were written. Fixed by making `make_gcn` a
  module-level `functools.partial`; the model and hyperparameters are unchanged. The GCN
  block is re-queued (`run_queue_gcn.sh`) to start after the limma-heavy controls finish.

## 2026-09-30 06:32 — Container restart
- The container restarted at about 06:31, killing both queues. All files survived: data,
  limma cache, raw cells, and the R and Python packages. Cells are written atomically
  (tmp + rename), so none are partial. Relaunched as one resumable queue
  (`scripts/run_all.sh`), which skips cached cells. At restart, C2b_V2 had 65/125 cells.

## 2026-09-30 07:05 — Second container restart
- The container is reclaimed a few minutes after the agent's turn goes idle, which kills
  detached jobs (C2b_V2 advanced only 65 → 81 cells between 06:33 and 07:05). From now on the
  agent stays attached to the queue with foreground log watches until the runs finish.
  The resumable cache means nothing already written is lost.

### V1 — reference GCN (variant #7 of 12, pre-registered PLAN §7) — `results/raw/V1_gcn_coexpr`
Macro-F1 **0.268** (chance 0.25), fold SD 0.066, **seed SD 0.060** (seed noise ≈ fold noise);
AUC 0.578. Versus C1: Δ = −0.042, NB corrected p = 0.43 (uncorrected 0.039). Recall
BD 0.18 / Ctrl 0.22 / MDD 0.24 / SCZ 0.56. **Pre-registered failure criterion met** (not
above C1). Its graph controls (C2, C2b, C3) are pending.
Also finished: C2b_V2 (degree-preserving random STRING graph): 0.293, fold SD 0.069,
seed SD 0.066.

### V1 graph controls
- C3 no-graph (A = I): 0.275 (fold SD 0.065, seed SD 0.059). V1 − C3 = −0.007, NB p = 0.82.
- C2 uniform random, same edge count: 0.240. V1 − C2 = +0.028, NB p = 0.49.
- C2b degree-preserving random: 0.276. V1 − C2b = −0.007, NB p = 0.77.
**The co-expression graph carries no detectable signal beyond its degree structure or no
graph at all.** This agrees with Brouard et al. 2024.

### V2 — STRING PPI GCN (variant #8) and controls
V2 0.277 (fold SD 0.067, seed SD 0.059; 224–869 edges per fold, mean 411). Versus C1:
−0.033 (NB p = 0.41). Versus C3 no-graph: +0.002 (p = 0.94). Versus C2 uniform (0.290):
−0.013 (p = 0.75). Versus C2b degree-preserving (0.293): −0.016 (p = 0.58).
**Failure criterion met.** The prior-knowledge graph adds nothing either.

### C8 — selection-bias demo (gen_prep-style global top table)
Macro-F1 **0.476** vs C1 0.311: **+0.165**, NB p = 0.002 (uncorrected 1.8e-9). Versus chance
p = 0.0001. Choosing genes once on all development samples before CV produces an apparently
"significant" classifier from a dataset that, analysed honestly, doesn't beat chance
(Ambroise & McLachlan 2002). It is the same size as the per-sample-split inflation the
protocol cites (+0.156).

### C5 — binary SCZ vs Control (signal check)
Macro-F1 **0.721** (chance 0.50), fold SD 0.147, MCC 0.49, AUC 0.82. Versus chance: **NB
corrected p = 0.010** (uncorrected 9e-8). **Signal exists for SCZ vs Control and survives
the corrected test.** Precision/recall: Control 0.76/0.81, SCZ 0.76/0.66.
**PVALB positive control passes:** lower in SCZ in 25/25 folds (mean logFC −0.62, median
P 0.009, median rank 1462 of about 15.6k genes).

### C4 — leakage demonstration (sample-level vs donor-grouped split)
Sample-level 5×5 CV: macro-F1 **0.653** (AUC 0.88) vs donor-grouped C1 0.311: **+0.343**,
NB p < 0.001. **Conflict with the protocol's figure (+0.156):** the inflation here is about
twice as large. A likely reason is that our pipeline standardises within region and uses
all 3 regions per donor, so a donor's other two regions in training are near-duplicates of
the held-out sample. Per protocol I report the observed value, not the cited one.

### V5 — covariate residualisation (variant #9)
0.301 (fold SD 0.093, AUC 0.616). Versus C1: −0.009, NB p = 0.85. **Failure criterion met.**
Removing pH/RIN/PMI/age/sex doesn't help (and doesn't hurt: the SCZ signal isn't a
pH/RIN artefact at this resolution).

### C6b — primary + 5 mixed-sex donors in training (same validation folds)
0.341; versus C1 +0.031, NB p = 0.40. The random-drop reference draws (C6c) are running.

## 2026-09-30 ~10:25 — FINALS FROZEN (rule 9), before any test-set access

All 9 pre-registered variants are complete (**9 / 12 used; the 3 reserve slots are left unused
on purpose**, because adding variants after seeing these results would be an outcome-driven
search). Development-CV macro-F1 ranking of eligible configurations (V4, the confound probe,
is ineligible per PLAN §7):
C1 0.3106 > V9 0.3017 > V5 0.3011 > V8 0.2922 > V3 0.2891 > V2 0.2771 > V1 0.2682 >
V7 0.2658 > V6 0.2656.

Applying the rule fixed in PLAN §7:
- **Best model = C1_logreg** (highest development-CV mean). No variant beat the baseline.
- **Best simple baseline = V5_logreg_residualised.** It is the highest among {C1, V3, V5, V8}
  excluding the best model.

Both are frozen as registered in `scripts/06_run.py`. The final evaluation
(`scripts/08_final_test.py C1_logreg V5_logreg_residualised`) refits each on all 143 clean
development samples and predicts the 49 locked test samples **once**. Both models are
deterministic (1 seed). CI: donor-level cluster bootstrap, 10,000 resamples.
