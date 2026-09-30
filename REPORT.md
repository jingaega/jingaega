# Classifying psychiatric diagnosis from postmortem brain expression (GSE53987): results

**Short answer.** On the 4-class task (SCZ / BD / MDD / Control), no model beats chance in a way
that survives the Nadeau–Bengio corrected test during development (best: logistic regression,
macro-F1 0.311 ± 0.088 vs chance 0.25, corrected p = 0.21). The locked test set gives
0.419 [0.249, 0.545] for that model. That is above chance on its face, but the interval
nearly reaches 0.25, because the test set is only 18 donors. Whatever signal exists is
**schizophrenia vs everything else**: SCZ vs Control alone is clearly learnable (macro-F1
0.721 vs chance 0.50, corrected p = 0.010), while BD and MDD are not separable from
controls or from each other. **Graph neural networks added nothing:** the reference GCN
(0.268) was not better than the same network with no graph (0.275) or with a
degree-preserving random graph (0.276). **As designed, this dataset cannot support 4-class
diagnosis; it supports, at best, a modest SCZ-vs-rest signal.**

Everything below is reproducible from `scripts/` (numbered in run order). Every number
comes from a JSON under `results/`. The decision trail is in `LOG.md`; the pre-registered
plan is in `PLAN.md`.

## 1. What was done

- **Data.** GEO GSE53987 raw CEL files (205 arrays, 75 donors, 3 regions), normalised with
  **frozen RMA** (fRMA; per-array, so no held-out array influences training features). Probe
  sets were collapsed to 20,794 genes, and sex-chromosome-specific genes were removed.
  QA: r ≥ 0.96 with the submitters' RMA values, and every array matches its own submitted
  profile.
- **Mixed-sex donors.** QA confirmed exactly BD_18, MDD_11, MDD_2, SCZ_4, SCZ_7. All 8
  discordant arrays are hippocampus arrays. These 5 donors were excluded from the primary
  analysis and from the test pool, per the pre-modelling decision. The QA also found
  **Control_7** (all arrays male by expression, metadata female). It was kept, with sex
  taken from expression, per the user's decision.
- **Locked test set** (created first): 18 donors / 49 samples, donor-grouped and
  class-stratified, scikit-learn 1.9.1, sha256 `c822ca51…`. Development: 52 donors / 143
  samples.
- **Pipeline, fitted inside every training fold:** expression filter, then limma moderated
  F with `duplicateCorrelation` blocking on donor (region in the design), then the top 500
  genes, then per-region standardisation using training statistics. Fixed hyperparameters
  for every model, with no tuning, so no model got a larger search budget than any other.
- **Evaluation:** donor-grouped stratified 5×5 CV, the Nadeau–Bengio corrected resampled
  t-test on 25 paired fold differences (assumes between-fold correlation
  ρ = n_test/(n_train + n_test) ≈ 0.2), with uncorrected p and Holm-adjusted p alongside.

## 2. Every configuration (winners and losers) — `results/variants_table.md`

Exploratory variants used: **9 / 12** (the 3 reserve slots were deliberately not used).
Macro-F1 chance = 0.25.

| Config | What | CV macro-F1 | fold SD | seed SD | Δ vs C1 | NB p |
|---|---|---|---|---|---|---|
| **C1** | logistic regression, top-500 limma genes (baseline) | **0.311** | 0.088 | 0 | – | vs chance 0.21 |
| V9 | MLP, same capacity as the GCN head | 0.302 | 0.077 | 0.040 | −0.009 | 0.65 |
| V5 | C1 after in-fold pH/RIN/PMI/age/sex residualisation | 0.301 | 0.093 | 0 | −0.009 | 0.85 |
| V8 | linear SVM | 0.292 | 0.101 | 0 | −0.018 | 0.51 |
| V3 | C1 + donor-level soft vote | 0.289 | 0.139 | 0 | −0.021 | 0.69 |
| V2 | GCN on STRING PPI graph | 0.277 | 0.067 | 0.059 | −0.033 | 0.41 |
| V1 | **reference GCN** (in-fold Pearson graph, 3×64 GCN, hybrid SAGPool + limma prior) | 0.268 | 0.066 | 0.060 | −0.042 | 0.43 |
| V7 | random forest | 0.266 | 0.082 | 0.040 | −0.045 | 0.28 |
| V6 | co-expression module eigengenes + LR | 0.266 | 0.092 | 0 | −0.045 | 0.43 |
| V4 | covariates only (confound probe) | 0.181 | 0.086 | 0 | −0.129 | 0.13 |

No variant beats the baseline (all Holm-adjusted p = 1.0).

**Graph controls** (same GCN architecture):

| Graph | V1 (co-expression) | V2 (STRING) |
|---|---|---|
| real graph | 0.268 | 0.277 |
| no graph (A = I), C3 | 0.275 (Δ −0.007, p 0.82) | 0.275 (Δ +0.002, p 0.94) |
| uniform random, same edge count, C2 | 0.240 (Δ +0.028, p 0.49) | 0.290 (Δ −0.013, p 0.75) |
| degree-preserving random, C2b | 0.276 (Δ −0.007, p 0.77) | 0.293 (Δ −0.016, p 0.58) |

Neither graph carries detectable signal beyond "no graph" or a random graph with the same
degrees. GCN **seed variance (≈ 0.06) is as large as fold variance (≈ 0.065)**, so single-seed
GCN comparisons in this regime are uninterpretable.

## 3. Final evaluation — both tracks side by side — `results/final_comparison.md`

Finals chosen by the rule fixed in PLAN §7 before any test access (LOG.md, "FINALS FROZEN"):
the best model is the configuration with the highest CV mean (**C1**), and the best simple
baseline is the next-highest linear model (**V5**). Each was evaluated **once** on the locked
test set.

| Role | Config | **Test** macro-F1 [95% CI, donor bootstrap] | Test MCC / AUC | **CV** macro-F1 mean ± SD (25 folds) |
|---|---|---|---|---|
| Best model | C1 logistic regression | **0.419** [0.249, 0.545] | 0.264 / 0.673 | **0.311 ± 0.088** |
| Best simple baseline | V5 residualised LR | **0.466** [0.270, 0.609] | 0.333 / 0.716 | **0.301 ± 0.093** |

Chance = 0.25. **A test set of 49 samples from 18 donors makes these intervals wide**
(≈ 0.3 macro-F1). Removing any single test donor moves the estimate by up to ±0.04. The
test set cannot distinguish C1 from V5, nor either from the CV estimates.

Test confusion matrices (rows = true BD, Control, MDD, SCZ; columns = predicted):

```
C1                    V5
BD       3  3  3  5   BD       3  0  5  6
Control  3  4  3  4   Control  3  6  1  4
MDD      3  1  6  2   MDD      4  0  5  3
SCZ      0  0  1  8   SCZ      0  0  0  9
```

SCZ is recalled 8/9 and 9/9, but at a cost: many BD and Control samples are *called* SCZ
(SCZ precision ≈ 0.4). BD recall is 3/14 for both models.

## 4. Controls

| Control | Result | Reading |
|---|---|---|
| C4 leakage: sample-level vs donor-grouped split | **0.653 vs 0.311 (+0.343)** | Per-sample splitting more than doubles the score. **This conflicts with the protocol's cited +0.156**: the inflation here is about 2× larger, plausibly because each donor's other two regions sit in training as near-duplicates after region-wise standardisation. |
| C8 selection bias (the supplied `gen_prep-v7.py` design: top table chosen once on all dev samples) | **0.476 (+0.165, NB p = 0.002)**, "significant" vs chance p = 0.0001 | Selecting genes before splitting manufactures an apparently significant 4-class classifier out of a dataset that, analysed properly, does not beat chance. |
| C5 SCZ vs Control | **0.721** (chance 0.50), AUC 0.82, NB p = 0.010 | Real signal, confined to SCZ. |
| PVALB positive control | lower in SCZ in **25/25** folds (mean logFC −0.62) | Reproduces the known parvalbumin deficit (Lanz et al. 2015). Preprocessing is sound. |
| C6 mixed-sex sensitivity | SEE_C6 | SEE_C6_READ |
| C7 label-permutation null (C1, rep 0) | SEE_C7 | SEE_C7_READ |
| V4 covariates only | 0.181, AUC 0.36 (below chance) | The donors are tetrad-matched on age, sex and PMI. Holding out a fold leaves the training set imbalanced in the opposite direction, so a covariate model *anti*-predicts. Matched designs can push CV scores below chance when there is no signal. |

## 5. Where the signal is (and isn't)

Across every model, per-class CV recall is about 0.55 for SCZ and about 0.2 for BD, MDD and
Control. The binary SCZ-vs-Control task is solidly above chance, and a known SCZ marker
(PVALB) behaves as expected. BD and MDD produce no separable expression signature in these
three regions at n ≈ 13 donors per class. The 4-class macro-F1 of about 0.3 is essentially
one learnable class plus three at chance.

## 6. Literature, weighed against what we observed

- **Pro-graph papers supplied** (Xing et al. 2022, MLA-GNN; Zhang et al. 2023, GAMB-GNN;
  Chereda et al. 2021): their reported gains come from cancer datasets with large effect
  sizes, and none of the passages read reports a density-matched random-graph control.
  Chereda's own Table 1 shows Graph-CNN vs random forest AUC 82.6 ± 1.3 vs 81.3 ± 1.7, with
  hyperparameters tuned on the reported folds. **On this data the graph adds nothing**,
  which conflicts with the direction of those papers' claims. It is not a contradiction
  of their specific datasets.
- **Brouard et al. 2024** ("Should we really use GNNs for transcriptomic prediction?")
  concluded that GNNs "rarely provide a real improvement". **Our results agree.**
- **Ambroise & McLachlan 2002** (selection bias): reproduced directly by C8 (+0.165).
- **Tomita et al. 2004** (pH/agonal effects): V5 residualisation neither helped nor hurt,
  so the SCZ signal isn't reducible to pH/RIN at this resolution.

## 7. Limitations and deviations (all logged)

- fRMA instead of plain RMA (rule-3 argument; PLAN §1, §10).
- Control_7 kept with expression-derived sex (user decision).
- Two container restarts during the runs; the resumable cache meant no cell was recomputed
  inconsistently.
- The NB test's ρ assumption may understate the correlation in donor-clustered folds, so
  even corrected p-values may be somewhat optimistic.
- C7 uses repetition 0 only (5 folds × 20 permutations) to stay within budget.

## 8. Versions

Linux 6.18 x86_64, Python 3.11.15, scikit-learn 1.9.1, numpy 2.4.6, scipy 1.17.1,
pandas 3.0.6, torch 2.14.0+cpu, R 4.3.3, limma 3.58.1, affy 1.80.0, frma 1.54.0,
hgu133plus2frmavecs 1.5.0, hgu133plus2.db 3.13.0. The same information is embedded in every
result JSON (`provenance`). Fold sizes are in every raw cell (`fold_sizes`) and in
`splits/cv_*.json`.
