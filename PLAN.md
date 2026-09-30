# PLAN — Psychiatric diagnosis from postmortem brain expression (GSE53987)

Status: **DRAFT, awaiting approval.** Nothing has been split, normalised, trained or evaluated.
Written 2026-09-30.

---

## 0. Environment (verified, not assumed)

| Item | Version | How verified |
|---|---|---|
| Platform | Linux-6.18.44-fc-v50-x86_64, glibc 2.39, 4 CPU, 15 GB RAM, no GPU | `platform`, `nproc`, `free` |
| Python | 3.11.15 | `python3 --version` |
| scikit-learn | **1.9.1** | import |
| numpy / scipy / pandas / statsmodels | 2.4.6 / 1.17.1 / 3.0.6 / 0.15.0 | import |
| PyTorch | 2.14.0+cpu | import |
| R | 4.3.3 (2024-02-29) | `Rscript` |
| limma | 3.58.1 | `packageVersion` |
| affy | 1.80.0 | `packageVersion` |
| hgu133plus2cdf | 2.18.0 (Bioconductor source tarball, md5 284fef2f…de3 matches bioconda recipe) | md5 |
| frma, hgu133plus2frmavecs, hgu133plus2.db | Bioconductor 3.18 — **installing now; plan depends on them, see §1 fallback** | pending |

PyTorch Geometric is not used; the GCN and pooling are written in plain PyTorch with dense
adjacency (≤ 1000 nodes), which removes one dependency whose version would affect results.

Every result JSON will carry a `provenance` block: platform, Python, sklearn, numpy, torch,
R, limma, affy, frma versions, fold sizes (train/val samples and donors, per class), seeds,
git commit hash.

## 1. Data starting point

**Raw CEL files, processed in R.** Downloaded from NCBI GEO on 2026-09-30:

| File | md5 |
|---|---|
| `GSE53987_RAW.tar` (205 CEL.gz) | 67576fd879ea2ef079aec41f723b7f09 |
| `GSE53987_series_matrix.txt.gz` | 26b164ffb27bb115709a06b6f1e9d635 |

Metadata (parsed from the series matrix): 205 samples, 75 donors (19 BD, 19 Control, 19 SCZ,
18 MDD); per class 52 BD / 55 Control / 50 MDD / 48 SCZ samples; 3 regions (hip, pfc, str).
Donors contribute 3 samples (58 donors), 2 (14) or 1 (3). Donor id = diagnosis + the number in
the sample title (e.g. `bipolar_hip_18`, `bipolar_pfc_18` → `BD_18`). Available covariates:
age, sex, race, PMI, pH, RIN.

**Normalisation: frozen RMA (fRMA), not plain RMA. This is a deliberate choice, argued here
(rule 12).** Plain RMA is not a per-sample transform: quantile normalisation fits a reference
distribution from all arrays, and median polish fits probe effects from all arrays. Running it
on all 205 arrays (which is what the GEO submitters did, "Brain regions were normalized
separately") lets the locked test arrays shape every training feature, which violates rule 3.
It carries no label information, so the leak is probably small, but I can't show it's zero
without running the comparison. fRMA (McCall, Bolstad & Irizarry, 2010) is RMA with the quantile reference
and probe effects *frozen* from thousands of external GEO arrays, so each array is normalised
independently and no data-dependent step touches held-out arrays. It is still background
correction + quantile normalisation + robust summarisation, i.e. the RMA family.
I have not read the fRMA paper in this session; the citation is for provenance only, not
for any claim.

- Probe sets → genes: `hgu133plus2.db` symbols; drop probe sets with no or multiple symbols;
  average probe sets per gene. This is a per-sample operation with a fixed annotation, so it
  involves no fitting.
- Remove sex-chromosome-specific genes (chrY genes and XIST) from the feature set. If the
  diagnosis groups differ in sex ratio, these genes are a shortcut that has nothing to do with
  psychiatric disease. The sex QA in §3 uses them, but no model sees them.
- **Fallback, decided now:** if `frma` fails to install or run on R 4.3, I will run RMA
  (`affy::rma`) *separately on the development arrays and on the test arrays*. That is still
  leak-free, at the cost of the test arrays being normalised as their own batch. I will log
  which path was taken. I will not silently use the submitters' all-array matrix.
- The submitters' series matrix is used only to cross-check sample identity: the per-array
  correlation between my fRMA values and theirs should be high. It is not a modelling input.

## 2. Locking the test set (done first, before anything else is computed on expression)

1. Build a donor table: donor id, diagnosis, sample GSM ids.
2. Remove the 5 mixed-sex donors from the pool (§3), leaving **70 donors**
   (18 BD, 19 Control, 16 MDD, 17 SCZ).
3. `sklearn.model_selection.train_test_split(donors, test_size=0.25, stratify=diagnosis,
   random_state=20260930)` at the **donor** level. That gives about 18 test donors
   (about 4–5 per class) and about 49 test samples. The exact counts per class (donors and
   samples) and the sklearn version go into `splits/test_lock.json`.
4. Write `splits/test_lock.json` with test GSM ids, donor ids, per-class counts, sklearn
   version, random_state, and the SHA-256 of the file itself (recorded in LOG.md and
   committed).
5. Every downstream loader asserts that no test GSM id is in its input. The test arrays are
   not even fRMA-processed until the single final evaluation. fRMA is per-array, so this
   costs nothing.

Development set = 52 clean donors (about 146 samples) + the 5 quarantined mixed-sex donors
(used only in the sensitivity run).

## 3. The five mixed-sex donors (decision written before modelling)

Donors BD_18, MDD_11, MDD_2, SCZ_4, SCZ_7 have samples of both sexes by Y-gene/XIST
expression. Metadata sex is consistent within every donor, so at least one array per donor
is not from the recorded person.

**Decision: exclude all samples of these five donors from the primary analysis and from the
test set.**

Reasons:
1. Keeping them risks **cross-fold donor leakage**. A foreign array really belongs to one of
   3–4 other donors, who may sit in a different fold. That is exactly the per-sample
   splitting leak that rule 5 forbids, and it can't be fixed because the owner is
   unidentifiable.
2. Keeping them adds **label noise**. The foreign array may come from a donor with a
   different diagnosis.
3. Dropping only the sex-discordant array isn't well defined. With 2 arrays of different
   sexes, it isn't clear which one is foreign.

**Pre-check (QA, not modelling):** before the split, compute a sex score (mean of RPS4Y1,
DDX3Y, KDM5D, USP9Y, EIF1AY minus XIST) per array on the submitters' matrix. Confirm that
exactly these five donors are discordant and that no other donor disagrees with its
metadata. If the list differs, I stop and report before splitting. This step uses only sex
genes and metadata sex, never diagnosis.

**Sensitivity check (required control):** on the development CV, compare
(a) primary (52 clean dev donors),
(b) primary + the 5 mixed-sex donors,
(c) 20 draws of (b) minus 5 random clean dev donors with the same diagnosis composition
(1 BD, 2 MDD, 2 SCZ).
If (a) − (b) sits inside the (c) − (b) distribution, the effect of the mixed-sex donors is
indistinguishable from losing any 5 donors.

## 4. Fixed pipeline shared by every model (equal treatment)

Inside each training fold only, in this order:
1. Low-expression filter: keep genes whose training-set mean is above the 25th percentile
   of training-set gene means.
2. **Region adjustment:** per region, centre and scale each gene with the training samples'
   mean and SD for that region. Held-out samples use the training statistics of their own
   region (region is metadata known at prediction time, not a label). Region is the largest
   variance component in multi-region data; without this, a classifier has to learn
   diagnosis within three shifted clouds.
3. **Differential expression:** limma, design `~ 0 + diagnosis + region`,
   `duplicateCorrelation(block = donor)` because donors contribute up to 3 correlated
   samples, then `eBayes` and moderated F across the 4 diagnoses. Keep the top **k = 500**
   genes. Called from Python through `Rscript` on the training fold; the result is cached per
   fold under `cache/`.
4. Model-specific fit on those k genes.

Every model sees identical folds, identical k, identical preprocessing. Hyperparameters are
**fixed a priori for all models; there is no hyperparameter search.** Nested tuning of a GCN
would cost about 12× per configuration, which is over the budget in rule 11. Tuning only the
cheap models would favour them unfairly. So nobody is tuned, and a model that needs tuning
to work has failed under these conditions. Neural models use a fixed epoch count, no early
stopping, and no checkpoint selection (rule 4).

Prediction unit: the **sample**. Macro-F1 is computed over held-out samples. Variant V3
explores donor-level pooling of predictions.

## 5. Evaluation machinery

- **Development CV:** `StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=r)`,
  r ∈ {0,…,4}, groups = donor, y = diagnosis. That gives 25 folds. The fold assignments are
  saved to `splits/cv_folds.json` once and reused by every configuration.
- **Seeds:** stochastic models use 5 seeds per fold. Deterministic models (logistic
  regression) are run once per fold, and seed variance is reported as 0 by construction.
- **Variance reporting:** fold variance = SD of the 25 per-fold seed-means; seed variance =
  mean over folds of the within-fold SD across seeds.
- **Comparisons (rule 8):** Nadeau–Bengio corrected resampled t-test on the 25 paired
  fold-level differences (seed-averaged within fold, never repetition-averaged):
  t = mean(d) / sqrt((1/25 + n_test/n_train) · var(d)), df = 24. It assumes the correlation
  between fold estimates is ρ = n_test/(n_train + n_test) ≈ 0.2, a heuristic from the
  train/test overlap, not a measured quantity. The uncorrected paired t p is reported
  beside it. Vs chance: the same corrected test on (macro-F1 − 0.25). Across the family of
  "variant vs baseline" comparisons I also report Holm-adjusted p, as a secondary number.
- **Metrics per configuration:** macro-F1 (primary; chance 0.25 stated beside each),
  MCC, per-class precision and recall, pooled confusion matrix (summed over folds per
  repetition, then averaged), one-vs-rest macro AUC.
- **Final test set:** evaluated once per final configuration. 95% CI by **donor-level**
  cluster bootstrap (10,000 resamples), because samples within a donor are correlated and the
  effective n is about 18 donors, not about 49 samples. The interval will be wide: at this
  size a CI of roughly ±0.15 on macro-F1 is expected, and the report will say so.

## 6. Required controls (outside the 12-variant budget)

| Id | Control | What it separates |
|---|---|---|
| C1 | **Baseline:** §4 pipeline → L2 multinomial logistic regression (C = 1.0, lbfgs, balanced class weights) | "simplest sensible classifier" |
| C2 | Density-matched **random graph**, same edge count, uniform random topology, rebuilt per fold and per seed, for every graph variant | this graph vs any graph of this density |
| C3 | **No-graph ablation**: identical GNN with A = I (self-loops only), for every graph architecture | graph vs the same network without message passing |
| C4 | **Leakage demo:** C1 with sample-level `StratifiedKFold` (5×5) vs donor-grouped | size of the per-sample inflation (claimed ≈ +0.156) |
| C5 | **Binary SCZ vs Control** with C1, and with the final best model | whether signal exists at all |
| C6 | **Mixed-sex sensitivity** (§3), run on C1 and on the final best model | composition vs real effect |
| C7 | **Label-permutation null** for C1 (labels shuffled at donor level inside training folds, 20 permutations) | empirical chance distribution, a check on the 0.25 figure |

C7 is an addition beyond the list; it is cheap and it doesn't count as a variant.

## 7. Exploratory variants (budget 12; 9 pre-registered now, 3 held in reserve)

Each variant is compared to C1 with the corrected test. A reserve variant gets written into
LOG.md with its full pre-registration before it runs.

| # | Variant | Why it could help (mechanism / source) | Counts as failure if |
|---|---|---|---|
| V1 | **Reference config (rule 10):** in-fold Pearson co-expression graph over the k=500 limma genes (keep the top 2,500 \|r\| edges, i.e. mean degree 10); 3-layer GCN, 64 hidden; hybrid pooling Q = λ·P + (1−λ)·Z with Z = SAGPool attention score and P = min–max-scaled −log10 limma p per gene, λ = 0.5, keep 50% of nodes; mean+max readout → linear; Adam lr 1e-3, wd 5e-4, 200 epochs, batch 16. With C2 and C3. | Co-expression modules may carry disease signal that per-gene models miss (argued by Xing et al. 2022, MLA-GNN); the pooling formula is Zhang et al. 2023 (GAMB-GNN) eq. 17. | macro-F1 not above C1 with corrected p < 0.05, **or** not above its own C2 random graph |
| V2 | Same GCN on a **prior-knowledge PPI graph** (STRING v12 human, combined score ≥ 700, restricted to the k=500 genes) with C2 and C3. C3 is shared with V1 because the architecture is identical. | A graph fitted on ~110 training samples is noisy; a prior graph isn't estimated from these samples. Chereda et al. 2021 use a PPI prior (HPRD). | same as V1 |
| V3 | C1 + **donor-level soft vote:** each held-out sample gets the mean predicted probability over all held-out samples of its donor | Diagnosis is a donor property. If region-specific noise is partly independent, averaging 2–3 arrays reduces variance. Uses no labels, and a donor's samples are always held out together. | Δ macro-F1 ≤ 0 or corrected p ≥ 0.05 |
| V4 | **Covariate-only** model: LR on age, sex, PMI, pH, RIN (in-fold scaled) | Confound probe, not an attempt to win. Brain pH and RIN are known to shape postmortem expression and may differ by diagnosis. | (interpretive) If V4 ≥ C1, the expression signal may be confounded. That gets reported whatever happens. |
| V5 | C1 after **in-fold covariate residualisation**: per gene, regress out pH, RIN, PMI, age, sex fitted on training samples; apply the training coefficients to held-out samples | Removes technical/agonal variance that isn't diagnosis, raising signal-to-noise for the linear model | Δ ≤ 0 or p ≥ 0.05 |
| V6 | **Co-expression module features:** in-fold hierarchical clustering (1 − \|r\|, average linkage) of the k genes into 30 modules; module eigengene (PC1) → LR | Tests the graph papers' key idea (module-level signal) without a GNN. If V6 ≈ V1, the benefit, if any, is aggregation rather than message passing. | Δ ≤ 0 or p ≥ 0.05 |
| V7 | **Random forest** (500 trees, max_features = sqrt, balanced) on the same k genes, 5 seeds | Standard non-linear baseline in Chereda 2021 (Table 1) and Xing 2022 (Table 2); captures interactions | Δ ≤ 0 or p ≥ 0.05 |
| V8 | **Linear SVM** (C = 1.0, balanced) on the same k genes | Margin-based linear model, often strongest at p ≫ n | Δ ≤ 0 or p ≥ 0.05 |
| V9 | **MLP** with the same capacity as the GCN head (2 × 64, dropout 0.3, same optimiser and epochs), 5 seeds | Isolates "neural net" from "graph": if V9 ≈ V1, the graph adds nothing | Δ ≤ 0 or p ≥ 0.05 |
| R1–R3 | reserve | pre-registered in LOG.md before running | — |

**Choosing the finals (rule 9), fixed now:**
- *Best model* = the exploratory variant or baseline with the highest development CV mean
  macro-F1. V4, the confound probe, is not eligible.
- *Best simple baseline* = the highest among C1, V3, V5, V8 (linear models without learned
  structure). If the best model is itself one of these, the baseline is the next highest
  among them.
- Ties are broken by the mean; no re-running.

## 8. Literature, as read (passages actually checked)

- **Chereda et al. 2021 (Genome Med 13:42):** Graph-CNN on a PPI prior, breast-cancer
  metastasis, 969 patients, 10-fold CV. Table 1: Graph-CNN AUC 82.57 ± 1.25 vs Random
  Forest (no prior) 81.27 ± 1.66, with overlapping standard errors. Hyperparameters were
  "tweaked manually on this 10-fold cross validation", the same folds that were reported.
  This is closer to evidence that the graph adds little than that it helps, and it comes
  from a sample about 5× larger than ours.
- **Xing et al. 2022 (Bioinformatics 38:2178, MLA-GNN):** WGCNA co-expression graph built on
  training data only (good practice). The edge threshold `adj_thresh` was "optimized by an
  automated machine learning algorithm", and I couldn't tell on which data. Reports gains
  over SVM/RF on TCGA glioma. In the sections I read, I found no random-graph control.
  Their graph comparison is against other *real* graphs (HumanNet/PPI).
- **Zhang et al. 2023 (Chemom Intell Lab Syst 232:104713, GAMB-GNN):** hybrid pooling
  Q = λP + (1−λ)Z (eq. 17), the source of the rule-10 reference pooling. Evaluation uses
  random 0.8/0.1/0.1 splits × 10 trials with a Wilcoxon rank-sum test across those
  overlapping trials. The trials are not independent, which is the misapplication rule 8
  guards against. Their datasets are cancer microarrays, where effect sizes are far larger
  than in psychiatric postmortem brain.
- **Parwati et al. 2024 (ICoDSA):** GCN for drug–drug interactions on warfarin from
  molecular graphs. It isn't about expression data, so I don't cite it for any claim here.

Expectation, stated in advance: psychiatric case–control differences in postmortem brain are
small, and n = 52 development donors is tiny. The most likely outcome is every model near
0.25–0.35 macro-F1 with overlapping intervals. SCZ vs Control (C5) is the best chance to see
any signal.

## 9. Runtime estimate (rule 11)

Anchor: 1.5 h per GNN configuration (25 folds × 5 seeds), from the protocol. I will time the
first V1 fold and re-estimate in LOG.md before continuing. If the estimate exceeds 20 h, I
will cut GNN seeds for *controls only* from 5 to 3, and say so in the log.

| Block | h |
|---|---|
| fRMA on 156 dev arrays (+ ~49 test arrays at the end) | 0.5 |
| limma per fold: 25 folds × ~4 sample configurations, cached | 0.5 |
| C1, C4, C7, V3–V9 (non-GNN; RF and MLP × 5 seeds) | 1.5 |
| V1 + C2(V1) + C3 | 4.5 |
| V2 + C2(V2) (C3 shared) | 3.0 |
| C5 binary SCZ/Ctrl: C1 + best model (worst case a GNN at ~½ the data) | 0.8 |
| C6 sensitivity on the best model if it's a GNN: (b) + 3 random-drop draws at 3 seeds | ~3.6 |
| Final test evaluation (2 configurations) | 0.2 |
| **Total (worst case, GNN is best)** | **~14.6 h** |

Caching: each (configuration, repetition, fold, seed) writes its own JSON under
`results/raw/`. A run skips cells that already exist, so an interrupt costs at most one
fold.

## 10. Escalations (rule 12). I will follow the rules as written; these are for the record.

1. **RMA vs rule 3.** Plain RMA on all arrays conflicts with rule 3; I use fRMA (§1). If you
   want plain whole-dataset RMA for comparability with other studies, it can be run as a
   reserve variant, and the difference would measure the transductive-normalisation effect.
2. **Sample-level macro-F1 with correlated samples.** About 49 test samples come from about
   18 donors. A sample-level CI would be too narrow, so I bootstrap by donor. The Nadeau–
   Bengio correction assumes ρ = n_test/n_total. With donor clustering, the true
   between-fold correlation may be higher, so even the corrected test may be somewhat
   anti-conservative.
3. **Mixed-sex sensitivity direction.** The protocol phrases C6 as "a run without the five".
   Since my primary analysis already excludes them, C6 compares with vs without, plus the
   random-drop draws. The comparison is the same, but the primary/sensitivity roles are
   swapped.

## 11. Deliverables

`scripts/` (numbered, each writes JSON), `splits/`, `results/raw/*.json`,
`results/summary.json`, `results/variants_table.md` (every variant, running count),
`results/confusion/*.json`, `results/final_comparison.md` (test-set evaluation and 5×5 CV
mean ± SD side by side, for the best model and the best simple baseline), `LOG.md`.
Raw CEL data is not committed (1 GB); the md5 values above identify it.
