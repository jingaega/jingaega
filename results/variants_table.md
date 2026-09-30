| # | Config | Kind | Macro-F1 mean (chance) | fold SD | seed SD | MCC | AUC | Δ vs C1 | NB p (corr.) | p (uncorr.) | Holm | Folds×seeds |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
|  | C1_logreg | control | 0.311 (0.25) | 0.088 | 0.000 | 0.096 | 0.588 | – | – | – | – | 25×1 |
|  | C3_gcn_nograph | control | 0.275 (0.25) | 0.065 | 0.059 | 0.082 | 0.589 | – | – | – | – | 25×5 |
|  | C2_V1_gcn_uniform | control | 0.240 (0.25) | 0.046 | 0.044 | 0.037 | 0.556 | – | – | – | – | 25×5 |
|  | C2b_V1_gcn_degree | control | 0.276 (0.25) | 0.063 | 0.057 | 0.085 | 0.580 | – | – | – | – | 25×5 |
|  | C2_V2_gcn_uniform | control | 0.290 (0.25) | 0.057 | 0.077 | 0.090 | 0.588 | – | – | – | – | 25×5 |
|  | C2b_V2_gcn_degree | control | 0.293 (0.25) | 0.069 | 0.066 | 0.098 | 0.582 | – | – | – | – | 25×5 |
|  | C4_logreg_sample_level | control | 0.653 (0.25) | 0.094 | 0.000 | 0.554 | 0.878 | 0.343 | 0.000 | 0.000 | – | 25×1 |
|  | C5_binary_logreg | control | 0.721 (0.50) | 0.147 | 0.000 | 0.490 | 0.821 | – | – | – | – | 25×1 |
|  | C6b_logreg_with_mixed | control | 0.341 (0.25) | 0.104 | 0.000 | 0.142 | 0.613 | 0.031 | 0.399 | 0.037 | – | 25×1 |
|  | C8_logreg_global_selection | control | 0.476 (0.25) | 0.092 | 0.000 | 0.323 | 0.747 | 0.165 | 0.002 | 0.000 | – | 25×1 |
| 1 | V1_gcn_coexpr | variant | 0.268 (0.25) | 0.066 | 0.060 | 0.069 | 0.578 | -0.042 | 0.425 | 0.039 | 1.000 | 25×5 |
| 2 | V2_gcn_string | variant | 0.277 (0.25) | 0.067 | 0.059 | 0.082 | 0.580 | -0.033 | 0.407 | 0.032 | 1.000 | 25×5 |
| 3 | V3_logreg_donor_vote | variant | 0.289 (0.25) | 0.139 | 0.000 | 0.098 | 0.599 | -0.021 | 0.685 | 0.280 | 1.000 | 25×1 |
| 4 | V4_covariates_only | variant | 0.181 (0.25) | 0.086 | 0.000 | -0.070 | 0.359 | -0.129 | 0.126 | 0.000 | 1.000 | 25×1 |
| 5 | V5_logreg_residualised | variant | 0.301 (0.25) | 0.093 | 0.000 | 0.087 | 0.616 | -0.009 | 0.847 | 0.604 | 1.000 | 25×1 |
| 6 | V6_modules_logreg | variant | 0.266 (0.25) | 0.092 | 0.000 | 0.029 | 0.549 | -0.045 | 0.429 | 0.040 | 1.000 | 25×1 |
| 7 | V7_rf | variant | 0.266 (0.25) | 0.082 | 0.040 | 0.048 | 0.575 | -0.045 | 0.281 | 0.007 | 1.000 | 25×5 |
| 8 | V8_svm | variant | 0.292 (0.25) | 0.101 | 0.000 | 0.074 | 0.577 | -0.018 | 0.514 | 0.087 | 1.000 | 25×1 |
| 9 | V9_mlp | variant | 0.302 (0.25) | 0.077 | 0.040 | 0.092 | 0.567 | -0.009 | 0.649 | 0.227 | 1.000 | 25×5 |

Running exploratory-variant count: **9 / 12**. Chance macro-F1 = 0.25 (4-class), 0.50 (binary).
NB p = Nadeau–Bengio corrected resampled t-test on 25 paired fold-level differences (seed-averaged within fold); assumes between-fold correlation ρ = n_test/(n_train+n_test).
