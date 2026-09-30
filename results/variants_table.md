| # | Config | Kind | Macro-F1 mean (chance) | fold SD | seed SD | MCC | AUC | Δ vs C1 | NB p (corr.) | p (uncorr.) | Holm | Folds×seeds |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
|  | C1_logreg | control | 0.311 (0.25) | 0.088 | 0.000 | 0.096 | 0.588 | – | – | – | – | 25×1 |
| 1 | V3_logreg_donor_vote | variant | 0.289 (0.25) | 0.139 | 0.000 | 0.098 | 0.599 | -0.021 | 0.685 | 0.280 | 1.000 | 25×1 |
| 2 | V4_covariates_only | variant | 0.181 (0.25) | 0.086 | 0.000 | -0.070 | 0.359 | -0.129 | 0.126 | 0.000 | 0.503 | 25×1 |
| 3 | V6_modules_logreg | variant | 0.266 (0.25) | 0.092 | 0.000 | 0.029 | 0.549 | -0.045 | 0.429 | 0.040 | 1.000 | 25×1 |
| 4 | V7_rf | variant | 0.248 (0.25) | 0.072 | 0.038 | 0.027 | 0.555 | – | – | – | – | 12×5 |
| 5 | V8_svm | variant | 0.292 (0.25) | 0.101 | 0.000 | 0.074 | 0.577 | -0.018 | 0.514 | 0.087 | 1.000 | 25×1 |

Running exploratory-variant count: **5 / 12**. Chance macro-F1 = 0.25 (4-class), 0.50 (binary).
NB p = Nadeau–Bengio corrected resampled t-test on 25 paired fold-level differences (seed-averaged within fold); assumes between-fold correlation ρ = n_test/(n_train+n_test).
