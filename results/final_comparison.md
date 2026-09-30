# Final comparison — both evaluations side by side

Chance macro-F1 = **0.25** (4 balanced classes). Test set = 49 samples from 18 donors; **a test set this small makes the confidence interval wide**: the donor-bootstrap 95% CI spans about 0.3 macro-F1, and removing any single test donor moves the estimate by up to ±0.04. CV = donor-grouped stratified 5 folds × 5 repeats on 143 development samples / 52 donors.

| Role | Config | **Test** macro-F1 [95% CI, donor bootstrap] | Test MCC | Test AUC | **CV** macro-F1 mean ± SD (25 folds) | CV seed SD | CV MCC | CV AUC |
|---|---|---|---|---|---|---|---|---|
| Best model | `C1_logreg` | **0.419** [0.249, 0.545] | 0.264 | 0.673 | **0.311 ± 0.088** | 0.000 | 0.096 | 0.588 |
| Best simple baseline | `V5_logreg_residualised` | **0.466** [0.270, 0.609] | 0.333 | 0.716 | **0.301 ± 0.093** | 0.000 | 0.087 | 0.616 |

Both models are deterministic, so seed SD = 0 by construction.

## Best model: `C1_logreg`

| Class | Test precision | Test recall | CV precision | CV recall |
|---|---|---|---|---|
| BD | 0.33 | 0.21 | 0.26 | 0.24 |
| Control | 0.50 | 0.29 | 0.24 | 0.29 |
| MDD | 0.46 | 0.50 | 0.24 | 0.21 |
| SCZ | 0.42 | 0.89 | 0.60 | 0.56 |

Test confusion matrix (rows = true, cols = predicted; BD, Control, MDD, SCZ):

```
BD          3    3    3    5
Control     3    4    3    4
MDD         3    1    6    2
SCZ         0    0    1    8
```
CV confusion (summed over the 5 folds of a repetition, averaged over 5 repetitions):

```
BD         8.4  11.8   8.4   7.4
Control   12.8  11.2  12.8   4.2
MDD        5.8  17.2   7.2   2.8
SCZ        6.4   3.6   4.4  18.6
```

CV vs chance: Nadeau–Bengio corrected p = 0.211 (uncorrected 0.00204).

## Best simple baseline: `V5_logreg_residualised`

| Class | Test precision | Test recall | CV precision | CV recall |
|---|---|---|---|---|
| BD | 0.30 | 0.21 | 0.23 | 0.22 |
| Control | 1.00 | 0.43 | 0.23 | 0.26 |
| MDD | 0.45 | 0.42 | 0.25 | 0.25 |
| SCZ | 0.41 | 1.00 | 0.57 | 0.56 |

Test confusion matrix (rows = true, cols = predicted; BD, Control, MDD, SCZ):

```
BD          3    0    5    6
Control     3    6    1    4
MDD         4    0    5    3
SCZ         0    0    0    9
```
CV confusion (summed over the 5 folds of a repetition, averaged over 5 repetitions):

```
BD         7.6  13.4   6.8   8.2
Control   10.8  10.0  15.8   4.4
MDD        7.0  14.6   8.2   3.2
SCZ        5.8   6.0   2.8  18.4
```

CV vs chance: Nadeau–Bengio corrected p = 0.315 (uncorrected 0.0108).

