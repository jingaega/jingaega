"""Development CV folds (rule 2, PLAN §5). Saved once, reused by every configuration."""
import json
import os

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedGroupKFold, StratifiedKFold

from common import CLASSES, PROC, SPLITS, provenance, write_json

meta = pd.read_csv(os.path.join(PROC, "samples.tsv"), sep="\t", index_col=0)
lock = json.load(open(os.path.join(SPLITS, "test_lock.json")))
dev = meta[~meta.index.isin(lock["test_gsm"])]
clean = dev[~dev.mixed_sex_donor]
mixed = dev[dev.mixed_sex_donor]


def fold_record(sub, tr_idx, va_idx):
    tr, va = sub.iloc[tr_idx], sub.iloc[va_idx]
    assert not set(tr.donor) & set(va.donor) or grouped is False
    return dict(train=list(tr.index), val=list(va.index),
                sizes=dict(train_samples=len(tr), val_samples=len(va),
                           train_donors=tr.donor.nunique(), val_donors=va.donor.nunique(),
                           val_per_class={c: int((va.diagnosis == c).sum()) for c in CLASSES},
                           train_per_class={c: int((tr.diagnosis == c).sum()) for c in CLASSES}))


def make(sub, grouped_split, name):
    global grouped
    grouped = grouped_split
    folds = []
    for rep in range(5):
        if grouped_split:
            cv = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=rep)
            it = cv.split(sub, sub.diagnosis, groups=sub.donor)
        else:
            cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=rep)
            it = cv.split(sub, sub.diagnosis)
        for k, (tr, va) in enumerate(it):
            r = fold_record(sub, tr, va)
            r.update(rep=rep, fold=k)
            folds.append(r)
    out = dict(name=name, grouped_by_donor=grouped_split, n_samples=len(sub), n_donors=sub.donor.nunique(),
               folds=folds, provenance=provenance())
    write_json(os.path.join(SPLITS, f"cv_{name}.json"), out)
    vs = [f["sizes"]["val_samples"] for f in folds]
    print(name, "samples", len(sub), "donors", sub.donor.nunique(), "val size range", min(vs), max(vs))


make(clean, True, "multiclass")                                       # primary
make(clean[clean.diagnosis.isin(["SCZ", "Control"])], True, "binary_scz_ctrl")   # C5
make(clean, False, "multiclass_sample_level")                          # C4 leakage demo (NOT for model comparison)

# C6: training-set modifications on the SAME validation folds as the primary.
rng = np.random.default_rng(20260930)
dons = clean.groupby("donor").diagnosis.first()
draws = []
for _ in range(20):
    pick = []
    for dx, n in (("BD", 1), ("MDD", 2), ("SCZ", 2)):
        pick += list(rng.choice(sorted(dons[dons == dx].index), size=n, replace=False))
    draws.append(sorted(pick))
write_json(os.path.join(SPLITS, "c6_sensitivity.json"), dict(
    design="Validation folds identical to cv_multiclass. (a) primary train; (b) train + all samples of the 5 "
           "mixed-sex donors; (c_i) (b) with 5 random clean donors (1 BD, 2 MDD, 2 SCZ) removed from TRAINING only.",
    mixed_sex_gsm=list(mixed.index), mixed_sex_donors=sorted(mixed.donor.unique()),
    random_drop_draws=draws, seed=20260930, provenance=provenance()))
print("mixed-sex dev samples", len(mixed), "; random-drop draws", len(draws))
