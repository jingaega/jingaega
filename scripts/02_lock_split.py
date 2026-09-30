"""Lock the final test set (rule 1, PLAN §2). Run once; refuses to overwrite."""
import json
import os
import sys

import pandas as pd
import sklearn
from sklearn.model_selection import train_test_split

from common import CLASSES, MIXED_SEX, PROC, RESULTS, SPLITS, load_metadata, provenance, sha256, write_json

OUT = os.path.join(SPLITS, "test_lock.json")
if os.path.exists(OUT):
    sys.exit(f"{OUT} exists; the test split is locked and must not be regenerated.")

meta = load_metadata()
qa = json.load(open(os.path.join(RESULTS, "qa_sex.json")))
assert qa["matches_expected"], "sex QA list differs from protocol list"

# Sample table with expression-derived sex (user decision 2026-09-30: Control_7 kept, sex from expression)
sex_expr = {r["gsm"]: r["expr_sex"] for r in qa["samples"]}
meta["sex_expr"] = [sex_expr.get(g, s) for g, s in zip(meta.index, meta.sex)]
meta["mixed_sex_donor"] = meta.donor.isin(MIXED_SEX)
meta.to_csv(os.path.join(PROC, "samples.tsv"), sep="\t")

donors = meta.groupby("donor").agg(diagnosis=("diagnosis", "first"), n=("region", "size")).reset_index()
pool = donors[~donors.donor.isin(MIXED_SEX)].reset_index(drop=True)
RANDOM_STATE = 20260930
dev_d, test_d = train_test_split(pool, test_size=0.25, stratify=pool.diagnosis, random_state=RANDOM_STATE)

test_s = meta[meta.donor.isin(test_d.donor)]
dev_s = meta[meta.donor.isin(dev_d.donor) | meta.mixed_sex_donor]
assert not set(test_s.index) & set(dev_s.index)
assert len(test_s) + len(dev_s) == len(meta)

lock = dict(
    created="2026-09-30", rule="protocol rule 1; PLAN.md §2",
    sklearn_version=sklearn.__version__, random_state=RANDOM_STATE,
    call="train_test_split(donors_without_mixed_sex, test_size=0.25, stratify=diagnosis)",
    pool_donors=int(len(pool)), excluded_mixed_sex_donors=MIXED_SEX,
    test_donors=sorted(test_d.donor), test_gsm=sorted(test_s.index),
    test_counts={c: dict(donors=int((test_d.diagnosis == c).sum()), samples=int((test_s.diagnosis == c).sum()))
                 for c in CLASSES},
    dev_counts_clean={c: dict(donors=int((dev_d.diagnosis == c).sum()),
                              samples=int(((dev_s.diagnosis == c) & ~dev_s.mixed_sex_donor).sum()))
                      for c in CLASSES},
    dev_mixed_sex_samples=int(dev_s.mixed_sex_donor.sum()),
    provenance=provenance(),
)
write_json(OUT, lock)
digest = sha256(OUT)
with open(OUT + ".sha256", "w") as f:
    f.write(digest + "  test_lock.json\n")
print(json.dumps({k: lock[k] for k in ("test_counts", "dev_counts_clean", "dev_mixed_sex_samples", "sklearn_version")}, indent=1))
print("test donors", len(test_d), "test samples", len(test_s), "sha256", digest)
