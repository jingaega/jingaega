"""Shared paths, metadata parsing and provenance for the GSE53987 study."""
import gzip
import hashlib
import json
import os
import platform
import subprocess
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw")
PROC = os.path.join(ROOT, "data", "processed")
SPLITS = os.path.join(ROOT, "splits")
RESULTS = os.path.join(ROOT, "results")
CACHE = os.path.join(ROOT, "cache")
SERIES = os.path.join(RAW, "GSE53987_series_matrix.txt.gz")

CLASSES = ["BD", "Control", "MDD", "SCZ"]
DX_MAP = {"bipolar disorder": "BD", "control": "Control",
          "major depressive disorder": "MDD", "schizophrenia": "SCZ"}
MIXED_SEX = ["BD_18", "MDD_11", "MDD_2", "SCZ_4", "SCZ_7"]
SEX_GENES_Y = ["RPS4Y1", "DDX3Y", "KDM5D", "USP9Y", "EIF1AY"]


def load_metadata():
    """Sample table from the series-matrix header (no expression values)."""
    rows = {}
    with gzip.open(SERIES, "rt") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Sample_"):
                key = line.split("\t")[0]
                vals = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
                rows.setdefault(key, []).append(vals)
    gsm = rows["!Sample_geo_accession"][0]
    titles = rows["!Sample_title"][0]
    cel = [os.path.basename(u) for u in rows["!Sample_supplementary_file"][0]]
    recs = []
    for i in range(len(gsm)):
        d = {c[i].split(": ", 1)[0]: c[i].split(": ", 1)[1]
             for c in rows["!Sample_characteristics_ch1"] if ": " in c[i]}
        dx = DX_MAP[d["disease state"]]
        num = titles[i].rsplit("_", 1)[1]
        region = titles[i].rsplit("_", 2)[1]
        recs.append(dict(gsm=gsm[i], title=titles[i], cel=cel[i], diagnosis=dx,
                         donor=f"{dx}_{num}" if dx != "Control" else f"Control_{num}",
                         region=region, sex=d["gender"], age=float(d["age"]),
                         race=d["race"], pmi=float(d["pmi"]), ph=float(d["ph"]),
                         rin=float(d["rin"])))
    meta = pd.DataFrame(recs).set_index("gsm")
    # The protocol names donors BD_18, MDD_11, SCZ_4 ...; our ids follow the same pattern.
    return meta


def load_series_values():
    """Submitters' expression matrix (probes x samples), used ONLY for QA."""
    return pd.read_csv(SERIES, sep="\t", comment="!", index_col=0)


def _cmd(args):
    try:
        return subprocess.run(args, capture_output=True, text=True, timeout=60).stdout.strip()
    except Exception as e:  # pragma: no cover
        return f"unavailable ({e})"


_PROV = None


def provenance():
    global _PROV
    if _PROV is None:
        import scipy
        import sklearn
        try:
            import torch
            tv = torch.__version__
        except ImportError:
            tv = None
        rv = _cmd(["Rscript", "-e", 'cat(R.version.string, "|limma", as.character(packageVersion("limma")), '
                   '"|affy", as.character(packageVersion("affy")), "|frma", as.character(packageVersion("frma")))'])
        _PROV = dict(platform=platform.platform(), python=platform.python_version(),
                     sklearn=sklearn.__version__, numpy=np.__version__, scipy=scipy.__version__,
                     pandas=pd.__version__, torch=tv, R=rv,
                     git_commit=_cmd(["git", "-C", ROOT, "rev-parse", "HEAD"]))
    return dict(_PROV)


def write_json(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(obj, f, indent=1, default=_default)
    os.replace(tmp, path)


def _default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(type(o))


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        h.update(f.read())
    return h.hexdigest()
