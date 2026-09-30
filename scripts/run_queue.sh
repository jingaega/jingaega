#!/bin/bash
# Runs every registered configuration in order; each is resumable (per-cell JSON cache).
cd "$(dirname "$0")"
run() { echo "=== $(date +%H:%M:%S) $1"; python3 06_run.py "$@" 2>&1 | grep -v Warning | tail -2; }
for c in V3_logreg_donor_vote V4_covariates_only V6_modules_logreg V8_svm; do run $c; done
run V7_rf --workers 4
run V9_mlp --workers 4
for c in V1_gcn_coexpr C3_gcn_nograph C2_V1_gcn_uniform C2b_V1_gcn_degree V2_gcn_string C2_V2_gcn_uniform C2b_V2_gcn_degree; do run $c --workers 4; done
for c in C8_logreg_global_selection C5_binary_logreg C4_logreg_sample_level V5_logreg_residualised C6b_logreg_with_mixed; do run $c; done
for i in $(seq -w 0 19); do run C6c_logreg_drop$i; done
for i in $(seq -w 0 19); do run C7_logreg_perm$i; done
echo "=== ALL DONE $(date +%H:%M:%S)"
