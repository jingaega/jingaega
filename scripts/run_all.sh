#!/bin/bash
# Single resumable queue (every cell cached as JSON; re-running skips finished cells).
cd "$(dirname "$0")"
run() { echo "=== $(date +%H:%M:%S) $1"; python3 06_run.py "$@" 2>&1 | grep -v Warning | tail -2; }
for c in C2b_V2_gcn_degree V1_gcn_coexpr C3_gcn_nograph C2_V1_gcn_uniform C2b_V1_gcn_degree V2_gcn_string C2_V2_gcn_uniform; do run $c --workers 4; done
for c in C8_logreg_global_selection C5_binary_logreg C4_logreg_sample_level V5_logreg_residualised C6b_logreg_with_mixed; do run $c; done
for i in $(seq -w 0 19); do run C6c_logreg_drop$i; done
for i in $(seq -w 0 19); do run C7_logreg_perm$i; done
echo "=== ALL DONE $(date +%H:%M:%S)"
