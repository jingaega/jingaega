#!/bin/bash
# GCN block (re-run after the pickling fix). Waits for the main queue to finish first.
cd "$(dirname "$0")"
while kill -0 16235 2>/dev/null; do sleep 30; done
run() { echo "=== $(date +%H:%M:%S) $1"; python3 06_run.py "$@" 2>&1 | grep -v Warning | tail -2; }
for c in V1_gcn_coexpr C3_gcn_nograph C2_V1_gcn_uniform C2b_V1_gcn_degree V2_gcn_string C2_V2_gcn_uniform C2b_V2_gcn_degree; do run $c --workers 4; done
echo "=== GCN ALL DONE $(date +%H:%M:%S)"
