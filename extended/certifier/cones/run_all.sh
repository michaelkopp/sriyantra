#!/bin/bash
# The whole configuration space, without the chart A cells.
# usage: ./run_all.sh [WORKERS]    (resumable: run it again after an interruption)
#   BUDGET=60000 REDO=1 ./run_all.sh 10      runs again the cones left open (larger budget per cone); with REDO=1
#                                            the open cones are shared out through a claim queue, so all workers stay busy
#   NOREGIONS=0 ...                          also lets a box close when it lies in the region D or in a chart B
#                                            region (certified separately by fanD.js and fanRBB.js); default 1 (not used)
# Cones closed in cloud_all_*.jsonl (the run from the other end) are skipped and counted in the summary.
# Needs a python3 with numpy and scipy; the first one found below is used (or set PYTHON=/path/to/python3).
N=${1:-10}
BUDGET=${BUDGET:-5000}
NOREG=${NOREGIONS:-1}
export SKIPFILES=$(ls cloud_all_*.jsonl 2>/dev/null | paste -sd, -)
PY=${PYTHON:-}
if [ -z "$PY" ]; then
  for c in python3 /opt/homebrew/bin/python3 /usr/local/bin/python3 \
           /Library/Frameworks/Python.framework/Versions/Current/bin/python3 \
           "$HOME/miniconda3/bin/python3" "$HOME/anaconda3/bin/python3" "$HOME/miniforge3/bin/python3" \
           /opt/miniconda3/bin/python3 /opt/anaconda3/bin/python3 /opt/homebrew/Caskroom/miniconda/base/bin/python3 \
           /usr/bin/python3; do
    p=$(command -v "$c" 2>/dev/null) || continue
    if "$p" -c 'import numpy, scipy' 2>/dev/null; then PY="$p"; break; fi
  done
fi
if [ -z "$PY" ]; then echo "No python3 with numpy and scipy found; run with PYTHON=/path/to/python3 ./run_all.sh"; exit 2; fi
echo "python: $PY   workers: $N   budget per cone: $BUDGET   NOREGIONS=$NOREG"
"$PY" test_cone.py || exit 1
"$PY" test_mvf.py 1 20 || exit 1
"$PY" test_contract.py 1 6 || exit 1
"$PY" test_split.py || exit 1
"$PY" test_grad.py 1 20 || exit 1
CLAIMS=""
if [ "${REDO:-0}" = "1" ]; then CLAIMS="claims_$(date +%Y%m%d_%H%M%S)"; fi
trap 'kill $(jobs -p) 2>/dev/null; exit 130' INT TERM      # Ctrl-C stops the workers too
for w in $(seq 0 $((N-1))); do
  CLAIMS=$CLAIMS HB=${HB:-600} CELLS=all NOREGIONS=$NOREG "$PY" runcones.py roots_all_surv.json $w $N all_$w.jsonl $BUDGET > all_$w.log 2>&1 &
done
echo "started $N workers; progress: tail -n 2 all_*.log"
wait
"$PY" summarize.py roots_all_surv.json all_*.jsonl cloud_all_*.jsonl
