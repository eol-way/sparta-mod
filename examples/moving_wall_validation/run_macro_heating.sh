#!/bin/bash
# long-time energy budget of a gas in a rotating container (in.macro_heating)
# not part of run_all.sh: two runs of 50000 steps, about 35 s each on 4 procs
# usage: ./run_macro_heating.sh [path/to/spa_executable] [nprocs] [nrun]
#   default executable = ../../src/spa_mpi, nprocs = 4, nrun = 50000
#   set MPIRUN to change the launcher, default "mpirun --oversubscribe -np N"
# output goes to ./out_macro, verify_macro_heating.py prints PASS/FAIL and
#   plots out_macro/macro_heating.png (if matplotlib is available)

set -u

HERE=$(cd "$(dirname "$0")" && pwd)
SPA=${1:-$HERE/../../src/spa_mpi}
NP=${2:-4}
NRUN=${3:-50000}
MPIRUN=${MPIRUN:-mpirun --oversubscribe -np $NP}
OUT=$HERE/out_macro

if [ ! -x "$SPA" ]; then
  echo "SPARTA executable not found: $SPA"
  exit 1
fi
SPA=$(cd "$(dirname "$SPA")" && pwd)/$(basename "$SPA")

rm -rf "$OUT"
mkdir -p "$OUT"
cp "$HERE"/data.square "$HERE"/air.* "$HERE"/in.macro_heating "$OUT"/
cd "$OUT"

for wv in yes no; do
  printf "running %-34s" "log.macro_heating.$wv"
  if $MPIRUN "$SPA" -log log.macro_heating.$wv -in in.macro_heating \
       -var wallvel $wv -var nrun $NRUN > log.macro_heating.$wv.screen 2>&1; then
    echo "done"
  else
    echo "ABORTED (see out_macro/log.macro_heating.$wv.screen)"
  fi
done

echo
PY=python3
[ -x "$HERE/.venv/bin/python" ] && PY=$HERE/.venv/bin/python
"$PY" "$HERE/verify_macro_heating.py" "$OUT" "$OUT"
