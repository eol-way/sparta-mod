#!/bin/bash
# GIF animations of the moving wall validation runs
# usage: ./generate_all_gifs.sh [path/to/spa_executable] [nprocs]
#   default executable = ../../src/spa_mpi, default nprocs = 4
#   set MPIRUN to change the launcher, default "mpirun --oversubscribe -np N"
#   set PYTHON to choose the interpreter, default ./.venv/bin/python if it
#     exists, else python3; it needs numpy, matplotlib and Pillow:
#       python3 -m venv .venv && .venv/bin/pip install matplotlib pillow
# runs the inputs with -var anim N in ./anim_out, then make_gifs.py writes
#   ./gifs/*.gif
# exit code 1 if a run aborts or no GIF is written

set -u

HERE=$(cd "$(dirname "$0")" && pwd)
SPA=${1:-$HERE/../../src/spa_mpi}
NP=${2:-4}
MPIRUN=${MPIRUN:-mpirun --oversubscribe -np $NP}
WORK=$HERE/anim_out
GIFS=$HERE/gifs

if [ -z "${PYTHON:-}" ]; then
  if [ -x "$HERE/.venv/bin/python" ]; then PYTHON=$HERE/.venv/bin/python
  else PYTHON=python3
  fi
fi
if ! "$PYTHON" -c "import numpy, matplotlib, PIL" 2> /dev/null; then
  echo "$PYTHON lacks numpy, matplotlib or Pillow, install them with"
  echo "  python3 -m venv $HERE/.venv && $HERE/.venv/bin/pip install matplotlib pillow"
  exit 1
fi

if [ ! -x "$SPA" ]; then
  echo "SPARTA executable not found: $SPA"
  exit 1
fi
SPA=$(cd "$(dirname "$SPA")" && pwd)/$(basename "$SPA")   # runs are in ./anim_out

rm -rf "$WORK"
mkdir -p "$WORK" "$GIFS"
cp "$HERE"/data.* "$HERE"/air.* "$HERE"/in.* "$WORK"/
cd "$WORK"

nabort=0
run() {
  local log=$1
  shift
  printf "running %-36s" "$log"
  if $MPIRUN "$SPA" -log "$log" "$@" > "$log.screen" 2>&1; then
    echo "done"
  else
    echo "ABORTED (see anim_out/$log.screen)"
    nabort=$((nabort+1))
  fi
}

# frame interval -var anim N: about 30-50 frames per run

run log.ellipse2d.delete        -in in.ellipse2d -var mode delete -var anim 20
run log.ellipse2d.project       -in in.ellipse2d -var mode project -var anim 20
run log.two_rotors              -in in.two_rotors -var anim 20
run log.ellipsoid3d.delete      -in in.ellipsoid3d -var mode delete -var anim 10
run log.ellipsoid3d.project     -in in.ellipsoid3d -var mode project -var anim 10
run log.piston.project_yes      -in in.piston -var engulf project -var wallvel yes -var anim 10
run log.piston.project_no       -in in.piston -var engulf project -var wallvel no -var anim 10
run log.piston.delete_no        -in in.piston -var engulf delete -var wallvel no -var anim 10
run log.rotating_box.project_yes -in in.rotating_box -var engulf project -var wallvel yes -var anim 100
run log.rotating_box.project_no  -in in.rotating_box -var engulf project -var wallvel no -var anim 100

echo
"$PYTHON" "$HERE/make_gifs.py" "$WORK" "$GIFS"
status=$?

if [ $nabort -gt 0 ]; then
  echo "$nabort runs aborted"
  exit 1
fi
if [ $status -ne 0 ] || ! ls "$GIFS"/*.gif > /dev/null 2>&1; then
  echo "GIF generation failed"
  exit 1
fi
exit 0
