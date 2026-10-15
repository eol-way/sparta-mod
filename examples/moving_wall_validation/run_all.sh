#!/bin/bash
# regression tests for moving surfs: engulf project and wallvel
# usage: ./run_all.sh [path/to/spa_executable] [nprocs]
#   default executable = ../../src/spa_mpi, default nprocs = 4
#   set MPIRUN to change the launcher, default "mpirun --oversubscribe -np N"
# all output goes to ./out, check.py prints PASS/FAIL per check
# exit code 1 if any check fails or a run aborts

set -u

HERE=$(cd "$(dirname "$0")" && pwd)
SPA=${1:-$HERE/../../src/spa_mpi}
NP=${2:-4}
MPIRUN=${MPIRUN:-mpirun --oversubscribe -np $NP}
OUT=$HERE/out

if [ ! -x "$SPA" ]; then
  echo "SPARTA executable not found: $SPA"
  exit 1
fi
SPA=$(cd "$(dirname "$SPA")" && pwd)/$(basename "$SPA")   # runs are in ./out

rm -rf "$OUT"
mkdir -p "$OUT"
cp "$HERE"/data.* "$HERE"/air.* "$HERE"/in.* "$OUT"/
cd "$OUT"

nabort=0
run() {
  local log=$1
  shift
  printf "running %-34s" "$log"
  if $MPIRUN "$SPA" -log "$log" "$@" > "$log.screen" 2>&1; then
    echo "done"
  else
    echo "ABORTED (see out/$log.screen)"
    nabort=$((nabort+1))
  fi
}

# runs which must stop with an error, check.py looks for the message

run_expect_error() {
  local log=$1
  shift
  printf "running %-34s" "$log"
  if $MPIRUN "$SPA" -log "$log" "$@" > "$log.screen" 2>&1; then
    echo "done"
  else
    echo "stopped (error expected)"
  fi
}

run log.ellipse2d.delete       -in in.ellipse2d -var mode delete
run log.ellipse2d.project      -in in.ellipse2d -var mode project
run log.ellipsoid3d.delete     -in in.ellipsoid3d -var mode delete
run log.ellipsoid3d.project    -in in.ellipsoid3d -var mode project
run log.two_rotors             -in in.two_rotors
run log.move_surf_cmd.delete   -in in.move_surf_cmd -var mode delete
run log.move_surf_cmd.project  -in in.move_surf_cmd -var mode project
run log.piston.static          -in in.piston -var engulf project -var wallvel yes -var dx 0.0
run log.piston.project_yes     -in in.piston -var engulf project -var wallvel yes
run log.piston.project_no      -in in.piston -var engulf project -var wallvel no
run log.piston.delete_no       -in in.piston -var engulf delete -var wallvel no
run log.piston.startup         -in in.piston -var engulf project -var wallvel yes \
                               -var nevery 20 -var nrun 40 -var fnum 5.0e12
run log.rotating_box.project_yes -in in.rotating_box -var engulf project -var wallvel yes
run log.rotating_box.project_no  -in in.rotating_box -var engulf project -var wallvel no
run log.amr_static.none        -in in.amr_static -var amr none
run log.amr_static.surf        -in in.amr_static -var amr surf
run log.amr_static.region      -in in.amr_static -var amr region
run log.amr_dynamic.surf       -in in.amr_dynamic -var adapt surf
run log.amr_dynamic.particle   -in in.amr_dynamic -var adapt particle
run log.amr_3d                 -in in.amr_3d
run log.test_translation.diffuse_yes -in in.test_translation -var model diffuse -var wallvel yes
run log.test_translation.diffuse_no  -in in.test_translation -var model diffuse -var wallvel no
run log.test_translation.vanish_yes  -in in.test_translation -var model vanish -var wallvel yes
run log.test_rotation.diffuse_yes    -in in.test_rotation -var model diffuse -var wallvel yes
run log.test_rotation.diffuse_no     -in in.test_rotation -var model diffuse -var wallvel no
run log.test_rotation.vanish_yes     -in in.test_rotation -var model vanish -var wallvel yes
run log.guard_rotate.no        -in in.guard_rotate -var wallvel no
run_expect_error log.guard_rotate.yes   -in in.guard_rotate -var wallvel yes
run_expect_error log.guard_cll_translate -in in.guard_cll_translate
run_expect_error log.guard_piston       -in in.guard_piston

echo
python3 "$HERE/check.py" "$OUT" "$NP"
status=$?

# velocity distributions of diffusely reflected particles
# plots need matplotlib, ./.venv is used if it exists

echo
PY=python3
[ -x "$HERE/.venv/bin/python" ] && PY=$HERE/.venv/bin/python
"$PY" "$HERE/verify_distribution.py" "$OUT" "$OUT" || status=1

if [ $nabort -gt 0 ]; then
  echo "$nabort runs aborted"
  exit 1
fi
exit $status
