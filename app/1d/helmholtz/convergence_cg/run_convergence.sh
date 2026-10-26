#!/usr/bin/env bash
#===============================================================================
# 1D Helmholtz SEM convergence study with cg_helmholtz
#
# Series
#   p1  : po = 1,  ne = 2, 4, ..., 1024
#   p2  : po = 2,  ne = 2, 4, ..., 512
#   ne1 : ne = 1,  po = 2, 3, ..., 16
#
# Base parameters are taken from ../cg_helmholtz.prm; only po and ne are
# varied and SVV is switched off (svv = F). Each run is performed in its own
# temporary directory, so the original .prm and .dat files remain untouched.
#
# Output (in the directory of this script)
#   series_p1.dat, series_p2.dat, series_ne1.dat  with columns
#   po  ne  e_max  e_rms  t_pre  t_sol
#   logs/<series>_po<po>_ne<ne>.log  solver output of each run
#
# The environment (compiler runtime etc.) must be set before calling this script.
#===============================================================================

set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
solver_dir=$(cd "$here/.." && pwd)
solver="$solver_dir/cg_helmholtz"
base_prm="$solver_dir/cg_helmholtz.prm"

[[ -x $solver   ]] || { echo "solver not found: $solver"   >&2; exit 1; }
[[ -f $base_prm ]] || { echo "parameters not found: $base_prm" >&2; exit 1; }

mkdir -p "$here/logs"

#-------------------------------------------------------------------------------
# write parameter file with given po and ne to stdout

make_prm() {
  local po=$1 ne=$2
  sed -E \
    -e "s/^([[:space:]]*po[[:space:]]*=)[[:space:]]*[0-9]+/\1 $po/" \
    -e "s/^([[:space:]]*ne[[:space:]]*=)[[:space:]]*[0-9]+/\1 $ne/" \
    -e "s/^([[:space:]]*svv[[:space:]]*=)[[:space:]]*[TtFf][A-Za-z.]*/\1 F/" \
    "$base_prm"
}

#-------------------------------------------------------------------------------
# run one case and append result line to the series file

run_case() {
  local series=$1 po=$2 ne=$3
  local out="$here/series_$series.dat"
  local log="$here/logs/${series}_po${po}_ne${ne}.log"
  local tmp
  tmp=$(mktemp -d "${TMPDIR:-/tmp}/cg_helmholtz.XXXXXX")

  make_prm "$po" "$ne" > "$tmp/cg_helmholtz.prm"
  if ! (cd "$tmp" && "$solver") > "$log" 2>&1; then
    echo "  run failed: po=$po ne=$ne (see $log)" >&2
    rm -rf "$tmp"
    return 1
  fi
  rm -rf "$tmp"

  awk -v po="$po" -v ne="$ne" '
    $1 == "e_max" { e_max = $3 }
    $1 == "e_rms" { e_rms = $3 }
    $1 == "t_pre" { t_pre = $3 }
    $1 == "t_sol" { t_sol = $3 }
    END { printf "%4d %6d %13s %13s %13s %13s\n", po, ne, e_max, e_rms, t_pre, t_sol }
  ' "$log" | tee -a "$out"
}

#-------------------------------------------------------------------------------
# start a series file with header

new_series() {
  local series=$1 descr=$2
  echo
  echo "series $series: $descr"
  printf "#%3s %6s %13s %13s %13s %13s\n" po ne e_max e_rms t_pre t_sol \
    > "$here/series_$series.dat"
}

#===============================================================================
# convergence series

new_series p1 "po = 1, ne = 2 ... 1024"
for (( ne = 2; ne <= 1024; ne *= 2 )); do run_case p1 1 $ne; done

new_series p2 "po = 2, ne = 2 ... 512"
for (( ne = 2; ne <= 512; ne *= 2 )); do run_case p2 2 $ne; done

new_series ne1 "ne = 1, po = 2 ... 16"
for (( po = 2; po <= 16; po++ )); do run_case ne1 $po 1; done

echo
echo "results written to $here/series_{p1,p2,ne1}.dat"
