#!/usr/bin/env bash
#
# Reproduce every computational claim in NOTES.md for Erdos problem #1162.
#
# Stages:
#   1. fetch primary sources        (refs/fetch.sh; nothing is redistributed)
#   2. build the lattice enumerator (cc -O3; NEON on arm64, portable elsewhere)
#   3. full subgroup lattice        n = 1..8, validated against A005432/A000638
#   4. brute-force cross-check      verify_elemab.py, which shares no code with elemab.py
#   5. exact elementary abelians    n up to 512: rank law, Galois constant
#   6. the orbit model              second-order term, incl. the held-out extrapolation
#   7. final_check.py               re-derives every claim from the artifacts; exit 0 = pass
#
# Needs: cc (C99), curl, a python3 with numpy and mpmath, and poppler's
# pdftotext for the quotation checks in stage 7.  $PY overrides the
# interpreter; otherwise ~/.venvs/main/bin/python is used when present and
# python3 from $PATH otherwise.  Without the PDF those checks report "skip",
# not "pass" -- final_check.py distinguishes them.
#
# Runtime on an otherwise-idle M-series laptop.  FAST figures are the range
# over two runs; the full-only figures are from one run each, so treat them as
# indicative rather than a mean.
#
#   stage            FAST=1                    full
#   3 lattice        8.8-11.6 s  (n<=7)        ~15 min  (n=8 alone: 883 s)
#   4 brute force    43.9-50.8 s (n<=9)        913 s    (n=10)
#   5 elementary ab. 0.8-0.9 s   (n<=128)      ~4 min   (n=512, needs ~2 GB)
#   6 orbit model    0.2-0.3 s                 ~90 s
#   6b sd precision  ~2 s        (n<=128)      115 s    (n=512, 60 dp)
#   7 final check    1.1-1.2 s                 ~3 s
#
# FAST=1 totalled ~60 s end to end excluding the downloads; full is ~35 min.
# The n=10 brute force, not the lattice, is what FAST is avoiding.
#
# FAST=1 reproduces every claim in NOTES.md except the n=8 lattice row, the
# n=10 brute-force row, the 512-column of the rank table and the 256/512 rows
# of the high-precision sd table; final_check.py then reports the held-out
# extrapolation as "skip" (not enough points) and prints fewer checks.  Both
# were run: full gives 118 checks and 0 failed.  A skip is not a pass and the
# summary line says so.
#
# Everything is redirected to results/*.log.  Do not pipe these to `tail`:
# stdout is block-buffered through a pipe and you will see nothing until the
# end of a 15-minute run.
set -euo pipefail
cd "$(dirname "$0")/.."

PY="${PY:-}"
if [ -z "$PY" ]; then
  if [ -x "$HOME/.venvs/main/bin/python" ]; then
    PY="$HOME/.venvs/main/bin/python"
  else
    PY="$(command -v python3 || true)"
  fi
fi
CC="${CC:-cc}"
FAST="${FAST:-0}"

if [ -z "$PY" ] || ! "$PY" -c 'import numpy, mpmath' >/dev/null 2>&1; then
  echo "need a python3 with numpy and mpmath -- set PY=/path/to/python" >&2
  exit 1
fi

mkdir -p data results

if [ "$FAST" = "1" ]; then
  LAT_MAX=7; N_MAX=128; BF_MAX=9
else
  LAT_MAX=8; N_MAX=512; BF_MAX=10
fi

echo "=== 1. primary sources ==============================================="
# Not fatal: every stage below runs without them, and final_check.py reports
# the quotation checks as skipped rather than passed.
sh refs/fetch.sh > results/fetch.log 2>&1 || \
  echo "  fetch failed (offline?) -- quotation checks will report 'skip'"

echo "=== 2. build ========================================================="
"$CC" -O3 -std=c99 -Wall -o src/subgroups src/subgroups.c
echo "  src/subgroups built"

echo "=== 3. full subgroup lattice of S_n, n <= $LAT_MAX ========================="
# Counts every subgroup, not just up to conjugacy, by expanding each class.
# n=8 is the slow one; it is also the only new row relative to the literature
# that this stage produces by an independent route.
for n in $(seq 1 "$LAT_MAX"); do
  ./src/subgroups "$n" "data/lattice_n$n.json" > "results/lattice_n$n.log" 2>&1
  "$PY" - "$n" <<'EOF'
import json, sys
n = int(sys.argv[1])
d = json.load(open(f"data/lattice_n{n}.json"))
print(f"  n={n:2d}  subgroups={d['total_subgroups']:>8d}  classes={d['total_classes']:>4d}")
EOF
done

echo "=== 4. brute-force cross-check of a(n,k), n <= $BF_MAX ================="
# Builds every elementary abelian 2-subgroup by multiplying permutations.
# Shares no code with elemab.py, which gets the same numbers from a formula.
"$PY" src/verify_elemab.py --n-max "$BF_MAX" --out results/verify_elemab.json \
  > results/verify_elemab.log 2>&1
tail -3 results/verify_elemab.log

echo "=== 5. elementary abelian 2-subgroups, n <= $N_MAX ====================="
"$PY" src/asymptotics.py --n-max "$N_MAX" --dense-to 64 --step 8 --outdir data \
  > "results/asymptotics_$N_MAX.log" 2>&1
tail -6 "results/asymptotics_$N_MAX.log"

echo "=== 6. orbit model for the n log n term =============================="
if [ "$FAST" = "1" ]; then
  "$PY" src/orbit_model.py --profile-n 64 128 --pure-n 64 128 --outdir data \
    > results/orbit_model.log 2>&1
else
  "$PY" src/orbit_model.py --outdir data > results/orbit_model.log 2>&1
fi
tail -8 results/orbit_model.log

echo "=== 6b. rank sd against its closed form, exact at 60 dp =============="
# The doubles in elemab_stats.json truncate at ~1e-16 relative, which is
# coarser than the actual convergence from n=256 on.  This recomputes the
# variance from the exact integer counts in mpmath.
if [ "$FAST" = "1" ]; then
  "$PY" src/sd_precision.py --n 64 128 --out data/sd_precision.json \
    > results/sd_precision.log 2>&1
else
  "$PY" src/sd_precision.py --out data/sd_precision.json \
    > results/sd_precision.log 2>&1
fi
grep -E "^  n=" results/sd_precision.log

echo "=== 7. final check ==================================================="
# The gate.  Exit 0 means every number in NOTES.md was re-derived from the
# artifacts above.  It was itself validated by planting three defects (see
# NOTES.md section 9); a clean pass of an unchallenged checker proves nothing.
"$PY" src/final_check.py 2>&1 | tee results/final_check.log
