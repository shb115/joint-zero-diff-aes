#!/bin/sh
# Run joint_exp on several cores by splitting the runs into contiguous ranges.
#
# Usage: sh scripts/run_parallel.sh [-j JOBS] [-k KEYS] [-s SEED] [-o OUTDIR] [-r REF] [-- joint_exp options]
#
#   -j JOBS    parallel processes (default: number of logical CPUs)
#   -k KEYS    total number of runs (default: 100)
#   -s SEED    master seed (default: current Unix time; printed and recorded)
#   -o OUTDIR  output directory, relative to the repository root (default: out/parallel)
#   -r REF     also compare the merged output with the reference CSV REF
#              (relative to the repository root, like OUTDIR)
#   --         further options for joint_exp, e.g. --log2-trials 28 or --mu-trials 0
#
# With the same SEED the result is identical to a single-process run
#   ./joint_exp --seed SEED --keys KEYS [the same joint_exp options]
# because each run is seeded from (SEED, run index) only.
# Writes OUTDIR/part_XX.csv, OUTDIR/part_XX.log and OUTDIR/merged.csv
# (existing files with these names in OUTDIR are removed first).
# Progress is in OUTDIR/part_XX.log; Ctrl-C or kill PID stops all (exit 130).
# Exit status: 2 on a usage error, 1 if a process, Python or the binary fails,
# otherwise the status of scripts/summarize.py (0 = all checks pass).

set -e
cd "$(dirname "$0")/.."
SELF=$PWD/scripts/$(basename "$0")

JOBS=$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)
KEYS=100
SEED=$(date +%s)
OUTDIR=out/parallel
REF=

usage() { sed -n '2,22p' "$SELF" | sed 's/^# \{0,1\}//'; }

while [ $# -gt 0 ]; do
    case "$1" in
        -j|-k|-s|-o|-r) [ $# -ge 2 ] || { echo "missing value for $1" >&2; exit 2; } ;;
    esac
    case "$1" in
        -j) JOBS=$2; shift 2 ;;
        -k) KEYS=$2; shift 2 ;;
        -s) SEED=$2; shift 2 ;;
        -o) OUTDIR=$2; shift 2 ;;
        -r) REF=$2; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        --) shift; break ;;
        *) echo "unknown option: $1 (see --help)" >&2; exit 2 ;;
    esac
done

positive() { case "$1" in ''|*[!0-9]*) return 1 ;; *) [ "$1" -ge 1 ] ;; esac; }
positive "$JOBS" || { echo "-j must be a positive integer" >&2; exit 2; }
positive "$KEYS" || { echo "-k must be a positive integer" >&2; exit 2; }
case "$SEED" in ''|*[!0-9]*) echo "-s must be a non-negative integer" >&2; exit 2 ;; esac
[ -z "$REF" ] || [ -f "$REF" ] || { echo "reference $REF not found (paths are relative to the repository root)" >&2; exit 2; }
for a in "$@"; do
    case "$a" in
        -k|--keys|-f|--first-run|-s|--seed|--legacy-seed|-o|--output|-q|--quick|--force|-h|--help|--self-test)
            echo "option $a is set by this script or conflicts with it; it cannot be passed to joint_exp" >&2; exit 2 ;;
    esac
done

BIN=./joint_exp
[ -x "$BIN" ] || [ -x "$BIN.exe" ] || { echo "joint_exp not found: run make first" >&2; exit 1; }
if [ -z "$PYTHON" ]; then
    for p in python3 python; do
        if "$p" -c 'import sys' >/dev/null 2>&1; then PYTHON=$p; break; fi
    done
fi
"${PYTHON:-python3}" -c 'import sys; sys.exit(sys.version_info < (3, 6))' >/dev/null 2>&1 \
    || { echo "Python >= 3.6 not found; set PYTHON=..." >&2; exit 1; }

if [ "$JOBS" -gt "$KEYS" ]; then JOBS=$KEYS; fi
mkdir -p "$OUTDIR"
rm -f "$OUTDIR"/part_*.csv "$OUTDIR"/part_*.log "$OUTDIR/merged.csv"
echo "runs 1..$KEYS on $JOBS processes, seed $SEED, output $OUTDIR"

pids=""
trap 'kill $pids 2>/dev/null || true; echo "interrupted" >&2; exit 130' INT TERM
first=1
i=0
while [ $i -lt "$JOBS" ]; do
    count=$(( (KEYS - first + 1) / (JOBS - i) ))     # contiguous, nearly equal ranges
    tag=$(printf "%02d" $i)
    "$BIN" "$@" --seed "$SEED" --first-run "$first" --keys "$count" --force \
           -o "$OUTDIR/part_$tag.csv" > "$OUTDIR/part_$tag.log" 2>&1 &
    pids="$pids $!"
    first=$(( first + count ))
    i=$(( i + 1 ))
done

failed=0
for p in $pids; do
    # exit status 3 = a mu-relation failed; still summarized (and reported as FAIL)
    if wait "$p"; then :; else rc=$?; [ "$rc" -eq 3 ] || failed=1; fi
done
trap - INT TERM
[ "$failed" -eq 0 ] || { echo "error: a process failed, see $OUTDIR/part_*.log" >&2; exit 1; }
echo "all processes finished"

if [ -n "$REF" ]; then
    "$PYTHON" scripts/summarize.py "$OUTDIR"/part_*.csv --expect-runs "$KEYS" --merge "$OUTDIR/merged.csv" --reference "$REF"
else
    "$PYTHON" scripts/summarize.py "$OUTDIR"/part_*.csv --expect-runs "$KEYS" --merge "$OUTDIR/merged.csv"
fi
