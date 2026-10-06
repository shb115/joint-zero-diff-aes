#!/usr/bin/env python3
"""Summarize joint_exp CSV output into the statistics of Table 3 of the paper.

Usage:
    python3 scripts/summarize.py Results/results.csv
    python3 scripts/summarize.py out/full/part_*.csv --merge out/full/merged.csv
    python3 scripts/summarize.py out/full/merged.csv --reference Results/results_v2.csv

Several CSV files (e.g. from scripts/run_parallel.sh) are merged by run index.
Both the original CSV format (8 columns) and the current one are accepted.
Exit status: 0 if all checks pass, 1 otherwise, 2 on bad input or if the
--merge file cannot be written.
Standard library only (Python >= 3.6).
"""

import argparse
import csv
import glob
import io
import math
import os
import sys
from collections import Counter

DEFAULT_TRIALS = 2 ** 30  # the original CSV format has no `trials` column
ALPHA = 0.00135           # per-tail level of the Poisson tests (as for 3 sigma)
REQUIRED = ("run", "mu_verified", "mu_total", "joint_aes", "joint_rand")
INT_COLS = REQUIRED + ("single_xp", "single_xxp", "elapsed_sec", "mu_nontrivial",
                       "half_aes", "half_rand", "trials")


def fail(msg):
    sys.stderr.write("error: %s\n" % msg)
    sys.exit(2)


def load(patterns):
    """Read all files matching the patterns; return rows sorted by run."""
    files = []
    for pattern in patterns:
        matches = [pattern] if os.path.exists(pattern) else sorted(glob.glob(pattern))
        if not matches:
            fail("no file matches %s" % pattern)
        files.extend(matches)
    header, rows = None, {}
    for path in files:
        try:
            with open(path, newline="", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                fields = reader.fieldnames or []
                file_rows = list(reader)
        except (OSError, UnicodeDecodeError, csv.Error) as e:
            fail("cannot read %s as CSV (%s)" % (path, e))
        if header is None:
            header = fields
            missing = [c for c in REQUIRED if c not in header]
            if missing:
                fail("%s is not joint_exp output (missing columns: %s)" % (path, ", ".join(missing)))
        elif fields != header:
            fail("%s has different columns than %s" % (path, files[0]))
        for row in file_rows:
            try:
                run = int(row["run"])
                for c in INT_COLS:
                    if c in header:
                        int(row[c])
            except (TypeError, ValueError):
                fail("%s: malformed row %s" % (path, row))
            if run in rows:
                same = all(rows[run][c] == row[c] for c in header if c != "elapsed_sec")
                if not same:
                    fail("run %d appears twice with different data (%s)" % (run, path))
            rows[run] = row
    if not rows:
        fail("no rows found")
    for col in ("trials", "seed_mode", "seed"):
        if col in header and len({r[col] for r in rows.values()}) != 1:
            fail("rows with different %s cannot be pooled" % col)
    return header, [rows[k] for k in sorted(rows)]


def ints(rows, col):
    return [int(r[col]) for r in rows] if col in rows[0] else None


def mean(xs):
    return sum(xs) / len(xs)


def pvar(xs):
    m = mean(xs)
    return sum((x - m) ** 2 for x in xs) / len(xs)


def svar(xs):
    if len(xs) < 2:
        return float("nan")
    m = mean(xs)
    return sum((x - m) ** 2 for x in xs) / (len(xs) - 1)


def poisson_pmf(k, lam):
    if lam == 0:
        return 1.0 if k == 0 else 0.0
    return math.exp(-lam + k * math.log(lam) - math.lgamma(k + 1))


def poisson_test(total, mu):
    """Two-sided exact test of a Poisson(mu) total; True if not rejected."""
    cdf = sum(poisson_pmf(k, mu) for k in range(total + 1))         # P(X <= total)
    upper = 1.0 - cdf + poisson_pmf(total, mu)                        # P(X >= total)
    return cdf > ALPHA and upper > ALPHA


def fmt_lambda(lam):
    if lam >= 1:
        return "%.1f" % lam
    return "%.2f" % lam if lam >= 0.1 else "%.3g" % lam


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", nargs="+", help="CSV file(s) or glob pattern(s) written by joint_exp")
    ap.add_argument("--merge", metavar="FILE", help="also write the merged rows (sorted by run) to FILE")
    ap.add_argument("--reference", metavar="FILE",
                    help="check that every column except elapsed_sec equals FILE (same seed and options)")
    ap.add_argument("--expect-runs", metavar="N", type=int, help="check that exactly N runs are present")
    ap.add_argument("--checks-only", action="store_true", help="print only the checks")
    args = ap.parse_args()

    header, rows = load(args.csv)
    ref = load([args.reference]) if args.reference else None   # fail early if missing
    n = len(rows)
    trials_col = ints(rows, "trials")
    trials = trials_col[0] if trials_col else DEFAULT_TRIALS
    runs = [int(r["run"]) for r in rows]

    mu_ok, mu_tot = ints(rows, "mu_verified"), ints(rows, "mu_total")
    aes, rnd = ints(rows, "joint_aes"), ints(rows, "joint_rand")
    nontriv = ints(rows, "mu_nontrivial")
    half_aes, half_rnd = ints(rows, "half_aes"), ints(rows, "half_rand")

    # Expected hits per run. lam_aes depends on dx' and on joint_hit; lam_ctl
    # (random function, control condition) and lam_rnd (random baseline, joint
    # condition) depend only on joint_hit. For the provided experiment
    # lam_aes == lam_ctl.
    lam_aes = 12 * 2.0 ** -32 * trials
    lam_ctl = 12 * 2.0 ** -32 * trials
    lam_rnd = 12 * 2.0 ** -64 * trials

    real_stdout = sys.stdout
    if args.checks_only:
        sys.stdout = io.StringIO()
    print("Runs: %d (run %d..%d), quartets per run: 2^%g" % (n, runs[0], runs[-1], math.log2(trials)))
    if trials_col is None:
        print("(original CSV format: 2^30 quartets per run assumed)")
    print()

    # ---- Table 3, left ------------------------------------------------------
    print("Table 3 (left): summary statistics")
    print("  %-24s %22s %10s" % ("", "Small AES", "Random"))
    print("  %-24s %22s %10s" % ("mu-relation", "%d/%d" % (sum(mu_ok), sum(mu_tot)), "---"))
    print("  %-24s %22.2f %10.2f" % ("Hits (mean)", mean(aes), mean(rnd)))
    print("  %-24s %22.2f %10.2f" % ("Hits (var., population)", pvar(aes), pvar(rnd)))
    if n > 1:
        print("  %-24s %22.2f %10.2f" % ("Hits (var., sample)", svar(aes), svar(rnd)))
    else:
        print("  %-24s %22s %10s" % ("Hits (var., sample)", "n/a", "n/a"))
    print("  %-24s %22s %10.1e" % ("Poisson lambda", fmt_lambda(lam_aes), lam_rnd))
    if mean(aes) > 0 and n > 1:
        print("  (index of dispersion N*var/mean = %.1f on %d df; Poisson: mean %d, sd %.1f)"
              % (n * pvar(aes) / mean(aes), n - 1, n - 1, math.sqrt(2 * (n - 1))))
    if nontriv:
        print("  (mu-relation: %d of the checked quartets are non-trivial, i.e. some" % sum(nontriv))
        print("   difference has a zero inverse diagonal; ~1.8e-4 for the provided experiment)")
    print()

    # ---- Table 3, right -----------------------------------------------------
    hist = Counter(aes)
    kmax = max(7, max(aes))
    print("Table 3 (right): distribution of joint hits (small AES), %d run%s" % (n, "" if n == 1 else "s"))
    print("  %6s %9s %9s" % ("Hits", "Observed", "Poisson"))
    for k in range(kmax + 1):
        print("  %6d %9d %9.1f" % (k, hist.get(k, 0), n * poisson_pmf(k, lam_aes)))
    tail = max(0.0, n * (1 - sum(poisson_pmf(k, lam_aes) for k in range(kmax + 1))))
    print("  %6s %9d %9.1f" % (">=%d" % (kmax + 1), 0, tail))
    print()

    # ---- Control: only C0^C2 and C0^C3 ---------------------------------------
    if half_aes:
        print("Control (only C0^C2 and C0^C3 tested; expected per run: small AES %s, random %s):"
              % (fmt_lambda(lam_aes), fmt_lambda(lam_ctl)))
        print("  mean small AES %.2f, mean random %.2f" % (mean(half_aes), mean(half_rnd)))
        print()

    # ---- Checks ---------------------------------------------------------------
    aes_total_expected = n * lam_aes
    checks = [
        ("mu_verified == mu_total in every run", all(a == b for a, b in zip(mu_ok, mu_tot))),
        ("total joint_rand = %d consistent with Poisson(%.2g)" % (sum(rnd), n * lam_rnd),
         poisson_test(sum(rnd), n * lam_rnd)),
        ("total joint_aes = %d consistent with Poisson(%.3g) (exact two-sided test)"
         % (sum(aes), aes_total_expected), poisson_test(sum(aes), aes_total_expected)),
    ]
    if half_aes:
        checks.append(("joint_aes == half_aes in every run", aes == half_aes))
        checks.append(("total half_rand = %d consistent with Poisson(%.3g)" % (sum(half_rnd), n * lam_ctl),
                       poisson_test(sum(half_rnd), n * lam_ctl)))
    checks.append(("runs are consecutive (%d..%d)" % (runs[0], runs[-1]), runs == list(range(runs[0], runs[0] + n))))
    if args.expect_runs is not None:
        checks.append(("exactly %d runs present" % args.expect_runs, n == args.expect_runs))
    if ref:
        ref_header, ref_rows = ref
        cols = [c for c in ref_header if c != "elapsed_sec"]
        mismatch = None
        if len(ref_rows) != n or ref_header != header:
            mismatch = "different number of runs or columns"
        else:
            for a, b in zip(rows, ref_rows):
                bad = [c for c in cols if a[c] != b[c]]
                if bad:
                    mismatch = "run %s, column %s: %s vs %s" % (a["run"], bad[0], a[bad[0]], b[bad[0]])
                    break
        label = "identical to %s in all %d columns except elapsed_sec" % (args.reference, len(cols))
        if mismatch:
            label += " (first difference: %s)" % mismatch
        checks.append((label, mismatch is None))

    sys.stdout = real_stdout
    print("Checks:")
    for text, ok in checks:
        print("  [%s] %s" % ("PASS" if ok else "FAIL", text))
    all_ok = all(ok for _, ok in checks)

    if half_aes and all_ok and not args.checks_only:
        print()
        if aes_total_expected == 0:
            print("lambda for small AES is 0: the joint condition cannot occur in this setting.")
        elif sum(half_aes) == 0 or sum(half_rnd) == 0:
            print("Too few control hits (small AES %d, random function %d; expected %.3g and %.3g)"
                  % (sum(half_aes), sum(half_rnd), aes_total_expected, n * lam_ctl))
            print("to show the joint property; use more runs or more quartets per run")
            print("(for the provided experiment: make medium or larger).")
        else:
            print("Under the control condition both small AES and the random function give hits;")
            print("the joint condition keeps every small-AES hit (joint_aes == half_aes) and")
            print("removes the random ones (joint_rand == 0), as predicted by Theorem 3.")

    if args.merge:
        try:
            with open(args.merge, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=header, lineterminator="\n")
                w.writeheader()
                w.writerows(rows)
        except OSError as e:
            fail("cannot write %s (%s)" % (args.merge, e))
        if not args.checks_only:
            print("\nMerged rows written to %s" % args.merge)

    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
