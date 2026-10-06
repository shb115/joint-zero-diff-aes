# Joint Generalized Zero-Difference Property: Experiments on Small-Scale AES

This repository contains the code and data for the experiments in Section 3.3 of

> H. Shin, S. Kim, B. Seok, D. Hong, J. Sung, S. Hong, S. Lee, D. Lee.
> **Key-Independent Secret-Key Distinguisher for 7-Round AES based on the Joint Generalized Zero-Difference Property.**
> ASIACRYPT 2026.

The program verifies the Joint Generalized Zero-Difference Property (Theorem 3 of the paper) on 4-round small-scale AES (64-bit block, 4-bit S-box). The repository also contains:
- the raw data behind Table 3, and a full re-run of the experiment;
- scripts that run the experiment in parallel and regenerate Table 3.

## Quick start

All commands are run from the repository root in a POSIX shell. They need a C99 compiler and Python ≥ 3.6; GNU make is optional (see *Requirements* for tested versions and for Windows).

```bash
make           # build ./joint_exp
make test      # check the cipher implementation
make table3    # print Table 3 of the paper from the stored data
make quick     # a small experiment, compared with a stored reference output
```

| Command | Time | What it does |
|---|---|---|
| `make` | seconds | Builds `./joint_exp`. |
| `make test` | seconds | Checks the cipher: the S-box is the small-scale AES S-box (affine map of the GF(2^4) inverse), and the table-based encryption agrees with a straightforward implementation. |
| `make table3` | seconds | Prints the values of Table 3 from the paper's data `Results/results.csv`. |
| `make quick` | < 1 min | 4 runs × 2^24 quartets. The μ-relations hold, and the output is compared with `Results/quick_seed20261002.csv` in every column except `elapsed_sec`. |
| `make medium` | ~2–3 min on 16 cores (~35 core-minutes) | 16 runs × 2^28 quartets in parallel, enough to see the joint property (see *Reproducing the paper's results*). The output is compared with `Results/medium_seed20261002.csv`. |
| `make full` | ~1 h on 16 cores (~15 core-hours) | The paper's experiment: 100 runs × 2^30 quartets. The output is compared with `Results/results_v2.csv`. Its statistics agree with the Poisson model of Table 3 (λ = 3.0 for small AES, ≈ 0 for random) within the tolerances in *Interpreting the results*. |

Make variables:
- `JOBS=N`: parallel processes. The default is all logical CPUs. The results do not depend on it.
- `SEED=S`: default 20261002. Reference comparisons are made only for the default.
- `PYTHON`, `CC`, `CFLAGS`, `CPPFLAGS` (e.g. `CPPFLAGS=-DNUM_ROUNDS=5`, see *Building and running*).

Example: `make full JOBS=8`. Each target exits with a non-zero status if a check fails. `make clean` removes the binaries (including `joint_exp5`, see *Building and running*) and `out/quick.csv`; `make distclean` also removes `out/`.

## Reproducing the paper's results

**Table 3.** `make table3` prints the values of Table 3 from `Results/results.csv`, the data of the paper's experiment. Its random seed was not recorded (see *Data*), so a repetition of the experiment gives a different sample of Poisson(3) counts, not exactly the mean 2.86 and variance 2.38 of Table 3.

**The full experiment.** `make full` re-runs the experiment with the default seed 20261002. Its output is then identical to `Results/results_v2.csv` except for `elapsed_sec`. The statistics of that re-run (mean 3.10, variance 2.41, no random hits) are consistent with the Poisson model of Table 3 (see *Interpreting the results*).

**What `make medium` shows.**
- Expected totals: 12 joint hits for small AES, ≈0 for the random baseline, and 12 for a random function under the weaker control condition, which tests only C0⊕C2 and C0⊕C3 (columns `half_aes`, `half_rand`; see *Interpreting the results*).
- With the default seed:
  - small AES gives 4 hits, a low but accepted Poisson draw (P[X ≤ 4] ≈ 0.008);
  - the random baseline gives 0, and the control gives 12;
  - every small-AES hit satisfies the full condition (`joint_aes == half_aes`).
- These 16 runs are the first 2^28 quartets of runs 1–16 of `Results/results_v2.csv`; over their full 2^30 quartets, those runs have 46 hits against 48 expected.

**Checking a single full-size run** (about 5–10 minutes on one core; any K in 1..100, here K = 12):
```bash
mkdir -p out && ./joint_exp --seed 20261002 -f 12 -k 1 -m 0 --force -o out/run12.csv
{ head -n 1 Results/results_v2.csv; grep "^12," Results/results_v2.csv; } > out/ref12.csv
python3 scripts/summarize.py out/run12.csv --reference out/ref12.csv
```
(Use `python` instead of `python3` if that is the interpreter's name.)
The run is identical to row 12 of `Results/results_v2.csv` except for `elapsed_sec`.

## Paper-to-code mapping

Section, theorem and table numbers follow the ASIACRYPT 2026 proceedings version of the paper.

| Paper | Claim | Code / data | How to check |
|---|---|---|---|
| Sec. 3.3, *Setup* | 4-round small-scale AES (a 2-round SPN in the superbox representation), 64-bit block, 4-bit S-box. Δx has no zero nibble; Δx' is zero on diagonals 0 and 2 and equals Δx elsewhere. A fresh random key and fresh differences per run; 2^30 quartets per run; 100 runs. | `Code/joint_property_experiment.c`: cipher (`SBOX`, `gf_mul`, `round_mc`, `round_last`, `encrypt`) and difference sampling (`run_experiment`). Implementation choices: the last round has no MixColumns (see *Notes*); the key is five independent uniform round keys (no key schedule). | `make quick` |
| Theorem 3 (Sec. 3.1); Sec. 3.3, *Experiment 1* | The three μ-relations hold for every quartet: 10^8/10^8 (the first 10^6 quartets of each of 100 runs). | μ-check in `run_experiment`; CSV columns `mu_verified`, `mu_total` | `mu_verified == mu_total` in every run |
| Sec. 3.3, *Experiment 2* | Joint zero inverse-diagonal hits. Per quartet the probability is 12·2^-32 for small-scale AES (λ = 3.0 per run) and 12·2^-64 for a random permutation (λ ≈ 7.0·10^-10 per run). | CSV columns `joint_aes`, `joint_rand`. The random permutation is modeled by a random function (see *Notes*). | `joint_rand == 0`; total of `joint_aes` consistent with Poisson(100·3) |
| Table 3 (left) | μ 10^8/10^8; hits mean 2.86 and variance 2.38 (small AES) vs 0.00 and 0.00 (random); λ = 3.0 vs ≈ 0. | `Results/results.csv` | `make table3` |
| Table 3 (right) | Histogram of joint hits vs Poisson(3): 4, 16, 26, 20, 18, 11, 4, 1 runs with 0..7 hits. | `Results/results.csv` | `make table3` |
| *Not covered* | Theorem 4. The 7-round characteristic (Sec. 3.2). The distinguisher of Sec. 4 (Algorithm 1), with complexity 2^126.2 and success probability ≥ 75.7%. The transfer of the μ-relations to full AES (Sec. 3.3, in the paragraph beginning "The observed mean of 2.86"). | — | These are analytical results. As explained at the end of Sec. 3.3, an end-to-end experiment would need about 2^62.2 chosen plaintexts even for small-scale AES. |

## Requirements

| Component | Needed for | Tested versions |
|---|---|---|
| C99 compiler with the standard C library and libm | `joint_exp` | GCC 11.2.0 (MinGW-w64, Windows 11); GCC 9.4.0 (Ubuntu 20.04) |
| Python ≥ 3.6, standard library only | `scripts/summarize.py` | 3.10.2 (Windows), 3.8.10 (Ubuntu) |
| POSIX `sh` | `scripts/run_parallel.sh` | `dash` (Ubuntu 20.04); bash 5.2 from Git for Windows 2.53 |
| GNU Make (optional) | `Makefile` | 4.3 (MinGW, Windows); 4.2.1 (Ubuntu 20.04) |

There are no other dependencies. The CSV output contains no floating-point values (only integers, the seeding mode, and Δx in hex). For a given seed it is identical on the platforms above, except for `elapsed_sec`.

**Windows.** Use Git Bash, MSYS2 or WSL; PowerShell and cmd are not supported. Git Bash does not ship a compiler or make: install them separately, plus Python from python.org. Options are MinGW-w64 GCC and make (if it provides only `mingw32-make`, use that name instead of `make`) or MSYS2 (`pacman -S mingw-w64-x86_64-gcc make`).

**Without make** (on any platform), the commands behind each target are listed in the `Makefile`. For example, `make quick` is:
`mkdir -p out && ./joint_exp --quick --seed 20261002 --force -o out/quick.csv && python3 scripts/summarize.py out/quick.csv --reference Results/quick_seed20261002.csv` (use `python` if that is the interpreter's name).

## Building and running

```bash
gcc -O2 -std=c99 -o joint_exp Code/joint_property_experiment.c -lm
./joint_exp --help
```

| Option | Default | Meaning |
|---|---|---|
| `-o, --output FILE` | `results.csv` | CSV output, also accepted as a positional argument. |
| `--force` | — | Overwrite an existing output file (otherwise the program refuses). |
| `-k, --keys N` | 100 | Number of runs (independent keys). |
| `-f, --first-run K` | 1 | Index of the first run; the program executes runs K..K+N-1. |
| `-t, --log2-trials T` | 30 | Quartets per run = 2^T, 1 ≤ T ≤ 40. |
| `-m, --mu-trials M` | 10^6 | Quartets per run checked for the μ-relations (0 = all). |
| `-s, --seed S` | current time | Master seed, recorded in the CSV. |
| `--legacy-seed T` | — | Seeding scheme of the original version (see *Data*). |
| `-q, --quick` | — | Smoke test: 4 runs of 2^24 quartets. An explicit `-k` or `-t` takes precedence. |
| `--self-test` | — | Check the cipher implementation and exit (as `make test`). Acts as soon as it is read: later options and the range checks are skipped; earlier options must still be valid. |
| `-h, --help` | — | Usage, including exit codes. Acts as soon as it is read, like `--self-test`. |

Exit codes:
- 0: success
- 1: I/O error, the output file exists (without `--force`), or the self-test failed
- 2: usage error (including an option given twice, or `--seed` together with `--legacy-seed`)
- 3: a quartet violated a μ-relation

Run k is seeded from (S, k) only, so a set of runs gives the same result whether it runs in one process or is split over several. `scripts/run_parallel.sh` uses this to split the runs over all cores. A small example, which reproduces `make quick` in parallel in a few seconds:

```bash
sh scripts/run_parallel.sh -j 4 -k 4 -s 20261002 -o out/qp -r Results/quick_seed20261002.csv -- --log2-trials 24
```

The full experiment (what `make full JOBS=16` runs):

```bash
sh scripts/run_parallel.sh -j 16 -k 100 -s 20261002 -o out/full -r Results/results_v2.csv -- --mu-trials 0
```

The script writes `part_XX.csv`, `part_XX.log` (progress) and `merged.csv` into the `-o` directory and prints the summary.
- `-r FILE` adds the comparison with a reference output.
- Options after `--` are passed to `joint_exp`, except those the script sets itself or that conflict with it (`-k`/`--keys`, `-f`/`--first-run`, `-s`/`--seed`, `--legacy-seed`, `-o`/`--output`, `-q`/`--quick`, `--force`, `-h`/`--help`, `--self-test`), which it rejects with exit status 2.
- The `-o` and `-r` paths are relative to the repository root.
- Existing `part_*` files and `merged.csv` in the `-o` directory are removed first.
- Ctrl-C (or `kill PID`) stops all processes; the script then exits with status 130.
- Exit status: 2 on a usage error, 1 if a process, Python or the binary fails, and otherwise that of `summarize.py`.

`python3 scripts/summarize.py --help` lists its options (`--merge`, `--reference`, `--expect-runs`, `--checks-only`). It exits with 0 if all checks pass, 1 if a check fails, and 2 on bad input or if the `--merge` file cannot be written.

**Running time.** One run of 2^30 quartets takes about 5–10 minutes on one core. On an i7-12700K (12 cores, 20 threads) single runs took 5.3–7.4 minutes depending on the core and the load, and 6–9 minutes per run with 16 parallel processes. All 100 runs take about 15 core-hours, i.e. about 1 hour with 16 processes on that machine. The original run behind `Results/results.csv` took about 11.6 minutes per run in its own environment (see *Data*); built with the same compiler on the same machine, the original program runs at about the same speed as this version.

**Changing the number of rounds.** `NUM_ROUNDS` is a compile-time constant (`make -B CPPFLAGS=-DNUM_ROUNDS=5`, or `-DNUM_ROUNDS=5` on the gcc command line). Building through make overwrites `./joint_exp`; run `make -B` afterwards to restore the 4-round build. The experiment matches the setting of Theorem 3 only for the default value of 4.
- With 5 or 6 rounds the μ-check fails and `joint_exp` exits with status 3. Under `make quick`, make stops at that point; under `make medium` and `make full`, the summary still runs and reports FAIL.
- With fewer than 4 rounds the μ-check passes and `joint_exp` exits with 0. The outputs then differ from the reference files, so the reference comparison fails.

The 5-round case doubles as a negative control showing that the μ-check is not vacuous. It uses a separate binary, so `./joint_exp` stays at 4 rounds:
- Build: `gcc -O2 -std=c99 -DNUM_ROUNDS=5 -o joint_exp5 Code/joint_property_experiment.c -lm`
- Run: `mkdir -p out && ./joint_exp5 --seed 1 -k 2 -t 22 -m 0 --force -o out/r5.csv`
- Result: 8,385,492 of 8,388,608 quartets verified. All 3,116 non-trivial quartets violate the μ-relations, and the program exits with status 3. The 4-round build with the same options verifies all quartets.

## Output format (CSV, one row per run)

| Column | Meaning |
|---|---|
| `run` | Run index (1-based) |
| `mu_verified` | Checked quartets satisfying all three μ-relations |
| `mu_total` | Quartets checked for the μ-relations |
| `single_xp` | Quartets where C0⊕C2 has a zero inverse diagonal (expected 4·2^-16 per quartet) |
| `single_xxp` | Same for C0⊕C3 |
| `joint_aes` | Joint hits for small-scale AES: C0⊕C2 and C1⊕C3 share a zero inverse diagonal, and C0⊕C3 and C1⊕C2 share a zero inverse diagonal at a disjoint position |
| `joint_rand` | Same condition for the random baseline |
| `elapsed_sec` | Running time of the run |
| `mu_nontrivial` | Checked quartets in which at least one of the six differences has a zero inverse diagonal. The other quartets satisfy the μ-relations trivially. For the provided experiment this is about 1.8·10^-4 of the checked quartets. |
| `half_aes`, `half_rand` | Control: hits when only C0⊕C2 and C0⊕C3 are tested |
| `trials` | Quartets per run |
| `seed_mode`, `seed` | Seeding scheme (`run` or `legacy`) and master seed |
| `dx` | Δx in hex: 16 nibbles in the code's state order (index 4r+c, i.e. row by row; see *Notes*) |

`Results/results.csv` uses the original format, which has only the first eight columns.

## Interpreting the results

`scripts/summarize.py FILE...` prints both halves of Table 3 for any set of runs and checks the following:

- `mu_verified == mu_total` in every run. This is Experiment 1 / Theorem 3, which is deterministic.
- The total of `joint_rand` is consistent with Poisson(N·λ_rnd), where λ_rnd = 12·2^-64·2^T ≈ 7·10^-10 per run for 2^30 quartets. In practice this means no hit at all.
- The total of `joint_aes` is consistent with Poisson(N·λ), where λ = 12·2^-32·2^T and N is the number of runs. This is an exact two-sided test, rejecting at 0.135% per tail as for 3σ.
  - For 100 runs of 2^30 quartets the mean must lie in [2.49, 3.53].
  - An independent repetition gives a different sample of Poisson(3) counts. The exact values of Table 3 (mean 2.86, variance 2.38) therefore come only from the stored data.
- For CSV files with the control columns (not `Results/results.csv`): `joint_aes == half_aes` in every run, and the total of `half_rand` is consistent with Poisson(N·λ_ctl). Under the weaker control condition (only C0⊕C2 and C0⊕C3), a random function gives λ_ctl = 12·2^-32·2^T hits per run, the same as λ for the provided experiment. The joint condition removes the random hits (`joint_rand == 0`) but keeps every small-AES hit. This is the joint property predicted by Theorem 3.
- The run indices are consecutive, and with `--expect-runs N` (used by `run_parallel.sh`) exactly N runs are present.
- With `--reference FILE`: the output is identical to FILE in every column except `elapsed_sec`.

Expected values for the provided modes:

| Mode | Runs × quartets | λ per run (`joint_aes`, `half_aes`, `half_rand`) | `single_xp` per run |
|---|---|---|---|
| `make quick` | 4 × 2^24 | 0.047 (too small to show hits) | ≈ 1,024 |
| `make medium` | 16 × 2^28 | 0.75 (12 in total) | ≈ 16,384 |
| `make full` | 100 × 2^30 | 3.0 | ≈ 65,536 |

## Data

**`Results/results.csv`** holds the data of Table 3.
- It was produced by the version used for the paper (commit `542ca32`, see *Version history*). The configuration was 100 runs of 2^30 quartets, with the μ-relations checked on the first 10^6 quartets of each run. It ran on the first author's desktop (Intel Core i7-12700K, Windows 11), built with `gcc -O2` (compiler version not recorded) and run as one sequential process: about 694 s per run, 19.3 h in total. The counts do not depend on the compiler. For a fixed seed, GCC 9.4.0 on Ubuntu and GCC 11.2.0 on Windows give identical outputs.
- That version seeded the PRNG with the start time of the program and did not record it. Candidate seeds derived from the file's timestamp did not reproduce run 1, so the start time could not be recovered. This file can therefore be reproduced statistically but not bit-exactly. `--legacy-seed T` re-implements that seeding scheme for a given Unix time T: it reproduces the original eight columns (except `elapsed_sec`) with the default `-t` and `-m`.
- That version counted the AES hits by testing only C0⊕C2 and C0⊕C3. By Theorem 3 this gives the same counts as the four-difference condition (see *Version history*).

**`Results/results_v2.csv`** is a full re-run with this version.
- Command: `sh scripts/run_parallel.sh -j 16 -k 100 -s 20261002 -o out/full -- --mu-trials 0`, which is what `make full JOBS=16` runs apart from the comparison with this file (`-r Results/results_v2.csv`): seed 20261002, 100 runs of 2^30 quartets, μ-relations checked on all quartets.
- Platform: Intel Core i7-12700K, Windows 11, GCC 11.2.0, 16 processes, 60 minutes wall time (357–557 s per run).
- Summary: μ-relations hold for 107,374,182,400 / 107,374,182,400 quartets (19,658,705 of them non-trivial); joint hits mean 3.10, variance 2.41 (population; sample 2.43); histogram 1, 16, 21, 25, 18, 11, 6, 2 runs with 0..7 hits; random baseline 0 in every run; control condition: small AES 3.10, random function 3.27 hits per run.

**`Results/quick_seed20261002.csv`** and **`Results/medium_seed20261002.csv`** are the reference outputs of `make quick` and `make medium`.

## Version history

The paper's footnote cites this repository (https://github.com/shb115/joint-zero-diff-aes). The version used for the paper, which existed when the paper was written and produced the data of Table 3, is commit [`542ca32`](https://github.com/shb115/joint-zero-diff-aes/tree/542ca32ff12516e3eb92ffbe06686c4f342cc165) (`542ca32ff12516e3eb92ffbe06686c4f342cc165`). Its program and data are unchanged since commit `610c2b2`, which added `Results/results.csv`.

Changes since that version:
- **Experiment 2.** The AES hits are now tested on all four ciphertext differences, as described in the paper and as already done for the random baseline. The previous version tested only C0⊕C2 and C0⊕C3 and relied on the μ-relations for C1⊕C3 and C1⊕C2.
  - By Theorem 3 the two conditions give the same counts.
  - In `Results/results_v2.csv`, `joint_aes == half_aes` in all 100 runs, and the μ-relations held for all 100·2^30 quartets.
  - For `Results/results.csv` this equality follows from Theorem 3; it cannot be re-checked, since that run's seed is unknown.
- **New CSV columns:** `mu_nontrivial`, `half_aes`, `half_rand`, `trials`, `seed_mode`, `seed`, `dx`.
- **Seeding.** Per-run seeding makes runs reproducible and parallelizable. `--legacy-seed` keeps the original scheme.
- **Command line.** Options are validated, and an existing output file is not overwritten without `--force`.
- **CSV output.** Each row is written as soon as its run finishes, with LF line endings on every platform.
- **Exit codes.** These are listed under *Building and running*.
- **Standard output.** The per-run table and the "theoretical advantage" line were removed; the summary prints the control counts. The comment in the source header now refers to Theorem 3 (it said Theorem 4).
- **New files:** `Makefile`, `scripts/summarize.py`, `scripts/run_parallel.sh`, the reference outputs in `Results/`, `.gitattributes` and `.gitignore`.

The cipher, the sampling of keys, differences and plaintexts, and the order in which the PRNG is used are unchanged.

## Notes

- **State layout.** The code stores cell (row r, column c) at index 4r+c, while the paper uses r+4c. `DIAG[j]` and `INV_DIAG[j]` are exactly the paper's j-th diagonal and inverse diagonal.
- **MC⁻¹.** The paper writes the conditions on MC⁻¹(C_i ⊕ C_j). Its μ is ν∘L⁻¹ (Sec. 2.3), applied to 4-round AES in the superbox model R⁴ = L'∘S'∘L'∘S' with L' = AK∘MC (Sec. 2.3), whose 4th round has MixColumns. The cipher here omits that last MixColumns, as AES does (Sec. 2.1). Its C_i ⊕ C_j therefore already equals the difference before that MixColumns (the last round key cancels in differences), and the zero inverse diagonals are tested on it directly. Applying MC⁻¹ to these ciphertexts would be wrong.
- **Random baseline.** For each quartet the baseline draws four independent uniform 64-bit values, i.e. an ideal random function. This is the model used in the paper's analysis of a random permutation. For four distinct inputs, a random permutation and a random function have output distributions within statistical distance about 6·2^-64 ≈ 2^-61.
- **Joint condition.** A hit requires the two sets of zero inverse diagonals (one per pair of differences) to be non-empty and disjoint. This differs from "a disjoint pair exists" only when the sets overlap and one of them has a further element.
  - For small-scale AES with the provided Δx' this never happens. Before the MixColumns between the two superbox layers, the differences for Δx' and Δx⊕Δx' have at most two active nibbles per column. MixColumns is MDS, so each non-zero column afterwards has at most one zero nibble. Hence at most one inverse diagonal of each of C0⊕C2, C1⊕C3, C0⊕C3 and C1⊕C2 can be zero, and the two definitions coincide (also for `half_aes`).
  - For the random baseline (`joint_rand`) they differ with probability about 2^-75 per quartet.
  - For the control count `half_rand` they differ with probability about 2^-43 per quartet (< 0.01 expected over 100·2^30 quartets).
- **Variance.** The variance in Table 3 is the population variance (divided by N). `summarize.py` prints both the population and the sample variance.
- **Erratum (Sec. 3.3, paragraph beginning "The observed mean of 2.86").** The expected number of random-permutation hits is 12 × 2^-64 × 2^30 ≈ 7.0 × 10^-10 per run, not 2.0 × 10^-8 as printed. No conclusion is affected.

## Source code organization

```
.
├── Code/joint_property_experiment.c   experiment program (single file)
├── scripts/summarize.py               Table 3 and checks from one or more CSV files
├── scripts/run_parallel.sh            split runs over several processes
├── Results/results.csv                data of Table 3 (commit 542ca32, used for the paper)
├── Results/results_v2.csv             full re-run with this version (seed 20261002)
├── Results/quick_seed20261002.csv     reference output of `make quick`
├── Results/medium_seed20261002.csv    reference output of `make medium`
├── Makefile
├── LICENSE
├── README.md
├── .gitattributes                     LF line endings on every platform
└── .gitignore
```

Sections of `joint_property_experiment.c`:

| Section | Contents |
|---|---|
| Small-scale AES primitives | `SBOX`, `gf_mul`, `MDS`, `init_tables` (combined SubBytes/MixColumns table) |
| State layout | `DIAG`, `INV_DIAG` |
| Encryption | `round_mc`, `round_last`, `encrypt` |
| Zero inverse-diagonal patterns | `zmask`, `joint_hit` |
| PRNG | xoroshiro128++, splitmix64 seeding, legacy seeding, `rng_skip_run` |
| One run | `run_experiment`: key and difference sampling, Experiments 1 and 2, random baseline |
| Self-test | `sbox_affine`, `encrypt_reference`, `self_test` (`--self-test`) |
| Command line | `usage`, `parse_u64` |
| Main | option handling, CSV output, summary |

**Extending the experiment.** The code is written for the paper's experiment; the following are starting points for extensions.
- To change the choice of Δx', edit the diagonals copied from Δx in `run_experiment`.
  - Δx' must stay related to Δx at the diagonal level, as in Theorem 3: each diagonal of Δx' is either entirely zero or equal to that diagonal of Δx. Otherwise, for example if Δx' takes only some nibbles of a diagonal, the μ-relations do not hold and `joint_exp` exits with status 3.
  - Zero diagonals {1, 3} give the same experiment as {0, 2} with P2 and P3 swapped.
- To change the hit condition, edit `joint_hit`.
- To change the cipher, edit the *Small-scale AES primitives* (`SBOX`, `gf_mul`, `MDS`) and *Encryption* sections. The self-test (`sbox_affine` and the GF(2^4) inverse check in `self_test`) is specific to the small-scale AES S-box: adapt or remove that check, otherwise `make test` fails. 4-bit cells and a 64-bit state are also assumed in:
  - `rand_nibbles` (16 nibbles from one 64-bit PRNG output);
  - `sb_mc`/`init_tables` (4-bit fields packed into a 16-bit column);
  - the Δx sampling (`& 0xF`);
  - the hex output of `dx`.
- Update the expected values: `lam_aes`, `lam_ctl` and `lam_rnd` in `scripts/summarize.py`, and `lambda_aes` and `lambda_rand` in `main()` of the program (used only in its printed summary):
  - `lam_aes` depends on Δx' and on the hit condition. It must be re-derived and can be 0; for example, it is 0 when Δx' or Δx⊕Δx' has only one active diagonal, by the branch number of MixColumns.
  - `lam_ctl` (control condition) and `lam_rnd` (random baseline) depend only on the hit condition.
- A change that alters what the experiment computes generally changes the outputs, so the reference comparisons of `make quick`/`medium`/`full` with the default seed fail by design. Use `SEED=<another seed>` (or `run_parallel.sh` without `-r`) and rely on the statistical checks.
- If you change how `run_experiment` uses the PRNG, update `rng_skip_run` as well. It is used only by `--legacy-seed` with `--first-run` > 1.

## Third-party material

- **PRNG.** `rng_next` is xoroshiro128++ 1.0 by D. Blackman and S. Vigna, and `splitmix64` is by S. Vigna. Both were dedicated to the public domain by their authors: https://prng.di.unimi.it/. Apart from their names, they are modified only to take their state by pointer instead of from a global variable: `rng_next` takes a struct (`rng_t`), `splitmix64` a `uint64_t`.
- **Cipher components.** The S-box and the MixColumns matrix are those of small-scale AES with a 4×4 state of 4-bit words (the cipher is the variant SR*(4,4,4,4), without MixColumns in the last round, with independent round keys): C. Cid, S. Murphy, M. Robshaw, *Small Scale Variants of the AES*, FSE 2005.
  - As in AES, the last round has no MixColumns.
  - The paper's "fresh random key" is implemented as independent uniform round keys. Theorem 3 holds for any round keys.

## Use of AI tools

This statement follows the ASIACRYPT 2026 AI tool policy.

- **Tool.** Claude Code (Anthropic).
- **Assistance received.**
  - Writing the experiment program: the small-scale AES implementation, Experiments 1 and 2, the random baseline, the command-line interface, seeding, and the self-test.
  - Writing `scripts/`, the `Makefile` and this README.
  - Checking the code and the documentation against the paper and the ASIACRYPT 2026 artifact requirements, and running the experiments for `Results/results_v2.csv` and the reference outputs.
- The authors reviewed and validated the results and take full responsibility for the contents of this repository.

## License

MIT License, see [LICENSE](LICENSE).

## Contact

Hanbeom Shin (newonetiger@korea.ac.kr) and Dongjae Lee (dongjae.lee@kangwon.ac.kr).
