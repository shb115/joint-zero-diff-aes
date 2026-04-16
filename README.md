# Experimental Verification of the Joint Generalized Zero-Difference Property

This repository contains the source code and experimental results for verifying the **Joint Generalized Zero-Difference Property** on small-scale AES.

## Overview

The Joint Generalized Zero-Difference Property states that a quartet constructed from related differences satisfies three distinct generalized zero-difference relations (mu-relations) simultaneously through a 2-round SPN. This repository provides experimental verification of this property using small-scale AES (4-bit S-box, 64-bit block size, 4 rounds).

Two experiments are performed:

- **Experiment 1 (Deterministic mu-relation):** Verifies that all three mu-relations hold simultaneously for every quartet formed by related differences. Across 10^8 quartets over 100 independent keys, a 100% success rate is expected.

- **Experiment 2 (Joint zero inverse-diagonal hits):** Counts quartets whose ciphertext differences exhibit joint zero inverse-diagonal patterns at disjoint positions. Under small-scale AES, the expected hit count per run follows a Poisson distribution with lambda = 12 * 2^{-32} * 2^{30} = 3.0. Under a random permutation, the expected count is negligible (~2.0 * 10^{-8}).

## Repository Structure

```
.
├── README.md
├── Code/
│   └── joint_property_experiment.c    # Main experiment (100 independent keys)
└── Results/
    └── results.csv                    # Experimental results (100 runs)
```

## Build and Run

### Requirements

- C compiler (GCC recommended)
- Standard C math library (`-lm`)

### Compilation

```bash
gcc -O2 -o joint_exp Code/joint_property_experiment.c -lm
```

### Execution

```bash
./joint_exp [output.csv]
```

- Default output file: `results.csv`
- Each run encrypts 2^{30} quartets with a fresh random key and fresh random related differences.
- Total: 100 runs (configurable via `NUM_KEYS` in the source code).
- Approximate runtime: \~12 minutes per run (\~20 hours total for 100 runs).

## Configuration

The following parameters can be adjusted in `Code/joint_property_experiment.c`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `NUM_ROUNDS` | 4 | Number of AES rounds (4 = 2-round SPN in superbox representation) |
| `TRIALS` | 2^{30} | Number of quartets per key |
| `MU_CHECK_TRIALS` | 10^6 | Number of quartets for mu-relation verification |
| `NUM_KEYS` | 100 | Number of independent key experiments |

## Output Format (CSV)

| Column | Description |
|--------|-------------|
| `run` | Run index (1-based) |
| `mu_verified` | Number of quartets passing all three mu-relations |
| `mu_total` | Total quartets tested for mu-relation |
| `single_xp` | Quartets with at least one zero inverse diagonal for Delta x' |
| `single_xxp` | Quartets with at least one zero inverse diagonal for Delta x + Delta x' |
| `joint_aes` | Joint hits under small-scale AES |
| `joint_rand` | Joint hits under random permutation |
| `elapsed_sec` | Runtime in seconds |

## Summary of Results

Over 100 independent experiments (2^{30} quartets each):

| Metric | Small-scale AES | Random Permutation |
|--------|----------------:|-------------------:|
| mu-relation verified | 10^8 / 10^8 (100%) | --- |
| Joint hits (mean) | 2.86 | 0.00 |
| Joint hits (variance) | 2.38 | 0.00 |
| Expected (Poisson lambda) | 3.0 | ~0 |

The observed distribution of joint hit counts closely matches the Poisson distribution with lambda = 3.0, confirming the theoretical predictions.

## Small-Scale AES Specification

- **Block size:** 64 bits (4x4 state of 4-bit nibbles)
- **S-box:** 4-bit affine S-box over GF(2^4) with irreducible polynomial x^4 + x + 1
- **MixColumns:** Same MDS matrix as full AES, operating over GF(2^4)
- **Key schedule:** Independent random round keys
- **Rounds:** 4 (last round without MixColumns)

## License

This code is provided for academic and research purposes.
