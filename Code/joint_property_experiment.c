/**
 * @file joint_property_experiment.c
 * @brief Experimental verification of the Joint Generalized Zero-Difference
 *        Property on 4-round small-scale AES (Section 3.3 and Table 3 of the
 *        paper).
 *
 * Paper: H. Shin, S. Kim, B. Seok, D. Hong, J. Sung, S. Hong, S. Lee, D. Lee,
 *   "Key-Independent Secret-Key Distinguisher for 7-Round AES based on the
 *    Joint Generalized Zero-Difference Property", ASIACRYPT 2026.
 *
 * Each run samples fresh independent round keys and fresh related
 * differences (dx, dx'), and encrypts quartets
 *   (P0, P1, P2, P3) = (a, a ^ dx, a ^ dx', a ^ dx ^ dx')
 * for uniformly random a.
 *
 *  - Experiment 1 (Theorem 3): the three mu-relations
 *      (C0^C1 vs C2^C3), (C0^C2 vs C1^C3), (C0^C3 vs C1^C2)
 *    hold for every checked quartet.
 *  - Experiment 2: count quartets for which C0^C2 and C1^C3 share a zero
 *    inverse diagonal and C0^C3 and C1^C2 share a zero inverse diagonal at a
 *    disjoint position (precisely: the two sets of shared zero inverse
 *    diagonals are non-empty and disjoint, see joint_hit). Expected
 *    12 * 2^-32 per quartet for small-scale AES, 12 * 2^-64 for an ideal
 *    random function.
 *
 * Build: gcc -O2 -std=c99 -o joint_exp Code/joint_property_experiment.c -lm
 * Usage: ./joint_exp --help
 */

#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <math.h>

// ===========================================================================
// Compile-time defaults (all but NUM_ROUNDS can also be set on the command line)
// ===========================================================================
#ifndef NUM_ROUNDS
#define NUM_ROUNDS          4         /* 4 rounds = 2-round SPN (superbox view) */
#endif
#ifndef DEFAULT_LOG2_TRIALS
#define DEFAULT_LOG2_TRIALS 30        /* quartets per run (paper: 2^30)         */
#endif
#ifndef DEFAULT_MU_TRIALS
#define DEFAULT_MU_TRIALS   1000000ULL /* quartets per run checked for mu (10^6) */
#endif
#ifndef DEFAULT_NUM_KEYS
#define DEFAULT_NUM_KEYS    100       /* independent runs (paper: 100)          */
#endif

// ===========================================================================
// Small-scale AES primitives (4-bit S-box, GF(2^4) with x^4 + x + 1)
// ===========================================================================

/** S-box of small-scale AES with 4-bit words [Cid, Murphy, Robshaw, FSE 2005]:
 *  inversion in GF(2^4) followed by an affine map. */
static const uint8_t SBOX[16] = {
    0x6, 0xB, 0x5, 0x4, 0x2, 0xE, 0x7, 0xA,
    0x9, 0xD, 0xF, 0xC, 0x3, 0x1, 0x0, 0x8
};

/** Multiplication in GF(2^4) = GF(2)[x]/(x^4 + x + 1). */
static uint8_t gf_mul(uint8_t a, uint8_t b) {
    uint8_t p = 0;
    for (int i = 0; i < 4; i++) {
        if (b & 1) p ^= a;
        int hi = a & 8;
        a = (uint8_t)((a << 1) & 0xF);
        if (hi) a ^= 0x3;
        b >>= 1;
    }
    return p;
}

/** MixColumns matrix circ(2, 3, 1, 1), as in AES but over GF(2^4). */
static const uint8_t MDS[4][4] = {
    {2, 3, 1, 1}, {1, 2, 3, 1}, {1, 1, 2, 3}, {3, 1, 1, 2}
};

/** sb_mc[row][v]: MixColumns output column (4 packed nibbles, row 0 in the
 *  most significant nibble) for S-box input v placed in row `row`. */
static uint16_t sb_mc[4][16];

static void init_tables(void) {
    for (int row = 0; row < 4; row++)
        for (int v = 0; v < 16; v++) {
            uint8_t s = SBOX[v];
            uint16_t p = 0;
            for (int r = 0; r < 4; r++)
                p |= (uint16_t)(gf_mul(MDS[r][row], s) << (4 * (3 - r)));
            sb_mc[row][v] = p;
        }
}

// ===========================================================================
// State layout
// ===========================================================================
// state[4*r + c] holds the nibble in row r, column c (row-major). The paper
// indexes cell (r, c) as r + 4c; the cell sets below are nevertheless exactly
// the paper's j-th diagonal {(r, r+j mod 4)} and j-th inverse diagonal
// {(r, j-r mod 4)}. ShiftRows maps diagonal j to column j, and column j to
// inverse diagonal j.

static const int DIAG[4][4] = {
    {0, 5, 10, 15}, {1, 6, 11, 12}, {2, 7, 8, 13}, {3, 4, 9, 14}
};

static const int INV_DIAG[4][4] = {
    {0, 7, 10, 13}, {1, 4, 11, 14}, {2, 5, 8, 15}, {3, 6, 9, 12}
};

// ===========================================================================
// Encryption
// ===========================================================================

/** Full round: SubBytes, ShiftRows, MixColumns, AddRoundKey. */
static void round_mc(uint8_t *state, const uint8_t *rk) {
    uint8_t tmp[16];
    for (int c = 0; c < 4; c++) {
        uint16_t col = sb_mc[0][state[DIAG[c][0]]]
                     ^ sb_mc[1][state[DIAG[c][1]]]
                     ^ sb_mc[2][state[DIAG[c][2]]]
                     ^ sb_mc[3][state[DIAG[c][3]]];
        tmp[c]      = (uint8_t)(((col >> 12) & 0xF) ^ rk[c]);
        tmp[4 + c]  = (uint8_t)(((col >> 8)  & 0xF) ^ rk[4 + c]);
        tmp[8 + c]  = (uint8_t)(((col >> 4)  & 0xF) ^ rk[8 + c]);
        tmp[12 + c] = (uint8_t)((col         & 0xF) ^ rk[12 + c]);
    }
    memcpy(state, tmp, 16);
}

/** Last round: SubBytes, ShiftRows, AddRoundKey (no MixColumns, as in AES). */
static void round_last(uint8_t *state, const uint8_t *rk) {
    uint8_t tmp[16];
    for (int i = 0; i < 16; i++) tmp[i] = SBOX[state[i]];
    uint8_t tmp2[16];
    for (int r = 0; r < 4; r++)
        for (int c = 0; c < 4; c++)
            tmp2[r * 4 + c] = tmp[r * 4 + (c + r) % 4];
    for (int i = 0; i < 16; i++) state[i] = tmp2[i] ^ rk[i];
}

/** NUM_ROUNDS-round small-scale AES with independent round keys rk[0..NUM_ROUNDS].
 *
 *  The paper writes the conditions on MC^-1(C_i ^ C_j) for its superbox model
 *  R^4 = L'S'L'S' with L' = AK o MC (Sec. 2.3), whose 4th round has MixColumns.
 *  The last round here omits that MixColumns (as in AES), so C_i ^ C_j already
 *  equals the difference before it (the last round key cancels in differences)
 *  and the zero inverse diagonals are tested on it directly. Applying MC^-1 to
 *  these ciphertexts would be wrong. */
static void encrypt(uint8_t *ct, const uint8_t *pt, uint8_t rk[][16]) {
    memcpy(ct, pt, 16);
    for (int i = 0; i < 16; i++) ct[i] ^= rk[0][i];
    for (int r = 1; r < NUM_ROUNDS; r++)
        round_mc(ct, rk[r]);
    round_last(ct, rk[NUM_ROUNDS]);
}

// ===========================================================================
// Zero inverse-diagonal patterns
// ===========================================================================

/** Bit c of the result is set iff inverse diagonal c of s is zero. */
static inline unsigned zmask(const uint8_t *s) {
    unsigned m = 0;
    for (int c = 0; c < 4; c++)
        if ((s[INV_DIAG[c][0]] | s[INV_DIAG[c][1]]
           | s[INV_DIAG[c][2]] | s[INV_DIAG[c][3]]) == 0)
            m |= 1u << c;
    return m;
}

/** Joint condition of Experiment 2: both zero sets are non-empty and disjoint.
 *  This differs from "some i != j with i in z1, j in z3" only when the sets
 *  overlap and one of them has a further element: never for small-scale AES
 *  with the provided dx' (each set has at most one element, by the MDS
 *  property of MixColumns), about 2^-75 per quartet for the random baseline
 *  and about 2^-43 for the control count half_rand. */
static inline int joint_hit(unsigned z1, unsigned z3) {
    return z1 && z3 && !(z1 & z3);
}

// ===========================================================================
// PRNG
// ===========================================================================

typedef struct { uint64_t s[2]; } rng_t;

static inline uint64_t rotl64(uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

/** xoroshiro128++ 1.0 by D. Blackman and S. Vigna (public domain,
 *  https://prng.di.unimi.it/xoroshiro128plusplus.c). */
static inline uint64_t rng_next(rng_t *g) {
    uint64_t s0 = g->s[0], s1 = g->s[1];
    uint64_t res = rotl64(s0 + s1, 17) + s0;
    s1 ^= s0;
    g->s[0] = rotl64(s0, 49) ^ s1 ^ (s1 << 21);
    g->s[1] = rotl64(s1, 28);
    return res;
}

/** 16 uniform nibbles from one PRNG output. */
static void rand_nibbles(rng_t *g, uint8_t *s) {
    uint64_t r = rng_next(g);
    for (int i = 0; i < 16; i++)
        s[i] = (uint8_t)((r >> (i * 4)) & 0xF);
}

/** splitmix64 (S. Vigna, public domain), used only to seed xoroshiro128++. */
static uint64_t splitmix64(uint64_t *x) {
    uint64_t z = (*x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

/** Default seeding: run k (1-based) gets its own stream, taken from outputs
 *  2k-1 and 2k of splitmix64 seeded with `seed`. A run therefore depends only
 *  on (seed, k), so runs can be split across processes or machines. */
static void rng_seed_run(rng_t *g, uint64_t seed, uint64_t run) {
    uint64_t x = seed + 2 * (run - 1) * 0x9E3779B97F4A7C15ULL;
    g->s[0] = splitmix64(&x);
    g->s[1] = splitmix64(&x);
    if ((g->s[0] | g->s[1]) == 0) g->s[1] = 1;
}

/** Legacy seeding of the original code (commit 542ca32), which produced
 *  Results/results.csv: one stream for all runs, seeded from a Unix time t. */
static void rng_seed_legacy(rng_t *g, uint64_t t) {
    g->s[0] = t * 6364136223846793005ULL + 1;
    g->s[1] = g->s[0] ^ 0xdeadbeefcafebabeULL;
    for (int i = 0; i < 20; i++) rng_next(g);
}

/** Advance a legacy stream past one whole run without computing it
 *  (same PRNG consumption as run_experiment). */
static void rng_skip_run(rng_t *g, uint64_t trials) {
    for (int r = 0; r <= NUM_ROUNDS; r++) rng_next(g);
    for (int i = 0; i < 16; i++) {
        while (!(rng_next(g) & 0xF)) { }
    }
    for (uint64_t t = 0; t < 5 * trials; t++) rng_next(g);
}

// ===========================================================================
// One run (one key)
// ===========================================================================

typedef struct {
    uint64_t  run;            /* 1-based run index                              */
    uint64_t  mu_ok;          /* quartets satisfying all three mu-relations     */
    uint64_t  mu_total;       /* quartets checked for the mu-relations          */
    uint64_t  mu_nontrivial;  /* checked quartets with any zero inv. diagonal   */
    uint64_t  single_xp;      /* C0^C2 has a zero inverse diagonal              */
    uint64_t  single_xxp;     /* C0^C3 has a zero inverse diagonal              */
    uint64_t  hit_aes;        /* joint hits, all four differences (paper)       */
    uint64_t  hit_rand;       /* same for the random-function baseline          */
    uint64_t  half_aes;       /* control: only C0^C2 and C0^C3 are tested       */
    uint64_t  half_rand;      /* same control for the random-function baseline  */
    long long elapsed;        /* seconds                                        */
    uint8_t   dx[16];
} RunResult;

/**
 * Runs one experiment with fresh round keys and related differences.
 *
 * The PRNG is consumed in the same order as in the original code: 5 outputs
 * for the round keys, rejection sampling for dx, then per quartet one output
 * for `a` followed by four outputs for the random baseline.
 */
static void run_experiment(rng_t *g, uint64_t run, uint64_t trials,
                           uint64_t mu_trials, RunResult *res) {
    uint8_t rk[NUM_ROUNDS + 1][16];
    for (int r = 0; r <= NUM_ROUNDS; r++) rand_nibbles(g, rk[r]);

    // Related differences: dx has no zero nibble; dx' equals dx on diagonals
    // 1 and 3 and is zero on diagonals 0 and 2 (paper, Sec. 3.3 Setup).
    uint8_t dx[16], dxp[16], dxxp[16];
    for (int i = 0; i < 16; i++) {
        do { dx[i] = (uint8_t)(rng_next(g) & 0xF); } while (!dx[i]);
    }
    memset(dxp, 0, 16);
    for (int k = 0; k < 4; k++) {
        dxp[DIAG[1][k]] = dx[DIAG[1][k]];
        dxp[DIAG[3][k]] = dx[DIAG[3][k]];
    }
    for (int i = 0; i < 16; i++) dxxp[i] = dx[i] ^ dxp[i];

    memset(res, 0, sizeof *res);
    res->run = run;
    memcpy(res->dx, dx, 16);

    printf("  Run %3llu: dx=", (unsigned long long)run);
    for (int i = 0; i < 16; i++) printf("%X", dx[i]);
    printf("\n");
    fflush(stdout);

    time_t t0 = time(NULL);
    uint64_t quarter = trials / 4;

    for (uint64_t t = 0; t < trials; t++) {
        uint8_t a[16];
        rand_nibbles(g, a);

        uint8_t P0[16], P1[16], P2[16], P3[16];
        for (int i = 0; i < 16; i++) {
            P0[i] = a[i];
            P1[i] = a[i] ^ dx[i];
            P2[i] = a[i] ^ dxp[i];
            P3[i] = a[i] ^ dxxp[i];
        }

        uint8_t C0[16], C1[16], C2[16], C3[16];
        encrypt(C0, P0, rk);
        encrypt(C1, P1, rk);
        encrypt(C2, P2, rk);
        encrypt(C3, P3, rk);

        uint8_t d1[16], d2[16], d3[16], d4[16];
        for (int i = 0; i < 16; i++) {
            d1[i] = C0[i] ^ C2[i];   /* dx'       */
            d2[i] = C1[i] ^ C3[i];   /* dx'       */
            d3[i] = C0[i] ^ C3[i];   /* dx ^ dx'  */
            d4[i] = C1[i] ^ C2[i];   /* dx ^ dx'  */
        }
        unsigned m1 = zmask(d1), m2 = zmask(d2), m3 = zmask(d3), m4 = zmask(d4);

        // Experiment 1: mu-relations (Theorem 3)
        if (t < mu_trials) {
            uint8_t d5[16], d6[16];
            for (int i = 0; i < 16; i++) {
                d5[i] = C0[i] ^ C1[i];   /* dx */
                d6[i] = C2[i] ^ C3[i];   /* dx */
            }
            unsigned m5 = zmask(d5), m6 = zmask(d6);
            res->mu_total++;
            if (m1 == m2 && m3 == m4 && m5 == m6) res->mu_ok++;
            if (m1 | m2 | m3 | m4 | m5 | m6) res->mu_nontrivial++;
        }

        // Experiment 2: joint zero inverse-diagonal hits (small-scale AES)
        if (m1) res->single_xp++;
        if (m3) res->single_xxp++;
        if (joint_hit(m1 & m2, m3 & m4)) res->hit_aes++;
        if (joint_hit(m1, m3))           res->half_aes++;

        // Random-function baseline: four independent uniform 64-bit outputs
        uint8_t R0[16], R1[16], R2[16], R3[16];
        rand_nibbles(g, R0); rand_nibbles(g, R1);
        rand_nibbles(g, R2); rand_nibbles(g, R3);

        uint8_t rd1[16], rd2[16], rd3[16], rd4[16];
        for (int i = 0; i < 16; i++) {
            rd1[i] = R0[i] ^ R2[i];
            rd2[i] = R1[i] ^ R3[i];
            rd3[i] = R0[i] ^ R3[i];
            rd4[i] = R1[i] ^ R2[i];
        }
        unsigned r1 = zmask(rd1), r2 = zmask(rd2), r3 = zmask(rd3), r4 = zmask(rd4);
        if (joint_hit(r1 & r2, r3 & r4)) res->hit_rand++;
        if (joint_hit(r1, r3))           res->half_rand++;

        if (quarter && (t + 1) % quarter == 0 && t + 1 < trials) {
            printf("  Run %3llu: %3d%%  joint AES=%llu  random=%llu  (%.0f s)\n",
                   (unsigned long long)run, (int)(100 * ((t + 1) / quarter) / 4),
                   (unsigned long long)res->hit_aes,
                   (unsigned long long)res->hit_rand,
                   difftime(time(NULL), t0));
            fflush(stdout);
        }
    }

    res->elapsed = (long long)difftime(time(NULL), t0);
    printf("  Run %3llu: done  joint AES=%llu  random=%llu  mu %llu/%llu  (%lld s)\n",
           (unsigned long long)run,
           (unsigned long long)res->hit_aes,
           (unsigned long long)res->hit_rand,
           (unsigned long long)res->mu_ok,
           (unsigned long long)res->mu_total,
           res->elapsed);
    fflush(stdout);
}

// ===========================================================================
// Self-test of the cipher implementation (--self-test)
// ===========================================================================

/** Affine map of the small-scale AES S-box [Cid, Murphy, Robshaw, FSE 2005]:
 *  matrix rows 1011, 1101, 1110, 0111 (from the most significant output bit),
 *  constant 0x6. ROW[j] is the input mask of output bit j. */
static uint8_t sbox_affine(uint8_t y) {
    static const uint8_t ROW[4] = {0x7, 0xE, 0xD, 0xB};
    uint8_t out = 0;
    for (int j = 0; j < 4; j++) {
        uint8_t b = (uint8_t)(ROW[j] & y);
        b ^= b >> 2; b ^= b >> 1;
        out |= (uint8_t)((b & 1) << j);
    }
    return out ^ 0x6;
}

/** Straightforward implementation of the same cipher (SubBytes, ShiftRows,
 *  MixColumns with gf_mul, AddRoundKey on a 4x4 array), used to check the
 *  table-based encrypt(). */
static void encrypt_reference(uint8_t *ct, const uint8_t *pt, uint8_t rk[][16]) {
    uint8_t s[4][4], t[4][4];
    for (int r = 0; r < 4; r++)
        for (int c = 0; c < 4; c++) s[r][c] = pt[4 * r + c] ^ rk[0][4 * r + c];
    for (int round = 1; round <= NUM_ROUNDS; round++) {
        for (int r = 0; r < 4; r++)                          /* SubBytes, ShiftRows */
            for (int c = 0; c < 4; c++) t[r][c] = SBOX[s[r][(c + r) % 4]];
        for (int c = 0; c < 4; c++)                          /* MixColumns, AddRoundKey */
            for (int r = 0; r < 4; r++) {
                uint8_t v = t[r][c];
                if (round < NUM_ROUNDS) {
                    v = 0;
                    for (int k = 0; k < 4; k++) v ^= gf_mul(MDS[r][k], t[k][c]);
                }
                s[r][c] = v ^ rk[round][4 * r + c];
            }
    }
    for (int r = 0; r < 4; r++)
        for (int c = 0; c < 4; c++) ct[4 * r + c] = s[r][c];
}

/** Returns 0 if (1) SBOX equals the affine map applied to the inverse in
 *  GF(2^4) and (2) encrypt() agrees with encrypt_reference() on random
 *  inputs; 1 otherwise. Uses its own PRNG stream. */
static int self_test(void) {
    int fail = 0;
    for (int x = 0; x < 16; x++) {
        uint8_t inv = 0;
        for (int y = 1; y < 16; y++)
            if (gf_mul((uint8_t)x, (uint8_t)y) == 1) inv = (uint8_t)y;
        if (SBOX[x] != sbox_affine(inv)) fail = 1;
    }
    printf("S-box = affine map of the GF(2^4) inverse (small-scale AES): %s\n",
           fail ? "FAIL" : "PASS");

    const int n = 100000;
    int bad = 0;
    rng_t g;
    rng_seed_run(&g, 0x5E1F7E57ULL, 1);
    for (int i = 0; i < n; i++) {
        uint8_t rk[NUM_ROUNDS + 1][16], pt[16], c1[16], c2[16];
        for (int r = 0; r <= NUM_ROUNDS; r++) rand_nibbles(&g, rk[r]);
        rand_nibbles(&g, pt);
        encrypt(c1, pt, rk);
        encrypt_reference(c2, pt, rk);
        if (memcmp(c1, c2, 16)) bad++;
    }
    printf("Table-based encryption = reference implementation on %d random inputs: %s\n",
           n, bad ? "FAIL" : "PASS");
    return fail || bad;
}

// ===========================================================================
// Command line
// ===========================================================================

static void usage(const char *prog) {
    printf(
"Usage: %s [options] [output.csv]\n"
"\n"
"  -o, --output FILE     CSV output file (default: results.csv)\n"
"  -k, --keys N          number of runs, i.e. independent keys (default: %d)\n"
"  -f, --first-run K     index of the first run (default: 1); runs K..K+N-1\n"
"  -t, --log2-trials T   quartets per run = 2^T, 1 <= T <= 40 (default: %d)\n"
"  -m, --mu-trials M     quartets per run checked for the mu-relations\n"
"                        (default: %llu; 0 = all quartets)\n"
"  -s, --seed S          master seed (default: current time; always recorded)\n"
"      --legacy-seed T   original single-stream seeding with Unix time T\n"
"                        (the scheme that produced Results/results.csv)\n"
"  -q, --quick           smoke test: 4 runs of 2^24 quartets (explicit -k/-t win)\n"
"      --force           overwrite an existing output file\n"
"      --self-test       check the cipher implementation and exit\n"
"  -h, --help            show this help\n"
"  (--self-test and --help act as soon as they are read: later options and the\n"
"   range checks are skipped; earlier options must still be valid)\n"
"\n"
"Exit status: 0 on success, 1 on an I/O error, an existing output file or a\n"
"failed self-test, 2 on a usage error, 3 if some quartet violated a\n"
"mu-relation. Each option may be given at most once.\n"
"\n"
"Runs with the same --seed are identical whether executed in one process or\n"
"split over several (see scripts/run_parallel.sh).\n",
        prog, DEFAULT_NUM_KEYS, DEFAULT_LOG2_TRIALS,
        (unsigned long long)DEFAULT_MU_TRIALS);
}

/** Parses a non-negative decimal integer; rejects signs, spaces and overflow. */
static int parse_u64(const char *s, uint64_t *out) {
    char *end = NULL;
    if (!s || !isdigit((unsigned char)*s)) return 0;
    errno = 0;
    unsigned long long v = strtoull(s, &end, 10);
    if (errno == ERANGE || *end != '\0') return 0;
    *out = (uint64_t)v;
    return 1;
}

// ===========================================================================
// Main
// ===========================================================================

int main(int argc, char *argv[]) {
    const char *csv_path   = "results.csv";
    uint64_t    num_keys   = DEFAULT_NUM_KEYS;
    uint64_t    first_run  = 1;
    uint64_t    log2_trials = DEFAULT_LOG2_TRIALS;
    uint64_t    mu_trials  = DEFAULT_MU_TRIALS;
    uint64_t    seed       = 0;
    int         have_seed  = 0, legacy = 0, force = 0, quick = 0;
    int         set_keys   = 0, set_trials = 0, have_path = 0;
    int         seen_first = 0, seen_mu = 0, seen_seed = 0;

    for (int i = 1; i < argc; i++) {
        const char *o = argv[i];
        const char *v = (i + 1 < argc) ? argv[i + 1] : NULL;
        int ok = 1;
        if (!strcmp(o, "-h") || !strcmp(o, "--help")) { usage("joint_exp"); return 0; }
        else if (!strcmp(o, "-o") || !strcmp(o, "--output"))      { ok = v && v[0] != '-' && !have_path; csv_path = v; have_path = 1; i++; }
        else if (!strcmp(o, "-k") || !strcmp(o, "--keys"))        { ok = !set_keys && parse_u64(v, &num_keys); set_keys = 1; i++; }
        else if (!strcmp(o, "-f") || !strcmp(o, "--first-run"))   { ok = !seen_first && parse_u64(v, &first_run); seen_first = 1; i++; }
        else if (!strcmp(o, "-t") || !strcmp(o, "--log2-trials")) { ok = !set_trials && parse_u64(v, &log2_trials); set_trials = 1; i++; }
        else if (!strcmp(o, "-m") || !strcmp(o, "--mu-trials"))   { ok = !seen_mu && parse_u64(v, &mu_trials); seen_mu = 1; i++; }
        else if (!strcmp(o, "-s") || !strcmp(o, "--seed"))        { ok = !seen_seed && parse_u64(v, &seed); seen_seed = have_seed = 1; i++; }
        else if (!strcmp(o, "--legacy-seed"))                     { ok = !seen_seed && parse_u64(v, &seed); seen_seed = have_seed = legacy = 1; i++; }
        else if (!strcmp(o, "--self-test"))                       { init_tables(); return self_test(); }
        else if (!strcmp(o, "-q") || !strcmp(o, "--quick"))       { ok = !quick; quick = 1; }
        else if (!strcmp(o, "--force"))                           { ok = !force; force = 1; }
        else if (o[0] != '-')                                     { ok = !have_path; csv_path = o; have_path = 1; }
        else { fprintf(stderr, "Unknown option: %s (see --help)\n", o); return 2; }
        if (!ok) { fprintf(stderr, "Missing, invalid or repeated option or value: %s\n", o); return 2; }
    }
    if (quick) {
        if (!set_keys)   num_keys = 4;
        if (!set_trials) log2_trials = 24;
    }
    if (num_keys == 0 || first_run == 0 || first_run + num_keys - 1 < first_run
        || log2_trials < 1 || log2_trials > 40) {
        fprintf(stderr, "Invalid --keys, --first-run or --log2-trials (1..40)\n");
        return 2;
    }
    const uint64_t trials = 1ULL << log2_trials;
    if (mu_trials == 0 || mu_trials > trials) mu_trials = trials;
    if (!have_seed) seed = (uint64_t)time(NULL);

    // Open the output first so that a bad path fails immediately.
    FILE *chk = fopen(csv_path, "r");
    if (chk) {
        fclose(chk);
        if (!force) {
            fprintf(stderr, "Output file %s exists (use --force to overwrite)\n", csv_path);
            return 1;
        }
    }
    FILE *fp = fopen(csv_path, "wb");   /* binary: LF line endings on every platform */
    if (!fp) { perror(csv_path); return 1; }
    fprintf(fp, "run,mu_verified,mu_total,single_xp,single_xxp,joint_aes,joint_rand,"
                "elapsed_sec,mu_nontrivial,half_aes,half_rand,trials,seed_mode,seed,dx\n");
    fflush(fp);

    init_tables();

    printf("=============================================================\n");
    printf(" Joint Generalized Zero-Difference Property - small-scale AES\n");
    printf("=============================================================\n");
    printf(" Block: 64-bit (4-bit S-box), rounds: %d (last without MC)\n", NUM_ROUNDS);
    printf(" Runs: %llu..%llu, quartets per run: 2^%llu, mu-checked: %llu\n",
           (unsigned long long)first_run,
           (unsigned long long)(first_run + num_keys - 1),
           (unsigned long long)log2_trials, (unsigned long long)mu_trials);
    printf(" Seed: %llu (%s)\n", (unsigned long long)seed, legacy ? "legacy" : "per-run");
    printf(" Output: %s\n", csv_path);
    printf("-------------------------------------------------------------\n");
    fflush(stdout);

    rng_t g;
    if (legacy) {
        rng_seed_legacy(&g, seed);
        for (uint64_t k = 1; k < first_run; k++) {
            printf("  (legacy stream: skipping run %llu)\n", (unsigned long long)k);
            fflush(stdout);
            rng_skip_run(&g, trials);
        }
    }

    uint64_t sum_mu_ok = 0, sum_mu_total = 0, sum_nontriv = 0;
    double s_aes = 0, s_rand = 0, ss_aes = 0, ss_rand = 0, s_half_aes = 0, s_half_rand = 0;
    time_t total_t0 = time(NULL);

    for (uint64_t k = 0; k < num_keys; k++) {
        uint64_t run = first_run + k;
        RunResult r;
        if (!legacy) rng_seed_run(&g, seed, run);
        run_experiment(&g, run, trials, mu_trials, &r);

        fprintf(fp, "%llu,%llu,%llu,%llu,%llu,%llu,%llu,%lld,%llu,%llu,%llu,%llu,%s,%llu,",
                (unsigned long long)r.run,
                (unsigned long long)r.mu_ok, (unsigned long long)r.mu_total,
                (unsigned long long)r.single_xp, (unsigned long long)r.single_xxp,
                (unsigned long long)r.hit_aes, (unsigned long long)r.hit_rand,
                r.elapsed,
                (unsigned long long)r.mu_nontrivial,
                (unsigned long long)r.half_aes, (unsigned long long)r.half_rand,
                (unsigned long long)trials, legacy ? "legacy" : "run",
                (unsigned long long)seed);
        for (int i = 0; i < 16; i++) fprintf(fp, "%X", r.dx[i]);
        fprintf(fp, "\n");
        fflush(fp);

        sum_mu_ok += r.mu_ok; sum_mu_total += r.mu_total; sum_nontriv += r.mu_nontrivial;
        s_aes += (double)r.hit_aes;   ss_aes += (double)r.hit_aes * (double)r.hit_aes;
        s_rand += (double)r.hit_rand; ss_rand += (double)r.hit_rand * (double)r.hit_rand;
        s_half_aes += (double)r.half_aes; s_half_rand += (double)r.half_rand;
    }
    if (ferror(fp) | fclose(fp)) { perror(csv_path); return 1; }

    // =======================================================================
    // Summary (Table 3 in full, including the histogram: scripts/summarize.py)
    // =======================================================================
    double n = (double)num_keys;
    double mean_aes = s_aes / n, mean_rand = s_rand / n;
    double var_aes = ss_aes / n - mean_aes * mean_aes;
    double var_rand = ss_rand / n - mean_rand * mean_rand;
    double lambda_aes = (double)trials * 12.0 * pow(2, -32);
    double lambda_rand = (double)trials * 12.0 * pow(2, -64);

    printf("\n=============================================================\n");
    printf(" SUMMARY (%llu runs x 2^%llu quartets)\n",
           (unsigned long long)num_keys, (unsigned long long)log2_trials);
    printf("=============================================================\n");
    printf("[1] mu-relations: %llu / %llu quartets verified (%s), %llu non-trivial\n",
           (unsigned long long)sum_mu_ok, (unsigned long long)sum_mu_total,
           sum_mu_ok == sum_mu_total ? "ALL PASS" : "FAILURES",
           (unsigned long long)sum_nontriv);
    printf("[2] Joint hits per run           small AES      random\n");
    printf("    mean                         %-12.2f   %-12.2f\n", mean_aes, mean_rand);
    printf("    variance (population)        %-12.2f   %-12.2f\n", var_aes, var_rand);
    printf("    expected (Poisson lambda)    %-12.3g   %-12.3e\n", lambda_aes, lambda_rand);
    printf("    control, C0^C2 and C0^C3 only (mean): small AES %.2f, random %.2f\n",
           s_half_aes / n, s_half_rand / n);
    printf("\n Results saved to: %s\n", csv_path);
    printf(" Total time: %.0f seconds\n", difftime(time(NULL), total_t0));
    printf(" Table 3: python3 (or python) scripts/summarize.py %s\n", csv_path);
    printf("=============================================================\n");
    return sum_mu_ok == sum_mu_total ? 0 : 3;
}
