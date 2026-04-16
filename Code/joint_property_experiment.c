/*
 * joint_property_experiment.c
 * 
 * Verifies the Joint Generalized Zero-Difference Property (Theorem 4)
 * on 4-round small-scale AES (4-bit S-box, 64-bit block).
 *
 * Runs NUM_KEYS independent experiments (each with a fresh random key
 * and fresh random related differences), and saves per-run results
 * to a CSV file. Also prints a summary with mean/variance.
 *
 * Compile: gcc -O2 -o joint_exp joint_property_experiment.c -lm
 * Usage:   ./joint_exp [output.csv]
 *          (default output: results.csv)
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <math.h>

// ===========================================================================
// Configuration
// ===========================================================================
#define NUM_ROUNDS       4
#define TRIALS           (1ULL << 30)
#define MU_CHECK_TRIALS  1000000ULL
#define NUM_KEYS         100

// ===========================================================================
// Small AES Primitives (4-bit S-box, GF(2^4) with x^4 + x + 1)
// ===========================================================================
static const uint8_t SBOX[16] = {
    0x6, 0xB, 0x5, 0x4, 0x2, 0xE, 0x7, 0xA,
    0x9, 0xD, 0xF, 0xC, 0x3, 0x1, 0x0, 0x8
};

static uint8_t gf_mul(uint8_t a, uint8_t b) {
    uint8_t p = 0;
    for (int i = 0; i < 4; i++) {
        if (b & 1) p ^= a;
        int hi = a & 8;
        a = (a << 1) & 0xF;
        if (hi) a ^= 0x3;
        b >>= 1;
    }
    return p;
}

static uint16_t sb_mc[4][16];
static const uint8_t MDS[4][4] = {
    {2, 3, 1, 1}, {1, 2, 3, 1}, {1, 1, 2, 3}, {3, 1, 1, 2}
};

static void init_tables(void) {
    for (int row = 0; row < 4; row++)
        for (int v = 0; v < 16; v++) {
            uint8_t s = SBOX[v];
            uint16_t p = 0;
            for (int r = 0; r < 4; r++)
                p |= (uint16_t)gf_mul(MDS[r][row], s) << (4 * (3 - r));
            sb_mc[row][v] = p;
        }
}

// ===========================================================================
// State Layout
// ===========================================================================
static const int DIAG[4][4] = {
    {0, 5, 10, 15}, {1, 6, 11, 12}, {2, 7, 8, 13}, {3, 4, 9, 14}
};

static const int INV_DIAG[4][4] = {
    {0, 7, 10, 13}, {1, 4, 11, 14}, {2, 5, 8, 15}, {3, 6, 9, 12}
};

// ===========================================================================
// Encryption
// ===========================================================================
static void round_mc(uint8_t *state, const uint8_t *rk) {
    uint8_t tmp[16];
    for (int c = 0; c < 4; c++) {
        uint16_t col = sb_mc[0][state[DIAG[c][0]]]
                     ^ sb_mc[1][state[DIAG[c][1]]]
                     ^ sb_mc[2][state[DIAG[c][2]]]
                     ^ sb_mc[3][state[DIAG[c][3]]];
        tmp[c]    = ((col >> 12) & 0xF) ^ rk[c];
        tmp[4+c]  = ((col >> 8)  & 0xF) ^ rk[4+c];
        tmp[8+c]  = ((col >> 4)  & 0xF) ^ rk[8+c];
        tmp[12+c] = (col         & 0xF) ^ rk[12+c];
    }
    memcpy(state, tmp, 16);
}

static void round_last(uint8_t *state, const uint8_t *rk) {
    uint8_t tmp[16];
    for (int i = 0; i < 16; i++) tmp[i] = SBOX[state[i]];
    uint8_t tmp2[16];
    for (int r = 0; r < 4; r++)
        for (int c = 0; c < 4; c++)
            tmp2[r * 4 + c] = tmp[r * 4 + (c + r) % 4];
    for (int i = 0; i < 16; i++) state[i] = tmp2[i] ^ rk[i];
}

static void encrypt(uint8_t *ct, const uint8_t *pt, const uint8_t rk[][16]) {
    memcpy(ct, pt, 16);
    for (int i = 0; i < 16; i++) ct[i] ^= rk[0][i];
    for (int r = 1; r < NUM_ROUNDS; r++)
        round_mc(ct, rk[r]);
    round_last(ct, rk[NUM_ROUNDS]);
}

// ===========================================================================
// Helpers
// ===========================================================================
static inline int inv_diag_zero(const uint8_t *s, int d) {
    return (s[INV_DIAG[d][0]] | s[INV_DIAG[d][1]]
          | s[INV_DIAG[d][2]] | s[INV_DIAG[d][3]]) == 0;
}

// ===========================================================================
// PRNG: xoroshiro128++
// ===========================================================================
static uint64_t rng[2];

static inline uint64_t rotl64(uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

static inline uint64_t next_rng(void) {
    uint64_t s0 = rng[0], s1 = rng[1];
    uint64_t res = rotl64(s0 + s1, 17) + s0;
    s1 ^= s0;
    rng[0] = rotl64(s0, 49) ^ s1 ^ (s1 << 21);
    rng[1] = rotl64(s1, 28);
    return res;
}

static void rand_nibbles(uint8_t *s) {
    uint64_t r = next_rng();
    for (int i = 0; i < 16; i++)
        s[i] = (r >> (i * 4)) & 0xF;
}

// ===========================================================================
// Single-key experiment
// ===========================================================================
typedef struct {
    int      run_id;
    uint64_t mu_ok;
    uint64_t mu_total;
    uint64_t single_xp;
    uint64_t single_xxp;
    uint64_t hit_aes;
    uint64_t hit_rand;
    long     elapsed;
} RunResult;

static void run_experiment(int run_id, RunResult *res) {
    // Generate random round keys
    uint8_t rk[NUM_ROUNDS + 1][16];
    for (int r = 0; r <= NUM_ROUNDS; r++) rand_nibbles(rk[r]);

    // Generate related differences
    uint8_t dx[16], dxp[16], dxxp[16];
    for (int i = 0; i < 16; i++) {
        do { dx[i] = next_rng() & 0xF; } while (!dx[i]);
    }
    memset(dxp, 0, 16);
    for (int k = 0; k < 4; k++) {
        dxp[DIAG[1][k]] = dx[DIAG[1][k]];
        dxp[DIAG[3][k]] = dx[DIAG[3][k]];
    }
    for (int i = 0; i < 16; i++) dxxp[i] = dx[i] ^ dxp[i];

    printf("  Run %2d: dx=", run_id + 1);
    for (int i = 0; i < 16; i++) printf("%X", dx[i]);
    printf("\n");

    res->run_id = run_id + 1;
    res->mu_ok = 0; res->mu_total = 0;
    res->single_xp = 0; res->single_xxp = 0;
    res->hit_aes = 0; res->hit_rand = 0;

    time_t t0 = time(NULL);

    for (uint64_t t = 0; t < TRIALS; t++) {
        uint8_t a[16];
        rand_nibbles(a);

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
            d1[i] = C0[i] ^ C2[i];
            d2[i] = C1[i] ^ C3[i];
            d3[i] = C0[i] ^ C3[i];
            d4[i] = C1[i] ^ C2[i];
        }

        // mu-relation verification (first MU_CHECK_TRIALS)
        if (t < MU_CHECK_TRIALS) {
            int ok = 1;
            uint8_t dx_d1[16], dx_d2[16];
            for (int i = 0; i < 16; i++) {
                dx_d1[i] = C0[i] ^ C1[i];
                dx_d2[i] = C2[i] ^ C3[i];
            }
            for (int c = 0; c < 4 && ok; c++) {
                if (inv_diag_zero(d1, c) != inv_diag_zero(d2, c)) ok = 0;
                if (inv_diag_zero(d3, c) != inv_diag_zero(d4, c)) ok = 0;
                if (inv_diag_zero(dx_d1, c) != inv_diag_zero(dx_d2, c)) ok = 0;
            }
            res->mu_total++;
            if (ok) res->mu_ok++;
        }

        // Joint hit counting (AES)
        uint8_t zm1 = 0, zm3 = 0;
        for (int c = 0; c < 4; c++) {
            if (inv_diag_zero(d1, c)) zm1 |= (1 << c);
            if (inv_diag_zero(d3, c)) zm3 |= (1 << c);
        }
        if (zm1) res->single_xp++;
        if (zm3) res->single_xxp++;
        if (zm1 && zm3 && !(zm1 & zm3))
            res->hit_aes++;

        // Random baseline
        uint8_t R0[16], R1[16], R2[16], R3[16];
        rand_nibbles(R0); rand_nibbles(R1);
        rand_nibbles(R2); rand_nibbles(R3);

        uint8_t rd1[16], rd2[16], rd3[16], rd4[16];
        for (int i = 0; i < 16; i++) {
            rd1[i] = R0[i] ^ R2[i];
            rd2[i] = R1[i] ^ R3[i];
            rd3[i] = R0[i] ^ R3[i];
            rd4[i] = R1[i] ^ R2[i];
        }
        uint8_t rzm1 = 0, rzm3 = 0;
        for (int c = 0; c < 4; c++) {
            if (inv_diag_zero(rd1, c) && inv_diag_zero(rd2, c))
                rzm1 |= (1 << c);
            if (inv_diag_zero(rd3, c) && inv_diag_zero(rd4, c))
                rzm3 |= (1 << c);
        }
        if (rzm1 && rzm3 && !(rzm1 & rzm3))
            res->hit_rand++;

        // Progress
        if ((t & 0xFFFFFFFULL) == 0 && t > 0) {
            time_t now = time(NULL);
            double pct = 100.0 * t / TRIALS;
            printf("\r  Run %2d: [%5.1f%%] AES=%llu Rand=%llu (%lds)",
                   run_id + 1, pct,
                   (unsigned long long)res->hit_aes,
                   (unsigned long long)res->hit_rand,
                   now - t0);
            fflush(stdout);
        }
    }

    res->elapsed = time(NULL) - t0;
    printf("\r  Run %2d: [100.0%%] AES=%llu Rand=%llu (%lds)        \n",
           run_id + 1,
           (unsigned long long)res->hit_aes,
           (unsigned long long)res->hit_rand,
           res->elapsed);
}

// ===========================================================================
// Main
// ===========================================================================
int main(int argc, char *argv[]) {
    const char *csv_path = (argc > 1) ? argv[1] : "results.csv";

    init_tables();

    // Seed PRNG
    rng[0] = (uint64_t)time(NULL) * 6364136223846793005ULL + 1;
    rng[1] = rng[0] ^ 0xdeadbeefcafebabeULL;
    for (int i = 0; i < 20; i++) next_rng();

    printf("=============================================================\n");
    printf(" Joint Zero-Difference Property — Multi-Key Experiment\n");
    printf("=============================================================\n");
    printf(" Block: 64-bit (4-bit S-box), Rounds: %d (last w/o MC)\n", NUM_ROUNDS);
    printf(" Trials per key: 2^%.0f = %llu\n", log2((double)TRIALS), (unsigned long long)TRIALS);
    printf(" Number of keys: %d\n", NUM_KEYS);
    printf(" Output: %s\n", csv_path);
    printf("-------------------------------------------------------------\n\n");

    RunResult results[NUM_KEYS];
    time_t total_t0 = time(NULL);

    for (int k = 0; k < NUM_KEYS; k++) {
        run_experiment(k, &results[k]);
    }

    time_t total_t1 = time(NULL);

    // ===========================================================================
    // Write CSV
    // ===========================================================================
    FILE *fp = fopen(csv_path, "w");
    if (!fp) { perror("fopen"); return 1; }
    fprintf(fp, "run,mu_verified,mu_total,single_xp,single_xxp,joint_aes,joint_rand,elapsed_sec\n");
    for (int k = 0; k < NUM_KEYS; k++) {
        fprintf(fp, "%d,%llu,%llu,%llu,%llu,%llu,%llu,%ld\n",
                results[k].run_id,
                (unsigned long long)results[k].mu_ok,
                (unsigned long long)results[k].mu_total,
                (unsigned long long)results[k].single_xp,
                (unsigned long long)results[k].single_xxp,
                (unsigned long long)results[k].hit_aes,
                (unsigned long long)results[k].hit_rand,
                results[k].elapsed);
    }
    fclose(fp);

    // ===========================================================================
    // Summary
    // ===========================================================================
    printf("\n=============================================================\n");
    printf(" SUMMARY (%d runs x 2^%.0f trials)\n", NUM_KEYS, log2((double)TRIALS));
    printf("=============================================================\n\n");

    // mu-relation
    int mu_all_pass = 1;
    for (int k = 0; k < NUM_KEYS; k++)
        if (results[k].mu_ok != results[k].mu_total) mu_all_pass = 0;
    printf("[1] mu-relation: %s (%llu verified per run)\n",
           mu_all_pass ? "ALL PASS (100%%)" : "SOME FAILURES",
           (unsigned long long)MU_CHECK_TRIALS);

    // Joint hits statistics
    double sum_aes = 0, sum_rand = 0, sumsq_aes = 0, sumsq_rand = 0;
    printf("\n[2] Joint hits per run:\n");
    printf("    %-6s  %-12s  %-12s\n", "Run", "AES", "Random");
    printf("    %-6s  %-12s  %-12s\n", "---", "---", "------");
    for (int k = 0; k < NUM_KEYS; k++) {
        printf("    %-6d  %-12llu  %-12llu\n",
               results[k].run_id,
               (unsigned long long)results[k].hit_aes,
               (unsigned long long)results[k].hit_rand);
        sum_aes += results[k].hit_aes;
        sum_rand += results[k].hit_rand;
        sumsq_aes += results[k].hit_aes * results[k].hit_aes;
        sumsq_rand += results[k].hit_rand * results[k].hit_rand;
    }

    double mean_aes = sum_aes / NUM_KEYS;
    double mean_rand = sum_rand / NUM_KEYS;
    double var_aes = sumsq_aes / NUM_KEYS - mean_aes * mean_aes;
    double var_rand = sumsq_rand / NUM_KEYS - mean_rand * mean_rand;

    double expected_aes = (double)TRIALS * 12.0 * pow(2, -32);
    double expected_rand = (double)TRIALS * 12.0 * pow(2, -64);

    printf("\n    %-6s  %-12.2f  %-12.6f\n", "Mean", mean_aes, mean_rand);
    printf("    %-6s  %-12.2f  %-12.6f\n", "Var", var_aes, var_rand);
    printf("    %-6s  %-12.2f  %-12.6f\n", "E[X]", expected_aes, expected_rand);
    printf("    (Poisson: E[X]=Var=lambda)\n");

    printf("\n    Theoretical advantage: 2^%.1f\n",
           log2(12.0 * pow(2, -32)) - log2(12.0 * pow(2, -64)));

    printf("\n Results saved to: %s\n", csv_path);
    printf(" Total time: %ld seconds (%.1f minutes)\n",
           total_t1 - total_t0, (total_t1 - total_t0) / 60.0);
    printf("=============================================================\n");

    return 0;
}
