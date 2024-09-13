#include <string.h>
#include "prng.h"
#include "mirath_parameters.h"
#include "hash.h"
#include "mirath_matrix_ff.h"
#include "mirath_tcith.h"

void unparse_public_key(uint8_t *pk, const seed_t seed_pk, const ff_t y[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1)]) {
    memcpy(pk, seed_pk, MIRATH_SECURITY_BYTES);

    memcpy(pk + MIRATH_SECURITY_BYTES, y, mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1));
}

void parse_public_key(seed_t seed_pk, ff_t y[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1)], const uint8_t *pk) {
    memcpy(seed_pk, pk, MIRATH_SECURITY_BYTES);

    memcpy(y, pk + MIRATH_SECURITY_BYTES, mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1));
}

void unparse_secret_key(uint8_t *sk, const seed_t seed_sk, const seed_t seed_pk) {
    memcpy(sk, seed_sk, MIRATH_SECURITY_BYTES);

    memcpy(sk + MIRATH_SECURITY_BYTES, seed_pk, MIRATH_SECURITY_BYTES);
}

void parse_secret_key(ff_t S[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_R)],
                      ff_t C[mirath_matrix_ff_bytes_size(MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R)],
                      ff_t H[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, MIRATH_PARAM_K)],
                      const uint8_t *sk) {
    seed_t seed;
    mirath_prng_t prng;

    // secret seed
    memcpy(seed, sk, MIRATH_SECURITY_BYTES);

    mirath_prng_init(&prng, NULL, seed, MIRATH_SECURITY_BYTES);
    // TODO: use pointers in the optimized version
    ff_t T[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_R) + mirath_matrix_ff_bytes_size(MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R)];
    mirath_prng(&prng, T, mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_R) + mirath_matrix_ff_bytes_size(MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R));
    memcpy(S, T, mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_R));
    memcpy(C, T + mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_R), mirath_matrix_ff_bytes_size(MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R));

    mirath_matrix_set_to_ff(S, MIRATH_PARAM_M, MIRATH_PARAM_R);
    mirath_matrix_set_to_ff(C, MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R);

    // public seed
    memcpy(seed, sk + MIRATH_SECURITY_BYTES, MIRATH_SECURITY_BYTES);
    mirath_prng_init(&prng, NULL, seed, MIRATH_SECURITY_BYTES);
    mirath_prng(&prng, H, mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, MIRATH_PARAM_K));

    mirath_matrix_set_to_ff(H, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, MIRATH_PARAM_K);
}

void unparse_signature(uint8_t *signature, const uint8_t salt[MIRATH_PARAM_SALT_BYTES], const uint64_t ctr,
                       const hash_t hash2, const uint8_t *path, const uint64_t path_length,
                       mirath_tcith_commit_t *commits[MIRATH_PARAM_TAU],
                       const ff_t aux[MIRATH_PARAM_TAU][mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_R + MIRATH_PARAM_R * (MIRATH_PARAM_N - MIRATH_PARAM_R), 1)],
                       const ff_mu_t mid_alpha[MIRATH_PARAM_TAU][MIRATH_PARAM_RHO],
                       const mirath_tcith_challenge_t i_star) {

    uint8_t *ptr;

    ptr = signature;

    memcpy(ptr, salt, MIRATH_PARAM_SALT_BYTES);
    ptr += MIRATH_PARAM_SALT_BYTES;

    memcpy(ptr, &ctr, sizeof(uint64_t));
    ptr += sizeof(ctr);

    memcpy(ptr, hash2, 2 * MIRATH_SECURITY_BYTES);
    ptr += 2 * MIRATH_SECURITY_BYTES;

    memcpy(ptr, path, path_length * MIRATH_SECURITY_BYTES);
    ptr += path_length * MIRATH_SECURITY_BYTES;

    for (uint32_t e = 0; e < MIRATH_PARAM_TAU; e++) {
        memcpy(ptr, commits[e][i_star[e]], 2 * MIRATH_SECURITY_BYTES);
        ptr += 2 * MIRATH_SECURITY_BYTES;
    }

    const uint32_t n_bytes = mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_R + MIRATH_PARAM_R * (MIRATH_PARAM_N - MIRATH_PARAM_R), 1);
    for (uint32_t e = 0; e < MIRATH_PARAM_TAU; e++) {
        memcpy(ptr, aux[e], n_bytes);
        ptr += n_bytes;

        memcpy(ptr, mid_alpha[e], sizeof(ff_mu_t) * MIRATH_PARAM_RHO);
        ptr += sizeof(ff_mu_t) * MIRATH_PARAM_RHO;
    }
}

void parse_signature(uint8_t salt[MIRATH_PARAM_SALT_BYTES], uint64_t *ctr, hash_t hash2,
                     uint8_t *path, mirath_tcith_commit_t commits_i_star[MIRATH_PARAM_TAU],
                     ff_t aux[MIRATH_PARAM_TAU][mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_R + MIRATH_PARAM_R * (MIRATH_PARAM_N - MIRATH_PARAM_R), 1)],
                     ff_mu_t mid_alpha[MIRATH_PARAM_TAU][MIRATH_PARAM_RHO], const uint8_t *signature) {

    uint8_t *ptr = (uint8_t *)signature;

    memcpy(salt, ptr, MIRATH_PARAM_SALT_BYTES);
    ptr += MIRATH_PARAM_SALT_BYTES;

    memcpy(ctr, ptr, sizeof(uint64_t));
    ptr += sizeof(uint64_t);

    memcpy(hash2, ptr, 2 * MIRATH_SECURITY_BYTES);
    ptr += 2 * MIRATH_SECURITY_BYTES;

    memcpy(ptr, path, MIRATH_SECURITY_BYTES * MIRATH_PARAM_T_OPEN);
    ptr += MIRATH_SECURITY_BYTES * MIRATH_PARAM_T_OPEN;

    for (uint32_t e = 0; e < MIRATH_PARAM_TAU; e++) {
        memcpy(commits_i_star[e], ptr, 2 * MIRATH_SECURITY_BYTES);
        ptr += 2 * MIRATH_SECURITY_BYTES;
    }

    const uint32_t n_bytes = mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_R + MIRATH_PARAM_R * (MIRATH_PARAM_N - MIRATH_PARAM_R), 1);
    for (uint32_t e = 0; e < MIRATH_PARAM_TAU; e++) {
        memcpy(aux[e], ptr, n_bytes);
        ptr += n_bytes;

        memcpy(mid_alpha[e], ptr, MIRATH_PARAM_RHO);
        ptr += MIRATH_PARAM_RHO;
    }
}