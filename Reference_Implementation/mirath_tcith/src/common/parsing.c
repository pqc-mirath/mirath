#include <string.h>
#include "prng.h"
#include "mirath_parameters.h"
#include "hash.h"
#include "mirath_matrix_ff.h"

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