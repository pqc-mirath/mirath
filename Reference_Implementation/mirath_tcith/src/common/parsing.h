#ifndef MIRATH_TCITH_PARSING_H
#define MIRATH_TCITH_PARSING_H

#include "mirath_parameters.h"
#include "mirath_matrix_ff.h"
#include "mirath_tcith.h"

void unparse_public_key(uint8_t *pk, const seed_t seed_pk, const ff_t y[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1)]);

void parse_public_key(seed_t seed_pk, ff_t y[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1)], const uint8_t *pk);

void unparse_secret_key(uint8_t *sk, const seed_t seed_sk, const seed_t seed_pk);

void parse_secret_key(ff_t S[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_R)],
                      ff_t C[mirath_matrix_ff_bytes_size(MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R)],
                      ff_t H[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, MIRATH_PARAM_K)],
                      const uint8_t *sk);

void unparse_signature(uint8_t *signature, const uint8_t salt[MIRATH_PARAM_SALT_BYTES], uint64_t ctr,
                       const hash_t hash2, const uint8_t *path, uint64_t path_length,
                       mirath_tcith_commit_t *commits[MIRATH_PARAM_TAU],
                       const ff_t aux[MIRATH_PARAM_TAU][MIRATH_VAR_FF_AUX_BYTES],
                       const ff_mu_t mid_alpha[MIRATH_PARAM_TAU][MIRATH_PARAM_RHO],
                       const mirath_tcith_challenge_t i_star);

void parse_signature(uint8_t salt[MIRATH_PARAM_SALT_BYTES], uint64_t *ctr, hash_t hash2,
                     uint8_t *path, mirath_tcith_commit_t commits_i_star[MIRATH_PARAM_TAU],
                     ff_t aux[MIRATH_PARAM_TAU][MIRATH_VAR_FF_AUX_BYTES],
                     ff_mu_t mid_alpha[MIRATH_PARAM_TAU][MIRATH_PARAM_RHO], const uint8_t *signature);

#endif //MIRATH_TCITH_PARSING_H
