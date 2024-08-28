#ifndef MIRATH_TCITH_PARSING_H
#define MIRATH_TCITH_PARSING_H

#include "mirath_parameters.h"
#include "mirath_matrix_ff.h"

void unparse_public_key(uint8_t *pk, const seed_t seed_pk, const ff_t y[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1)]);

void parse_public_key(seed_t seed_pk, ff_t y[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1)], const uint8_t *pk);

void unparse_secret_key(uint8_t *sk, const seed_t seed_sk, const seed_t seed_pk);

void parse_secret_key(ff_t S[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_R)],
                      ff_t C[mirath_matrix_ff_bytes_size(MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R)],
                      ff_t H[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, MIRATH_PARAM_K)],
                      const uint8_t *sk);

#endif //MIRATH_TCITH_PARSING_H
