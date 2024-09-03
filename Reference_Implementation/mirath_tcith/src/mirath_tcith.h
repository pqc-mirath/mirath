#ifndef MIRATH_TCITH_H
#define MIRATH_TCITH_H

#include <stdint.h>
#include <stdlib.h>
#include "mirath_parameters.h"
#include "mirath_tree.h"

typedef uint8_t mirath_tcith_commit_t[2 * MIRATH_SECURITY_BYTES];
typedef mirath_tcith_commit_t mirath_tcith_commit_1_t[MIRATH_PARAM_TAU_1][MIRATH_PARAM_N_1];
typedef mirath_tcith_commit_t mirath_tcith_commit_2_t[MIRATH_PARAM_TAU_2][MIRATH_PARAM_N_2];

void mirath_tcith_commit_set_as_a_grid_list(mirath_tcith_commit_t *seeds[MIRATH_PARAM_TAU],
                                            mirath_tcith_commit_1_t *input_1,
                                            mirath_tcith_commit_2_t *input_2);

void mirath_tcith_commit(mirath_tcith_commit_t commit, const uint8_t *salt, uint16_t e, uint16_t i, const uint8_t *seed);
size_t mirath_tcith_psi(size_t i, size_t e);

void build_sharing_N(ff_t aux[MIRATH_PARAM_TAU][mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_R + MIRATH_PARAM_R * (MIRATH_PARAM_N - MIRATH_PARAM_R), 1)],
                     ff_mu_t rnd_S[MIRATH_PARAM_TAU][MIRATH_PARAM_M * MIRATH_PARAM_R],
                     ff_mu_t rnd_C[MIRATH_PARAM_TAU][MIRATH_PARAM_R * (MIRATH_PARAM_N - MIRATH_PARAM_R)],
                     ff_mu_t rnd_v[MIRATH_PARAM_TAU][MIRATH_PARAM_RHO], hash_ctx_t hash_ctx,
                     const mirath_tree_leaves_t seeds,
                     const ff_t S[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_R)],
                     const ff_t C[mirath_matrix_ff_bytes_size(MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R)],
                     const uint8_t salt[MIRATH_PARAM_SALT_BYTES]);

void split_codeword_ff_mu(ff_mu_t e_A[MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K], ff_mu_t e_B[MIRATH_PARAM_K],
                          const ff_mu_t in_X[MIRATH_PARAM_M * MIRATH_PARAM_R],
                          const ff_mu_t in_Y[MIRATH_PARAM_M * (MIRATH_PARAM_N - MIRATH_PARAM_R)]);

void emulateMPC_mu(ff_mu_t base_alpha[MIRATH_PARAM_RHO],
                   ff_mu_t mid_alpha[MIRATH_PARAM_RHO],
                   const ff_t S[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_R)],
                   const ff_mu_t rnd_S[MIRATH_PARAM_M * MIRATH_PARAM_R],
                   const ff_t C[mirath_matrix_ff_bytes_size(MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R)],
                   const ff_mu_t rnd_C[MIRATH_PARAM_R * (MIRATH_PARAM_N - MIRATH_PARAM_R)],
                   const ff_mu_t v[MIRATH_PARAM_RHO], ff_mu_t rnd_v[MIRATH_PARAM_RHO],
                   const ff_mu_t gamma[MIRATH_PARAM_RHO * (MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K)],
                   const ff_t H[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, MIRATH_PARAM_K)]);

void emulateparty_mu(ff_mu_t base_alpha[MIRATH_PARAM_RHO], const ff_mu_t p,
                     const ff_mu_t share_S[MIRATH_PARAM_M * MIRATH_PARAM_R],
                     const ff_mu_t share_C[MIRATH_PARAM_R * (MIRATH_PARAM_N - MIRATH_PARAM_R)],
                     const ff_mu_t share_v[MIRATH_PARAM_RHO],
                     const ff_mu_t gamma[MIRATH_PARAM_RHO * (MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K)],
                     const ff_t H[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, MIRATH_PARAM_K)],
                     const ff_t y[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1)],
                     const ff_mu_t mid_alpha[MIRATH_PARAM_RHO]);

#endif //MIRATH_TCITH_H
