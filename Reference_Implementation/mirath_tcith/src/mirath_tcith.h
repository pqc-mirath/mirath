#ifndef MIRATH_TCITH_H
#define MIRATH_TCITH_H

#include <stdint.h>
#include <stdlib.h>
#include "mirath_parameters.h"
#include "mirath_tree.h"

typedef uint8_t mirath_tcith_commit_t[2 * MIRATH_SECURITY_BYTES];
typedef mirath_tcith_commit_t mirath_tcith_commit_1_t[MIRATH_PARAM_TAU_1][MIRATH_PARAM_N_1];
typedef mirath_tcith_commit_t mirath_tcith_commit_2_t[MIRATH_PARAM_TAU_2][MIRATH_PARAM_N_2];

typedef uint16_t mirath_tcith_challenge_t[MIRATH_PARAM_TAU];

void mirath_tcith_internal_steps_pk(ff_t y[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1)],
                                    const ff_t S[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_R)],
                                    const ff_t C[mirath_matrix_ff_bytes_size(MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R)],
                                    const ff_t H[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, MIRATH_PARAM_K)]);

void mirath_tciht_compute_public_key(uint8_t *pk, const uint8_t *sk,
                                     const ff_t S[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_R)],
                                     const ff_t C[mirath_matrix_ff_bytes_size(MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R)],
                                     const ff_t H[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, MIRATH_PARAM_K)]);

void mirath_tcith_commit_set_as_a_grid_list(mirath_tcith_commit_t *seeds[MIRATH_PARAM_TAU],
                                            mirath_tcith_commit_1_t *input_1,
                                            mirath_tcith_commit_2_t *input_2);

void mirath_tcith_commit(mirath_tcith_commit_t commit, const uint8_t *salt, uint16_t e, uint16_t i, const uint8_t *seed);
size_t mirath_tcith_psi(size_t i, size_t e);

void build_sharing_N(ff_t aux[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_R + MIRATH_PARAM_R * (MIRATH_PARAM_N - MIRATH_PARAM_R), 1)],
                     ff_mu_t rnd_S[MIRATH_PARAM_M * MIRATH_PARAM_R],
                     ff_mu_t rnd_C[MIRATH_PARAM_R * (MIRATH_PARAM_N - MIRATH_PARAM_R)],
                     ff_mu_t rnd_v[MIRATH_PARAM_RHO], mirath_tcith_commit_t *commits[MIRATH_PARAM_TAU],
                     const mirath_tree_leaves_t seeds,
                     const ff_t S[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_R)],
                     const ff_t C[mirath_matrix_ff_bytes_size(MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R)],
                     const uint8_t salt[MIRATH_PARAM_SALT_BYTES], const uint32_t e);

void compute_share(ff_mu_t share_S[MIRATH_PARAM_M * MIRATH_PARAM_R],
                   ff_mu_t share_C[MIRATH_PARAM_R * (MIRATH_PARAM_N - MIRATH_PARAM_R)],
                   ff_mu_t share_v[MIRATH_PARAM_RHO], mirath_tcith_commit_t *commits[MIRATH_PARAM_TAU],
                   const uint32_t i_star, const mirath_tree_leaves_t seeds, const uint32_t e,
                   const uint8_t salt[MIRATH_PARAM_SALT_BYTES], mirath_tcith_commit_t commits_i_star[MIRATH_PARAM_TAU],
                   const ff_t aux[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_R + MIRATH_PARAM_R * (MIRATH_PARAM_N - MIRATH_PARAM_R), 1)]);

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


void mirath_tcith_compute_challenge_2(mirath_tcith_challenge_t challenge, const uint8_t *seed_input, const uint8_t *salt);
uint8_t mirath_tcith_discard_input_challenge_2(const uint8_t *seed_input);

#define mirath_tcith_shift_to_right(shiftOut, highIn, lowIn, shift, DigitSize)  \
    (shiftOut) = ((lowIn) >> (shift)) ^ ((highIn) << ((DigitSize) - (shift)));
void mirath_tcith_shift_to_right_array(uint8_t *string, size_t length);

#endif //MIRATH_TCITH_H
