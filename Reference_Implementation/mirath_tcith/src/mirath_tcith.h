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

void mirath_tcith_internal_steps_pk(ff_t y[MIRATH_VAR_FF_Y_BYTES],
                                    const ff_t S[MIRATH_VAR_FF_S_BYTES], const ff_t C[MIRATH_VAR_FF_C_BYTES],
                                    const ff_t H[MIRATH_VAR_FF_H_BYTES]);

void mirath_tciht_compute_public_key(uint8_t *pk, const uint8_t *sk,
                                     const ff_t S[MIRATH_VAR_FF_S_BYTES], const ff_t C[MIRATH_VAR_FF_C_BYTES],
                                     const ff_t H[MIRATH_VAR_FF_H_BYTES]);

void mirath_tcith_commit_set_as_a_grid_list(mirath_tcith_commit_t *seeds[MIRATH_PARAM_TAU],
                                            mirath_tcith_commit_1_t *input_1,
                                            mirath_tcith_commit_2_t *input_2);

void mirath_tcith_commit(mirath_tcith_commit_t commit, const uint8_t *salt, uint16_t e, uint16_t i, const uint8_t *seed);
size_t mirath_tcith_psi(size_t i, size_t e);

void build_sharing_N(ff_t aux[MIRATH_VAR_FF_AUX_BYTES],
                     ff_mu_t rnd_S[MIRATH_VAR_S],
                     ff_mu_t rnd_C[MIRATH_VAR_C],
                     ff_mu_t rnd_v[MIRATH_PARAM_RHO], ff_mu_t v[MIRATH_PARAM_RHO],
                     mirath_tcith_commit_t *commits[MIRATH_PARAM_TAU],
                     const mirath_tree_leaves_t seeds,
                     const ff_t S[MIRATH_VAR_FF_S_BYTES],
                     const ff_t C[MIRATH_VAR_FF_C_BYTES],
                     const uint8_t salt[MIRATH_PARAM_SALT_BYTES], uint32_t e);

void compute_share(ff_mu_t share_S[MIRATH_VAR_S],
                   ff_mu_t share_C[MIRATH_VAR_C],
                   ff_mu_t share_v[MIRATH_PARAM_RHO], mirath_tcith_commit_t *commits[MIRATH_PARAM_TAU],
                   uint32_t i_star, const mirath_tree_leaves_t seeds, uint32_t e,
                   const uint8_t salt[MIRATH_PARAM_SALT_BYTES], mirath_tcith_commit_t commits_i_star[MIRATH_PARAM_TAU],
                   const ff_t aux[MIRATH_VAR_FF_AUX_BYTES]);

void split_codeword_ff_mu(ff_mu_t e_A[MIRATH_VAR_E_A], ff_mu_t e_B[MIRATH_PARAM_K],
                          const ff_mu_t in_X[MIRATH_VAR_S],
                          const ff_mu_t in_Y[MIRATH_VAR_BASE_MID]);

void emulateMPC_mu(ff_mu_t base_alpha[MIRATH_PARAM_RHO], ff_mu_t mid_alpha[MIRATH_PARAM_RHO],
                   const ff_t S[MIRATH_VAR_FF_S_BYTES], const ff_mu_t rnd_S[MIRATH_VAR_S],
                   const ff_t C[MIRATH_VAR_FF_C_BYTES], const ff_mu_t rnd_C[MIRATH_VAR_C],
                   const ff_mu_t v[MIRATH_PARAM_RHO], ff_mu_t rnd_v[MIRATH_PARAM_RHO],
                   const ff_mu_t gamma[MIRATH_VAR_GAMMA], const ff_t H[MIRATH_VAR_FF_H_BYTES]);

void emulateparty_mu(ff_mu_t base_alpha[MIRATH_PARAM_RHO], ff_mu_t p,
                     const ff_mu_t share_S[MIRATH_VAR_S], const ff_mu_t share_C[MIRATH_VAR_C],
                     const ff_mu_t share_v[MIRATH_PARAM_RHO], const ff_mu_t gamma[MIRATH_VAR_GAMMA],
                     const ff_t H[MIRATH_VAR_FF_H_BYTES], const ff_t y[MIRATH_VAR_FF_Y_BYTES],
                     const ff_mu_t mid_alpha[MIRATH_PARAM_RHO]);


void mirath_tcith_compute_challenge_2(mirath_tcith_challenge_t challenge, const uint8_t *seed_input, const uint8_t *salt);
uint8_t mirath_tcith_discard_input_challenge_2(const uint8_t *seed_input);

#define mirath_tcith_shift_to_right(shiftOut, highIn, lowIn, shift, DigitSize)  \
    (shiftOut) = ((lowIn) >> (shift)) ^ ((highIn) << ((DigitSize) - (shift)));

#endif //MIRATH_TCITH_H
