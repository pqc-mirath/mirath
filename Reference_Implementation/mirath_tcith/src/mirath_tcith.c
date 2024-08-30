/**
 * @file mirath_tcith.c
 * @brief Implementation of TCitH related functions
 */

#include "prng.h"
#include "random.h"
#include "matrix_ff_mu.h"
#include "vector_ff_mu.h"
#include "mirath_matrix_ff.h"
#include "mirath_tcith.h"

void mirath_tcith_commit_set_as_a_grid_list(mirath_tcith_commit_t *seeds[MIRATH_PARAM_TAU],
                                            mirath_tcith_commit_1_t *input_1,
                                            mirath_tcith_commit_2_t *input_2) {
    // The below lines represent the leaves as a tau-dimensional grid
    for (size_t i = 0; i < MIRATH_PARAM_TAU_1; i++) {
        seeds[i] = (*input_1)[i];
    }
    for (size_t i = 0; i < MIRATH_PARAM_TAU_2; i++) {
        seeds[i + MIRATH_PARAM_TAU_1] = (*input_2)[i];
    }
}

void mirath_tcith_commit(mirath_tcith_commit_t commit, const uint8_t *salt, uint8_t e, size_t i, const uint8_t *seed) {
    uint8_t domain_separator = DOMAIN_SEPARATOR_COMMITMENT;

    hash_t ctx;
    hash_init(&ctx);
    hash_update(&ctx, &domain_separator, sizeof(uint8_t));
    hash_update(&ctx, salt, 2 * MIRATH_SECURITY_BYTES);
    hash_update(&ctx, &e, sizeof(uint8_t));
    hash_update(&ctx, (uint8_t *) &i, sizeof(uint16_t));
    hash_update(&ctx, seed, MIRATH_SECURITY_BYTES);
    hash_finalize(&ctx, commit);
}

size_t mirath_tcith_psi(size_t i, size_t e) {
    if (i < MIRATH_PARAM_N_2) {
        return i * MIRATH_PARAM_TAU + e;
    }
    else {
        return MIRATH_PARAM_N_2 * MIRATH_PARAM_TAU + (i - MIRATH_PARAM_N_2) * MIRATH_PARAM_TAU_1 + e;
    }
}

// this is just for the reference implementation
// use pointers in the optimized implementation
void split_codeword_ff_mu(ff_mu_t e_A[MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K], ff_mu_t e_B[MIRATH_PARAM_K],
                          const ff_mu_t in_X[MIRATH_PARAM_M * MIRATH_PARAM_R],
                          const ff_mu_t in_Y[MIRATH_PARAM_M * (MIRATH_PARAM_N - MIRATH_PARAM_R)]) {
    ff_mu_t T[MIRATH_PARAM_M * MIRATH_PARAM_N];

    uint32_t n_bytes1 = mirath_matrix_ff_mu_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_R);
    uint32_t n_bytes2 = mirath_matrix_ff_mu_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_N - MIRATH_PARAM_R);

    memcpy((uint8_t *)T, (uint8_t *)in_X, n_bytes1);
    memcpy((uint8_t *)T + n_bytes1, (uint8_t *)in_Y, n_bytes2);

    n_bytes1 = mirath_matrix_ff_mu_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1);
    n_bytes2 = mirath_matrix_ff_mu_bytes_size(MIRATH_PARAM_K, 1);

    memcpy((uint8_t *)e_A, (uint8_t *)T, n_bytes1);
    memcpy((uint8_t *)e_B, (uint8_t *)T + n_bytes1, n_bytes2);
}

void emulateMPC_mu(ff_mu_t base_alpha[MIRATH_PARAM_RHO],
                   ff_mu_t mid_alpha[MIRATH_PARAM_RHO],
                   const ff_t S[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_R)],
                   const ff_mu_t rnd_S[MIRATH_PARAM_M * MIRATH_PARAM_R],
                   const ff_t C[mirath_matrix_ff_bytes_size(MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R)],
                   const ff_mu_t rnd_C[MIRATH_PARAM_R * (MIRATH_PARAM_N - MIRATH_PARAM_R)],
                   const ff_mu_t v[MIRATH_PARAM_RHO], ff_mu_t rnd_v[MIRATH_PARAM_RHO],
                   const ff_mu_t gamma[MIRATH_PARAM_RHO * (MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K)],
                   const ff_t H[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, MIRATH_PARAM_K)]) {

    ff_mu_t aux_E[MIRATH_PARAM_M * (MIRATH_PARAM_N - MIRATH_PARAM_R)];
    ff_mu_t e_A[MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K];
    ff_mu_t e_B[MIRATH_PARAM_K];

    // rnd_S * rnd_C
    mirath_matrix_ff_mu_product(aux_E, rnd_S, rnd_C, MIRATH_PARAM_M, MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R);

    ff_mu_t zero[MIRATH_PARAM_M * MIRATH_PARAM_R] = {0};

    split_codeword_ff_mu(e_A, e_B, zero, aux_E);

    ff_mu_t tmp[MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K];

    // H * e_B
    mirath_matrix_ff_mu_product_ff1mu(tmp, H, e_B, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, MIRATH_PARAM_K, 1);
    // e_A + (H * e_B)
    mirath_vector_ff_mu_add(tmp, tmp, e_A, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K);
    // gamma * [e_A + (H * e_B)]
    mirath_matrix_ff_mu_product(base_alpha, gamma, tmp, MIRATH_PARAM_RHO, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1);
    // gamma * [e_A + (H * e_B)] + rnf_V
    mirath_vector_ff_mu_add(base_alpha, base_alpha, rnd_v, MIRATH_PARAM_RHO);

    ff_mu_t aux_s[MIRATH_PARAM_M * MIRATH_PARAM_R];
    ff_mu_t aux_c[MIRATH_PARAM_R * MIRATH_PARAM_N - MIRATH_PARAM_R];
    ff_mu_t aux_sc[MIRATH_PARAM_M * (MIRATH_PARAM_N - MIRATH_PARAM_R)];

    ff_t sc[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_N - MIRATH_PARAM_R)];

    // S + rnd_S
    mirath_matrix_ff_mu_add_mu1ff(aux_s, rnd_S, S, MIRATH_PARAM_M, MIRATH_PARAM_R);
    // C + rnd_C
    mirath_matrix_ff_mu_add_mu1ff(aux_c, rnd_C, C, MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R);
    // (S + rnd_S)(C + rnd_C)
    mirath_matrix_ff_mu_product(aux_sc, aux_s, aux_c, MIRATH_PARAM_M, MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R);
    // (S + rnd_S)(C + rnd_C) - base_E
    mirath_matrix_ff_mu_add(aux_E, aux_E, aux_sc, MIRATH_PARAM_M, MIRATH_PARAM_N - MIRATH_PARAM_R);
    // S * C
    mirath_matrix_ff_product(sc, S, C, MIRATH_PARAM_M, MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R);
    // (S + rnd_S)(C + rnd_C) - base_E - (S * C)
    mirath_matrix_ff_mu_add_mu1ff(aux_E, aux_E, sc, MIRATH_PARAM_M, MIRATH_PARAM_N - MIRATH_PARAM_R);

    split_codeword_ff_mu(e_A, e_B, rnd_S, aux_E);

    // H * e'_B
    mirath_matrix_ff_mu_product_ff1mu(tmp, H, e_B, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, MIRATH_PARAM_K, 1);
    // e'_A + (H * e'_B)
    mirath_vector_ff_mu_add(tmp, tmp, e_A, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K);
    // gamma * [e'_A + (H * e'_B)]
    mirath_matrix_ff_mu_product(mid_alpha, gamma, tmp, MIRATH_PARAM_RHO, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1);
    // gamma * [e_A + (H * e_B)] + v
    mirath_vector_ff_mu_add(mid_alpha, mid_alpha, v, MIRATH_PARAM_RHO);
}

void emulateparty_mu(ff_mu_t base_alpha[MIRATH_PARAM_RHO], const ff_mu_t p,
                     const ff_mu_t share_S[MIRATH_PARAM_M * MIRATH_PARAM_R],
                     const ff_mu_t share_C[MIRATH_PARAM_R * (MIRATH_PARAM_N - MIRATH_PARAM_R)],
                     const ff_mu_t share_v[MIRATH_PARAM_RHO],
                     const ff_mu_t gamma[MIRATH_PARAM_RHO * (MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K)],
                     const ff_t H[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, MIRATH_PARAM_K)],
                     const ff_t y[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1)],
                     const ff_mu_t mid_alpha[MIRATH_PARAM_RHO]) {

    ff_mu_t share_E[MIRATH_PARAM_M * MIRATH_PARAM_N];
    ff_mu_t e_A[MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K];
    ff_mu_t e_B[MIRATH_PARAM_K];

    ff_mu_t aux[MIRATH_PARAM_M * (MIRATH_PARAM_N - MIRATH_PARAM_R)];
    ff_mu_t Ts[MIRATH_PARAM_M * MIRATH_PARAM_R] = {0};

    // p * share_S
    mirath_matrix_ff_mu_add_multiple_2(Ts, p, share_S, MIRATH_PARAM_M, MIRATH_PARAM_R);
    // share_S * share_C
    mirath_matrix_ff_mu_product(aux, share_S, share_C, MIRATH_PARAM_M, MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R);

    split_codeword_ff_mu(e_A, e_B, Ts, aux);

    ff_mu_t tmp[MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K];

    // H * e_B
    mirath_matrix_ff_mu_product_ff1mu(tmp, H, e_B, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, MIRATH_PARAM_K, 1);
    // e_A + (H * e_B)
    mirath_vector_ff_mu_add(tmp, tmp, e_A, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K);
    // (e_A + (H * e_B)) - y * p^2
    mirath_vector_ff_mu_add_multiple_ff(tmp, tmp, mirath_ff_mu_mult(p, p), y, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K);
    // gamma * [(e_A + (H * e_B)) - y * p^2]
    mirath_matrix_ff_mu_product(base_alpha, gamma, tmp, MIRATH_PARAM_RHO, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1);
    // gamma * [(e_A + (H * e_B)) - y * p^2] + share_v
    mirath_vector_ff_mu_add(base_alpha, base_alpha, share_v, MIRATH_PARAM_RHO);

    // share_alpha - mid_alpha * p
    mirath_vector_ff_mu_add_multiple(base_alpha, base_alpha, p, mid_alpha, MIRATH_PARAM_RHO);
}
