/**
 * @file mirath_tcith.c
 * @brief Implementation of TCitH related functions
 */

#include "prng.h"
#include "random.h"
#include "q_16/matrix_ff_mu.h"
#include "q_16/vector_ff_mu.h"
#include "mirath_matrix_ff.h"
#include "mirath_tcith.h"
#include "mirath_tree.h"
#include "parsing.h"

/**
* \fn void mirath_tcith_shift_to_left_array(uint8_t *inout_a, size_t length)
* \brief This function performs a shift to right of the input.
*
* \param[in/out] inout_a uint8_t* Representation of a byte string
* \param[in] length size_t Representation of the byte string length
*/
static inline void mirath_tcith_shift_to_right_array(uint8_t *string, const size_t length) {
    for(size_t i = 0; i < length - 1; i++) {
        mirath_tcith_shift_to_right(string[i], string[i + 1], string[i], 1, 8);
    }
    string[length - 1] >>= 1;
}

void mirath_tcith_internal_steps_pk(ff_t y[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1)],
                                     const ff_t S[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_R)],
                                     const ff_t C[mirath_matrix_ff_bytes_size(MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R)],
                                     const ff_t H[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, MIRATH_PARAM_K)]) {
    ff_t e_A[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1)] = {0};
    ff_t e_B[mirath_matrix_ff_bytes_size(MIRATH_PARAM_K, 1)] = {0};

    ff_t T[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, (MIRATH_PARAM_N - MIRATH_PARAM_R))];
    ff_t E[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_N)];

    mirath_matrix_ff_product(T, S, C, MIRATH_PARAM_M, MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R);
    mirath_matrix_ff_horizontal_concat(E, S, T, MIRATH_PARAM_M, MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R);

    const uint32_t n_bytes_A = mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1);

    memcpy(e_A, E, n_bytes_A);

    if ((MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K) & 1u) {
        const uint32_t n_bytes_B = mirath_matrix_ff_bytes_size(MIRATH_PARAM_K, 1);

        e_A[n_bytes_A - 1] &= 0x0F;
        memcpy(e_B, E + (n_bytes_A - 1), mirath_matrix_ff_bytes_size(MIRATH_PARAM_K, 1));

        for (uint32_t i = 0; i < n_bytes_B - 1; i++) {
            mirath_tcith_shift_to_right(e_B[i], e_B[i + 1], e_B[i], 4, 8);
        }
        e_B[n_bytes_B - 1] >>= 4;
    } else {
        memcpy(e_B, E + n_bytes_A, mirath_matrix_ff_bytes_size(MIRATH_PARAM_K, 1));
    }

    mirath_matrix_ff_product(y, H, e_B, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, MIRATH_PARAM_K, 1);

    mirath_vec_ff_add_arith(y, y, e_A, n_bytes_A);
}

void mirath_tciht_compute_public_key(uint8_t *pk, const uint8_t *sk,
                                     const ff_t S[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_R)],
                                     const ff_t C[mirath_matrix_ff_bytes_size(MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R)],
                                     const ff_t H[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, MIRATH_PARAM_K)]) {

    ff_t y[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1)];

    mirath_tcith_internal_steps_pk(y, S, C, H);

    unparse_public_key(pk, sk + MIRATH_SECURITY_BYTES, y);
}

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

void mirath_tcith_commit(mirath_tcith_commit_t commit, const uint8_t *salt, uint16_t e, uint16_t i, const uint8_t *seed) {
    const uint8_t domain_separator = DOMAIN_SEPARATOR_COMMITMENT;

    hash_ctx_t ctx;
    hash_init(&ctx);
    hash_update(ctx, &domain_separator, sizeof(uint8_t));
    hash_update(ctx, salt, MIRATH_PARAM_SALT_BYTES);
    hash_update(ctx, (uint8_t *)&e, sizeof(uint16_t));
    hash_update(ctx, (uint8_t *)&i, sizeof(uint16_t));
    hash_update(ctx, seed, MIRATH_SECURITY_BYTES);
    hash_finalize(commit, ctx);
}

size_t mirath_tcith_psi(size_t i, size_t e) {
    if (i < MIRATH_PARAM_N_2) {
        return i * MIRATH_PARAM_TAU + e;
    }
    else {
        return MIRATH_PARAM_N_2 * MIRATH_PARAM_TAU + (i - MIRATH_PARAM_N_2) * MIRATH_PARAM_TAU_1 + e;
    }
}

void build_sharing_N(ff_t aux[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_R + MIRATH_PARAM_R * (MIRATH_PARAM_N - MIRATH_PARAM_R), 1)],
                   ff_mu_t rnd_S[MIRATH_PARAM_M * MIRATH_PARAM_R],
                   ff_mu_t rnd_C[MIRATH_PARAM_R * (MIRATH_PARAM_N - MIRATH_PARAM_R)],
                   ff_mu_t rnd_v[MIRATH_PARAM_RHO], ff_mu_t v[MIRATH_PARAM_RHO],
                   mirath_tcith_commit_t *commits[MIRATH_PARAM_TAU],
                   const mirath_tree_leaves_t seeds,
                   const ff_t S[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_R)],
                   const ff_t C[mirath_matrix_ff_bytes_size(MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R)],
                   const uint8_t salt[MIRATH_PARAM_SALT_BYTES], const uint32_t e) {

    const uint16_t N = e < MIRATH_PARAM_TAU_1 ? MIRATH_PARAM_N_1 : MIRATH_PARAM_N_2;

    memset(aux, 0, mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_R + MIRATH_PARAM_R * (MIRATH_PARAM_N - MIRATH_PARAM_R), 1));
    memset(rnd_S, 0, sizeof(ff_mu_t) * (MIRATH_PARAM_M * MIRATH_PARAM_R));
    memset(rnd_C, 0, sizeof(ff_mu_t) * (MIRATH_PARAM_R * (MIRATH_PARAM_N - MIRATH_PARAM_R)));
    memset(rnd_v, 0, sizeof(ff_mu_t) * MIRATH_PARAM_RHO);
    memset(v, 0, sizeof(ff_mu_t) * MIRATH_PARAM_RHO);

    ff_t acc_share_S[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_R)] = {0};
    ff_t acc_share_C[mirath_matrix_ff_bytes_size(MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R)] = {0};

    for (uint16_t i = 0; i < N; i++) {
        mirath_prng_t prng;
        ff_t Si[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_R)];
        ff_t Ci[mirath_matrix_ff_bytes_size(MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R)];
        ff_mu_t vi[MIRATH_PARAM_RHO];

        const uint32_t idx = mirath_tcith_psi((size_t)i, (size_t)e);

        mirath_tcith_commit(commits[e][i], salt, e, i, seeds[idx]);

        mirath_prng_init(&prng, NULL, seeds[idx], MIRATH_SECURITY_BYTES);

        mirath_matrix_ff_init_random(Si, MIRATH_PARAM_M, MIRATH_PARAM_R, &prng);
        mirath_matrix_ff_init_random(Ci, MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R, &prng);
        mirath_prng(&prng, vi, sizeof(ff_mu_t) * MIRATH_PARAM_RHO);

        mirath_matrix_ff_add(acc_share_S, acc_share_S, Si, MIRATH_PARAM_M, MIRATH_PARAM_R);
        mirath_matrix_ff_add(acc_share_C, acc_share_C, Ci, MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R);
        mirath_vector_ff_mu_add(v, v, vi, MIRATH_PARAM_RHO);

        const ff_mu_t phi_i = (ff_mu_t)i;
        mirath_matrix_ff_mu_add_multiple_ff(rnd_S, phi_i, Si, MIRATH_PARAM_M, MIRATH_PARAM_R);
        mirath_matrix_ff_mu_add_multiple_ff(rnd_C, phi_i, Ci, MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R);
        mirath_vector_ff_mu_add_multiple(rnd_v, rnd_v, phi_i, vi, MIRATH_PARAM_RHO);
    }

    // S - acc_S
    mirath_matrix_ff_add(aux, S, acc_share_S, MIRATH_PARAM_M, MIRATH_PARAM_R);
    const uint32_t n_bytes = mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_R);
    // C - acc_C
    mirath_matrix_ff_add(aux + n_bytes, C, acc_share_C, MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R);
}

void compute_share(ff_mu_t share_S[MIRATH_PARAM_M * MIRATH_PARAM_R],
                   ff_mu_t share_C[MIRATH_PARAM_R * (MIRATH_PARAM_N - MIRATH_PARAM_R)],
                   ff_mu_t share_v[MIRATH_PARAM_RHO], mirath_tcith_commit_t *commits[MIRATH_PARAM_TAU],
                   const uint32_t i_star, const mirath_tree_leaves_t seeds, const uint32_t e,
                   const uint8_t salt[MIRATH_PARAM_SALT_BYTES], mirath_tcith_commit_t commits_i_star[MIRATH_PARAM_TAU],
                   const ff_t aux[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_R + MIRATH_PARAM_R * (MIRATH_PARAM_N - MIRATH_PARAM_R), 1)]) {

    const ff_t *aux_S = aux;
    const uint32_t n_bytes = mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_R);
    // TODO: here we could have the same situation when nrows is odd
    const ff_t *aux_C = aux + n_bytes;

    memcpy(commits[e][i_star], commits_i_star[e], 2 * MIRATH_SECURITY_BYTES);

    const uint16_t N = e < MIRATH_PARAM_TAU_1 ? MIRATH_PARAM_N_1 : MIRATH_PARAM_N_2;
    for (uint16_t i = 0; i < N; i++) {
        if ((uint32_t)i != i_star) {
            mirath_prng_t prng;
            ff_t Si[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_R)];
            ff_t Ci[mirath_matrix_ff_bytes_size(MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R)];
            ff_mu_t vi[MIRATH_PARAM_RHO];

            const uint32_t idx = mirath_tcith_psi((size_t) i, (size_t) e);

            mirath_tcith_commit(commits[e][i], salt, e, i, seeds[idx]);

            mirath_prng_init(&prng, NULL, seeds[idx], MIRATH_SECURITY_BYTES);

            mirath_matrix_ff_init_random(Si, MIRATH_PARAM_M, MIRATH_PARAM_R, &prng);
            mirath_matrix_ff_init_random(Ci, MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R, &prng);
            mirath_prng(&prng, vi, sizeof(ff_mu_t) * MIRATH_PARAM_RHO);

            const ff_mu_t sc = (ff_mu_t) ((uint16_t) i_star ^ i);

            mirath_matrix_ff_mu_add_multiple_ff(share_S, sc, Si, MIRATH_PARAM_M, MIRATH_PARAM_R);
            mirath_matrix_ff_mu_add_multiple_ff(share_C, sc, Ci, MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R);
            mirath_vector_ff_mu_add_multiple(share_v, share_v, sc, vi, MIRATH_PARAM_RHO);
        }
    }

    const ff_mu_t phi_i = i_star;
    mirath_matrix_ff_mu_add_multiple_ff(share_S, phi_i, aux_S, MIRATH_PARAM_M, MIRATH_PARAM_R);
    mirath_matrix_ff_mu_add_multiple_ff(share_C, phi_i, aux_C, MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R);
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
    // gamma * [e_A + (H * e_B)] + rnd_V
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

/**
* \fn void mirath_tcith_compute_challenge_2(mirath_tcith_challenge_t challenge, const uint8_t *seed_input, const uint8_t *salt)
* \brief This function generates the second challenges from an input seed
*
* \param[out] challenge mirath_tcith_challenge_t Representation of challenge
* \param[in] seed_input String containing the input seed
* \param[in] salt String containing the salt
*/
void mirath_tcith_compute_challenge_2(mirath_tcith_challenge_t challenge, const uint8_t *seed_input, const uint8_t *salt) {
    uint8_t random[MIRATH_PARAM_CHALLENGE_2_BYTES] = {0}, mask;

    mirath_prng_t prng;
    mirath_prng_init(&prng, salt, seed_input, 2 * MIRATH_SECURITY_BYTES); // seed is the double size

    memset(challenge, 0, sizeof(mirath_tcith_challenge_t));

    // generate MIRATH_PARAM_CHALLENGE_2_BYTES random bytes and store them in random
    mirath_prng(&prng, random, MIRATH_PARAM_CHALLENGE_2_BYTES);

    // Challenges concerning N_1
    mask = MIRATH_PARAM_N_1_MASK;
    for(size_t i = 0; i < MIRATH_PARAM_TAU_1; i++) {
        uint8_t block[MIRATH_PARAM_N_1_BYTES] = {0};
        memcpy(block, random, MIRATH_PARAM_N_1_BYTES);
        block[MIRATH_PARAM_N_1_BYTES - 1] &= mask;
        memcpy((uint8_t *)&challenge[i], block, MIRATH_PARAM_N_1_BYTES);
        // Left shift bits: starts
        for (size_t j = 0; j < MIRATH_PARAM_N_1_BITS; j++) {
            mirath_tcith_shift_to_right_array(random, MIRATH_PARAM_CHALLENGE_2_BYTES);
        }
        // Left shift bits: ends
    }

    // Challenges concerning N_2
    mask = MIRATH_PARAM_N_2_MASK;
    for(size_t i = 0; i < MIRATH_PARAM_TAU_2; i++) {
        uint8_t block[MIRATH_PARAM_N_2_BYTES] = {0};
        memcpy(block, random, MIRATH_PARAM_N_2_BYTES);
        block[MIRATH_PARAM_N_2_BYTES - 1] &= mask;
        memcpy((uint8_t *)&challenge[i + MIRATH_PARAM_TAU_1], block, MIRATH_PARAM_N_2_BYTES);
        // Left shift bits: starts
        for (size_t j = 0; j < MIRATH_PARAM_N_2_BITS; j++) {
            mirath_tcith_shift_to_right_array(random, MIRATH_PARAM_CHALLENGE_2_BYTES);
        }
        // Left shift bits: ends
    }
}

/**
* \fn uint8_t mirath_tcith_discard_input_challenge_2(const uint8_t *seed_input)
* \brief This function determines if the w most significant bits of the input are zero.
*
* \param[in] seed_input String containing the input seed
*/
uint8_t mirath_tcith_discard_input_challenge_2(const uint8_t *seed_input) {
    uint8_t output = 0x00;
    uint8_t mask = MIRATH_PARAM_HASH_2_MASK;
    for(size_t i = 0; i < MIRATH_PARAM_HASH_2_MASK_BYTES; i++) {
        output |= (uint8_t)((seed_input[i] & mask) != 0);
        mask = 0xFF;
    }

    return output;
}