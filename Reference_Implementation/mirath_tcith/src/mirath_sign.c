#include <stdint.h>
#include <string.h>

#include "mirath_parameters.h"
#include "prng.h"
#include "random.h"
#include "matrix_ff_mu.h"
#include "vector_ff_mu.h"
#include "mirath_matrix_ff.h"
#include "parsing.h"
#include "mirath_tree.h"
#include "mirath_tcith.h"

const ff_mu_t phi_table_ff_mu[MIRATH_N_PARTIES] = {
        0, 2, 4, 8, 16, 32, 64, 128, 27, 54, 108, 216, 171, 77, 154, 47, 94, 188, 99, 198, 151, 53, 106, 212, 179, 125,
        250, 239, 197, 145, 57, 114, 228, 211, 189, 97, 194, 159, 37, 74, 148, 51, 102, 204, 131, 29, 58, 116, 232, 203,
        141, 1, 2, 4, 8, 16, 32, 64, 128, 27, 54, 108, 216, 171, 77, 154, 47, 94, 188, 99, 198, 151, 53, 106, 212, 179,
        125, 250, 239, 197, 145, 57, 114, 228, 211, 189, 97, 194, 159, 37, 74, 148, 51, 102, 204, 131, 29, 58, 116, 232,
        203, 141, 1, 2, 4, 8, 16, 32, 64, 128, 27, 54, 108, 216, 171, 77, 154, 47, 94, 188, 99, 198, 151, 53, 106, 212,
        179, 125, 250, 239, 197, 145, 57, 114, 228, 211, 189, 97, 194, 159, 37, 74, 148, 51, 102, 204, 131, 29, 58, 116,
        232, 203, 141, 1, 2, 4, 8, 16, 32, 64, 128, 27, 54, 108, 216, 171, 77, 154, 47, 94, 188, 99, 198, 151, 53, 106,
        212, 179, 125, 250, 239, 197, 145, 57, 114, 228, 211, 189, 97, 194, 159, 37, 74, 148, 51, 102, 204, 131, 29, 58,
        116, 232, 203, 141, 1, 2, 4, 8, 16, 32, 64, 128, 27, 54, 108, 216, 171, 77, 154, 47, 94, 188, 99, 198, 151, 53,
        106, 212, 179, 125, 250, 239, 197, 145, 57, 114, 228, 211, 189, 97, 194, 159, 37, 74, 148, 51, 102, 204, 131,
        29, 58, 116, 232, 203, 141, 1
};

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

int mirath_keypair(uint8_t *pk, uint8_t *sk) {
    seed_t seed_sk;
    seed_t seed_pk;

    mirath_prng_t prng;

    ff_t *S;
    ff_t *C;
    ff_t SC[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_R) + mirath_matrix_ff_bytes_size(MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R)];
    ff_t H[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, MIRATH_PARAM_K)];
    ff_t *e_A;
    ff_t *e_B;

    ff_t T[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, (MIRATH_PARAM_N - MIRATH_PARAM_R))];
    ff_t E[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_N)];

    // step 1
    randombytes(seed_sk, MIRATH_SECURITY_BYTES);
    //step 2
    randombytes(seed_pk, MIRATH_SECURITY_BYTES);

    // step 3
    mirath_prng_init(&prng, NULL, seed_sk);
    mirath_prng(&prng, SC, mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_R) + mirath_matrix_ff_bytes_size(MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R));
    S = SC;
    C = SC + mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_R);
    mirath_matrix_set_to_ff(S, MIRATH_PARAM_M, MIRATH_PARAM_R);
    mirath_matrix_set_to_ff(C, MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R);

    // step 4
    mirath_prng_init(&prng, NULL, seed_pk);
    mirath_matrix_ff_init_random(H, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, MIRATH_PARAM_K, &prng);

    // step 5
    mirath_matrix_ff_product(T, S, C, MIRATH_PARAM_M, MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R);
    mirath_matrix_ff_horizontal_concat(E, S, T, MIRATH_PARAM_M, MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R);

    // step 6
    // SplitCodeword(Flatten(E))
    e_A = E;
    e_B = E + mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1);

    // step 7
    ff_t y[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1)];
    mirath_matrix_ff_product(y, H, e_B, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, MIRATH_PARAM_K, 1);

    mirath_vec_ff_add_arith(y, y, e_A, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K);

    // step 8
    unparse_public_key(pk, seed_pk, y);

    unparse_secret_key(sk, seed_sk, seed_pk);

    return 0;
}

int mirath_sign(uint8_t *sig_msg, size_t *sig_msg_len, uint8_t *msg, size_t msg_len, uint8_t *sk) {
    uint8_t salt[MIRATH_PARAM_SALT_BYTES];
    seed_t rseed;
    uint8_t salt_and_rseed[MIRATH_PARAM_SALT_BYTES + MIRATH_SECURITY_BYTES] = {0};
    mirath_tree_leaves_t seeds;
    mirath_tree_t tree;

    mirath_tcith_commit_t *commits[MIRATH_PARAM_TAU];
    mirath_tcith_commit_1_t commits_1 = {0};
    mirath_tcith_commit_2_t commits_2 = {0};
    mirath_tcith_commit_set_as_a_grid_list(commits, &commits_1, &commits_2);

    hash_t hash1;
    hash_t hash2_partial;
    hash_t hash2;
    hash_t h_com;

    hash_ctx_t hash_ctx;

    mirath_prng_t prng;

    ff_t S[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_R)];
    ff_t C[mirath_matrix_ff_bytes_size(MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R)];
    ff_t H[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, MIRATH_PARAM_K)];

    ff_mu_t gamma[MIRATH_PARAM_RHO * (MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K)];
    ff_mu_t v[MIRATH_PARAM_TAU][MIRATH_PARAM_RHO];
    ff_mu_t rnd_S[MIRATH_PARAM_TAU][MIRATH_PARAM_M * MIRATH_PARAM_R];
    ff_mu_t rnd_C[MIRATH_PARAM_TAU][MIRATH_PARAM_R * (MIRATH_PARAM_N - MIRATH_PARAM_R)];
    ff_mu_t rnd_v[MIRATH_PARAM_TAU][MIRATH_PARAM_RHO];

    ff_t aux[MIRATH_PARAM_TAU][mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_R + MIRATH_PARAM_R * (MIRATH_PARAM_N - MIRATH_PARAM_R), 1)];
    ff_mu_t mid_alpha[MIRATH_PARAM_TAU][MIRATH_PARAM_RHO];

    // Phase 0: Initialization
    // step 1 and 2
    parse_secret_key(S, C, H, sk);

    // step 3
    randombytes(salt, MIRATH_PARAM_SALT_BYTES);

    // Phase 1: Sharing and commitments
    // step 4
    randombytes(rseed, MIRATH_SECURITY_BYTES);

    // step 5
    mirath_tree_init(&tree);
    memcpy(&salt, salt_and_rseed, MIRATH_PARAM_SALT_BYTES);                                 // salt
    memcpy(tree.nodes[0], &salt_and_rseed[MIRATH_PARAM_SALT_BYTES], MIRATH_SECURITY_BYTES); // root seed
    tree.nonempty[0] = 1;
    mirath_tree_prg(&tree, salt, 0);
    mirath_tree_get_leaves(seeds, &tree);

    // step 6 and 7
    hash_init(&hash_ctx);
    hash_update(hash_ctx, salt, MIRATH_PARAM_SALT_BYTES);
    hash_update(hash_ctx, h_com, 2 * MIRATH_SECURITY_BYTES);

    for (uint32_t e = 0; e < MIRATH_PARAM_TAU; e++) {
//        // TODO: buildsharing
//        // NOTE: Below is what TCitH.BuildSharing_N performs (we need some computed data in MPC steps and signature
//        size_t N = e < MIRATH_PARAM_TAU_1? MIRATH_PARAM_N_1 : MIRATH_PARAM_N_2;
//        for(size_t i = 0; i < N; i++) {
//            size_t idx = mirath_tcith_psi(i, e);
//
//            // Compute commit and add it to ctx_h1
//            mirath_tcith_commit(commits[e][i], salt, e, i, seeds[idx]);
//            hash_update(&hash_ctx, commits[e][i], 2 * MIRATH_SECURITY_BYTES);
//            mirath_tcith_share_sample(si, Ci, vi, seeds[idx], salt); // Sampling bytes and parsing to the vec/mat (Fq)
//
//            vec_add(ss[e], ss[e], si, MIRATH_PARAM_R - 1);
//            mat_fq_add(Cs[e], Cs[e], Ci, MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R);
//            vec_add(vs[e], vs[e], vi, MIRATH_PARAM_RHO);
//
//            // Compute random shares
//            ff_mu_t phi_i;
//            // Apply phi with index i;
//
//            vec_scalar_mul(si, si, phi_i, MIRATH_PARAM_R - 1);
//            mat_fq_mul_by_constant(Di, Ci, phi_i, MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R);
//            vec_scalar_mul(vi, vi, phi_i, MIRATH_PARAM_RHO);
//
//            vec_add(shares.s[e], shares.s[e], si, MIRATH_PARAM_R - 1);
//            mat_add(shares.C[e], shares.C[e], Di, MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R);
//            vec_add(shares.v[e], shares.v[e], vi, MIRATH_PARAM_RHO);
//        }
//
//        // Recall, add = sub in fields of characteristic two
//
//        // Compute (s - ss[e])
//        vec_add(ss[e], (rbc_vec)&s[1], ss[e], MIRATH_PARAM_R - 1);
//        // Compute (C - Cs[e])
//        mat_fq_add(Cs[e], C, Cs[e], MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R);
//
//        // Compute aux from s[e] and C[e]

        hash_update(hash_ctx, aux[e], mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_R + MIRATH_PARAM_R * (MIRATH_PARAM_N - MIRATH_PARAM_R), 1));
    }

    // Phase 2: First challenge (MPC challenge)
    // step 8
    hash_finalize(hash_ctx, hash1);

    // step 9
    mirath_prng_init(&prng, hash1, NULL);
    mirath_prng(&prng, gamma, sizeof(ff_mu_t) * (MIRATH_PARAM_RHO * (MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K)));

    // Phase 3: MPC simulation
    // step 10 adn 11
    hash_init(&hash_ctx);
    hash_update(hash_ctx, msg, msg_len);
    hash_update(hash_ctx, salt, MIRATH_PARAM_SALT_BYTES);
    hash_update(hash_ctx, hash1, 2 * MIRATH_SECURITY_BYTES);

    for (uint32_t e = 0; e < MIRATH_PARAM_TAU; e++) {
        ff_mu_t base_alpha[MIRATH_PARAM_M * (MIRATH_PARAM_N - MIRATH_PARAM_R)];

        emulateMPC_mu(base_alpha, mid_alpha[e], S, rnd_S[e], C, rnd_C[e], v[e], rnd_v[e], gamma, H);

        hash_update(hash_ctx, base_alpha, sizeof(ff_mu_t) * (MIRATH_PARAM_M * (MIRATH_PARAM_N - MIRATH_PARAM_R)));
        hash_update(hash_ctx, mid_alpha[e], sizeof(ff_mu_t) * MIRATH_PARAM_RHO);
    }

    // step 12
    hash_finalize(hash_ctx, hash2_partial);

    // Phase 4: Second challenge (view opening)
    // step from 13 to 20
    uint8_t ctr;
    ctr = 0;

    // Phase 5: Signature
    // step 21

    return 0;
}

int mirath_verify(uint8_t *msg, size_t *msg_len, uint8_t *sig_msg, size_t sig_msg_len, uint8_t *pk) {
    seed_t seed_pk;
    uint8_t salt[MIRATH_PARAM_SALT_BYTES];

    hash_t hash1;
    hash_t hash2_partial;
    hash_t hash2;
    hash_t hash2_computed;
    hash_t h_com;

    hash_ctx_t hash_ctx;
    mirath_prng_t prng;

    ff_mu_t share_S[MIRATH_PARAM_TAU][MIRATH_PARAM_M * MIRATH_PARAM_N];
    ff_mu_t share_C[MIRATH_PARAM_TAU][MIRATH_PARAM_R * (MIRATH_PARAM_N - MIRATH_PARAM_R)];
    ff_mu_t share_v[MIRATH_PARAM_TAU][MIRATH_PARAM_R * (MIRATH_PARAM_N - MIRATH_PARAM_R)];
    ff_mu_t gamma[MIRATH_PARAM_RHO * (MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K)];
    ff_t H[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, MIRATH_PARAM_K)];
    ff_t y[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1)];

    ff_t aux[MIRATH_PARAM_TAU][mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_R + MIRATH_PARAM_R * (MIRATH_PARAM_N - MIRATH_PARAM_R), 1)];
    ff_mu_t mid_alpha[MIRATH_PARAM_TAU][MIRATH_PARAM_RHO];

    uint32_t i_star[MIRATH_PARAM_TAU];

    uint8_t ctr;

    // step 1
    parse_public_key(seed_pk, y, pk);

    // Phase 0: Parsing and expansion
    // step 2

    // step 3
    mirath_prng_init(&prng, NULL, seed_pk);
    mirath_matrix_ff_init_random(H, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, MIRATH_PARAM_K, &prng);

    // step 4

    // Phase 1: Recomputing shares and commitments
    // step 5

    // step 6, and 7
    hash_init(&hash_ctx);
    hash_update(hash_ctx, salt, MIRATH_PARAM_SALT_BYTES);
    hash_update(hash_ctx, h_com, 2 * MIRATH_SECURITY_BYTES);

    for (uint32_t e = 0; e < MIRATH_PARAM_TAU; e++) {
        // TODO: computeshare

        hash_update(hash_ctx, aux[e], mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_R + MIRATH_PARAM_R * (MIRATH_PARAM_N - MIRATH_PARAM_R), 1));
    }

    // step 8
    hash_finalize(hash_ctx, hash1);

    // step 9
    mirath_prng_init(&prng, hash1, NULL);
    mirath_prng(&prng, gamma, sizeof(ff_mu_t) * (MIRATH_PARAM_RHO * (MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K)));

    // Phase 2: MPC simulation
    // step 10 and 11
    hash_init(&hash_ctx);
    hash_update(hash_ctx, msg, *msg_len);
    hash_update(hash_ctx, salt, MIRATH_PARAM_SALT_BYTES);
    hash_update(hash_ctx, hash1, 2 * MIRATH_SECURITY_BYTES);

    for (uint32_t e = 0; e < MIRATH_PARAM_TAU; e++) {
        ff_mu_t base_alpha[MIRATH_PARAM_RHO];

        emulateparty_mu(base_alpha, phi_table_ff_mu[i_star[e]], share_S[e], share_C[e], share_v[e], gamma, H, y, mid_alpha[e]);

        hash_update(hash_ctx, base_alpha, sizeof(ff_mu_t) * (MIRATH_PARAM_M * (MIRATH_PARAM_N - MIRATH_PARAM_R)));
        hash_update(hash_ctx, mid_alpha[e], sizeof(ff_mu_t) * MIRATH_PARAM_RHO);
    }

    // Phase 3: Verification
    // step 12
    hash_finalize(hash_ctx, hash2_partial);

    // step 12
    hash_init(&hash_ctx);
    hash_update(hash_ctx, &ctr, sizeof(ctr));
    hash_update(hash_ctx, hash2_partial, 2 * MIRATH_SECURITY_BYTES);
    hash_finalize(hash_ctx, hash2_computed);

    if (hash_equal(hash2, hash2_computed)) {
        return 0;
    }

    return -1;
}