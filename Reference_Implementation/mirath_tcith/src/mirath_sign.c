#include <stdint.h>
#include <string.h>

#include "prng.h"
#include "random.h"
#include "vector_ff_mu.h"
#include "mirath_matrix_ff.h"
#include "parsing.h"
#include "mirath_tree.h"
#include "mirath_tcith.h"

int mirath_keypair(uint8_t *pk, uint8_t *sk) {
    seed_t seed_sk;
    seed_t seed_pk;

    mirath_prng_t prng;

    ff_t *S;
    ff_t *C;
    ff_t SC[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_R) + mirath_matrix_ff_bytes_size(MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R)];
    ff_t H[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, MIRATH_PARAM_K)];

    ff_t y[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1)];

    // step 1
    randombytes(seed_sk, MIRATH_SECURITY_BYTES);
    //step 2
    randombytes(seed_pk, MIRATH_SECURITY_BYTES);

    // step 3
    mirath_prng_init(&prng, NULL, seed_sk, MIRATH_SECURITY_BYTES);
    mirath_prng(&prng, SC, mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_R) + mirath_matrix_ff_bytes_size(MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R));
    S = SC;
    C = SC + mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_R);
    mirath_matrix_set_to_ff(S, MIRATH_PARAM_M, MIRATH_PARAM_R);
    mirath_matrix_set_to_ff(C, MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R);

    // step 4
    mirath_prng_init(&prng, NULL, seed_pk, MIRATH_SECURITY_BYTES);
    mirath_matrix_ff_init_random(H, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, MIRATH_PARAM_K, &prng);

    // step 5, 6 and 7
    mirath_tcith_internal_steps_pk(y, S, C, H);

    // step 8
    unparse_public_key(pk, seed_pk, y);

    unparse_secret_key(sk, seed_sk, seed_pk);

    return 0;
}

int mirath_sign(uint8_t *sig_msg, size_t *sig_msg_len, uint8_t *msg, size_t msg_len, uint8_t *sk) {
    uint8_t salt[MIRATH_PARAM_SALT_BYTES];
    seed_t rseed;
    mirath_tree_leaves_t seeds;
    mirath_tree_t tree;

    mirath_tcith_commit_t *commits[MIRATH_PARAM_TAU];
    mirath_tcith_commit_1_t commits_1 = {0};
    mirath_tcith_commit_2_t commits_2 = {0};
    mirath_tcith_commit_set_as_a_grid_list(commits, &commits_1, &commits_2);
    
    uint8_t path[MIRATH_PARAM_TREE_LEAVES * MIRATH_SECURITY_BYTES] = {0};

    hash_t hash1;
    hash_t hash2_partial;
    hash_t hash2;

    hash_ctx_t hash_ctx;

    uint8_t domain_separator;

    mirath_prng_t prng;

    ff_t S[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_R)];
    ff_t C[mirath_matrix_ff_bytes_size(MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R)];
    ff_t H[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, MIRATH_PARAM_K)];

    uint8_t pk[MIRATH_PUBLIC_KEY_BYTES];

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
    mirath_tciht_compute_public_key(pk, sk, S, C, H);

    // step 3
    randombytes(salt, MIRATH_PARAM_SALT_BYTES);

    // Phase 1: Sharing and commitments
    // step 4
    randombytes(rseed, MIRATH_SECURITY_BYTES);

    // step 5, step 6 and 7
    mirath_tree_init(&tree);
    memcpy(tree.nodes[0], rseed, MIRATH_SECURITY_BYTES); // root seed
    tree.nonempty[0] = 1;
    mirath_tree_prg(&tree, salt, 0);
    mirath_tree_get_leaves(seeds, &tree);

    domain_separator = DOMAIN_SEPARATOR_HASH1;

    hash_init(&hash_ctx);
    hash_update(hash_ctx, &domain_separator, sizeof(uint8_t));
    hash_update(hash_ctx, salt, MIRATH_PARAM_SALT_BYTES);

    build_sharing_N(aux, rnd_S, rnd_C, rnd_v, hash_ctx, seeds, S, C, salt);

    // Phase 2: First challenge (MPC challenge)
    // step 8
    hash_finalize(hash1, hash_ctx);

    // step 9
    mirath_prng_init(&prng, hash1, NULL, 0);
    mirath_prng(&prng, gamma, sizeof(ff_mu_t) * (MIRATH_PARAM_RHO * (MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K)));

    domain_separator = DOMAIN_SEPARATOR_HASH2_PARTIAL;

    // Phase 3: MPC simulation
    // step 10 adn 11
    hash_init(&hash_ctx);
    hash_update(hash_ctx, &domain_separator, sizeof(uint8_t));
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
    hash_finalize(hash2_partial, hash_ctx);

    mirath_tcith_challenge_t i_star;
    uint16_t psi_i_star[MIRATH_PARAM_TAU];
    // Phase 4: Second challenge (view opening)
    uint16_t ctr = 0;
 retry:
     // Initialize hash_ctx
     hash_init(&hash_ctx);
     domain_separator = DOMAIN_SEPARATOR_HASH2;
     hash_update(hash_ctx, &domain_separator, sizeof(uint8_t));
     hash_update(hash_ctx, (uint8_t *)&ctr, sizeof(uint16_t));
     hash_update(hash_ctx, pk, MIRATH_PUBLIC_KEY_BYTES);
     hash_finalize(hash2, hash_ctx);

     if (mirath_tcith_discard_input_challenge_2(hash2)) {
         ctr += 1;
         goto retry;
     }

     mirath_tcith_compute_challenge_2(i_star, hash2, salt);
     // Next we map the challenges to the leaves position (GGM Tree optimization)
     for(size_t e = 0; e < MIRATH_PARAM_TAU; e++){
         size_t i = i_star[e];
         psi_i_star[e] = (uint16_t)mirath_tcith_psi(i, e); // store their respectively image under psi
     }

     size_t path_length = mirath_tree_get_sibling_path(path, &tree, psi_i_star, MIRATH_PARAM_TAU);
     if (path_length > MIRATH_PARAM_T_OPEN) {
         ctr += 1;
         memset(path, 0, sizeof(path_length) * MIRATH_SECURITY_BYTES);
         goto retry;
     }

    // Phase 5: Signature
    // step 21

    // --------------------------------------------------------------------------------------------------------- Signature
    // // se_str[] is the list of vectors aux_s
    // // Ce_str[] is the list of matrix aux_C
    // // mid_a_str[] is the list of vectors mid_alpha
    // memcpy(&sig_msg[0], salt, 2 * MIRATH_SECURITY_BYTES);
    // memcpy(&sig_msg[2 * MIRATH_SECURITY_BYTES], &ctr, sizeof(uint64_t));
    // memcpy(&sig_msg[2 * MIRATH_SECURITY_BYTES + sizeof(uint64_t)], h2, 2 * MIRATH_SECURITY_BYTES);
    // memcpy(&sig_msg[4 * MIRATH_SECURITY_BYTES + sizeof(uint64_t)], path, MIRATH_SECURITY_BYTES * MIRATH_PARAM_T_OPEN);

    // for(size_t e = 0; e < MIRATH_PARAM_TAU; e++) {
    //     // Vector s
    //     memcpy(&sig_msg[
    //             4 * MIRATH_SECURITY_BYTES + sizeof(uint64_t) + MIRATH_PARAM_T_OPEN * MIRATH_SECURITY_BYTES +
    //             e * (MIRATH_VEC_R_BYTES + MIRATH_MAT_BYTES + MIRATH_VEC_RHO_BYTES + 2 * MIRATH_SECURITY_BYTES)
    //             ],
    //             se_str[e],
    //             MIRATH_VEC_R_BYTES);

    //     // Matrix C
    //     memcpy(&sig_msg[
    //             4 * MIRATH_SECURITY_BYTES + sizeof(uint64_t) + MIRATH_PARAM_T_OPEN * MIRATH_SECURITY_BYTES +
    //             e * (MIRATH_VEC_R_BYTES + MIRATH_MAT_FQ_BYTES + MIRATH_VEC_RHO_BYTES + 2 * MIRATH_SECURITY_BYTES) +
    //             MIRATH_VEC_R_BYTES
    //             ],
    //             Ce_str[e],
    //             MIRATH_MAT_FQ_BYTES);

    //     // mid_a
    //     memcpy(&sig_msg[
    //             4 * MIRATH_SECURITY_BYTES + sizeof(uint64_t) + MIRATH_PARAM_T_OPEN * MIRATH_SECURITY_BYTES +
    //             e * (MIRATH_VEC_R_BYTES + MIRATH_MAT_FQ_BYTES + MIRATH_VEC_RHO_BYTES + 2 * MIRATH_SECURITY_BYTES) +
    //             MIRATH_VEC_R_BYTES +
    //             MIRATH_MAT_FQ_BYTES
    //             ],
    //             mid_a_str[e],
    //             MIRATH_VEC_RHO_BYTES);

    //     // Commitment concerning the hidden seed
    //     memcpy(&sig_msg[
    //             4 * MIRATH_SECURITY_BYTES + sizeof(uint64_t) + MIRATH_PARAM_T_OPEN * MIRATH_SECURITY_BYTES +
    //             e * (MIRATH_VEC_R_BYTES + MIRATH_MAT_FQ_BYTES + MIRATH_VEC_RHO_BYTES + 2 * MIRATH_SECURITY_BYTES) +
    //             MIRATH_VEC_R_BYTES +
    //             MIRATH_MAT_FQ_BYTES +
    //             MIRATH_VEC_RHO_BYTES
    //             ],
    //             commits[e][i_star[e]],
    //             2 * MIRATH_SECURITY_BYTES);
    // }

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
    mirath_prng_init(&prng, NULL, seed_pk, MIRATH_SECURITY_BYTES);
    mirath_matrix_ff_init_random(H, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, MIRATH_PARAM_K, &prng);

    // step 4

    // // Parse signature data
    // memcpy(salt, &signature[0], 2 * MIRATH_SECURITY_BYTES);
    // memcpy(&ctr, &signature[2 * MIRATH_SECURITY_BYTES], sizeof(uint64_t));
    // memcpy(h2,   &signature[2 * MIRATH_SECURITY_BYTES + sizeof(uint64_t)], 2 * MIRATH_SECURITY_BYTES);
    // memcpy(path, &signature[4 * MIRATH_SECURITY_BYTES + sizeof(uint64_t)], MIRATH_SECURITY_BYTES * MIRATH_PARAM_T_OPEN);

    // mirath_tcith_compute_challenge_2(i_star, h2, salt);
    // // Next we map the challenges to the leaves position (GGM Tree optimization)
    // for(size_t e = 0; e < MIRATH_PARAM_TAU; e++){
    //     size_t i = i_star[e];
    //     psi_i_star[e] = mirath_tcith_psi(i, e);
    // }

    // // Get sibling path length: starts
    // size_t path_length = 0;
    // for(size_t i = 0; i < MIRATH_PARAM_T_OPEN; i++) {
    //     uint8_t zero[MIRATH_SECURITY_BYTES] = {0};
    //     if (memcmp(zero, &path[i * MIRATH_SECURITY_BYTES], MIRATH_SECURITY_BYTES) == 0) { continue; }
    //     path_length += 1;
    // }
    // // Get sibling path length: ends

    // mirath_tree_init(&tree);
    // mirath_tree_get_seeds_from_path(&tree, psi_i_star, MIRATH_PARAM_TAU, path, path_length, salt, 0);
    // mirath_tree_get_leaves(seeds, &tree);

    // ++++++++
    // Phase 1: Recomputing shares and commitments
    // for(size_t e = 0; e < MIRATH_PARAM_TAU; e++) {
    //     // Set to zero the accumulator vector and matrix
    //     vec_set_zero(shares.s[e], MIRATH_PARAM_R - 1);
    //     mat_set_zero(shares.C[e], MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R);
    //     vec_set_zero(shares.v[e], MIRATH_PARAM_RHO);

    //     vec_set_zero(ss[e], MIRATH_PARAM_R - 1);
    //     mat_fq_set_zero(Cs[e], MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R);
    //     vec_set_zero(vs[e], MIRATH_PARAM_RHO);

    //     rbc_elt phi_i_star;
    //     mirath_tcith_phi(phi_i_star, i_star[e]);

    //     size_t N = e < MIRATH_PARAM_TAU_1? MIRATH_PARAM_N_1 : MIRATH_PARAM_N_2;
    //     for(size_t i = 0; i < N; i++) {
    //         // Set to zero the accumulator vector and matrix
    //         vec_set_zero(si, MIRATH_PARAM_R - 1);
    //         mat_fq_set_zero(Ci, MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R);
    //         mat_set_zero(Di, MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R);
    //         vec_set_zero(vi, MIRATH_PARAM_RHO);

    //         size_t idx = mirath_tcith_psi(i, e);
    //         if (i == (size_t)i_star[e]) {
    //             memcpy(commits[e][i],
    //                     &signature[
    //                     4 * MIRATH_SECURITY_BYTES + sizeof(uint64_t) + MIRATH_PARAM_T_OPEN * MIRATH_SECURITY_BYTES +
    //                     e * (MIRATH_VEC_R_BYTES + MIRATH_MAT_FQ_BYTES + MIRATH_VEC_RHO_BYTES + 2 * MIRATH_SECURITY_BYTES) +
    //                     MIRATH_VEC_R_BYTES +
    //                     MIRATH_MAT_FQ_BYTES +
    //                     MIRATH_VEC_RHO_BYTES
    //                     ],
    //                     2 * MIRATH_SECURITY_BYTES);
    //             hash_update(ctx_h1, commits[e][i], 2 * MIRATH_SECURITY_BYTES);
    //         } else {
    //             // Compute commit and add it to ctx_h1
    //             mirath_tcith_commit(commits[e][i], salt, e, i, seeds[idx]);
    //             hash_update(ctx_h1, commits[e][i], 2 * MIRATH_SECURITY_BYTES);
    //             mirath_tcith_share_sample(si, Ci, vi, seeds[idx], salt);

    //             // Compute shares
    //             // field element (ff_mu) phi_i;
    //             mirath_tcith_phi(phi_i, i);
    //             ff_mu_add(phi_i, phi_i_star, phi_i);

    //             vec_scalar_mul(si, si, phi_i, MIRATH_PARAM_R - 1);
    //             mat_fq_mul_by_constant(Di, Ci, phi_i, MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R);
    //             vec_scalar_mul(vi, vi, phi_i, MIRATH_PARAM_RHO);

    //             vec_add(shares.s[e], shares.s[e], si, MIRATH_PARAM_R - 1);
    //             mat_add(shares.C[e], shares.C[e], Di, MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R);
    //             vec_add(shares.v[e], shares.v[e], vi, MIRATH_PARAM_RHO);
    //         }
    //     }

    //     // Operations concerning vector se
    //     memcpy(se_str[e],
    //             &signature[
    //             4 * MIRATH_SECURITY_BYTES + sizeof(uint64_t) + MIRATH_PARAM_T_OPEN * MIRATH_SECURITY_BYTES +
    //             e * (MIRATH_VEC_R_BYTES + MIRATH_MAT_FQ_BYTES + MIRATH_VEC_RHO_BYTES + 2 * MIRATH_SECURITY_BYTES)
    //             ],
    //             MIRATH_VEC_R_BYTES);
    //     vec_from_string(si, MIRATH_PARAM_R - 1, se_str[e]);
    //     vec_scalar_mul(si, si, phi_i_star, MIRATH_PARAM_R - 1);
    //     vec_add(shares.s[e], shares.s[e], si, MIRATH_PARAM_R - 1);

    //     // Operations concerning matrix Ce
    //     memcpy(Ce_str[e],
    //             &signature[
    //             4 * MIRATH_SECURITY_BYTES + sizeof(uint64_t) + MIRATH_PARAM_T_OPEN * MIRATH_SECURITY_BYTES +
    //             e * (MIRATH_VEC_R_BYTES + MIRATH_MAT_FQ_BYTES + MIRATH_VEC_RHO_BYTES + 2 * MIRATH_SECURITY_BYTES) +
    //             MIRATH_VEC_R_BYTES
    //             ],
    //             MIRATH_MAT_FQ_BYTES);
    //     mat_fq_from_string(Ci, MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R, Ce_str[e]);
    //     mat_fq_mul_by_constant(Di, Ci, phi_i_star, MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R);
    //     mat_add(shares.C[e], shares.C[e], Di, MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R);
    // }
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
    hash_finalize(hash1, hash_ctx);

    // step 9
    mirath_prng_init(&prng, hash1, NULL, 0);
    mirath_prng(&prng, gamma, sizeof(ff_mu_t) * (MIRATH_PARAM_RHO * (MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K)));

    // Phase 2: MPC simulation
    // step 10 and 11
    hash_init(&hash_ctx);
    hash_update(hash_ctx, msg, *msg_len);
    hash_update(hash_ctx, salt, MIRATH_PARAM_SALT_BYTES);
    hash_update(hash_ctx, hash1, 2 * MIRATH_SECURITY_BYTES);

    for (uint32_t e = 0; e < MIRATH_PARAM_TAU; e++) {
        ff_mu_t base_alpha[MIRATH_PARAM_RHO];

        emulateparty_mu(base_alpha, i_star[e], share_S[e], share_C[e], share_v[e], gamma, H, y, mid_alpha[e]);

        hash_update(hash_ctx, base_alpha, sizeof(ff_mu_t) * (MIRATH_PARAM_M * (MIRATH_PARAM_N - MIRATH_PARAM_R)));
        hash_update(hash_ctx, mid_alpha[e], sizeof(ff_mu_t) * MIRATH_PARAM_RHO);
    }

    // Phase 3: Verification
    // step 12
    hash_finalize(hash2_partial, hash_ctx);

    // step 12
    hash_init(&hash_ctx);
    hash_update(hash_ctx, &ctr, sizeof(ctr));
    hash_update(hash_ctx, hash2_partial, 2 * MIRATH_SECURITY_BYTES);
    hash_finalize(hash2_computed, hash_ctx);

    if (hash_equal(hash2, hash2_computed)) {
        return 0;
    }

    return -1;
}