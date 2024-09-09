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

int mirath_sign(uint8_t *sig_msg, uint8_t *msg, size_t msg_len, uint8_t *sk) {
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

    domain_separator = DOMAIN_SEPARATOR_COMMITMENT;
    hash_t h_com;
    hash_ctx_t hash_commits;
    hash_init(&hash_commits);
    hash_update(hash_ctx, &domain_separator, sizeof(uint8_t));

    for (uint32_t e = 0; e < MIRATH_PARAM_TAU; e++) {
        const uint16_t N = e < MIRATH_PARAM_TAU_1 ? MIRATH_PARAM_N_1 : MIRATH_PARAM_N_2;

        build_sharing_N(aux[e], rnd_S[e], rnd_C[e], rnd_v[e], commits, seeds, S, C, salt, e);

        hash_update(hash_commits, (uint8_t *)commits[e], sizeof(mirath_tcith_commit_t) * N);
    }

    hash_finalize(h_com, hash_commits);

    hash_update(hash_ctx, h_com, 2 * MIRATH_SECURITY_BYTES);

    const uint32_t aux_n_bytes = mirath_matrix_ff_bytes_size(
            MIRATH_PARAM_M * MIRATH_PARAM_R + MIRATH_PARAM_R * (MIRATH_PARAM_N - MIRATH_PARAM_R), 1);
    for (uint32_t e = 0; e < MIRATH_PARAM_TAU; e++) {
        hash_update(hash_ctx, aux[e], aux_n_bytes);
    }

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
    uint64_t ctr = 0;
    uint64_t path_length;
 retry:
     // Initialize hash_ctx
     hash_init(&hash_ctx);
     domain_separator = DOMAIN_SEPARATOR_HASH2;
     hash_update(hash_ctx, &domain_separator, sizeof(uint8_t));
     hash_update(hash_ctx, (uint8_t *)&ctr, sizeof(uint64_t));
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

     path_length = mirath_tree_get_sibling_path(path, &tree, psi_i_star, MIRATH_PARAM_TAU);
     if (path_length > MIRATH_PARAM_T_OPEN) {
         ctr += 1;
         memset(path, 0, sizeof(uint8_t) * (MIRATH_PARAM_TREE_LEAVES * MIRATH_SECURITY_BYTES));
         goto retry;
     }

    // Phase 5: Signature
    // step 21
    unparse_signature(sig_msg, salt, ctr, hash2, path, path_length, commits, aux, mid_alpha, i_star);

    return 0;
}

int mirath_verify(uint8_t *msg, size_t *msg_len, uint8_t *sig_msg, size_t sig_msg_len, uint8_t *pk) {
    seed_t seed_pk;
    uint8_t salt[MIRATH_PARAM_SALT_BYTES];

    mirath_tree_t tree;
    mirath_tree_leaves_t seeds;

    mirath_tcith_commit_t *commits[MIRATH_PARAM_TAU];
    mirath_tcith_commit_1_t commits_1 = {0};
    mirath_tcith_commit_2_t commits_2 = {0};
    mirath_tcith_commit_set_as_a_grid_list(commits, &commits_1, &commits_2);

    hash_t hash1;
    hash_t hash2_partial;
    hash_t hash2;
    hash_t hash2_computed;

    hash_ctx_t hash_ctx;
    mirath_prng_t prng;

    uint8_t path[MIRATH_PARAM_TREE_LEAVES * MIRATH_SECURITY_BYTES] = {0};

    uint8_t domain_separator;

    ff_mu_t share_S[MIRATH_PARAM_TAU][MIRATH_PARAM_M * MIRATH_PARAM_N];
    ff_mu_t share_C[MIRATH_PARAM_TAU][MIRATH_PARAM_R * (MIRATH_PARAM_N - MIRATH_PARAM_R)];
    ff_mu_t share_v[MIRATH_PARAM_TAU][MIRATH_PARAM_RHO];
    ff_mu_t gamma[MIRATH_PARAM_RHO * (MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K)];
    ff_t H[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, MIRATH_PARAM_K)];
    ff_t y[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1)];

    ff_t aux[MIRATH_PARAM_TAU][mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_R + MIRATH_PARAM_R * (MIRATH_PARAM_N - MIRATH_PARAM_R), 1)];
    ff_mu_t mid_alpha[MIRATH_PARAM_TAU][MIRATH_PARAM_RHO];

    mirath_tcith_commit_t commits_i_star[MIRATH_PARAM_TAU];

    mirath_tcith_challenge_t i_star;
    uint16_t psi_i_star[MIRATH_PARAM_TAU];

    uint64_t ctr;

    // step 1
    parse_public_key(seed_pk, y, pk);

    // Phase 0: Parsing and expansion
    // step 2
    parse_signature(salt, &ctr, hash2, path, commits_i_star, aux, mid_alpha, sig_msg);

    // step 3
    mirath_prng_init(&prng, NULL, seed_pk, MIRATH_SECURITY_BYTES);
    mirath_matrix_ff_init_random(H, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, MIRATH_PARAM_K, &prng);

    // step 4
    mirath_tcith_compute_challenge_2(i_star, hash2, salt);
    // Next we map the challenges to the leaves position (GGM Tree optimization)
    for(size_t e = 0; e < MIRATH_PARAM_TAU; e++){
        size_t i = i_star[e];
        psi_i_star[e] = (uint16_t)mirath_tcith_psi(i, e);
    }

    size_t path_length = 0;
    for(uint32_t i = 0; i < MIRATH_PARAM_T_OPEN; i++) {
        const uint8_t zero[MIRATH_SECURITY_BYTES] = {0};
        if (memcmp(zero, &path[i * MIRATH_SECURITY_BYTES], MIRATH_SECURITY_BYTES) == 0) { continue; }
        path_length += 1;
    }

    // step 5, step 6, and 7
    mirath_tree_init(&tree);
    mirath_tree_get_seeds_from_path(&tree, psi_i_star, MIRATH_PARAM_TAU, path, path_length, salt, 0);
    mirath_tree_get_leaves(seeds, &tree);

    domain_separator = DOMAIN_SEPARATOR_HASH1;
    hash_init(&hash_ctx);
    hash_update(hash_ctx, &domain_separator, sizeof(uint8_t));
    hash_update(hash_ctx, salt, MIRATH_PARAM_SALT_BYTES);

    domain_separator = DOMAIN_SEPARATOR_COMMITMENT;
    hash_t h_com;
    hash_ctx_t hash_commits;
    hash_init(&hash_commits);
    hash_update(hash_ctx, &domain_separator, sizeof(uint8_t));

    for (uint32_t e = 0; e < MIRATH_PARAM_TAU; e++) {
        const uint16_t N = e < MIRATH_PARAM_TAU_1 ? MIRATH_PARAM_N_1 : MIRATH_PARAM_N_2;

        memset(share_S[e], 0, sizeof(ff_mu_t) * MIRATH_PARAM_M * MIRATH_PARAM_N);
        memset(share_C[e], 0, sizeof(ff_mu_t) * MIRATH_PARAM_R * (MIRATH_PARAM_N - MIRATH_PARAM_R));
        memset(share_v[e], 0, sizeof(ff_mu_t) * MIRATH_PARAM_RHO);

        compute_share(share_S[e], share_C[e], share_v[e], commits, i_star[e], seeds, e, salt, commits_i_star, aux[e]);

        hash_update(hash_commits, (uint8_t *)commits[e], sizeof(mirath_tcith_commit_t) * N);
    }

    hash_finalize(h_com, hash_commits);

    hash_update(hash_ctx, h_com, 2 * MIRATH_SECURITY_BYTES);

    const uint32_t aux_n_bytes = mirath_matrix_ff_bytes_size(
            MIRATH_PARAM_M * MIRATH_PARAM_R + MIRATH_PARAM_R * (MIRATH_PARAM_N - MIRATH_PARAM_R), 1);
    for (uint32_t e = 0; e < MIRATH_PARAM_TAU; e++) {
        hash_update(hash_ctx, aux[e], aux_n_bytes);
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
    hash_update(hash_ctx, (uint8_t *)ctr, sizeof(uint64_t));
    hash_update(hash_ctx, hash2_partial, 2 * MIRATH_SECURITY_BYTES);
    hash_finalize(hash2_computed, hash_ctx);

    if (hash_equal(hash2, hash2_computed)) {
        return 0;
    }

    return -1;
}