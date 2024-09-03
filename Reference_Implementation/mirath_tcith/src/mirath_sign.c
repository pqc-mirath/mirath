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
    mirath_tree_leaves_t seeds;
    mirath_tree_t tree;

    mirath_tcith_commit_t *commits[MIRATH_PARAM_TAU];
    mirath_tcith_commit_1_t commits_1 = {0};
    mirath_tcith_commit_2_t commits_2 = {0};
    mirath_tcith_commit_set_as_a_grid_list(commits, &commits_1, &commits_2);

    hash_t hash1;
    hash_t hash2_partial;
    hash_t hash2;

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

    // step 5, step 6 and 7
    mirath_tree_init(&tree);
    memcpy(tree.nodes[0], rseed, MIRATH_SECURITY_BYTES); // root seed
    tree.nonempty[0] = 1;
    mirath_tree_prg(&tree, salt, 0);
    mirath_tree_get_leaves(seeds, &tree);

    hash_init(&hash_ctx);
    hash_update(hash_ctx, salt, MIRATH_PARAM_SALT_BYTES);

    build_sharing_N(aux, rnd_S, rnd_C, rnd_v, hash_ctx, seeds, S, C, salt);

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

        emulateparty_mu(base_alpha, i_star[e], share_S[e], share_C[e], share_v[e], gamma, H, y, mid_alpha[e]);

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