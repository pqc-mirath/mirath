/**
 * @file mirath_tcith.c
 * @brief Implementation of TCitH related functions
 */

#include "prng.h"
#include "random.h"
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
