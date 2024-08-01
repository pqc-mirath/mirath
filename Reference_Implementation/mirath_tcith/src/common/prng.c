#include <string.h>
#include "prng.h"

// NOTE: this is only a skeleton

void mirath_prng_init(mirath_prng_t *prng, const hash_t salt, const seed_t seed) {
    uint8_t input[HASH_SIZE + SEED_SIZE];
    memset(input, 0, HASH_SIZE + SEED_SIZE);
	
	// NOTE: generalize this in the case of CAT 3 and 5
    mirath_fips202_ref_shake128_init(prng);

    /* Set 'buffer = salt | seed'. */
    if (salt != NULL) {
        memcpy(input, salt, HASH_SIZE);
    }

    if (seed != NULL){
        memcpy(input + HASH_SIZE, seed, SEED_SIZE);
    }

    mirath_fips202_ref_shake128_absorb(prng, input, HASH_SIZE + SEED_SIZE);
    shake128_finalize(prng);
}

void mirath_prng(mirath_prng_t *prng, void *target, size_t length) {
    mirath_fips202_ref_shake128_squeeze((uint8_t*) target, length, prng);
}
