#include <string.h>
#include "prng.h"

// NOTE: this is only a skeleton

// Agregar el tamaño de la semilla 2 * lambda como parametro de entrada
void mirath_prng_init(mirath_prng_t *prng, const hash_t salt, const seed_t seed, const uint32_t seed_size_bytes) {
    uint8_t input[2 * MIRATH_SECURITY_BYTES + seed_size_bytes];
    memset(input, 0, 2 * MIRATH_SECURITY_BYTES + seed_size_bytes);
	
	// NOTE: generalize this in the case of CAT 3 and 5
    mirath_fips202_ref_shake128_init(prng);

    /* Set 'buffer = salt | seed'. */
    if (salt != NULL) {
        memcpy(input, salt, 2 * MIRATH_SECURITY_BYTES);
    }

    if (seed != NULL){
        memcpy(input + 2 * MIRATH_SECURITY_BYTES, seed, seed_size_bytes);
    }

    mirath_fips202_ref_shake128_absorb(prng, input, 2 * MIRATH_SECURITY_BYTES + seed_size_bytes);
    shake128_finalize(prng);
}

void mirath_prng(mirath_prng_t *prng, void *target, size_t length) {
    mirath_fips202_ref_shake128_squeeze((uint8_t*) target, length, prng);
}
