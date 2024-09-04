
#ifndef PRNG_H
#define PRNG_H

#include <stdint.h>
#include "fips202.h"
#include "../mirath_parameters.h"
#include "hash.h"

typedef keccak_state mirath_prng_t;

/* Initialize 'prng' from 'salt' and 'seed'.
 * If 'salt == NULL' then 'salt' is ignored.
 * If 'seed == NULL' then 'seed' is ignored. */
void mirath_prng_init(mirath_prng_t *prng, const hash_t salt, const seed_t seed, const uint32_t seed_size_bytes);

/* Write 'length' pseudorandom bytes over 'target',
 * update the internal state of 'prng'. */
void mirath_prng(mirath_prng_t *prng, void *target, size_t length);

#endif
