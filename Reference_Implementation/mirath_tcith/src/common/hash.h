#ifndef HASH_H
#define HASH_H

#include <stdint.h>
#include "KeccakHash.h"
#include "../params.h"

typedef uint8_t seed_t[SEED_SIZE];
typedef uint8_t hash_t[HASH_SIZE];

/* Type for a hashing context. */
typedef void *hash_ctx_t;

/* Initialize the hashing context 'ctx'. */
void hash_init(hash_ctx_t *ctx);

/* Update the hashing context 'ctx' with 'data'. */
void hash_update(hash_ctx_t ctx, const uint8_t *data, size_t length);

/* Finalize the hashing, write the digest over 'hash'. */
void hash_finalize(hash_ctx_t ctx, hash_t hash);

/* Return 'True' if 'hash1' and 'hash2' are equal, and 'False' otherwise. */
int hash_equal(hash_t hash1, hash_t hash2);

/* Write over 'hash' the hash digest of 'salt', 'l', 'i', 'seed'. */
void hash_digest0(hash_t hash, const hash_t salt, uint32_t l, uint32_t i, const seed_t seed);
void hash_digest0x4(hash_t hash, const hash_t salt, uint32_t l, uint32_t i, const seed_t seed);
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#endif