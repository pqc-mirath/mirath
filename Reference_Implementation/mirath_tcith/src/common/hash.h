#ifndef HASH_H
#define HASH_H

#include <stdint.h>
#include "KeccakHash.h"
#include "../mirath_parameters.h"

typedef uint8_t seed_t[MIRATH_SECURITY_BYTES];
typedef uint8_t hash_t[2 * MIRATH_SECURITY_BYTES];

/* Type for a hashing context. */
typedef void *hash_ctx_t;

/* Initialize the hashing context 'ctx'. */
void hash_init(hash_ctx_t *ctx);

/* Update the hashing context 'ctx' with 'data'. */
void hash_update(hash_ctx_t ctx, const uint8_t *data, size_t length);

/* Finalize the hashing, write the digest over 'hash'. */
void hash_finalize(hash_t hash, hash_ctx_t ctx);

/* Return 'True' if 'hash1' and 'hash2' are equal, and 'False' otherwise. */
int hash_equal(hash_t hash1, hash_t hash2);

void hash_tree_digest(hash_t hash, const uint8_t *domain_separator, const uint8_t *salt, const uint8_t *e, const uint8_t* node);

#endif