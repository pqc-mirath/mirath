#include "KeccakHash.h"
#include <string.h>
#include <stdlib.h>
#include "hash.h"

void hash_init(hash_ctx_t *ctx)
{
    *ctx = (hash_ctx_t)malloc(sizeof(Keccak_HashInstance));
#if MIRATH_SECURITY_BYTES == 16
    Keccak_HashInitialize_SHA3_256(*ctx);
#elif MIRATH_SECURITY_BYTES == 24
    Keccak_HashInitialize_SHA3_384(*ctx);
#elif MIRATH_SECURITY_BYTES == 32
	Keccak_HashInitialize_SHA3_512(*ctx);
#else
#error "HASH_SIZE not implemented!"
#endif
}

/// length in byte
void hash_update(hash_ctx_t ctx, const uint8_t *data, size_t length)
{
    Keccak_HashUpdate(ctx, data, length*8);
}

void hash_finalize(hash_t hash, hash_ctx_t ctx)
{
    Keccak_HashFinal(ctx, hash);
    free(ctx);
}

int hash_equal(hash_t hash1, hash_t hash2)
{
    return memcmp(hash1, hash2, 2 * MIRATH_SECURITY_BYTES) == 0;
}

void hash_tree_digest(hash_t hash, const uint8_t *domain_separator, const uint8_t *salt, const uint8_t *e, const uint8_t *node) {
    hash_ctx_t ctx;

    hash_init(&ctx);
    hash_update(ctx, domain_separator, sizeof(uint8_t));
    hash_update(ctx, salt, 2 * MIRATH_SECURITY_BYTES);
    hash_update(ctx, e, sizeof(uint8_t));
    hash_update(ctx, node, MIRATH_SECURITY_BYTES);
    hash_finalize(hash, ctx);
}
