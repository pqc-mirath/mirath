#ifndef VECTOR_FF_MU_H
#define VECTOR_FF_MU_H

#include <stdint.h>
#include "prng.h"
#include "matrix_ff_arith.h"
#include "ff_mu.h"
#include "vector_ff_arith.h"

/// vector1 = vector2 + vector3
static inline void mirath_vector_ff_mu_add(ff_mu_t *vector1, const ff_mu_t *vector2, const ff_mu_t *vector3, const uint32_t d) {
    for (uint32_t i = 0; i < d; i++) {
        vector1[i] = vector2[i] ^ vector3[i];
    }
}

/// vector1 = vector2 * value
/// NOTE: vector2 is in base field
static inline void mirath_vector_ff_mu_mult_multiple_ff(ff_mu_t *vector1, const ff_mu_t scalar, const ff_t *vector2, const uint32_t d) {
    for (uint32_t i = 0; i < d; i++) {
        const ff_t t = mirath_matrix_ff_get_entry(vector2, 1, i, 0);
        vector1[i] = mirath_ff_mu_mult(scalar, mirath_map_ff_to_ff_mu[t]);
    }
}

/// vector1 = vector2 + scalar*vector3
/// NOTE: vector3 is in base field
static inline void mirath_vector_ff_mu_add_multiple_ff(ff_mu_t *vector1, const ff_mu_t *vector2, const ff_mu_t scalar, const ff_t *vector3, const uint32_t d) {
    for (uint32_t i = 0; i < d; i++) {
        const ff_t t = mirath_matrix_ff_get_entry(vector3, 1, i, 0);
        vector1[i] = vector2[i] ^ mirath_ff_mu_mult(scalar, mirath_map_ff_to_ff_mu[t]);
    }
}

/// vector1 = vector2 * value
static inline void mirath_vector_ff_mu_mult_multiple(ff_mu_t *vector1, const ff_mu_t *vector2, const ff_mu_t value, const uint32_t d) {
    for (uint32_t i = 0; i < d; i++) {
        vector1[i] = mirath_ff_mu_mult(vector2[i], value);
    }
}

///  vector1 = vector2 + vector3 * value
static inline void mirath_vector_ff_mu_add_multiple(ff_mu_t *vector1, const ff_mu_t *vector2, const ff_mu_t scalar, const ff_mu_t *vector3, const uint32_t d) {
    for (uint32_t i = 0; i < d; i++) {
        vector1[i] = vector2[i] ^ mirath_ff_mu_mult(scalar, vector3[i]);
    }
}

///
static inline ff_mu_t mirath_vector_ff_mu_mul_acc(const ff_mu_t *a, const ff_mu_t *b, const uint32_t n) {
    ff_mu_t acc = 0;
    for (uint32_t i = 0; i < n; ++i) {
        acc ^= mirath_ff_mu_mult(a[i], b[i]);
    }

    return acc;
}

#endif
