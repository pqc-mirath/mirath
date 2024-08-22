#ifndef VECTOR_FF_MU_H
#define VECTOR_FF_MU_H

#include <stdint.h>
#include "prng.h"
#include "matrix_ff_arith.h"
#include "ff_mu.h"
#include "vector_ff_arith.h"

/**
 * \fn static inline void mirath_vector_ff_mu_add(ff_mu_t *vector1, const ff_mu_t *vector2, const ff_mu_t *vector3, const uint32_t ncols)
 * \brief vector1 = vector2 + vector3
 *
 * \param[out] vector1 Vector over ff_mu
 * \param[in] vector2 Vector over ff_mu
 * \param[in] vector3 Vector over ff_mu
 * \param[in] ncols number of columns
 */
static inline void mirath_vector_ff_mu_add(ff_mu_t *vector1, const ff_mu_t *vector2, const ff_mu_t *vector3, const uint32_t ncols) {
    for (uint32_t i = 0; i < ncols; i++) {
        vector1[i] = vector2[i] ^ vector3[i];
    }
}

/**
 * \fn static inline void mirath_vector_ff_mu_add_ff(ff_mu_t *vector1, const ff_mu_t *vector2, const ff_t *vector3, const uint32_t ncols)
 * \brief vector1 = vector2 + vector3
 *
 * \param[out] vector1 Vector over ff_mu
 * \param[in] vector2 Vector over ff_mu
 * \param[in] vector3 Vector over ff
 * \param[in] ncols number of columns
 */
static inline void mirath_vector_ff_mu_add_ff(ff_mu_t *vector1, const ff_mu_t *vector2, const ff_t *vector3, const uint32_t ncols) {
    for (uint32_t i = 0; i < ncols; i++) {
        const ff_t t = mirath_matrix_ff_get_entry(vector3, 1, i, 0);
        vector1[i] = vector2[i] ^ mirath_map_ff_to_ff_mu[t];
    }
}

/**
 * \fn static inline void mirath_vector_ff_mu_mult_multiple_ff(ff_mu_t *vector1, const ff_mu_t scalar, const ff_t *vector2, const uint32_t ncols)
 * \brief vector1 = vector2 * scalar
 *
 * \param[out] vector1 Vector over ff_mu
 * \param[in] scalar Scalar over ff_mu
 * \param[in] vector2 Vector over ff
 * \param[in] ncols number of columns
 */
static inline void mirath_vector_ff_mu_mult_multiple_ff(ff_mu_t *vector1, const ff_mu_t scalar, const ff_t *vector2, const uint32_t ncols) {
    for (uint32_t i = 0; i < ncols; i++) {
        const ff_t t = mirath_matrix_ff_get_entry(vector2, 1, i, 0);
        vector1[i] = mirath_ff_mu_mult(scalar, mirath_map_ff_to_ff_mu[t]);
    }
}

/**
 * \fn static inline void mirath_vector_ff_mu_add_multiple_ff(ff_mu_t *vector1, const ff_mu_t *vector2, const ff_mu_t scalar, const ff_t *vector3, const uint32_t ncols)
 * \brief vector1 = vector2 + scalar * vector3
 *
 * \param[out] vector1 Vector over ff_mu
 * \param[in] vector2 Vector over ff_mu
 * \param[in] scalar Scalar over ff_mu
 * \param[in] vector3 Vector over ff
 * \param[in] ncols number of columns
 */
static inline void mirath_vector_ff_mu_add_multiple_ff(ff_mu_t *vector1, const ff_mu_t *vector2, const ff_mu_t scalar, const ff_t *vector3, const uint32_t ncols) {
    for (uint32_t i = 0; i < ncols; i++) {
        const ff_t t = mirath_matrix_ff_get_entry(vector3, 1, i, 0);
        vector1[i] = vector2[i] ^ mirath_ff_mu_mult(scalar, mirath_map_ff_to_ff_mu[t]);
    }
}

/**
 * \fn static inline void mirath_vector_ff_mu_mult_multiple(ff_mu_t *vector1, const ff_mu_t scalar, const ff_mu_t *vector2, const uint32_t ncols)
 * \brief vector1 = vector2 * scalar
 *
 * \param[out] vector1 Vector over ff_mu
 * \param[in] scalar Scalar over ff_mu
 * \param[in] vector2 Vector over ff_mu
 * \param[in] ncols number of columns
 */
static inline void mirath_vector_ff_mu_mult_multiple(ff_mu_t *vector1, const ff_mu_t scalar, const ff_mu_t *vector2, const uint32_t ncols) {
    for (uint32_t i = 0; i < ncols; i++) {
        vector1[i] = mirath_ff_mu_mult(vector2[i], scalar);
    }
}

/**
 * \fn static inline void mirath_vector_ff_mu_add_multiple(ff_mu_t *vector1, const ff_mu_t scalar, const ff_mu_t *vector2, const ff_mu_t *vector3, const uint32_t ncols)
 * \brief vector1 = vector2 + vector3 * value
 *
 * \param[out] vector1 Vector over ff_mu
 * \param[in] vector2 Vector over ff_mu
 * \param[in] scalar Scalar over ff_mu
 * \param[in] vector3 Vector over ff_mu
 * \param[in] ncols number of columns
 */
static inline void mirath_vector_ff_mu_add_multiple(ff_mu_t *vector1, const ff_mu_t *vector2, const ff_mu_t scalar, const ff_mu_t *vector3, const uint32_t ncols) {
    for (uint32_t i = 0; i < ncols; i++) {
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
