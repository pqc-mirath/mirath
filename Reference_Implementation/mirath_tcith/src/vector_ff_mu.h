#ifndef VECTOR_H
#define VECTOR_H

#include <stdint.h>
#include "common/prng.h"

#include "arith/vector_arith_ff.h"
#include "arith/vector_arith_ff_mu.h"

/// vector = 0
void mirath_vec_ff_mu_init(ff_mu_t *vector, const uint32_t d);
/// vector = 0
void mirath_vec_ff_mu_init_zero(ff_mu_t *vector, const uint32_t d);
/// returns an extension field entry
ff_mu_t mirath_vec_ff_mu_get_entry(const ff_mu_t *vector, const uint32_t i);
/// sets an entry
void mirath_vec_ff_mu_set_entry(ff_mu_t *vector, const uint32_t i, const ff_mu_t scalar);
/// inits the vector with random data
void mirath_vec_ff_mu_init_random(ff_mu_t *vector, const uint32_t d, prng_t *prng);
void mirath_vec_ff_mu_copy(ff_mu_t *p1, const ff_mu_t *p2, const uint32_t d);
void mirath_vec_ff_mu_negate(ff_mu_t *vector, const uint32_t d);
/// vector1 = vector2 + vector3
void mirath_vec_ff_mu_add(ff_mu_t *vector1, const ff_mu_t *vector2, const ff_mu_t *vector3, const uint32_t d);
/// vector1 = vector2*scalar
void mirath_vec_ff_mu_scalar_mult(ff_mu_t *vector1, ff_mu_t *vector2, const ff_mu_t value, const uint32_t d);
/// vector1 = vector2 + scalar * vector3
void mirath_vec_ff_mu_add_multiple(ff_mu_t *vector1, const ff_mu_t *vector2, const ff_mu_t scalar, const ff_mu_t *vector3, const uint32_t d);
/// vector1 = vector2 - vector3
void mirath_vec_ff_mu_sub(ff_mu_t *vector1, const ff_mu_t *vector2, const ff_mu_t *vector3, const uint32_t d);
/// vector1 = vector2 - scalar * vector3
void mirath_vec_ff_mu_sub_multiple(ff_mu_t *vector1, const ff_mu_t *vector2, const ff_mu_t scalar, const ff_mu_t *vector3, const uint32_t d);
/// \return sum_i=0,...,n (a_i * b_i)
ff_mu_t mirath_vec_ff_mu_mul_acc(const ff_mu_t *a, const ff_mu_t *b, const uint32_t n);
#endif
