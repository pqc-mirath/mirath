#ifndef VECTOR_ARITH_FF_MU_H
#define VECTOR_ARITH_FF_MU_H

#include <stdint.h>
#include <stdlib.h>
#include "ff_mu.h"

#define mirath_vec_ff_mu_bytes_size(d) (sizeof(ff_mu_t) * d)

/// arg1 = arg2 + arg3
static inline void mirath_vec_ff_mu_add_arith(ff_mu_t *arg1, const ff_mu_t *arg2, const ff_mu_t *arg3, const uint32_t d) {
	for (uint32_t i = 0; i < d; i++) {
		arg1[i] = arg2[i] ^ arg3[i];
	}
}

///  arg1 = arg2 + arg3 * value
static inline void mirath_vec_ff_mu_add_multiple_arith(ff_mu_t *arg1, const ff_mu_t *arg2, const ff_mu_t scalar, const ff_mu_t *arg3, const uint32_t d) {
	for (uint32_t i = 0; i < d; i++) {
		arg1[i] = arg2[i] ^ mirath_ff_mu_mult(scalar, arg3[i]);
	}
}

/// arg1 = arg2*value
static inline void mirath_vec_ff_mu_scalar_mult_arith(ff_mu_t *arg1, const ff_mu_t *arg2, const ff_mu_t value, const uint32_t d) {
	for (uint32_t i = 0; i < d; i++) {
		arg1[i] = mirath_ff_mu_mult(arg2[i], value);
	}
}

/// \return arg1(value)
static inline ff_mu_t mirath_vec_ff_mu_eval_arith(const ff_mu_t *arg1, const ff_mu_t value, const uint32_t d) {
    ff_mu_t ret = arg1[0]; // constant value
    ff_mu_t v = value;
	for (uint32_t i = 1; i < d; i++) {
		ret ^= mirath_ff_mu_mult(v, arg1[i]);
		v = mirath_ff_mu_mult(v, value);
	}

	return ret;
}

/// return arg_1(scalar_1), ..., arg_n(scalar_1)
static inline void mirath_vec_ff_mu_eval_parallel_arith(ff_mu_t *out,
                                                const ff_mu_t *arg1, const ff_mu_t scalar,
                                                const uint32_t d, const uint32_t n) {
    for (uint32_t i = 0; i < n; ++i) {
        out[i] = mirath_vec_ff_mu_eval_arith(arg1 + d*i, scalar, d);
    }
}

/// \return sum_i=0,...,n (a_i * b_i)
static inline ff_mu_t mirath_vec_ff_mu_mul_acc_arith(const ff_mu_t *a, const ff_mu_t *b, const uint32_t n) {
    ff_mu_t acc = 0;
    for (uint32_t i = 0; i < n; ++i) {
        acc ^= mirath_ff_mu_mult(a[i], b[i]);
    }

    return acc;
}

#endif
