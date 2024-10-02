#ifndef ARITH_FF_MU_H
#define ARITH_FF_MU_H

#include <stdint.h>
#include "../data_type.h"

/// AES modulus
#define MODULUS 0x1B

/// \return a+b
static inline ff_mu_t mirath_ff_mu_add(const ff_mu_t a, const ff_mu_t b) {
    return a^b;
}

/// \return a*b
static inline ff_mu_t mirath_ff_mu_mult(const ff_mu_t a, const ff_mu_t b) {
    ff_mu_t result = -(a & 1) & b;
    ff_mu_t tmp = b;
    for(int i=1 ; i<8 ; i++) {
        tmp = (tmp << 1) ^ (-(tmp >> 7) & MODULUS);
        result = result ^ (-(a >> i & 1) & tmp);
    }
    return result;
}

/// \return a^-1
static inline ff_mu_t mirath_ff_mu_inv(const ff_mu_t a) {
    ff_mu_t result = a;
    for(int i=0 ; i<6 ; i++) {
        result = mirath_ff_mu_mult(result, result);
        result = mirath_ff_mu_mult(result, a);
    }
    result = mirath_ff_mu_mult(result, result);
    return result;
}

#endif
