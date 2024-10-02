#ifndef ARITH_FF_H
#define ARITH_FF_H

#include <stdint.h>
#include "data_type.h"

/// NOTE: assumes a mod 2
/// \param a input element
/// \return a**{-1}
static inline ff_t mirath_ff_inv(const ff_t a) {
   return a;
}

/// NOTE: assumes that a and b are % 2.
/// \param a input element
/// \param b input element
/// \return a*b % 2
static inline ff_t mirath_ff_product(const ff_t a, const ff_t b) {
   return a & b;
}
#endif
