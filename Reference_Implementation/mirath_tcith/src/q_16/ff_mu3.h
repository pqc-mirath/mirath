#ifndef ARITH_FF_MU3_H
#define ARITH_FF_MU3_H

#include <stdint.h>
#include <stdint.h>

#include "ff.h"

typedef uint16_t ff_mu3_t;

/// \return a+b
static inline ff_mu3_t mirath_ff_mu_add(const ff_mu3_t a,
                                        const ff_mu3_t b) {
    return a^b;
}

/// \return a*b
static inline ff_mu3_t mirath_ff_mu_mult(const ff_mu3_t a,
                                         const ff_mu3_t b) {
    ff_mu3_t r;

    const ff_t a0 = a&0xF, a1 = (a>>4)&0XF, a2 = (a>>8)&0xF;
    const ff_t b0 = b&0xF, b1 = (b>>4)&0XF, b2 = (b>>8)&0xF;

    const ff_t p0 = mirath_ff_product(a0, b0);
    const ff_t p1 = mirath_ff_product(a1, b1);
    const ff_t p2 = mirath_ff_product(a2, b2);

    const ff_t a01 = mirath_ff_add(a0, a1);
    const ff_t a12 = mirath_ff_add(a1, a2);
    const ff_t b01 = mirath_ff_add(b0, b1);
    const ff_t b12 = mirath_ff_add(b1, b2);
    const ff_t p01 = mirath_ff_product(a01, b01);
    const ff_t p12 = mirath_ff_product(a12, b12);

    const ff_t a012 = mirath_ff_add(a01, a2);
    const ff_t b012 = mirath_ff_add(b01, b2);
    const ff_t p012 = mirath_ff_product(a012, b012);

    // compute lowest limb
    r = mirath_ff_add(p1, p2);
    r = mirath_ff_add(r, p12);
    r = mirath_ff_add(r, p0);

    r^= p2 << 4;
    r^= p01 << 4;
    r^= p0 << 4;
    r^= p1 << 4;

    r^= p012 << 8;
    r^= p01 << 8;
    r^= p12 << 8;
    return r;
}


#endif
