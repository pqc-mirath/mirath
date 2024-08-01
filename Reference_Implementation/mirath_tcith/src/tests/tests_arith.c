#include <stdio.h>
#include <assert.h>
#include "../arith/ff.h"
#include "../arith/ff_mu.h"
#include "../common/prng.h"

#define NTESTS 100

void tests_ff() {
    mirath_prng_t prng;
    mirath_prng_init(&prng, NULL, NULL);

    printf("Testing mirath_ff_mult(): ");
    for (uint32_t i = 0; i < NTESTS; i++) {
        ff_t a;

        mirath_prng(&prng, &a, sizeof(a));

        // NOTE: use the correct finite field
        a = a & 0x0F;

        assert(((mirath_ff_mult(a, 1) == a) ? 0 : 1) == 0);
        assert(((mirath_ff_mult(1, a) == a) ? 0 : 1) == 0);

        assert(((mirath_ff_mult(a, 0) == 0) ? 0 : 1) == 0);
        assert(((mirath_ff_mult(0, a) == 0) ? 0 : 1) == 0);
    }
    printf("OK!\n");

    printf("Testing mirath_ff_inv(): ");
    for (uint32_t i = 0; i < NTESTS; i++) {
        ff_t a, b;

        mirath_prng(&prng, &a, sizeof(a));
        mirath_prng(&prng, &b, sizeof(b));
        // NOTE: use the correct finite field
        // in order to avoid zero
        a = (a | b) & 0x0F;

        b = mirath_ff_inv(a);
        b = mirath_ff_mult(a, b);
        assert(((b == 1) ? 0 : 1) == 0);
    }
    printf("OK!\n");
}

void tests_ff_mu() {
    mirath_prng_t prng;
    mirath_prng_init(&prng, NULL, NULL);

    printf("Testing mirath_ff_mu_mult(): ");
    for (uint32_t i = 0; i < NTESTS; i++) {
        ff_mu_t a;

        mirath_prng(&prng, &a, sizeof(a));

        assert(((mirath_ff_mu_mult(a, 1) == a) ? 0 : 1) == 0);
        assert(((mirath_ff_mu_mult(1, a) == a) ? 0 : 1) == 0);

        assert(((mirath_ff_mu_mult(a, 0) == 0) ? 0 : 1) == 0);
        assert(((mirath_ff_mu_mult(0, a) == 0) ? 0 : 1) == 0);
    }
    printf("OK!\n");

    printf("Testing mirath_ff_mu_inv(): ");
    for (uint32_t i = 0; i < NTESTS; i++) {
        ff_mu_t a, b;

        mirath_prng(&prng, &a, sizeof(a));
        mirath_prng(&prng, &b, sizeof(b));
        // in order to avoid zero
        a = (a | b) & 0x0F;

        b = mirath_ff_mu_inv(a);
        b = mirath_ff_mu_mult(a, b);
        assert(((b == 1) ? 0 : 1) == 0);
    }
    printf("OK!\n");
}

int main(int argc, char **argv) {
    tests_ff();

    tests_ff_mu();

    return 0;
}