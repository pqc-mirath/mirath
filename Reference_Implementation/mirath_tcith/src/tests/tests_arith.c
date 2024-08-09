#include <stdio.h>
#include <assert.h>
#include "ff.h"
#include "ff_mu.h"
#include "prng.h"
#include "vector_ff_mu.h"
#include "matrix_ff_mu.h"

#define NTESTS 100

void tests_ff(mirath_prng_t prng) {

    printf("Testing mirath_ff_mult(): ");
    for (uint32_t i = 0; i < NTESTS; i++) {
        ff_t a;

        mirath_prng(&prng, &a, sizeof(a));

        // NOTE: use the correct finite field
        a = a & 0x0F;

        assert(((mirath_ff_product(a, 1) == a) ? 0 : 1) == 0);
        assert(((mirath_ff_product(1, a) == a) ? 0 : 1) == 0);

        assert(((mirath_ff_product(a, 0) == 0) ? 0 : 1) == 0);
        assert(((mirath_ff_product(0, a) == 0) ? 0 : 1) == 0);
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
        b = mirath_ff_product(a, b);
        assert(((b == 1) ? 0 : 1) == 0);
    }
    printf("OK!\n");
}

void tests_ff_mu(mirath_prng_t prng) {

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

int test_matrix_ff_mu(mirath_prng_t prng) {

    printf("Testing mirath_matrix_ff_mu_add_multiple_ff(): ");
    for (uint32_t i = 0; i < NTESTS; i++) {
        // TODO: use the correct parameters
        ff_mu_t A[MIRATH_PARAM_TAU * MIRATH_PARAM_RHO] = {0};
        ff_t B[mirath_matrix_ff_bytes_size(MIRATH_PARAM_TAU, MIRATH_PARAM_RHO)];
        ff_mu_t scalar;

        mirath_prng(&prng, B, sizeof(B));
        mirath_prng(&prng, &scalar, sizeof(scalar));

        mirath_matrix_ff_mu_add_multiple_ff(A, scalar, B, MIRATH_PARAM_TAU * MIRATH_PARAM_RHO, 1);

        ff_mu_t X[MIRATH_PARAM_TAU * MIRATH_PARAM_RHO] = {0};
        ff_mu_t Y[MIRATH_PARAM_TAU * MIRATH_PARAM_RHO];

        mirath_matrix_map_ff_to_ff_mu(Y, B, MIRATH_PARAM_TAU * MIRATH_PARAM_RHO, 1);

        mirath_matrix_ff_mu_add_multiple_2(X, scalar, Y, MIRATH_PARAM_TAU * MIRATH_PARAM_RHO, 1);

        for (uint32_t j = 0; j < MIRATH_PARAM_TAU * MIRATH_PARAM_RHO; j++) {
            if (A[j] != X[j]) {
                printf("FAIL (A != X)\n");
                return -1;
            }
        }
    }
    printf("OK!\n");

    printf("Testing mirath_matrix_ff_mu_product_ff1mu(): ");
    for (uint32_t i = 0; i < NTESTS; i++) {
        // TODO: use the correct parameters
        ff_mu_t R[MIRATH_PARAM_TAU * MIRATH_PARAM_RHO] = {0};
        ff_t A[mirath_matrix_ff_bytes_size(MIRATH_PARAM_TAU, MIRATH_PARAM_TAU)];
        ff_mu_t B[MIRATH_PARAM_TAU];

        mirath_prng(&prng, A, sizeof(A));
        mirath_prng(&prng, B, sizeof(B));

        mirath_matrix_ff_mu_product_ff1mu(R, A, B, MIRATH_PARAM_TAU, MIRATH_PARAM_RHO * 2, MIRATH_PARAM_RHO);

        ff_mu_t X[MIRATH_PARAM_TAU * MIRATH_PARAM_RHO] = {0};
        ff_mu_t Y[MIRATH_PARAM_TAU * MIRATH_PARAM_TAU];

        mirath_matrix_map_ff_to_ff_mu(Y, A, MIRATH_PARAM_TAU, MIRATH_PARAM_TAU);

        mirath_matrix_ff_mu_product(X, Y, B, MIRATH_PARAM_TAU, MIRATH_PARAM_RHO * 2, MIRATH_PARAM_RHO);

        for (uint32_t j = 0; j < MIRATH_PARAM_TAU * MIRATH_PARAM_RHO; j++) {
            if (R[j] != X[j]) {
                printf("FAIL (R != X)\n");
                return -1;
            }
        }
    }
    printf("OK!\n");

    return 0;
}

int main(int argc, char **argv) {
    mirath_prng_t prng;
    mirath_prng_init(&prng, NULL, NULL);

    tests_ff(prng);

    tests_ff_mu(prng);

    if (test_matrix_ff_mu(prng) != 0)
        return -1;

    return 0;
}