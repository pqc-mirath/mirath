#include <stdio.h>
#include <assert.h>
#include "q_16/ff.h"
#include "q_16/ff_mu.h"
#include "prng.h"
#include "q_16/vector_ff_mu.h"
#include "q_16/matrix_ff_mu.h"
#include "mirath_matrix_ff.h"

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
        a = (a | b) & 0x0F;
        // in order to avoid zero
        if (a == 0)
            a = 1;

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

        mirath_prng(&prng, &a, sizeof(ff_mu_t));
        mirath_prng(&prng, &b, sizeof(ff_mu_t));
        // in order to avoid zero
        a = a | b;

        b = mirath_ff_mu_inv(a);
        b = mirath_ff_mu_mult(a, b);
        assert(((b == 1) ? 0 : 1) == 0);
    }
    printf("OK!\n");
}

int test_matrix_ff_mu(mirath_prng_t prng) {

    printf("Testing: \n");
    printf("\tmirath_matrix_ff_mu_add()\n");
    printf("\tmirath_matrix_ff_mu_add_mu1ff()\n");
    printf("Status: ");
    for (uint32_t i = 0; i < NTESTS; i++) {
        ff_mu_t zero[MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K] = {0};

        ff_mu_t A[MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K] = {0};
        ff_mu_t B[MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K];
        ff_t C[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1)];

        mirath_prng(&prng, B, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K);
        mirath_matrix_ff_init_random(C, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1, &prng);

        mirath_matrix_ff_mu_add_mu1ff(A, B, C, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1);

        ff_mu_t X[MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K] = {0};
        ff_mu_t Y[MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K];

        mirath_matrix_map_ff_to_ff_mu(Y, C, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1);

        mirath_matrix_ff_mu_add(X, B, Y, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1);
        mirath_matrix_ff_mu_add(X, A, X, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1);

        if (memcmp(X, zero, mirath_matrix_ff_mu_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1)) != 0) {
            printf("FAIL (mirath_matrix_ff_mu_add != mirath_matrix_ff_mu_add_mu1ff)\n");
            return -1;
        }
    }
    printf("OK!\n");

    printf("Testing:\n");
    printf("\tmirath_matrix_ff_mu_add_multiple_ff()\n");
    printf("\tmirath_matrix_ff_mu_add_multiple_2()\n");
    printf("\tmirath_matrix_ff_mu_add_multiple_3()\n");
    printf("Status: ");
    for (uint32_t i = 0; i < NTESTS; i++) {
        // TODO: use the correct parameters
        ff_mu_t A[MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K] = {0};
        ff_t B[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1)];
        ff_mu_t scalar;

        mirath_matrix_ff_init_random(B, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1, &prng);
        mirath_prng(&prng, &scalar, sizeof(scalar));

        mirath_matrix_ff_mu_add_multiple_ff(A, scalar, B, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1);

        ff_mu_t X[MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K] = {0};
        ff_mu_t Y[MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K];

        mirath_matrix_map_ff_to_ff_mu(Y, B, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1);

        mirath_matrix_ff_mu_add_multiple_2(X, scalar, Y, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1);

        for (uint32_t j = 0; j < MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K; j++) {
            if (A[j] != X[j]) {
                printf("FAIL (mirath_matrix_ff_mu_add_multiple_ff != mirath_matrix_ff_mu_add_multiple_2)\n");
                return -1;
            }
        }

        memset(X, 0, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K);

        mirath_matrix_ff_mu_add_multiple_3(X, X, scalar, Y, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1);

        for (uint32_t j = 0; j < MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K; j++) {
            if (A[j] != X[j]) {
                printf("FAIL (mirath_matrix_ff_mu_add_multiple_ff != mirath_matrix_ff_mu_add_multiple_3)\n");
                return -1;
            }
        }
    }
    printf("OK!\n");

    printf("Testing:\n");
    printf("\tmirath_matrix_ff_mu_product_ff1mu()\n");
    printf("\tmirath_matrix_ff_mu_product()\n");
    printf("Status: ");
    for (uint32_t i = 0; i < NTESTS; i++) {
        // TODO: use the correct parameters
        ff_mu_t R[MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K] = {0};
        ff_t A[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, MIRATH_PARAM_K)];
        ff_mu_t B[MIRATH_PARAM_K];

        mirath_matrix_ff_init_random(A, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, MIRATH_PARAM_K, &prng);
        mirath_prng(&prng, B, MIRATH_PARAM_K);

        mirath_matrix_ff_mu_product_ff1mu(R, A, B, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, MIRATH_PARAM_K, 1);

        ff_mu_t X[MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K] = {0};
        ff_mu_t Y[mirath_matrix_ff_mu_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, MIRATH_PARAM_K)];

        mirath_matrix_map_ff_to_ff_mu(Y, A, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, MIRATH_PARAM_K);

        mirath_matrix_ff_mu_product(X, Y, B, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, MIRATH_PARAM_K, 1);

        for (uint32_t j = 0; j < MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K; j++) {
            if (R[j] != X[j]) {
                printf("FAIL (mirath_matrix_ff_mu_product_ff1mu != mirath_matrix_ff_mu_product)\n");
                return -1;
            }
        }
    }
    printf("OK!\n");

    printf("Testing mirath_matrix_ff_mu_product_mu1ff(): ");
    for (uint32_t i = 0; i < NTESTS; i++) {
        // TODO: use the correct parameters
        ff_mu_t R[MIRATH_PARAM_K] = {0};
        ff_t A[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1)];
        ff_mu_t B[MIRATH_PARAM_K * (MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K)];

        mirath_matrix_ff_init_random(A, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1, &prng);
        mirath_prng(&prng, B, MIRATH_PARAM_K * (MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K));

        mirath_matrix_ff_mu_product_mu1ff(R, B, A, MIRATH_PARAM_K, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1);

        ff_mu_t X[MIRATH_PARAM_K] = {0};
        ff_mu_t Y[MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K];

        mirath_matrix_map_ff_to_ff_mu(Y, A, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1);

        mirath_matrix_ff_mu_product(X, B, Y, MIRATH_PARAM_K, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1);

        for (uint32_t j = 0; j < MIRATH_PARAM_K; j++) {
            if (R[j] != X[j]) {
                printf("FAIL (mirath_matrix_ff_mu_product_mu1ff != mirath_matrix_ff_mu_product)\n");
                return -1;
            }
        }
    }
    printf("OK!\n");

    return 0;
}

int test_vector_ff_mu(mirath_prng_t prng) {
    printf("Testing: \n");
    printf("\tmirath_vector_ff_mu_add_ff()\n");
    printf("\tmirath_vector_ff_mu_add()\n");
    printf("Status: ");
    for (uint32_t i = 0; i < NTESTS; i++) {
        ff_mu_t zero[MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K] = {0};

        ff_mu_t A[MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K] = {0};
        ff_mu_t B[MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K];
        ff_t C[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1)];

        mirath_prng(&prng, B, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K);
        mirath_matrix_ff_init_random(C, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1, &prng);

        mirath_vector_ff_mu_add_ff(A, B, C, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K);

        ff_mu_t X[MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K] = {0};
        ff_mu_t Y[MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K];

        mirath_matrix_map_ff_to_ff_mu(Y, C, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1);

        mirath_vector_ff_mu_add(X, B, Y, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K);
        mirath_vector_ff_mu_add(X, A, X, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K);

        if (memcmp(X, zero, mirath_matrix_ff_mu_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1)) != 0) {
            printf("FAIL (mirath_vector_ff_mu_add_ff != mirath_vector_ff_mu_add)\n");
            return -1;
        }
    }
    printf("OK!\n");

    printf("Testing:\n");
    printf("\tmirath_vector_ff_mu_mult_multiple_ff()\n");
    printf("\tmirath_vector_ff_mu_add_multiple_ff()\n");
    printf("\tmirath_vector_ff_mu_mult_multiple()\n");
    printf("\tmirath_vector_ff_mu_add_multiple()\n");
    printf("Status: ");
    for (uint32_t i = 0; i < NTESTS; i++) {
        const ff_mu_t zero[MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K] = {0};
        // TODO: use the correct parameters
        ff_mu_t A[MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K] = {0};
        ff_t B[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1)];
        ff_mu_t C[MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K];
        ff_mu_t scalar;

        ff_mu_t X[MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K] = {0};
        ff_mu_t Y[MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K];
        ff_mu_t Z[MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K];

        mirath_matrix_ff_init_random(B, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1, &prng);
        mirath_prng(&prng, &scalar, sizeof(scalar));
        mirath_prng(&prng, C, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K);

        mirath_vector_ff_mu_mult_multiple_ff(A, scalar, B, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K);

        mirath_matrix_map_ff_to_ff_mu(Y, B, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1);

        mirath_vector_ff_mu_mult_multiple(X, scalar, Y, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K);

        mirath_vector_ff_mu_add(X, A, X, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K);

        if (memcmp(X, zero, mirath_matrix_ff_mu_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1)) != 0) {
            printf("FAIL (mirath_vector_ff_mu_mult_multiple_ff != mirath_vector_ff_mu_mult_multiple)\n");
            return -1;
        }

        mirath_vector_ff_mu_add_multiple_ff(Z, C, scalar, B, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K);
        mirath_vector_ff_mu_add(X, A, C, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K);
        mirath_vector_ff_mu_add(Z, X, Z, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K);

        if (memcmp(Z, zero, mirath_matrix_ff_mu_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1)) != 0) {
            printf("FAIL: mirath_vector_ff_mu_add_multiple_ff\n");
            return -1;
        }

        mirath_vector_ff_mu_add_multiple(Z, C, scalar, Y, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K);
        mirath_vector_ff_mu_add(Z, X, Z, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K);

        if (memcmp(Z, zero, mirath_matrix_ff_mu_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1)) != 0) {
            printf("FAIL: mirath_vector_ff_mu_add_multiple_ff\n");
            return -1;
        }
    }
    printf("OK!\n");

    return 0;
}

int main(int argc, char **argv) {
    mirath_prng_t prng;
    mirath_prng_init(&prng, NULL, NULL, 0);

    tests_ff(prng);

    tests_ff_mu(prng);

    if (test_matrix_ff_mu(prng) != 0)
        return -1;

    if (test_vector_ff_mu(prng) != 0)
        return -1;

    return 0;
}