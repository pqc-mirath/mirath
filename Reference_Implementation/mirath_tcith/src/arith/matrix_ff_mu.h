#ifndef MIRATH_MATRIX_FF_MU_H
#define MIRATH_MATRIX_FF_MU_H

#include <stdint.h>
#include "ff_mu.h"
#include "matrix_ff_arith.h"
#include "prng.h"

#define mirath_matrix_ff_mu_get_entry(m,n,i,j) m[j*n + i]
#define mirath_matrix_ff_mu_set_entry(m,n,i,j,v) m[j*n + i] = v

#define mirath_matrix_ff_mu_bytes_size(x, y) ((x) * (y) * sizeof(ff_mu_t))

#define MIRATH_VAR_FF_MU_S_BYTES (mirath_matrix_ff_mu_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_R))
#define MIRATH_VAR_FF_MU_T_BYTES (mirath_matrix_ff_mu_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_N - MIRATH_PARAM_R))
#define MIRATH_VAR_FF_MU_E_A_BYTES (mirath_matrix_ff_mu_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1))
#define MIRATH_VAR_FF_MU_K_BYTES (mirath_matrix_ff_mu_bytes_size(MIRATH_PARAM_K, 1))

static inline void mirath_matrix_ff_mu_random(ff_mu_t *matrix, const uint32_t n_rows, const uint32_t n_cols, mirath_prng_t *prng) {
    mirath_prng(prng, matrix, ((n_rows) * (n_cols) * sizeof(ff_mu_t)));
}

static inline void mirath_matrix_ff_mu_copy(ff_mu_t *matrix1, const ff_mu_t *matrix2, const uint32_t n_rows, const uint32_t n_cols) {
    memcpy(matrix1, matrix2, mirath_matrix_ff_mu_bytes_size(n_rows, n_cols));
}

/**
 * \fn static inline void mirath_matrix_ff_mu_add(ff_mu_t *matrix1, const ff_mu_t *matrix2,
 *                           const ff_mu_t *matrix3, const uint32_t n_rows,
 *                           const uint32_t n_cols)
 * \brief matrix1 = matrix2 + matrix3
 *
 * \param[out] matrix1 Matrix over ff_mu
 * \param[in] matrix2 Matrix over ff_mu
 * \param[in] matrix3 Matrix over ff_mu
 * \param[in] n_rows number of rows
 * \param[in] n_cols number of columns
 */
static inline void mirath_matrix_ff_mu_add(ff_mu_t *matrix1, const ff_mu_t *matrix2,
                             const ff_mu_t *matrix3, const uint32_t n_rows,
                             const uint32_t n_cols) {
    const uint32_t n_bytes = mirath_matrix_ff_mu_bytes_size(n_rows, n_cols);

    for (uint32_t i = 0; i < n_bytes; i++) {
        matrix1[i] = matrix2[i] ^ matrix3[i];
    }
}

/**
 * \fn static inline void mirath_matrix_ff_mu_add_mu1ff(ff_mu_t *matrix1, const ff_mu_t *matrix2, const ff_t *matrix3,
 *                                               const uint32_t n_rows, const uint32_t n_cols)
 * \brief matrix1 = matrix2 + matrix3
 *
 * \param[out] matrix1 Matrix over ff_mu
 * \param[in] matrix2 Matrix over ff_mu
 * \param[in] matrix3 Matrix over ff
 * \param[in] n_rows number of rows
 * \param[in] n_cols number of columns
 */
static inline void mirath_matrix_ff_mu_add_mu1ff(ff_mu_t *matrix1, const ff_mu_t *matrix2, const ff_t *matrix3,
                                                 const uint32_t n_rows, const uint32_t n_cols) {
    for (uint32_t i = 0; i < n_rows; i++) {
        for (uint32_t j = 0; j < n_cols; j++) {
            ff_mu_t entry2;

            const ff_mu_t entry1 = mirath_matrix_ff_mu_get_entry(matrix2, n_rows, i, j);
            entry2 = mirath_matrix_ff_get_entry(matrix3, n_rows, i, j);
            entry2 = mirath_map_ff_to_ff_mu[entry2];
            entry2 = entry1 ^ entry2;
            mirath_matrix_ff_mu_set_entry(matrix1, n_rows, i, j, entry2);
        }
    }
}

/**
 * \fn static inline void mirath_matrix_ff_mu_add_multiple_ff(ff_mu_t *matrix1, ff_mu_t scalar, const ff_t *matrix2,
 *                                       const uint32_t n_rows, const uint32_t n_cols)
 * \brief matrix1 += scalar * matrix2
 *
 * \param[out] matrix1 Matrix over ff_mu
 * \param[in] scalar scalar over ff_mu
 * \param[in] matrix2 Matrix over ff
 * \param[in] n_rows number of rows
 * \param[in] n_cols number of columns
 */
static inline void mirath_matrix_ff_mu_add_multiple_ff(ff_mu_t *matrix1, ff_mu_t scalar, const ff_t *matrix2,
                                         const uint32_t n_rows, const uint32_t n_cols) {
    for (uint32_t i = 0; i < n_rows; i++) {
        for (uint32_t j = 0; j < n_cols; j++) {
            ff_mu_t entry2, entry3;

            const ff_mu_t entry1 = mirath_matrix_ff_mu_get_entry(matrix1, n_rows, i, j);
            entry2 = mirath_matrix_ff_get_entry(matrix2, n_rows, i, j);
            entry2 = mirath_map_ff_to_ff_mu[entry2];
            entry3 = entry1 ^ mirath_ff_mu_mult(scalar, entry2);
            mirath_matrix_ff_mu_set_entry(matrix1, n_rows, i, j, entry3);
        }
    }
}

/**
 * \fn static void mirath_matrix_ff_mu_add_multiple_2(ff_mu_t *matrix1, ff_mu_t scalar, const ff_mu_t *matrix2,
 *                                            const uint32_t n_rows, const uint32_t n_cols)
 * \brief matrix1 += scalar * matrix2
 *
 * \param[out] matrix1 Matrix over ff_mu
 * \param[in] scalar scalar over ff_mu
 * \param[in] matrix2 Matrix over ff_mu
 * \param[in] n_rows number of rows
 * \param[in] n_cols number of columns
 */
static void mirath_matrix_ff_mu_add_multiple_2(ff_mu_t *matrix1, ff_mu_t scalar, const ff_mu_t *matrix2,
                                              const uint32_t n_rows, const uint32_t n_cols) {
    for (uint32_t i = 0; i < n_rows; i++) {
        for (uint32_t j = 0; j < n_cols; j++) {
            const ff_mu_t entry1 = mirath_matrix_ff_mu_get_entry(matrix1, n_rows, i, j);
            const ff_mu_t entry2 = mirath_matrix_ff_mu_get_entry(matrix2, n_rows, i, j);
            const ff_mu_t entry3 = entry1 ^ mirath_ff_mu_mult(scalar, entry2);

            mirath_matrix_ff_mu_set_entry(matrix1, n_rows, i, j, entry3);
        }
    }
}

/**
 * \fn static inline void mirath_matrix_ff_mu_add_multiple_3(ff_mu_t *matrix1, const ff_mu_t *matrix2,
 *                                    const ff_mu_t scalar, const ff_mu_t *matrix3,
 *                                    const uint32_t n_rows, const uint32_t n_cols)
 * \brief matrix1 = matrix2 + scalar * matrix3
 *
 * \param[out] matrix1 Matrix over ff_mu
 * \param[in] matrix2 Matrix over ff_mu
 * \param[in] scalar scalar over ff_mu
 * \param[in] matrix3 Matrix over ff_mu
 * \param[in] n_rows number of rows
 * \param[in] n_cols number of columns
 */
static inline void mirath_matrix_ff_mu_add_multiple_3(ff_mu_t *matrix1, const ff_mu_t *matrix2,
                                      const ff_mu_t scalar, const ff_mu_t *matrix3,
                                      const uint32_t n_rows, const uint32_t n_cols) {
    for (uint32_t i = 0; i < n_rows; i++) {
        for (uint32_t j = 0; j < n_cols; j++) {
            const ff_mu_t entry1 = mirath_matrix_ff_mu_get_entry(matrix2, n_rows, i, j);
            const ff_mu_t entry2 = mirath_matrix_ff_mu_get_entry(matrix3, n_rows, i, j);
            const ff_mu_t entry3 = entry1 ^ mirath_ff_mu_mult(scalar, entry2);

            mirath_matrix_ff_mu_set_entry(matrix1, n_rows, i, j, entry3);
        }
    }
}

/**
 * \fn static inline void mirath_matrix_ff_mu_product_ff1mu(ff_mu_t *result, const ff_t *matrix1,
 *                                     const ff_mu_t *matrix2, const uint32_t n_rows1,
 *                                     const uint32_t n_cols1, const uint32_t n_cols2)
 * \brief result = matrix1 * matrix2
 *
 * \param[out] result Matrix over ff_mu
 * \param[in] matrix1 Matrix over ff
 * \param[in] matrix2 Matrix over ff_mu
 * \param[in] n_rows1 number of rows in matrix1
 * \param[in] n_cols1 number of columns and rows in matrix1 and matrix2 respectively
 * \param[in] n_cols2 number of columns in matrix2
 */
static inline void mirath_matrix_ff_mu_product_ff1mu(ff_mu_t *result, const ff_t *matrix1,
                                       const ff_mu_t *matrix2, const uint32_t n_rows1,
                                       const uint32_t n_cols1, const uint32_t n_cols2) {
    ff_mu_t entry_i_k, entry_k_j, entry_i_j;

    for(uint32_t i = 0; i < n_rows1; i++) {
        for (uint32_t j = 0; j < n_cols2; j++) {
            entry_i_j = 0;

            for (uint32_t k = 0; k < n_cols1; k++) {
                entry_i_k = mirath_matrix_ff_get_entry(matrix1, n_rows1, i, k);
                entry_i_k = mirath_map_ff_to_ff_mu[entry_i_k];
                entry_k_j = mirath_matrix_ff_mu_get_entry(matrix2, n_cols1, k, j);
                entry_i_j ^= mirath_ff_mu_mult(entry_i_k, entry_k_j);
            }

            mirath_matrix_ff_mu_set_entry(result, n_rows1, i, j, entry_i_j);
        }
    }
}

/**
 * \fn static inline void mirath_matrix_ff_mu_product_mu1ff(ff_mu_t *result, const ff_mu_t *matrix1,
 *                                     const ff_t *matrix2, const uint32_t n_rows1,
 *                                     const uint32_t n_cols1, const uint32_t n_cols2)
 * \brief result = matrix1 * matrix2
 *
 * \param[out] result Matrix over ff_mu
 * \param[in] matrix1 Matrix over ff_mu
 * \param[in] matrix2 Matrix over ff
 * \param[in] n_rows1 number of rows in matrix1
 * \param[in] n_cols1 number of columns and rows in matrix1 and matrix2 respectively
 * \param[in] n_cols2 number of columns in matrix2
 */
static inline void mirath_matrix_ff_mu_product_mu1ff(ff_mu_t *result, const ff_mu_t *matrix1,
                                       const ff_t *matrix2, const uint32_t n_rows1,
                                       const uint32_t n_cols1, const uint32_t n_cols2) {
    ff_mu_t entry_i_k, entry_k_j, entry_i_j;

    for(uint32_t i = 0; i < n_rows1; i++) {
        for (uint32_t j = 0; j < n_cols2; j++) {
            entry_i_j = 0;

            for (uint32_t k = 0; k < n_cols1; k++) {
                entry_i_k = mirath_matrix_ff_mu_get_entry(matrix1, n_rows1, i, k);
                entry_k_j = mirath_matrix_ff_get_entry(matrix2, n_cols1, k, j);
                entry_k_j = mirath_map_ff_to_ff_mu[entry_k_j];
                entry_i_j ^= mirath_ff_mu_mult(entry_i_k, entry_k_j);
            }

            mirath_matrix_ff_mu_set_entry(result, n_rows1, i, j, entry_i_j);
        }
    }
}

/**
 * \fn static inline void mirath_matrix_ff_mu_product(ff_mu_t *result, const ff_mu_t *matrix1, const ff_mu_t *matrix2,
 *                               const uint32_t n_rows1, const uint32_t n_cols1,
 *                               const uint32_t n_cols2)
 * \brief result = matrix1 * matrix2
 *
 * \param[out] result Matrix over ff_mu
 * \param[in] matrix1 Matrix over ff_mu
 * \param[in] matrix2 Matrix over ff_mu
 * \param[in] n_rows1 number of rows in matrix1
 * \param[in] n_cols1 number of columns and rows in matrix1 and matrix2 respectively
 * \param[in] n_cols2 number of columns in matrix2
 */
static inline void mirath_matrix_ff_mu_product(ff_mu_t *result, const ff_mu_t *matrix1, const ff_mu_t *matrix2,
                                 const uint32_t n_rows1, const uint32_t n_cols1,
                                 const uint32_t n_cols2) {
    ff_mu_t entry_i_k, entry_k_j, entry_i_j;

    for(uint32_t i = 0; i < n_rows1; i++) {
        for (uint32_t j = 0; j < n_cols2; j++) {
            entry_i_j = 0;

            for (uint32_t k = 0; k < n_cols1; k++) {
                entry_i_k = mirath_matrix_ff_mu_get_entry(matrix1, n_rows1, i, k);
                entry_k_j = mirath_matrix_ff_mu_get_entry(matrix2, n_cols1, k, j);
                entry_i_j ^= mirath_ff_mu_mult(entry_i_k, entry_k_j);
            }

            mirath_matrix_ff_mu_set_entry(result, n_rows1, i, j, entry_i_j);
        }
    }
}

/**
 * \fn static inline void mirath_matrix_ff_mu_add_product(ff_mu_t *result, const ff_mu_t *matrix1,
 *                                   const ff_mu_t *matrix2, const uint32_t n_rows1,
 *                                   const uint32_t n_cols1, const uint32_t n_cols2)
 * \brief result += matrix1 * matrix2
 *
 * \param[out] result Matrix over ff_mu
 * \param[in] matrix1 Matrix over ff_mu
 * \param[in] matrix2 Matrix over ff_mu
 * \param[in] n_rows1 number of rows in matrix1
 * \param[in] n_cols1 number of columns and rows in matrix1 and matrix2 respectively
 * \param[in] n_cols2 number of columns in matrix2
 */
static inline void mirath_matrix_ff_mu_add_product(ff_mu_t *result, const ff_mu_t *matrix1,
                                     const ff_mu_t *matrix2, const uint32_t n_rows1,
                                     const uint32_t n_cols1, const uint32_t n_cols2) {
    ff_mu_t entry_i_k, entry_k_j, entry_i_j;

    for(uint32_t i = 0; i < n_rows1; i++) {
        for (uint32_t j = 0; j < n_cols2; j++) {
            entry_i_j = mirath_matrix_ff_mu_get_entry(result, n_rows1, i, j);

            for (uint32_t k = 0; k < n_cols1; k++) {
                entry_i_k = mirath_matrix_ff_mu_get_entry(matrix1, n_rows1, i, k);
                entry_k_j = mirath_matrix_ff_mu_get_entry(matrix2, n_cols1, k, j);
                entry_i_j ^= mirath_ff_mu_mult(entry_i_k, entry_k_j);
            }

            mirath_matrix_ff_mu_set_entry(result, n_rows1, i, j, entry_i_j);
        }
    }
}

/**
 * \fn static inline void mirath_matrix_map_ff_to_ff_mu(ff_mu_t *out, const uint8_t *input, const uint32_t nrows, const uint32_t ncols)
 * \brief mapping from ff to ff_mu
 *
 * \param[out] out Matrix over ff_mu
 * \param[in] input Matrix over ff
 * \param[in] nrows number of rows
 * \param[in] ncols number of columns
 */
static inline void mirath_matrix_map_ff_to_ff_mu(ff_mu_t *out, const uint8_t *input, const uint32_t nrows, const uint32_t ncols) {
    for (uint32_t i = 0; i < ncols; ++i) {
        for (uint32_t j = 0; j < nrows; ++j) {
            const uint8_t tmp = mirath_matrix_ff_get_entry(input, nrows, j, i);
            *out = mirath_map_ff_to_ff_mu[tmp];
            out += 1;
        }
    }
}

#endif
