#ifndef MIRATH_MATRIX_FF_H
#define MIRATH_MATRIX_FF_H

#include <stdint.h>
#include "arith/data_type_arith.h"
#include "prng.h"
#include "arith/matrix_ff_arith.h"

/* Return the number of bytes of a 'n_rows x n_cols' matrix. */
//int matrix_bytes_size(int n_rows, int n_cols);

/* Initialized 'matrix' with zero entries. */
void mirath_matrix_init_zero(ff_t *matrix, const uint32_t n_rows, const uint32_t n_cols);

/* Return the (i, j) entry of 'matrix'. */
ff_t mirath_matrix_ff_get_entry(const ff_t *matrix, const uint32_t n_rows, const uint32_t i, const uint32_t j);

/* Set the (i, j) entry of 'matrix' to be scalar.*/
void mirath_matrix_ff_set_entry(ff_t *matrix, const uint32_t n_rows, const uint32_t i, const uint32_t j, const ff_t scalar);

/* Initialized 'matrix' with random entries. */
void mirath_matrix_ff_init_random(ff_t *matrix, const uint32_t n_rows, const uint32_t n_cols, mirath_prng_t *prng);

/* Overwrite 'matrix1' with 'matrix2'. */
void mirath_matrix_ff_copy(ff_t *matrix1, const ff_t *matrix2, const uint32_t n_rows, const uint32_t n_cols);

/* Replace 'matrix' with '-matrix'. */
void mirath_matrix_ff_neg(ff_t *matrix, const uint32_t n_rows, const uint32_t n_cols);

/* set 'matrix1' with 'matrix2 + matrix3'. */
void mirath_matrix_ff_add(ff_t *matrix1, const ff_t *matrix2, const ff_t *matrix3, uint32_t n_rows, uint32_t n_cols);

/* Overwrite 'matrix1' with 'matrix1 + scalar * matrix2'. */
void mirath_matrix_ff_add_multiple(ff_t *matrix1, ff_t scalar, const ff_t *matrix2,
                                   const uint32_t n_rows, const uint32_t n_cols);

/* set 'matrix1' with 'matrix2 - matrix3'. */
void mirath_matrix_ff_sub(ff_t *matrix1, const ff_t *matrix2, const ff_t *matrix3, uint32_t n_rows, uint32_t n_cols);

/* Overwrite 'matrix1' with 'matrix1 - scalar * matrix2'. */
void mirath_matrix_ff_sub_multiple(ff_t *matrix1, ff_t scalar, const ff_t *matrix2,
                                   const uint32_t n_rows, const uint32_t n_cols);

/* Write 'matrix1 * matrix2' over 'result'. */
void mirath_matrix_ff_product(ff_t *result, const ff_t *matrix1, const ff_t *matrix2,
                              const uint32_t n_rows1, const uint32_t n_cols1, const uint32_t n_cols2);

/* Write '[matrix1 | matrix2]' over 'result'. */
void mirath_matrix_ff_horizontal_concat(ff_t *result, const ff_t *matrix1, const ff_t *matrix2,
                                        const uint32_t n_rows, const uint32_t n_cols1, const uint32_t n_cols2);

/* Split 'matrix' as 'matrix = [matrix1 | matrix2]. */
void mirath_matrix_ff_horizontal_split(ff_t *matrix1, ff_t *matrix2, const ff_t *matrix,
                                       const uint32_t n_rows, const uint32_t n_cols1, const uint32_t n_cols2);

/* Pack 'matrix' over 'dest (bytes) + bit_offset (bits)'.
 * Update 'dest' and 'bit_offset' for the next call of 'matrix_pack'.
 * If 'bit_offset == NULL', then an offset of 0 bits is used. */
//void matrix_pack(uint8_t **dest, uint32_t *bit_offset, const ff_t *matrix,
//    uint32_t n_rows, uint32_t n_cols);

/* Unpack 'matrix' from 'source (bytes) + bit_offset (bits)'.
 * Update 'source' and 'bit_offset' for the next call of 'matrix_unpack'.
 * If 'bit_offset == NULL', then an offset of 0 bits is used. */
//void matrix_unpack(ff_t *matrix, uint8_t **source, uint32_t *bit_offset,
//    uint32_t n_rows, uint32_t n_cols);

#endif
