
#include <assert.h>
#include <stdint.h>
#include <string.h>
#include "mirath_matrix_ff.h"

void mirath_matrix_init_zero(ff_t *matrix, const uint32_t n_rows, const uint32_t n_cols) {
    memset(matrix, 0, mirath_matrix_ff_bytes_size(n_rows, n_cols));
}

//
void mirath_matrix_ff_init_random(ff_t *matrix, const uint32_t n_rows, const uint32_t n_cols, mirath_prng_t *prng) {
    const uint32_t matrix_bytes = mirath_matrix_ff_bytes_size(n_rows, n_cols);
    mirath_prng(prng, matrix, matrix_bytes);

    mirath_matrix_set_to_ff(matrix, n_rows, n_cols);
}

void mirath_matrix_ff_copy(ff_t *matrix1, const ff_t *matrix2, const uint32_t n_rows, const uint32_t n_cols) {
    const uint32_t n_bytes = mirath_matrix_ff_bytes_size(n_rows, n_cols);
    memcpy(matrix1, matrix2, n_bytes);
}

void mirath_matrix_ff_neg(ff_t *matrix, const uint32_t n_rows, const uint32_t n_cols) {
    /* Nothing to do in characteristic 2. */

    /* Suppress 'unused parameter' warnings. */
    (void)(matrix); (void)(n_rows); (void)(n_cols);
}

void mirath_matrix_ff_add(ff_t *matrix1, const ff_t *matrix2, const ff_t *matrix3, uint32_t n_rows, uint32_t n_cols){
    mirath_matrix_ff_add_arith(matrix1, matrix2, matrix3, n_rows, n_cols);
}

void mirath_matrix_ff_add_multiple(ff_t *matrix1, ff_t scalar, const ff_t *matrix2,
    const uint32_t n_rows, const uint32_t n_cols) {
    mirath_matrix_ff_add_multiple_arith(matrix1, scalar, matrix2, n_rows, n_cols);
}

void mirath_matrix_ff_sub(ff_t *matrix1, const ff_t *matrix2, const ff_t *matrix3, uint32_t n_rows, uint32_t n_cols){
    mirath_matrix_ff_add_arith(matrix1, matrix2, matrix3, n_rows, n_cols);
}

void mirath_matrix_ff_sub_multiple(ff_t *matrix1, ff_t scalar, const ff_t *matrix2,
    const uint32_t n_rows, const uint32_t n_cols) {
    mirath_matrix_ff_add_multiple(matrix1, scalar, matrix2, n_rows, n_cols);
}

void mirath_matrix_ff_product(ff_t *result, const ff_t *matrix1, const ff_t *matrix2,
    const uint32_t n_rows1, const uint32_t n_cols1, const uint32_t n_cols2) {
    mirath_matrix_ff_product_arith(result, matrix1, matrix2, n_rows1, n_cols1, n_cols2);
}

void mirath_matrix_ff_horizontal_concat(ff_t *result, const ff_t *matrix1, const ff_t *matrix2,
    const uint32_t n_rows, const uint32_t n_cols1, const uint32_t n_cols2) {
    const uint32_t n_bytes1 = mirath_matrix_ff_bytes_size(n_rows, n_cols1);
    const uint32_t n_bytes2 = mirath_matrix_ff_bytes_size(n_rows, n_cols2);

    memcpy(result, matrix1, n_bytes1);
    memcpy(result + n_bytes1, matrix2, n_bytes2);
}

void mirath_matrix_ff_horizontal_split(ff_t *matrix1, ff_t *matrix2, const ff_t *matrix,
    const uint32_t n_rows, const uint32_t n_cols1, const uint32_t n_cols2) {
    const uint32_t n_bytes1 = mirath_matrix_ff_bytes_size(n_rows, n_cols1);
    const uint32_t n_bytes2 = mirath_matrix_ff_bytes_size(n_rows, n_cols2);

    if (matrix1 != NULL) {
        memcpy(matrix1, matrix, n_bytes1);
        memcpy(matrix2, matrix + n_bytes1, n_bytes2);
    }
    else {
        memcpy(matrix2, matrix + n_bytes1, n_bytes2);
    }
}

void _matrix_pack_nrows_even(uint8_t **dest, const uint32_t *bit_offset, const ff_t *matrix,
                             const uint32_t n_rows, const int n_cols) {

    /* the packing is done row-wise */
    uint32_t bo, n_bytes;

    if (bit_offset != NULL)
    {
        bo = *bit_offset;
    }
    else
    {
        bo = 0;
    }

    n_bytes = mirath_matrix_ff_bytes_size(n_rows, n_cols);

    if (bo)
    {
        /* Pack last entry of matrix (in column-major order) in the higher bits of dest[0]. */
        ((uint8_t *)*dest)[0] |= matrix[n_bytes - 1] & 0xf0;

        /* Pack all the bytes in matrix except the last one */
        memcpy(&(((uint8_t *)*dest)[1]), matrix, n_bytes - 1);

        /* Pack the second-to-last entry of matrix. */
        ((uint8_t *)*dest)[n_bytes] = matrix[n_bytes - 1] & 0x0f;
    }
    else
    {
        memcpy((uint8_t *)*dest, matrix, n_bytes);
    }

    *dest = &(((uint8_t *)*dest)[n_bytes]);
}

void _matrix_unpack_nrows_even(ff_t *matrix, uint8_t **source, const uint32_t *bit_offset,
                               const uint32_t n_rows, const int n_cols)
{
    uint32_t bo, n_bytes;

    if (bit_offset != NULL)
    {
        bo = *bit_offset;
    }
    else
    {
        bo = 0;
    }

    n_bytes = mirath_matrix_ff_bytes_size(n_rows, n_cols);

    if (bo)
    {
        /* Unpack all the bytes in matrix except the last one */
        memcpy(matrix, &(((uint8_t *)*source)[1]), n_bytes - 1); /* unpack all the bytes but the last one. */

        /* Unpack the last two entries of matrix. */
        matrix[n_bytes - 1] = (((uint8_t *)*source)[n_bytes] & 0x0f) | (((uint8_t *)*source)[0] & 0xf0);
    }
    else
    {
        memcpy(matrix, (uint8_t *)*source, n_bytes);
    }

    *source = &(((uint8_t *)*source)[n_bytes]);
}


/* Remove the last row of matrix and append it to matrix as additional column(s) */
void _matrix_pack_nrows_odd(uint8_t **dest, uint32_t *bit_offset, const ff_t *matrix,
                            const uint32_t n_rows, const uint32_t n_cols)
{
    assert((n_rows & 1) == 1);

    uint32_t j, bo, next_bo, matrix_height, matrix_height_x, n_bytes_not_in_last_row, n_bytes;
    uint32_t ad_bytes, jump_nbytes;
    uint8_t row_entry_j, row_entry_j_1;

    if (bit_offset != NULL)
    {
        bo = *bit_offset;
    }
    else
    {
        bo = 0;
    }

    matrix_height = (n_rows >> 1) + 1;
    matrix_height_x =  matrix_height - 1;
    n_bytes_not_in_last_row = matrix_height_x * n_cols;
    n_bytes = mirath_matrix_ff_bytes_size(n_rows, n_cols);

    /* Bytes that are not part of the last row. */
    for (j = 0u; j < n_cols; j++)
    {
        memcpy(&(((uint8_t *)*dest)[bo + j * matrix_height_x]), &matrix[matrix_height * j], matrix_height_x);
    }
    /* When n_cols is odd the maximum value of j is j_max = n_cols - 3, hence j_max + 1 = n_cols - 2.
     * Hence in the following loop wont add the entry n_cols - 1 (the last entry) of the last row.
     * When n_cols is even the maximum value of j is j_max = n_cols - 4, hence j_max + 1 = n_cols - 3.
     * Hence in the following loop wont add the entries n_cols - 2 and n_cols - 1 (the last entry) of the last row. */
    ad_bytes = bo;
    for (j = 0; (int)j < (int)n_cols - 2; j+=2)
    {
        row_entry_j = matrix[matrix_height * j + matrix_height_x] & 0x0f; /* j-th entry of the last row. */
        row_entry_j_1 = matrix[matrix_height * (j + 1) + matrix_height_x] & 0x0f; /* (j + 1)-th entry of the last row. */
        ((uint8_t *)*dest)[n_bytes_not_in_last_row + ad_bytes]  =  (row_entry_j_1 << 4) | row_entry_j;
        ad_bytes +=1;
    }
    /* When the is an odd number of columns and
     * bit_off_set = 1, we locate the last entry of matrix in higher bits
     * of the first byte of the source. Otherwise, if bit_off_set = 0, we locate ast entry of matrix
     * the next byte of dest. */
    if  (bo)
    {
        ((uint8_t *)*dest)[0]&= 0x0f;
        ((uint8_t *)*dest)[0] |= (matrix[n_bytes - 1] << 4);

        if ((n_cols & 1) == 0) /* case n_cols is even. */
        {
            /* Packing the second-last entry of the last row. */
            ((uint8_t *)*dest)[n_bytes_not_in_last_row + ad_bytes] = matrix[n_bytes - matrix_height - 1] & 0x0f;

        }
    }
    else
    {
        /* If n_cols is even and bo = 0,
         * we pack the entries n_cols - 2 and n_cols - 1
         * in the last row in the byte
         * of the current local buffer. */
        if ((n_cols & 1) == 0)
        {
            ((uint8_t *)*dest)[n_bytes_not_in_last_row + ad_bytes] = (matrix[n_bytes - 1] << 4) | matrix[n_bytes - matrix_height - 1];
        }
            /* Odd number of columns case. */
            /* In this case, we locate the last entry the next byte of dest. */
        else
        {
            ((uint8_t *)*dest)[n_bytes_not_in_last_row + ad_bytes] = matrix[n_bytes - 1] & 0x0f;
        }
    }

    jump_nbytes = matrix_height_x * n_cols + (n_cols >> 1);

    if (bo)
    {
        if ((n_cols & 1) == 0)
        {
            next_bo = 1;
        }
        else
        {
            next_bo = 0;
            jump_nbytes += (n_cols & 1);
        }
    }
    else
    {
        if ((n_cols & 1) == 0)
        {
            next_bo = 0;
        }
        else
        {
            next_bo = 1;
        }
    }

    if (bit_offset != NULL)
    {
        *bit_offset = next_bo;
    }

    *dest = &(((uint8_t *)*dest)[jump_nbytes]);
}

void _matrix_unpack_nrows_odd(ff_t *matrix, uint8_t **source, uint32_t *bit_offset,
                              const uint32_t n_rows, const uint32_t n_cols)
{
    assert((n_rows & 1) == 1);

    uint32_t j, bo, next_bo, matrix_height, matrix_height_x, n_bytes_not_in_last_row, n_bytes;
    uint32_t ad_bytes, jump_nbytes;
    uint8_t row_entries_j_and_j_1;

    if (bit_offset != NULL)
    {
        bo = *bit_offset;
    }
    else
    {
        bo = 0;
    }


    matrix_height = (n_rows >> 1) + 1;
    matrix_height_x =  matrix_height - 1;
    n_bytes_not_in_last_row = matrix_height_x * n_cols;
    n_bytes = mirath_matrix_ff_bytes_size(n_rows, n_cols);


    /* Bytes that are not part of the last row. */
    for (j = 0; j < n_cols; j++)
    {
        memcpy(&matrix[matrix_height * j], &(((uint8_t *)*source)[bo + j * matrix_height_x]), matrix_height_x);
    }

    /* When n_cols is odd the maximum value of j is j_max = n_cols - 3, hence j_max + 1 = n_cols - 2.
     * Hence in the following loop wont add the entry n_cols - 1 (the last entry) of the last row .
     * When n_cols is even the maximum value of j is j_max = n_cols - 4, hence j_max + 1 = n_cols - 3.
     * Hence in the following loop wont add the entries n_cols - 2 and n_cols - 1 (the last entry) of the last row. */
    ad_bytes = bo;
    for (j = 0; (int)j < (int)n_cols - 2; j+=2)
    {
        row_entries_j_and_j_1 = ((uint8_t *)*source)[n_bytes_not_in_last_row + ad_bytes];
        matrix[matrix_height * j + matrix_height_x] = row_entries_j_and_j_1 & 0x0f;
        matrix[matrix_height * (j + 1) + matrix_height_x] =  row_entries_j_and_j_1 >> 4;
        ad_bytes +=1;
    }
    /* When the is an odd number of columns and
     * bit_off_set = 1, we locate the last entry of matrix in higher bits
     * of the first byte of the source. Otherwise, if bit_off_set = 0, we locate ast entry of matrix
     * the next byte of dest. */
    if  (bo)
    {
        matrix[n_bytes - 1]  = ((uint8_t *)*source)[0] >> 4;

        if ((n_cols & 1) == 0) /* case n_rows is even. */
        {
            matrix[n_bytes - matrix_height - 1] = ((uint8_t *)*source)[n_bytes_not_in_last_row + ad_bytes] & 0x0f;
        }
    }
    else
    {
        /* If n_cols is even and bo = 0,
        * we unpack the last row in the byte of the current local buffer
        * into the entries n_cols - 2 and n_cols - 1
        * of the last row of matrix .*/
        if ((n_cols & 1) == 0)
        {
            row_entries_j_and_j_1 =  ((uint8_t *)*source)[n_bytes_not_in_last_row + ad_bytes];
            matrix[n_bytes - 1]  = row_entries_j_and_j_1 >> 4;
            matrix[n_bytes - matrix_height - 1] = row_entries_j_and_j_1 & 0x0f;
        }
            /* Odd number of columns case.
             * In this case, we locate the last entry the next byte of dest. */
        else
        {
            matrix[n_bytes -1 ] = ((uint8_t *)*source)[n_bytes_not_in_last_row + ad_bytes] & 0x0f;
        }

    }

    jump_nbytes = matrix_height_x * n_cols + (n_cols >> 1);

    if (bo)
    {
        if ((n_cols & 1) == 0)
        {
            next_bo = 1;
        }
        else
        {
            next_bo = 0;
            jump_nbytes +=  (n_cols & 1);
        }
    }
    else
    {
        if ((n_cols & 1) == 0)
        {
            next_bo = 0;
        }
        else
        {
            next_bo = 1;
        }
    }

    if (bit_offset != NULL)
    {
        *bit_offset = next_bo;
    }

    *source = &(((uint8_t *)*source)[jump_nbytes]);
}

void mirath_matrix_ff_unparse(uint8_t **dest, uint32_t *bit_offset, const ff_t *matrix,
                 const uint32_t n_rows, const uint32_t n_cols)
{
    /* An even number of rows. */
    if ((n_rows & 1) == 0)
    {
        _matrix_pack_nrows_even(dest, bit_offset, matrix, n_rows, n_cols);
    }
    else
    {
        _matrix_pack_nrows_odd(dest, bit_offset, matrix, n_rows, n_cols);
    }

}

void mirath_matrix_ff_parse(ff_t *matrix, uint8_t **source, uint32_t *bit_offset,
                   const uint32_t n_rows, const uint32_t n_cols)
{
    /* An even number of rows. */
    if ((n_rows & 1) == 0)
    {
        _matrix_unpack_nrows_even(matrix, source, bit_offset, n_rows, n_cols);
    }
    else
    {
        _matrix_unpack_nrows_odd(matrix, source, bit_offset, n_rows, n_cols);
    }
}


