#include "mvarmat.h"

#include <math.h>

void mvar_mat_copy(gsl_matrix *dst, gsl_matrix *src, size_t k1, size_t k2)
{
    size_t row, col;

    for (row = 0; row < src->size1; row++) {
        for (col = 0; col < src->size2; col++) {
            gsl_matrix_set(dst, k1 + row, k2 + col, gsl_matrix_get(src, row, col));
        }
    }
}

void mvar_mat_set_diag(gsl_matrix *mat, gsl_vector *vec)
{
    int i;

    for (i = 0; i < vec->size; i++)
        gsl_matrix_set(mat, i, i, gsl_vector_get(vec, i));
}

void mvar_mat_sum_sq_sqrt(gsl_matrix *mat, gsl_vector *vec)
{
    double acc;
    int i, j;

    for (i = 0; i < mat->size2; i++) {
        acc = 0.0;
        for (j = 0; j < mat->size1; j++) {
            acc += pow(gsl_matrix_get(mat, j, i), 2);
        }
        gsl_vector_set(vec, i, sqrt(acc));
    }
}

void mvar_mat_upper_tri(gsl_matrix *mat)
{
    size_t row, col;

    for (row = 0; row < mat->size1; row++) {
        for (col = 0; col < mat->size2; col++) {
            if (col >= row)
                gsl_matrix_set(mat, row, col, gsl_matrix_get(mat, row, col));
            else
                gsl_matrix_set(mat, row, col, 0);
        }
    }
}
