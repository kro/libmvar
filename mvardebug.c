#include "mvardebug.h"

#include <stdio.h>

void mvar_print_vector(gsl_vector *vec)
{
    int i;

    for (i = 0; i < vec->size; i++) {
        if (i > 0)
            printf(",");
        printf("%f", gsl_vector_get(vec, i));
    }

    printf("\n");
}

void mvar_print_matrix(gsl_matrix *mat)
{
    int i, j;

    for (i = 0; i < mat->size1; i++) {
        for (j = 0; j < mat->size2; j++) {
            if (j > 0)
                 printf(",");
            printf("%f", gsl_matrix_get(mat, i, j));
        }
        printf("\n");
    }
}

