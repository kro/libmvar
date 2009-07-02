#ifndef MVARMAT_H
#define MVARMAT_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/* Copies 'src' element by element to 'dst' starting from row 'k1' and column
 * 'k2'.  Unlike gsl_matrix_memcpy(), this function does not require that the
 * matrices are of same size as long as 'dst' is bigger than 'src'.
 */
void mvar_mat_copy(gsl_matrix *dst, gsl_matrix *src, size_t k1, size_t k2);

/* Set the elements of 'vec' as the diagonal of 'mat'. */
void mvar_mat_set_diag(gsl_matrix *mat, gsl_vector *vec);

/* Sets the elements of 'vec' as the square root of sum squares of elements of
 * each column of 'mat'.
 */
void mvar_mat_sum_sq_sqrt(gsl_matrix *mat, gsl_vector *vec);

/* Sets an upper triangular matrix by setting all element to zero except the
 * ones in the upper triangle.
 */
void mvar_mat_upper_tri(gsl_matrix *mat);

#endif /* MVARMAT_H */
