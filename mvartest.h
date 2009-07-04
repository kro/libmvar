#ifndef MVARTEST_H
#define MVARTEST_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#define assert_mat_equals(mat1, mat2, delta) \
    do_assert_mat_equals(__FILE__, __LINE__, mat1, mat2, delta);

#define assert_vec_equals(vec1, vec2, delta) \
    do_assert_vec_equals(__FILE__, __LINE__, vec1, vec2, delta);

void do_assert_mat_equals(const char *file, int line, gsl_matrix *mat1, gsl_matrix *mat2, double delta);
void do_assert_vec_equals(const char *file, int line, gsl_vector *vec1, gsl_vector *vec2, double delta);

#endif /* MVARTEST_H */
