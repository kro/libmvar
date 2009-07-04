#include "mvartest.h"

#include <stdarg.h>
#include <math.h>

static void fail(const char *file, int line, const char *format, ...)
{
    va_list ap;

    va_start(ap, format);
    fprintf(stderr, "Failed %s (%d): ", file, line);
    vfprintf(stderr, format, ap);
    va_end(ap);

    fputc('\n', stderr);

    abort();
}

static void
assert_mat_elem(const char *file, int line, gsl_matrix *mat1, gsl_matrix *mat2, size_t row, size_t col, double delta)
{
    double v1 = gsl_matrix_get(mat1, row, col);
    double v2 = gsl_matrix_get(mat1, row, col);

    if (fabs(v1 - v2) < delta)
        return;

    fail(file, line, "mat1(%lu,%lu)=%f, mat2(%lu,%lu)=%f\n", row + 1, col + 1, v1, row + 1, col + 1, v2);
}

void do_assert_mat_equals(const char *file, int line, gsl_matrix *mat1, gsl_matrix *mat2, double delta)
{
    size_t row, col;

    if (mat1->size1 != mat2->size1 || mat1->size2 != mat2->size2)
        fail(file, line, "matrices are of different dimensions");

    for (row = 0; row < mat1->size1; row++) {
        for (col = 0; col < mat1->size2; col++) {
            assert_mat_elem(file, line, mat1, mat2, row, col, delta);
        }
    }
}

static void
assert_vec_elem(const char *file, int line, gsl_vector *vec1, gsl_vector *vec2, size_t i, double delta)
{
    double val1, val2;

    val1 = gsl_vector_get(vec1, i);
    val2 = gsl_vector_get(vec2, i);

    if (fabs(val1 - val2) < delta)
        return;

    fail(file, line, "vec1(%lu)=%f, vec2(%lu)=%f\n", i + 1, val1, i + 1, val2);
}

void do_assert_vec_equals(const char *file, int line, gsl_vector *vec1, gsl_vector *vec2, double delta)
{
    size_t i;

    if (vec1->size != vec2->size)
        fail(file, line, "vectors are of different lengths");

    for (i = 0; i < vec1->size; i++)
        assert_vec_elem(file, line, vec1, vec2, i, delta);
}
