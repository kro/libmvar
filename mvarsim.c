#include "mvardebug.h"
#include "mvardie.h"
#include "mvar.h"

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <stdbool.h>
#include <math.h>

/* A very small number used for comparing floating points. */
#define EPSILON 1e-12

/* The default number of skipped time steps. */
#define DEFAULT_NR_STEPS_SKIPPED 1e+3

struct mvar_simulation {
    /* The number of skipped time steps. */
    size_t nr_steps_skipped;
    /* The number of time steps of the simulated process. */
    size_t nr_steps;
    /* The total number of time steps. */
    size_t nr_total_steps;
    /* The random number generator to use for generating samples. */
    gsl_rng *rng;
    /* The dimension of the state vectors. */
    size_t m;
    /* The order of the process. */
    size_t p;
};

/* Returns true if one or more of the components of the vector are non-zero,
 * otherwise returns false. 
 */
static bool is_non_zero(gsl_vector *vec)
{
    size_t i;

    for (i = 0; i < vec->size; i++) {
        if (gsl_vector_get(vec, i) > EPSILON)
            return true;
    }

    return false;
}

/* Aggregates each column by summing each row. The result is stored to the
 * given vector. 
 */
static void sum_mat_cols(gsl_vector *vec, gsl_matrix *mat)
{
    double acc;
    size_t i, j;

    for (i = 0; i < mat->size2; i++) {
        acc = 0.0;
        for (j = 0; j < mat->size1; j++) {
            acc += gsl_matrix_get(mat, j, i);
        }
        gsl_vector_set(vec, i, acc);
    }
}

static void complex_modulo(gsl_vector *dst, gsl_vector_complex *src)
{
    gsl_complex z;
    size_t i;

    for (i = 0; i < src->size; i++) {
        z = gsl_vector_complex_get(src, i);
        gsl_vector_set(dst, i, sqrt(GSL_REAL(z) * GSL_REAL(z) + GSL_IMAG(z) * GSL_IMAG(z)));
    }
}

static void copy_matrix(gsl_matrix *dst, gsl_matrix *src)
{
    size_t row, col;

    for (row = 0; row < src->size1; row++) {
        for (col = 0; col < src->size2; col++) {
            gsl_matrix_set(dst, row, col, gsl_matrix_get(src, row, col));
        }
    }
}

/* Sets a unit diagonal to the given matrix starting from the specified row
 * offset by assuming that the row offset is the first row of the matrix. Note
 * that the row offset is zero based.
 */
static void set_identity(gsl_matrix *mat, size_t row_offset)
{
    size_t row, col;

    for (row = row_offset; row < mat->size1; row++) {
        for (col = 0; col < mat->size2; col++) {
            if (row - row_offset == col)
                gsl_matrix_set(mat, row, col, 1);
            else
                gsl_matrix_set(mat, row, col, 0);
        }
    }
}

/* Tests whether the specified multivariate AR model is stable. */
static bool is_mvar_stable(struct mvar_model *model)
{
    gsl_eigen_nonsymmv_workspace *workspace;
    gsl_matrix_complex *eigenvectors;
    gsl_vector_complex *eigenvalues;
    gsl_vector *eigenvalues_mod;
    int nr_row, nr_col;
    gsl_matrix *A1;
    bool stable;
    int m, p;
    size_t i;

    stable = true;
    m = model->A->size1;
    p = model->A->size2 / m;
    nr_row = nr_col = m * p;

    A1 = gsl_matrix_alloc(nr_row, nr_col);

    copy_matrix(A1, model->A);
    set_identity(A1, m);

    eigenvalues = gsl_vector_complex_alloc(nr_row);
    eigenvectors = gsl_matrix_complex_alloc(nr_row, nr_col);
    eigenvalues_mod = gsl_vector_alloc(nr_col);

    workspace = gsl_eigen_nonsymmv_alloc(nr_row);
    gsl_eigen_nonsymmv(A1, eigenvalues, eigenvectors, workspace);
    gsl_eigen_nonsymmv_free(workspace);
    complex_modulo(eigenvalues_mod, eigenvalues);

    /* If a modulo of an eigenvalue is greater than zero, the multivariate AR
     * model is not stable. 
     */
    for (i = 0; i < eigenvalues_mod->size; i++) {
        if (stable > 1.0) {
            stable = false;
            break;
        }
    }

    gsl_matrix_free(A1);
    gsl_vector_complex_free(eigenvalues);
    gsl_matrix_complex_free(eigenvectors);
    gsl_vector_free(eigenvalues_mod);

    return stable;
}

static gsl_matrix *cor_rv_mat(struct mvar_model *model, struct mvar_simulation *sim)
{
    gsl_matrix *U, *mat;
    gsl_vector *vec;
    size_t i, j;

    U = gsl_matrix_alloc(model->C->size1, model->C->size2);
    gsl_matrix_memcpy(U, model->C);
    gsl_linalg_cholesky_decomp(U);

    vec = gsl_vector_alloc(sim->p);
    mat = gsl_matrix_alloc(sim->nr_total_steps, sim->p);

    for (i = 0; i < sim->nr_total_steps; i++) {
        for (j = 0; j < sim->p; j++) {
            gsl_vector_set(vec, j, gsl_ran_ugaussian(sim->rng));
            gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, U, vec);
            gsl_blas_daxpy(1.0, model->w, vec);
        }
        gsl_matrix_set_row(mat, i, vec);
    }

    gsl_matrix_free(U);
    gsl_vector_free(vec);

    return mat;
}

/* Returns the parameter matrix of a multivariate AR(p) process. */
static gsl_matrix *param_mat_B(struct mvar_model *model, struct mvar_simulation *sim)
{
    gsl_matrix *B;
    size_t i;

    B = gsl_matrix_alloc(sim->m, sim->p);
    gsl_matrix_set_identity(B);

    for (i = 1; i <= sim->p; i++) {
        gsl_matrix_view v = gsl_matrix_submatrix(model->A, 0, sim->m * i - sim->m, sim->m, sim->p);
        gsl_matrix_sub(B, &v.matrix);
    }

    return B;
}

/* Returns the optimal forecast of the next state given the previous states.
 * The optimal forecast of the next state is calculated from the coefficient
 * matrix 'A' and intercept terms 'w'.
 */
static gsl_matrix *proc_mean_vec(struct mvar_model *model, struct mvar_simulation *sim)
{
    gsl_matrix *col_mat, *B, *x;
    gsl_permutation *perm;
    gsl_vector *proc_mean;
    gsl_matrix_view view;
    int sign;

    x = gsl_matrix_alloc(sim->p, sim->m);
    gsl_matrix_set_all(x, 0);

    if (!is_non_zero(model->w))
        return x;

    B = param_mat_B(model, sim);
    proc_mean = gsl_vector_alloc(sim->p);

    perm = gsl_permutation_alloc(sim->p);
    gsl_linalg_LU_decomp(B, perm, &sign);
    gsl_linalg_LU_solve(B, perm, model->w, proc_mean);

    col_mat = gsl_matrix_alloc(sim->p, 1);
    gsl_matrix_set_all(col_mat, 1.0);

    view = gsl_matrix_view_vector(proc_mean, 1, proc_mean->size);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, col_mat, &view.matrix, 0.0, x);

    gsl_vector_free(proc_mean);
    gsl_matrix_free(col_mat);
    gsl_permutation_free(perm);
    gsl_matrix_free(B);

    return x;
}

/* Returns the initial state vectors of the specified multivariate AR(p)
 * process. 
 */
static gsl_matrix *proc_state_vecs(struct mvar_model *model, struct mvar_simulation *sim)
{
    gsl_matrix *u, *x;
    size_t n;

    n = sim->p + sim->nr_steps + sim->nr_steps_skipped;
    x = proc_mean_vec(model, sim);

    u = gsl_matrix_alloc(n, sim->m);
    copy_matrix(u, x);
    gsl_matrix_free(x);

    return u;
}

static void init_simulation(struct mvar_simulation *sim, struct mvar_model *model, size_t nr_steps, gsl_rng *rng)
{
    sim->nr_steps_skipped = DEFAULT_NR_STEPS_SKIPPED;
    sim->nr_steps = nr_steps;
    sim->nr_total_steps = DEFAULT_NR_STEPS_SKIPPED + nr_steps;
    sim->m = model->A->size1;
    sim->p = model->A->size2 / sim->m;
    sim->rng = rng;
}

gsl_matrix *mvar_sim(struct mvar_model *model, size_t nr_steps, gsl_rng *rng)
{
    gsl_matrix *rv_mat, *v, *u, *x;
    gsl_vector *x_row, *x_sum_vec;
    struct mvar_simulation sim;
    gsl_vector_view vec_view;
    gsl_matrix_view mat_view;
    size_t i, j;

    if (!is_mvar_stable(model))
        mvar_die("the specified multivariate AR model is unstable");

    init_simulation(&sim, model, nr_steps, rng);
    rv_mat = cor_rv_mat(model, &sim);

    u = proc_state_vecs(model, &sim);
    x = gsl_matrix_alloc(sim.p, sim.m);
    x_row = gsl_vector_alloc(sim.p);
    x_sum_vec = gsl_vector_alloc(sim.p);

    for (i = sim.p + 1; i <= sim.nr_total_steps + sim.p; i++) {
        for (j = 1; j <= sim.p; j++) {

            mat_view = gsl_matrix_submatrix(model->A, 0, sim.m * j - sim.m, sim.m, sim.p);
            vec_view = gsl_matrix_row(u, i - j - 1);

            gsl_blas_dgemv(CblasNoTrans, 1.0, &mat_view.matrix, &vec_view.vector, 0.0, x_row);
            gsl_matrix_set_row(x, j - 1, x_row);
        }

        sum_mat_cols(x_sum_vec, x);
        vec_view = gsl_matrix_row(rv_mat, i - sim.p - 1);
        gsl_vector_add(x_sum_vec, &vec_view.vector);
        gsl_matrix_set_row(u, i - 1, x_sum_vec);
    }

    v = gsl_matrix_alloc(sim.nr_steps, sim.m);
    mat_view = gsl_matrix_submatrix(u, sim.nr_steps_skipped + sim.p, 0, sim.nr_steps, sim.m);
    gsl_matrix_memcpy(v, &mat_view.matrix);

    gsl_matrix_free(x);
    gsl_matrix_free(u);
    gsl_vector_free(x_row);
    gsl_matrix_free(rv_mat);
    gsl_vector_free(x_sum_vec);

    return v;
}
