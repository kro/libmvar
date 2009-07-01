#include "mvarmat.h"
#include "mvardie.h"
#include "mvar.h"

#include <gsl/gsl_machine.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <math.h>

struct mvar_fit {
    /* The process following a multivariate AR(p) model. */
    gsl_matrix *v;
    /* The number of block equations. */
    size_t nr_eqs;
    /* The number of parameter vectors. */
    size_t nr_params;
    /* The number of state vectors. */
    size_t n;
    /* The dimension of the state vectors. */
    size_t m;
    /* The order of the process. */
    size_t p;
};

static void set_data_mat_K_intercept_col(struct mvar_fit *fit, gsl_matrix *K)
{
    gsl_vector *vec;

    vec = gsl_vector_alloc(fit->nr_eqs);
    gsl_vector_set_all(vec, 1.0);
    gsl_matrix_set_col(K, 0, vec);
    gsl_vector_free(vec);
}

static void set_data_mat_K_predictors(struct mvar_fit *fit, gsl_matrix *K)
{
    gsl_matrix_view mat_view;
    gsl_vector_view vec_view;
    size_t i, j, k;

    for (i = 1; i <= fit->p; i++) {
        mat_view = gsl_matrix_submatrix(fit->v, fit->p - i, 0, fit->n - fit->p, fit->m);
        for (j = 1 + fit->m * (i - 1), k = 0; j < 1 + fit->m * i; j++, k++) {
            vec_view = gsl_matrix_column(&mat_view.matrix, k);
            gsl_matrix_set_col(K, j, &vec_view.vector);
        }
    }
}

static void set_data_mat_K_observations(struct mvar_fit *fit, gsl_matrix *K)
{
    gsl_matrix_view mat_view;
    gsl_vector_view vec_view;
    size_t i, j;

    mat_view = gsl_matrix_submatrix(fit->v, fit->p, 0, fit->n - fit->p, fit->m);
    for (i = fit->nr_params, j = 0; i < K->size2; i++, j++) {
        vec_view = gsl_matrix_column(&mat_view.matrix, j);
        gsl_matrix_set_col(K, i, &vec_view.vector);
    }
}

/* Creates the data matrix K for QR decomposition. */
static gsl_matrix *data_mat_K(struct mvar_fit *fit)
{
    gsl_matrix *K;

    K =  gsl_matrix_alloc(fit->nr_eqs, fit->nr_params + fit->m);
    set_data_mat_K_intercept_col(fit, K);
    set_data_mat_K_predictors(fit, K);
    set_data_mat_K_observations(fit, K);

    return K;
}

/* Returns the R matrix of a QR factorization. */
static gsl_matrix *qr_fact(struct mvar_fit *fit, gsl_vector *scale)
{
    gsl_matrix *R, *K, *K_diag;
    gsl_vector *tau;
    double delta;

    K = data_mat_K(fit);
    delta = (pow(K->size2, 2) + K->size2 + 1) * GSL_DBL_EPSILON;
    mvar_mat_sum_sq_sqrt(K, scale);
    gsl_vector_scale(scale, sqrt(delta));

    K_diag = gsl_matrix_alloc(scale->size, scale->size);
    gsl_matrix_set_all(K_diag, 0.0);
    mvar_mat_set_diag(K_diag, scale);

    /* Combine the rows of K and K_diag into one big matrix R, which is then QR decomposed. */
    R = gsl_matrix_alloc(K->size1 + scale->size, K->size2);
    gsl_matrix_set_all(R, 0.0);
    mvar_mat_copy(R, K, 0, 0);
    mvar_mat_copy(R, K_diag, K->size1, 0);

    tau = gsl_vector_alloc(R->size2);
    gsl_linalg_QR_decomp(R, tau);
    mvar_mat_upper_tri(R);

    gsl_matrix_free(K);
    gsl_vector_free(tau);
    gsl_matrix_free(K_diag);

    return R;
}

static gsl_matrix *augmented_param_mat_A(gsl_matrix *R11, gsl_matrix *R12)
{
    gsl_matrix *t_aug_A, *aug_A, *LU, *R11_inv;
    gsl_permutation *p;
    int signum;

    p = gsl_permutation_alloc(R11->size2);
    LU = gsl_matrix_alloc(R11->size1, R11->size2);
    gsl_matrix_memcpy(LU, R11);

    R11_inv = gsl_matrix_alloc(R11->size1, R11->size2);
    gsl_linalg_LU_decomp(LU, p, &signum);
    gsl_linalg_LU_invert(LU, p, R11_inv);

    aug_A = gsl_matrix_alloc(R12->size1, R12->size2);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, R11_inv, R12, 0.0, aug_A);

    t_aug_A = gsl_matrix_alloc(aug_A->size2, aug_A->size1);
    gsl_matrix_transpose_memcpy(t_aug_A, aug_A);

    gsl_permutation_free(p);
    gsl_matrix_free(LU);
    gsl_matrix_free(aug_A);
    gsl_matrix_free(R11_inv);

    return t_aug_A;
}

static gsl_matrix_view mat_view_R11(struct mvar_fit *fit, gsl_matrix *R)
{
    return gsl_matrix_submatrix(R, 0, 0, fit->nr_params, fit->nr_params);
}

static gsl_matrix_view mat_view_R12(struct mvar_fit *fit, gsl_matrix *R)
{
    return gsl_matrix_submatrix(R, 0, fit->nr_params, fit->nr_params, fit->m);
}

static gsl_matrix_view mat_view_R22(struct mvar_fit *fit, gsl_matrix *R)
{
    return gsl_matrix_submatrix(R, fit->nr_params, fit->nr_params, fit->m, fit->m);
}

static double scale_factor(gsl_vector *scale)
{
    double scale_factor;
    size_t i;

    scale_factor = 0.0;
    for (i = 1; i < scale->size; i++)
        scale_factor = GSL_MAX(scale_factor, gsl_vector_get(scale, i));

    return scale_factor / gsl_vector_get(scale, 0);
}

static void rescale_R11(gsl_matrix *R11, gsl_vector *scale)
{
    gsl_vector_view vec_view = gsl_matrix_column(R11, 0);
    gsl_vector_scale(&vec_view.vector, scale_factor(scale));
}

static void init_mvar_fit(gsl_matrix *v, size_t p, struct mvar_fit *fit)
{
    fit->v = v;
    fit->p = p;

    fit->n = fit->v->size1;
    fit->m = fit->v->size2;

    fit->nr_eqs = fit->n - fit->p;
    fit->nr_params = fit->m * fit->p + 1;

    if (fit->nr_eqs <= fit->nr_params)
        mvar_die("time-series too short");
}

static void set_intercept_vec_w(struct mvar_model *model, gsl_matrix *aug_A, gsl_vector *scale)
{
    gsl_vector_view vec_view = gsl_matrix_column(aug_A, 0);

    gsl_vector_memcpy(model->w, &vec_view.vector);
    gsl_vector_scale(model->w, scale_factor(scale));
}

static void set_coef_mat_A(struct mvar_model *model, struct mvar_fit *fit, gsl_matrix *aug_A)
{
    gsl_matrix_view mat_view = gsl_matrix_submatrix(aug_A, 0, 1, aug_A->size1, aug_A->size2 - 1);
    gsl_matrix_memcpy(model->A, &mat_view.matrix);
}

static void set_noise_cov_mat_C(struct mvar_model *model, struct mvar_fit *fit, gsl_matrix *R22)
{
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, R22, R22, 0.0, model->C);
    gsl_matrix_scale(model->C, 1.0 / (fit->nr_eqs - fit->nr_params));
}

void mvar_fit(gsl_matrix *v, size_t p, struct mvar_model *model, double *sbc)
{
    gsl_matrix_view R11v, R12v, R22v;
    gsl_matrix *aug_A, *R, *R11;
    struct mvar_fit fit;
    gsl_vector *scale;

    init_mvar_fit(v, p, &fit);

    scale = gsl_vector_alloc(fit.nr_params + fit.m);
    R = qr_fact(&fit, scale);

    R11v = mat_view_R11(&fit, R);
    R12v = mat_view_R12(&fit, R);
    R22v = mat_view_R22(&fit, R);

    /* Make a copy of the R11 for scaling the first column. */
    R11 = gsl_matrix_alloc(R11v.matrix.size1, R11v.matrix.size2);
    gsl_matrix_memcpy(R11, &R11v.matrix);
    rescale_R11(R11, scale);
    aug_A = augmented_param_mat_A(R11, &R12v.matrix);

    set_intercept_vec_w(model, aug_A, scale);
    set_coef_mat_A(model, &fit, aug_A);
    set_noise_cov_mat_C(model, &fit, &R22v.matrix);

    gsl_matrix_free(aug_A);
    gsl_vector_free(scale);
    gsl_matrix_free(R11);
    gsl_matrix_free(R);
}
