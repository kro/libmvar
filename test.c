#include "mvar.h"
#include "mvardebug.h"
#include "mvartest.h"
#include "mvardie.h"

#include <math.h>

#define NR_STEPS 1e+4
#define MVAR_ORDER 2

static gsl_matrix *coef_mat_A(void)
{
    static double A[] = {
        0.40,  1.20,  0.30,  0.70,
        0.35, -0.30, -0.40, -0.50
    };

    static gsl_matrix_view view;

    view = gsl_matrix_view_array(A, MVAR_ORDER, MVAR_ORDER * 2);
    return &view.matrix;
}

static gsl_matrix *cov_noise_mat_C(void)
{
    static double C[] = {
        1.0, 0.50,
        0.5, 1.50
    };

    static gsl_matrix_view view;

    view = gsl_matrix_view_array(C, MVAR_ORDER, MVAR_ORDER);
    return &view.matrix;
}

static gsl_vector *intercept_terms(void)
{
    static double w[] = { 0.25, 0.1 };

    static gsl_vector_view view;

    view = gsl_vector_view_array(w, MVAR_ORDER);
    return &view.vector;
}

int main(int argc, char *argv[])
{
    struct mvar_model est_model;
    struct mvar_model model;
    gsl_matrix *v;
    gsl_rng *rng;
    double sbc;
    int p;

    model.A = coef_mat_A();
    model.C = cov_noise_mat_C();
    model.w = intercept_terms();

    /* Generate a process following a multivariate AR(2) model. */
    rng = gsl_rng_alloc(gsl_rng_ranlxs0);
    v = mvar_sim(&model, NR_STEPS, rng);

    /* Estimate the parameters of the generated process. */
    p = 2;
    est_model.A = gsl_matrix_alloc(MVAR_ORDER, MVAR_ORDER * 2);
    est_model.C = gsl_matrix_alloc(MVAR_ORDER, MVAR_ORDER);
    est_model.w = gsl_vector_alloc(MVAR_ORDER);

    mvar_fit(v, p, &est_model, &sbc);
    /* The intercept vector is proportional to the actual one. */
    gsl_vector_scale(est_model.w, 0.5);

    assert_vec_equals(est_model.w, model.w, 1e-1);
    assert_mat_equals(est_model.A, model.A, 1e-2);
    assert_mat_equals(est_model.C, model.C, 1e-2);

    gsl_matrix_free(est_model.A);
    gsl_matrix_free(est_model.C);
    gsl_vector_free(est_model.w);

    gsl_rng_free(rng);
    gsl_matrix_free(v);

    return 0;
}
