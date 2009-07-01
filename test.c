#include "mvar.h"
#include "mvardebug.h"

#define NR_STEPS 10
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
    struct mvar_model model;
    gsl_matrix *v;
    gsl_rng *rng;

    model.A = coef_mat_A();
    model.C = cov_noise_mat_C();
    model.w = intercept_terms();

    rng = gsl_rng_alloc(gsl_rng_ranlxs0);
    v = mvar_sim(&model, NR_STEPS, rng);

    mvar_print_matrix(v);

    gsl_rng_free(rng);
    gsl_matrix_free(v);

    return 0;
}
