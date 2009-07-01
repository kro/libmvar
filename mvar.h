#ifndef MVAR_H
#define MVAR_H

#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

struct mvar_model {
    /* Coefficient matrix A. */
    gsl_matrix *A;
    /* Noise covariance matrix C. */
    gsl_matrix *C;
    /* Vector of intercept terms. */
    gsl_vector *w;
};

/* Estimates a multivariate AR(p) model using stepwise least-squares estimation
 * as described in Section 3. of Neumaier, A. and Schneider, T. (2001). The
 * fields of 'model' are populated with the estimated parameters, an Schwarz's
 * Bayesian Criterion is set to 'sbc'. If an error occurs, the process is
 * terminated after outputting the reason of failure to the standard error
 * output.
 */
void mvar_fit(gsl_matrix *v, int p, struct mvar_model *model, double *sbc);

/* Simulates a process following a multivariate AR(p) model. If an error
 * occurs, the process is terminated after outputting the reason of failure to
 * the standard error output.
 */
gsl_matrix *mvar_sim(struct mvar_model *model, size_t nr_steps, gsl_rng *rng);

/* Extracts the eigenvalues of a fitted AR(p) process. If an error occurs, the
 * process is terminated after outputting the reason of failure to the standard
 * error output.
 */
void mvar_eig(gsl_matrix *v);

#endif /* MVAR_H */
