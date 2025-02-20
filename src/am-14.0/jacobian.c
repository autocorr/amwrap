/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* jacobian.c                      S. Paine rev. 2024 June 27
*
* Functions for computing Jacobians of designated model
* spectra differentiated with respect to designated model
* parameters.
************************************************************/

#include <float.h>
#include <math.h>
#include <stdlib.h>

#include "am_types.h"
#include "errlog.h"
#include "jacobian.h"
#include "mapping.h"
#include "model.h"
#include "output.h"
#include "simplex.h"

static int compare_numeric(const void*, const void*);

/*
 * Jacobian derivatives are computed using a finite difference
 * approximation.  Three-point and five-point approximations are
 * computed for the same outer step size h.  The five-point
 * approximation is used as the result, and the difference
 * between the three-point and five-point approximations is used
 * to estimate the truncation and roundoff error in the
 * approximation.
 *
 * The starting point for choosing the outer step size h is the
 * characteristic scale supplied by the user for each simplex
 * variable.  For the three-point approximation, the optimum step
 * size is
 *
 *    h_opt ~= eps^(1/3) * scale,
 *
 * where eps is the numerical fractional accuracy with which a
 * given spectral point S is computed.  Below, we assume that eps
 * ~ DBL_EPSILON.  The justification is that although computation
 * of spectra generally involves a large number Nops of
 * operations, three factors limit the growth of numerical errors
 * in the finite differences:
 *
 *  1. Expressions evaluated within processor registers will
 *  typically carry guard bits beyond double precision.
 *
 *  2. Many roundoff errors will be common-mode between model
 *  computations involving a single perturbed model parameter.
 *
 *  3. Roundoff errors will accumulate roughly as the square root
 *  of the number of operations subject to rounding.
 *
 * References:
 *
 *  Press, et al., Numerical Recipes in C, 2nd ed., section 5.7.
 *  (Cambridge University Press, 1992)
 *
 *  Abramowitz, M. A. and Stegun, I. A., 9th printing, ch. 25.
 *  (Dover, 1970)
 */


/*
 * Coefficients for forward, centered, and backward finite
 * difference approximations for the derivative.
 *
 *  a_...[] - step sizes in units of an outer step size h
 * w3_...[] - weights for 3-point difference approximation
 * w5_...[] - weights for 5-point difference approximation
 * 
 * Note that the steps a[] of the differentiation variable can be
 * carried out in any order.  For the forward and backward
 * differences, which include a model computation at the
 * reference state, it makes sense to do that computation last so
 * that any subsequent model update at the reference state will
 * be as fast as possible.
 */

enum {
    DIFF_CENTERED,
    DIFF_FORWARD,
    DIFF_BACKWARD
};

static double  a_forward[] = {     2.0,    1.5,    1.0,    0.5,    0.0};
static double w3_forward[] = {    -0.5,    0.0,    2.0,    0.0,   -1.5};
static double w5_forward[] = {    -0.5,  8./3.,   -6.0,    8.0,-25./6.};

static double  a_centered[] = {   -1.0,   -0.5,    0.5,    1.0};
static double w3_centered[] = {   -0.5,    0.0,    0.0,    0.5};
static double w5_centered[] = {  1./6., -4./3.,  4./3., -1./6.};

static double  a_backward[] = {   -2.0,   -1.5,   -1.0,   -0.5,    0.0};
static double w3_backward[] = {    0.5,    0.0,   -2.0,    0.0,    1.5};
static double w5_backward[] = {    0.5, -8./3.,    6.0,   -8.0, 25./6.};


/***********************************************************
* int compute_jacobians(
*         model_t   *model,
*         model_t   *lmodel,
*         simplex_t *simplex)
*
* Purpose:
*   For each model output spectrum with the OUTPUT_JACOBIAN
*   flag bit set in the global output[] table, this function
*   computes the Jacobian with respect to each model
*   parameter addressed in the variable pointer table
*   *simplex->varptr[], at a reference state corresponding
*   to the current values of the model parameters in
*   vertex[0] of the simplex.  For each parameter, there is
*   a corresponding user-defined characteristic scale
*   simplex->scale[].  This characteristic scale is used to
*   determine the step size to be used for numerical
*   differentiation with respect to that parameter.
*   
*   For mapped variables, the differentiation steps are
*   carried out in the unmapped, physical domain.
*
*   On exit, *lmodel will reflect the last model computation
*   performed by this function, which may have involved a
*   perturbed model variable, whereas *model will have been
*   restored to the reference state corresponding to
*   vertex[0] of the simplex, but not necessarily updated at
*   that state.  This is to avoid a possible redundant
*   upstream model update at the reference state.
*
* Arguments:
*   model_t *model     - pointer to a model structure.
*   model_t *lmodel    - pointer to a model structure
*                        containing scalar data from a prior
*                        computation, or initialized to
*                        MODEL_INIT if no prior computation
*                        has occurred.
*   simplex_t *simplex - pointer to a simplex structure
*                        defining differentiation variables
*                        and their characteristic scales.
*
* Return:
*   0 on success, 1 on error
************************************************************/

int compute_jacobians(
        model_t   *model,
        model_t   *lmodel,
        simplex_t *simplex)
{
    unsigned int j;

    /*
     * If there are no differentiation variables defined, or no
     * Jacobians being computed, just return.
     */
    if (!simplex->n || !(output[ALL_OUTPUTS].flags & OUTPUT_JACOBIAN))
        return 0;
    /*
     * Set the model variables to vertex[0] of the simplex.
     */
    for (j = 0; j < simplex->n; ++j)
        *simplex->varptr[j] = simplex->logarithmic ?
            exp(simplex->vertex[0][j]) : simplex->vertex[0][j];
    /*
     * Loop over differentiation variables
     */
    for (j = 0; j < simplex->n; ++j) {
        double *a, *w3, *w5;
        double mvar, var, h;
        volatile double vtemp;
        unsigned int k, m, p;
        int diff_type;
        /*
         * The reference model state for the Jacobians
         * corresponds to the variable values at the 0th simplex
         * vertex.  Note that we do not use simplex->init for the
         * reference state because we might be computing
         * Jacobians following an initial downhill-simplex fit.
         *
         * Differentiation is carried out by steps in the
         * physical domain.  For a given variable, there are
         * potentially two unmappings needed to get from simplex
         * coordinates to physical ones: (1) simplex vertices are
         * stored as the log of the corresponding model variable
         * unless this default setting is explicitly turned off;
         * and (2) to facilitate this mapping some model
         * variables (vmr's and signed quantities) are mapped
         * from the physical domain to [0, inf].
         *
         * The differentiation step size is computed from the
         * characteristic scale supplied by the user.
         */
        mvar = simplex->logarithmic ?
                exp(simplex->vertex[0][j]) : simplex->vertex[0][j];
        var  = unmap_variable(mvar, simplex->mapping[j]);
        h    = unmap_differential(
                mvar, simplex->scale[j], simplex->mapping[j]);
        /*
         * If the current differentiation variable is a mixing
         * ratio, and the characteristic scale exceeds 1.0, it
         * will be clamped to 1.0 with a warning to the user.
         */
        if ((simplex->mapping[j] == MAPPING_VMR) && (h > 1.0)) {
            errlog(107, 0);
            h = 1.0;
        }
        h *= pow(DBL_EPSILON, 1./3.);
        /*
         * If the current differentiation variable is the LO
         * frequency for an IF spectrum, and the step size h is
         * greater than or equal to 1/3 of the frequency grid
         * interval, then clamp h with a warning to the user.
         * This will prevent the LO grid size changing between
         * differentiation steps.
         * 
         * The factor (1.0 - FLT_EPSILON) ensures that h will
         * still be less than hmax after adjustment to get
         * exactly-representable double-precision differences.
         */
        if (simplex->varptr[j] == &(model->flo)) {
            double hmax = (1./3.) * model->df * (1.0 - FLT_EPSILON);
            if (h > hmax) {
                h = hmax;
                errlog(133, 0);
            }
        }
        /*
         * Do an add and subtract to find the value nearest to h
         * such that the differences (var + n * h/2) - var are
         * exactly-representable as doubles.
         */ 
        vtemp = var + h / 2.;
        h     = 2. * (vtemp - var);
        /*
         * Select forward, centered, or backward differences
         * based on h, the variable value, and its physical
         * domain.  The physical domain corresponds to the domain
         * mapping.
         *
         * A special case occurs when the LO frequency for an IF
         * spectrum is a differentiation variable.  In that case,
         * the difference formula has to be chosen so as to avoid
         * a change in the LO frequency grid size between
         * difference steps.
         */
        if (simplex->varptr[j] == &(model->flo)) { /* diff variable is flo */
            double f_offset = fmod(model->flo, model->df);
            if (f_offset <= h)
                diff_type = DIFF_FORWARD;
            else if (f_offset >= 1.0 - h)
                diff_type = DIFF_BACKWARD;
            else
                diff_type = DIFF_CENTERED;
        } else if (simplex->mapping[j] == MAPPING_NONE) { /* [0, inf]     */
            if (var <= h)
                diff_type = DIFF_FORWARD;
            else
                diff_type = DIFF_CENTERED;
        } else if (simplex->mapping[j] == MAPPING_VMR) {  /* [0, 1]       */
            if (var <= h)
                diff_type = DIFF_FORWARD;
            else if (var >= 1.0 - h)
                diff_type = DIFF_BACKWARD;        
            else
                diff_type = DIFF_CENTERED;
        } else if (simplex->mapping[j] == MAPPING_EXP) {  /* [-inf, inf] */
            diff_type = DIFF_CENTERED;
        } else { 
            errlog(78, simplex->mapping[j]);
            return 1;
        }
        switch (diff_type) {
            case DIFF_FORWARD:
                a  = a_forward;
                w3 = w3_forward;
                w5 = w5_forward;
                m  = sizeof(a_forward) / sizeof(double);
                break;
            case DIFF_CENTERED:
                a  = a_centered;
                w3 = w3_centered;
                w5 = w5_centered;
                m  = sizeof(a_centered) / sizeof(double);
                break;
            case DIFF_BACKWARD:
                a  = a_backward;
                w3 = w3_backward;
                w5 = w5_backward;
                m  = sizeof(a_backward) / sizeof(double);
                break;
            default:
                errlog(128, 0);
                return 1;
        }
        /*
         * Loop over steps of the current differentiation
         * variable.
         */
        for (p = 0; p < m; ++p) {
            double w_deriv = w5[p] / h;
            /*
             * The truncation error for the 3-point derivative
             * scales as h^2, whereas the truncation error for
             * the 5-point derivative, which uses half the step
             * size and scales as h^4, should be at least 4 times
             * smaller.   So, the estimated truncation error for
             * the 5-point derivative is taken as 1/3 of the
             * difference between the 3-point and 5-point
             * estimates.  The weights for this difference are
             * computed directly; the 3-point derivative is never
             * computed on its own.
             */
            double w_trunc_err = (w3[p] - w5[p]) / (3. * h);
            *simplex->varptr[j] =
                map_variable(var + a[p] * h, simplex->mapping[j]);
            if (compute_model(model, lmodel)) {
                errlog(72, 0);
                return 1;
            };
            /*
             * Sum the output arrays for which Jacobians are
             * being computed into the Jacobian arrays with
             * appropriate weights.
             */
            for (k = 1; k < OUTPUT_END_OF_TABLE; ++k) {
                gridsize_t i;
                gridsize_t npts = (model->ifmode ? model->nif : model->ngrid);
                if (!(output[k].flags & OUTPUT_JACOBIAN))
                    continue;
                if (p == 0) { /* on the first step, initialize   */
                    for (i = 0; i < npts; ++i) {
                        output[k].jacobian[j][i] =
                            w_deriv     * output[k].spectrum[i];
                        output[k].jacobian_trunc_err[j][i] =
                            w_trunc_err * output[k].spectrum[i];
                    }
                } else {      /* on subsequent steps, accumulate */
                    for (i = 0; i < npts; ++i) {
                        output[k].jacobian[j][i] +=
                            w_deriv     * output[k].spectrum[i];
                        output[k].jacobian_trunc_err[j][i] +=
                            w_trunc_err * output[k].spectrum[i];
                    }
                }
            }
        }
        /*
         * Take the absolute value of the estimated truncation
         * error, and estimate the rounding error.  Note that
         * output[k].spectrum[i] is from the most
         * recently-computed model, which corresponds to a = 1
         * for centered differences, and to a = 0 for the forward
         * and backward differences.  The difference is
         * negligible for the purpose of this estimate.
         *
         * The sums of the estimated absolute truncation and
         * rounding errors go into a third array which is sorted
         * to support later reporting of percentile statistics of
         * the estimated total errors.
         */
        for (k = 1; k < OUTPUT_END_OF_TABLE; ++k) {
            gridsize_t i;
            gridsize_t npts = (model->ifmode ? model->nif : model->ngrid);
            double w_round_err = DBL_EPSILON / h;
            if (!(output[k].flags & OUTPUT_JACOBIAN))
                continue;
            for (i = 0; i < npts; ++i) {
                double trunc_err, round_err, abs_err, denom;
                trunc_err = fabs(output[k].jacobian_trunc_err[j][i]);
                round_err = w_round_err * fabs(output[k].spectrum[i]); 
                abs_err   = trunc_err + round_err;
                denom     = fabs(output[k].jacobian[j][i]);
                output[k].jacobian_trunc_err[j][i]  = trunc_err;
                output[k].jacobian_round_err[j][i]  = round_err;
                output[k].jacobian_sorted_err[j][i] =
                    denom > DBL_MIN ? abs_err / denom : 0.0;
            }
            qsort(output[k].jacobian_sorted_err[j],
                npts, sizeof(double), compare_numeric);
        }
        /*
         * Restore the current differentiation variable to the
         * reference state.
         */
        *simplex->varptr[j] = mvar;
    } /* end of loop over differentiation variables */
    return 0;
}

/************************************************************
* static int compare_numeric(const void *u1, const void *u2)
*
* Purpose:
*   Comparison function for numeric qsort().
*
* Return:
*   -1 if u1 < u2
*    0 if u1 == u2
*   +1 if u1 > u2
************************************************************/

static int compare_numeric(const void *u1, const void *u2)
{
    double v1 = *(double*)u1;
    double v2 = *(double*)u2;
    return (v1 < v2) ? -1 : ((v1 > v2) ? 1 : 0);
}   /* compare_numeric() */
