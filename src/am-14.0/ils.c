/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* ils.c                             S. Paine rev. 2024 May 9
*
* Functions for computation of the instrumental line shape,
* and for convolution of the ILS with spectra.
************************************************************/

#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "am_types.h"
#include "errlog.h"
#include "ils.h"
#include "math_const.h"
#include "output.h"
#include "phys_const.h"
#include "specfunc.h"
#include "transform.h"

/*
 * Global array of instrumental line shape types.
 *
 * The Norton-Beer lineshapes are associated with a set of apodizing
 * functions presented in R.H Norton and R. Beer, "New apodizing
 * functions for  Fourier spectrometry," J. Opt. Soc. v. 66 p. 259
 * (March 1976).  See also the erratum in v. 67 p. 419 (March 1977).
 *
 * The second column in the table below is the FWHM of the ILS in
 * normal units.  In normal units, all of the ILS functions are
 * normalized to pi.
 */

struct ils_type_tabentry ils_type[] = {
{"none",           0.0},
{"sinc",           3.790988534067962},
{"sinc_squared",   5.566229513006039},
{"Bartlett",       5.566229513006039},   /* synonym for sinc_squared */
{"Hann",           6.283185307179586},
{"Hamming",        5.702697332538719},
{"Blackman",       7.221892673480022},
{"Norton-Beer_I1", 4.549284013281266},
{"Norton-Beer_I2", 5.307384112084631},
{"Norton-Beer_I3", 6.065555515678795},
{"Gaussian",       1.665109222315395},
{"rectangle",      1.0},
{"END",            0.0}
};

static int    ils_convolve(model_t*, double*);
static double ilsfunc(int, double, double, double);
static double impulse(double, double);
static double rect(double, double);
static double Q0(double);
static double Q1(double);
static double Q2(double);
static double Q4(double);


/***********************************************************
* int apply_instrumental_line_shape(
*   model_t *model, model_t *lmodel)
*
* Purpose:
*   Convolves the ILS with all active outputs to which it
*   to be applied.
*
* Arguments:
*   model_t *model - pointer to model structure
*   model_t *lmodel - pointer to a model data structure
*       containing scalar data from a prior computation.
*
* Return:
*   0 on success, 1 on failure.
************************************************************/

int apply_instrumental_line_shape(model_t *model, model_t *lmodel)
{
    int i;

    /*
     * Recompute the instrumental line shape, if ILS parameters
     * have changed.
     */
    if (lmodel != NULL) {
        if ((model->ils_typenum != lmodel->ils_typenum) ||
            (fabs(model->ils_fwhm - lmodel->ils_fwhm)
                > (DBL_EPSILON * model->ils_fwhm)) ||
            (fabs(model->ils_fif - lmodel->ils_fif)
                > (DBL_EPSILON * model->ils_fif)) ||
            (model->ilsmode != lmodel->ilsmode) ||
            (
                (model->ilsmode & ILSMODE_DSB) &&
                (fabs(model->dsb_utol_ratio - lmodel->dsb_utol_ratio)
                > (DBL_EPSILON * model->dsb_utol_ratio)))
            ) {
            initialize_ils(model);
        }
    }
    /*
     * Convolve active output arrays with the ILS.
     */
    for (i = 0; i < OUTPUT_END_OF_TABLE; ++i) {
        if ((output[i].flags & OUTPUT_ACTIVE)
            && (output[i].flags & ILS_APPLIED)) {
            if (ils_convolve(model, output[i].spectrum))
                return 1;
        }
    }
    return 0;
}   /* apply_instrumental_line_shape() */


/***********************************************************
* static int ils_convolve(model_t *model, double *s)
*
* Purpose:
*   Convolves a model spectrum s[] with the instrumental lineshape.
*   Convolution is done by multiplication of Hartley transforms,
*   taking advantage of the symmetry of the ILS.
*
* Arguments:
*   model_t *model - pointer to model data structure
*   double *s - spectrum to be convolved
*
* Return:
*   0 if OK, 1 on error.
************************************************************/

static int ils_convolve(model_t *model, double *s)
{
    gridsize_t i, mid;

    if (model->ilsworkspace == NULL) {
        errlog(47, 0);
        return 1;
    }
    /*
     * The spectrum s[] is padded, in wrap-around order, by extending its
     * terminal values to the midpoint of the padding area.  Since the padded
     * array is at least three times as long as the frequency grid, the
     * discontinuity at the midpoint of the padding area will not contaminate
     * the convolution.  See the comment in initialize_ils() for a diagram
     * of the array layouts.
     */
    for (i = 0; i < model->ngrid; ++i)
        model->ilsworkspace[i] = s[i];
    mid = (model->ngrid + model->npad) / 2;
    for (i = model->ngrid; i < mid; ++i)
        model->ilsworkspace[i] = s[model->ngrid - 1];
    for (i = mid; i < model->npad; ++i)
        model->ilsworkspace[i] = s[0];
    /*
     * Convolution by Hartley transform.  The Hartley convolution theorem is:
     *
     * H(a cnvl b) = 0.5 * [H(a)H(b) - _H(a)_H(b) + _H(a)H(b) + H(a)_H(b)],
     *
     * where _x is the reverse of sequence x, _x[0] = x[0], _x[1] = x[n-1], etc.
     * If the ILS is symmetric, this expression reduces to a straight element-
     * by-element multiplication in the transform domain.  Otherwise, it is
     * necessary to include all terms, taking into account the bit-reversed
     * ordering in the transform domain.
     */
    fht_dif(model->ilsworkspace, (unsigned long)model->npad);
    if ((model->ilsmode & ILSMODE_NORMAL) ||
        ((model->ilsmode & ILSMODE_DSB) && (fabs(model->dsb_utol_ratio - 1.0)
        <= DBL_EPSILON))
        ) {
        for (i = 0; i < model->npad; ++i)
            model->ilsworkspace[i] *= model->ils[i];
    } else {
        double wj, wk, sjk, djk;
        gridsize_t n, j, k;
        model->ilsworkspace[0] *= model->ils[0];
        if (model->npad > 1)
            model->ilsworkspace[1] *= model->ils[1];
        n = 2;
        while (n < model->npad) {
            for (j = n, k = 2*n-1; j < k; ++j, --k) {
                wj = model->ilsworkspace[j];
                wk = model->ilsworkspace[k];
                sjk = model->ils[j] + model->ils[k];
                djk = model->ils[j] - model->ils[k];
                model->ilsworkspace[j] = 0.5 * (wj * sjk + wk * djk);
                model->ilsworkspace[k] = 0.5 * (wk * sjk - wj * djk);
            }
            n *= 2;
        }
    }
    fht_dit(model->ilsworkspace, (unsigned long)model->npad);
    /*
     * Copy the spectrum back to s[].
     */
    for (i = 0; i < model->ngrid; ++i)
        s[i] = model->ilsworkspace[i];
    return 0;
}   /* ils_convolve() */


/***********************************************************
* static double ilsfunc(int num, double fwhm, double f, double df)
*
* Purpose:
*   Computes the amplitude of an unormalized instrumental
*   lineshape, given the function type number, num; the full
*   width at half maximum, fwhm; the frequency, f; and the
*   frequency grid interval df.
*
* Arguments:
*   int num - ils function type number
*   double fwhm - full width at half maximum, must be > 0
*   double f - frequency, in same units as fwhm
*   double df - frequency grid interval, in same units as fwhm
*
* Return:
*   value of the ils at f, as a double
************************************************************/

static double ilsfunc(int num, double fwhm, double f, double df)
{
    double fscale, a, da, q;

    fscale = ils_type[num].fwhm / fwhm;
    a = f * fscale;
    da = df * fscale;
    /*
     * If the ILS FWHM is less than twice the frequency grid spacing,
     * return the impulse function.  (The function initialize_ils()
     * will have logged a warning to the user.)
     */
    if (fwhm < 2.0 * df) {
        return impulse(a, da);
    }
    switch (num) {
    case ILS_SINC:
        return Q0(a);
    case ILS_SINC_SQUARED:
    case ILS_BARTLETT:
        q = Q0(0.5 * a);
        return 0.5 * q * q;
    case ILS_HANN:
        return 0.50 * Q0(a) + 0.25 * (Q0(a + PI) + Q0(a - PI));
    case ILS_HAMMING:
        return 0.54 * Q0(a) + 0.23 * (Q0(a + PI) + Q0(a - PI));
    case ILS_BLACKMAN:
        return 0.42 * Q0(a) + 0.25 * (Q0(a + PI) + Q0(a - PI))
            + 0.04 * (Q0(a + TWOPI) + Q0(a - TWOPI));
    case ILS_NORTON_BEER_I1:
        return 0.384093 * Q0(a) - 0.087577 * Q1(a) + 0.703484 * Q2(a);
    case ILS_NORTON_BEER_I2:
        return 0.152442 * Q0(a) - 0.136176 * Q1(a) + 0.983734 * Q2(a);
    case ILS_NORTON_BEER_I3:
        return 0.045335 * Q0(a) + 0.554883 * Q2(a) + 0.399782 * Q4(a);
    case ILS_GAUSSIAN:
        return SQRT_PI * am_exp(-a * a);
    case ILS_RECTANGLE:
        return rect(a, da);
    default: /* log an error and return zero */
        errlog(27, num);
        return 0.0;
    }
}   /* ilsfunc() */


/***********************************************************
* int initialize_ils(model_t *model)
*
* Purpose:
*   This function computes the instrumental line shape, and
*   takes its Hartley transform.  The result is stored in
*   bit-reversed order.
*
*   On the first call to initialize_ils(), memory is allocated
*   for the arrays model->ils[] and model->ilsworkspace[].
*
* Arguments:
*   model_t *model - pointer to model data structure
*
* Return:
*   1 on error, 0 otherwise
************************************************************/

int initialize_ils(model_t *model)
{
    gridsize_t i, j;
    double fwhm, ilsnorm;

    /*
     * If the ILS is too broad, or too narrow, log a warning.
     */
    if ((model->ils_fwhm + fabs(model->ils_fif)) > (model->fmax - model->fmin))
        errlog(77, 0);
    if (model->ils_fwhm < (2.0 * model->df))
        errlog(83,0);
    /*
     * To avoid division by zero, the FWHM is forced to be at least F_EPSILON.
     */
    if (model->ils_fwhm > F_EPSILON) {
        fwhm = model->ils_fwhm;
    } else {
        fwhm = F_EPSILON;
    }
    /*
     * If model->ils is NULL, then memory needs to be allocated.  Likewise
     * for model->ilsworkspace, which holds the padded spectrum to be convolved
     * with the ils in ils_convolve().  Memory allocated here will be freed at
     * program exit by a call to free_model_entities().
     */
    if (model->ils == NULL) {
        if ((model->ils = (double*)malloc(model->npad * sizeof(double)))
            == NULL) {
            errlog(30, 0);
            return 1;
        }
    }
    if (model->ilsworkspace == NULL) {
        if ((model->ilsworkspace
            = (double*)malloc(model->npad * sizeof(double))) == NULL) {
            errlog(30, 0);
            return 1;
        }
    }
    /*
     * For convolution, spectra are padded to at least three times the length
     * of the frequency grid by extending their terminal values.  The ILS is
     * computed over twice the length of the frequency grid, and zeroed
     * elsewhere.  The array layouts look like this:
     *
     * Spectrum: axxxxxxxxbbbbbbbbbbbaaaaaaaaaa
     *      ILS: xxxxxxxxxx0000000000xxxxxxxxxx
     *           ^        ^                   ^
     *    index: 0        ngrid-1             npad-1
     *
     * This guarantees that the cyclic convolution of these two arrays will be
     * free of truncation artifacts over the range 0..ngrid-1 of the model
     * frequency grid.
     *
     * The ILS is normalized such that its integral over frequency is 1/npad,
     * to account for the normalization of the discrete Hartley transforms in
     * the convolutions.
     */
    ilsnorm = model->df * (ils_type[model->ils_typenum].fwhm / fwhm);
    ilsnorm /= PI * model->npad;
    if (model->ilsmode & ILSMODE_NORMAL) {
        /*
         * The normal ILS is symmetrical about f = 0.  Fill it in two points
         * at a time.
         */
        double f;
        model->ils[0] =
            ilsnorm * ilsfunc(model->ils_typenum, fwhm, 0.0, model->df);
        for (i = 1, j = model->npad - 1 ; i < model->ngrid; ++i, --j) {
            f = (double)i * model->df;
            model->ils[i] =
                ilsnorm * ilsfunc(model->ils_typenum, fwhm, f, model->df);
            model->ils[j] = model->ils[i];
        }
    } else {
        /*
         * For LSB, USB, and DSB modes, start by computing the LSB-shifted
         * ILS.  (Note that, since spectra are to be convolved with the ILS
         * reflected about f=0, the LSB ILS is shifted _upwards_ in frequency.)
         * For DSB, sum the LSB with its reflection about f=0, weighted
         * according to the sideband ratio.  For USB, just reflect the LSB.
         */
        double f;
        model->ils[0] = ilsnorm *
            ilsfunc(model->ils_typenum, fwhm, -model->ils_fif, model->df);
        for (i = 1, j = model->npad - 1; i < model->ngrid; ++i, --j) {
            f = (double)i * model->df;
            model->ils[i] = ilsnorm * ilsfunc(
                model->ils_typenum, fwhm, f - model->ils_fif, model->df);
            model->ils[j] = ilsnorm * ilsfunc(
                model->ils_typenum, fwhm, -f - model->ils_fif, model->df);
        }
        if (model->ilsmode & ILSMODE_DSB) {
            double t, wtu, wtl;
            wtu = model->dsb_utol_ratio / (1.0 + model->dsb_utol_ratio);
            wtl = 1.0 - wtu;
            for (i = 1, j = model->npad - 1; i < model->ngrid; ++i, --j) {
                t = model->ils[i];
                model->ils[i] = wtl * model->ils[i] + wtu * model->ils[j];
                model->ils[j] = wtl * model->ils[j] + wtu * t;
            }
        } else if (model->ilsmode & ILSMODE_USB) {
            double t;
            for (i = 1, j = model->npad - 1; i < model->ngrid; ++i, --j) {
                t = model->ils[i];
                model->ils[i] = model->ils[j];
                model->ils[j] = t;
            }
        }
    }
    /*
     * Zero padding.
     */
    j = model->npad - model->ngrid;
    for (i = model->ngrid; i <= j; ++i)
        model->ils[i] = 0.0;
    /*
     * Compute the Hartley transform of the ILS.  This gets stored in
     * bit-reversed order to save a little time in the convolutions.
     */
    fht_dif(model->ils, (unsigned long)model->npad);
    return 0;
}   /* initialize_ils() */


/***********************************************************
* static double impulse(double a, double da)
*
* Purpose:
*   Computes an impulse function consisting of a triangle
*   with base width equal to twice the frequency grid
*   interval da, and area equal to pi.  If the impulse is
*   aligned to the frequency grid, the function has a
*   single non-zero point at a = 0, with amplitude pi/da.
*   Otherwise, there are two non-zero points bracketing
*   a = 0, whose amplitudes sum to pi/da.
*
* Arguments:
*   double a - frequency at which the impulse function is to
*    be evaluated.
*   double da - frequency grid interval.
*
* Return:
*   impulse(a) as a double
************************************************************/

static double impulse(double a, double da)
{
    double u, v;
    u = fabs(a);
    if (u >= da)
        return 0.0;
    v = 1.0 / da;
    return (PI * v * (1.0 - (u * v)));
}   /* impulse() */


/***********************************************************
* static double rect(double a, double da)
*
* Purpose:
*   Computes a rectangle function of unit width, centered on
*   a = 0, and with area equal to pi.  The edges are linear
*   ramps of width 2 * da, centered on a = +/- 0.5.
*
*   If da > 0.5, the  rect() function collapses to an
*   impulse() function.
*
* Arguments:
*   double a - frequency at which the rect function is to
*    be evaluated.
*   double da - frequency grid interval.
*
* Return:
*   rect(a) as a double
************************************************************/

static double rect(double a, double da)
{
    double u;
    u = fabs(a);
    if (da > 0.5)
        return impulse(a, da);
    if (u >= (0.5 + da))
        return 0.0;
    if (u <= (0.5 - da))
        return PI;
    return PI * (0.5 - ((u - 0.5) / (2.0 * da)));
}   /* rect() */


/***********************************************************
* static double Q0(double a)
* static double Q1(double a)
* static double Q2(double a)
* static double Q4(double a)
*
* Purpose:
*   Auxiliary functions for computing instrumental line shape
*   functions.  Qn is the Fourier transform of the time-domain
*   window function (1 - U^2)^n, where U = tau / Delta,
*   tau = lag, Delta = max. lag.
*
* Arguments:
*   double a - normalized frequency offset, a = 2*pi*Delta*f.
*
* Return:
*   Qn(a) as a double
************************************************************/

static double Q0(double a)
{
    if (fabs(a) < 0.7) {
        double a2 = a * a;
        return 1. +
            a2 * (-1.666666666666667e-1 +
            a2 * ( 8.333333333333333e-3 +
            a2 * (-1.984126984126984e-4 +
            a2 * ( 2.755731922398589e-6))));
    } else {
        return sin(a) / a;
    }
}   /* Q0() */


static double Q1(double a)
{
    double a2 = a * a;
    if (fabs(a) < 0.7) {
        return 6.666666666666667e-1 +
            a2 * (-6.666666666666667e-2 +
            a2 * ( 2.380952380952381e-3 +
            a2 * (-4.409171075837743e-5 +
            a2 * ( 5.010421677088344e-7))));
    } else {
        return 2. * ((sin(a) / a) - cos(a)) / a2;
    }
}   /* Q1() */


static double Q2(double a)
{
    double a2 = a * a;
    if (fabs(a) < 0.7) {
        return 5.333333333333333e-1 +
            a2 * (-3.809523809523809e-2 +
            a2 * ( 1.058201058201058e-3 +
            a2 * (-1.603334936668270e-5 +
            a2 * ( 1.541668208334875e-7))));
    } else {
        double r = 3. / a2;
        return -8. * ((1. - r) * (sin(a) / a) + r * cos(a)) / a2;
    }
}   /* Q2() */


static double Q4(double a)
{
    double a2 = a * a;
    if (fabs(a) < 0.7) {
        return 4.063492063492064e-1 +
            a2 * (-1.847041847041847e-2 +
            a2 * ( 3.552003552003552e-4 +
            a2 * (-3.946670613337280e-6 +
            a2 * ( 2.901963686277412e-8))));
    } else {
        double a4 = a2 * a2;
        double r;
        r = (1. - 45. / a2 + 105. / a4) * (sin(a) / a);
        r += (5. / a2) * (2. - 21. / a2) * cos(a);
        r *= 384. / a4;
        return r;
    }
}   /* Q4() */

#ifdef UNIT_TEST

/***********************************************************
* int main(int argc, char **argv)
*
* Purpose:
*   Tests ilsfunc(), evaluating the ILS functions in the
*   table ils_type for num_points points over a specified
*   range [amin, amax].  For each ILS function, the output
*   is written to a file ils_name.out, with three columns
*   for a, ils(a), and ils(a) / ils(0).
*
* Usage:
*   $ a.out amin amax num_points
*
* Example:
*   $ gcc -D UNIT_TEST ils.c
*   $ ./a.out 0 18.8496 1000
************************************************************/

#define BUFSIZE 32

int main(int argc, char **argv)
{
    int i, n;
    double amin, amax;
    char buf[BUFSIZE];
    FILE *fp;

    if (argc < 4) {
        fprintf(stderr, "ils.c unit test\n");
        fprintf(stderr, "usage:\n");
        fprintf(stderr, "  %s amin amax npoints\n", argv[0]);
        fprintf(stderr, "example:\n");
        fprintf(stderr, "  ./a.out 0 18.8496 1000\n");
        return 1;
    }
    amin = atof(argv[1]);
    amax = atof(argv[2]);
    n    = atoi(argv[3]);
    for (i = 1; i < ILS_END_OF_TABLE; ++i) {
        double da = (amax - amin) / (n - 1);
        double y0 = ilsfunc(i, ils_type[i].fwhm, 0.0, da);
        int j;
        snprintf(buf, sizeof(buf), "%s.out", ils_type[i].name);
        fp = fopen(buf, "w");

        for (j = 0; j < n - 1; ++j) {
            double a = amin + da * j;
            double y = ilsfunc(i, ils_type[i].fwhm, a, da);
            fprintf(fp, "%e %e %e\n", a, y, y / y0);
        }
        fclose(fp);
    }
    return 0;
} /* main() */


/***********************************************************
* int errlog(const int errnum, const int data)
*
* Purpose:
*   dummy errlog() function for unit test
************************************************************/

int errlog(const int errnum, const int data)
{
    printf("errlog(%d, %d)\n", errnum, data);
    return 0;
}   /* errlog() */


/***********************************************************
* void fht_dif(double *x, unsigned long n)
*
* Purpose:
*   dummy fht_dif() for unit test.
************************************************************/

void fht_dif(double *x, unsigned long n)
{
    return;
}   /* fht_dif() */


/***********************************************************
* void fht_dit(double *x, unsigned long n)
*
* Purpose:
*   dummy fht_dit() for unit test.
************************************************************/

void fht_dit(double *x, unsigned long n)
{
    return;
}   /* fht_dit() */


#endif /* UNIT_TEST */
