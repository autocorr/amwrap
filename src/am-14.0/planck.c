/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* planck.c                     S. Paine rev. 2023 December 1 
*
* Planck functions.
************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "am_sysdep.h"
#include "am_types.h"
#include "phys_const.h"
#include "planck.h"
#include "specfunc.h"


/***********************************************************
* void B_Planck(double *B, double *f, double df, gridsize_t ngrid, double T)
*
* Purpose:
*   Computes the Planck radiance spectrum B[](T) over the
*   frequency grid f[].  f[] must be in ascending frequency
*   order.
*
* Arguments:
*   double *B - computed radiance [watt / (cm^2 * GHz * sr)]
*   double *f - frequency [GHz]
*   double df - frequency grid interval [GHz]
*   ngrid - number of frequency grid points
*   double T - temperature [K]
************************************************************/

void B_Planck(double *B, double *f, double df, gridsize_t ngrid, double T)
{
    gridsize_t i;

    if (T < DBL_MIN) {
        for (i = 0; i < ngrid; ++i)
            B[i] = 0.0;
    } else {
        double r1, f_min;
        gridsize_t imin;
        r1 = H_ON_KB / T;
        /*
         * Below f[imin], a low-frequency limit of the Planck function
         * is used to avoid numeric trouble.
         */
        f_min = FLT_EPSILON / r1;
        if (f_min < f[0]) {
            imin = (gridsize_t)0;
        } else {
            imin = (gridsize_t)ceil((f_min - f[0]) / df);
            imin = (imin < ngrid) ? imin : ngrid;
        }
        for (i = 0; i < imin; ++i) {
            double r2;
            r2 = TWOKB_ON_CSQUARED * T;
            B[i] = f[i] * f[i] * r2;
        }
        #ifdef _OPENMP
        #pragma omp parallel for schedule(guided, 128) if (ngrid >= 16384)
        #endif
        for (i = imin; i < ngrid; ++i) {
            double r2;
            r2 = r1 * f[i];
            B[i] = TWOH_ON_CSQUARED * f[i] * f[i] * f[i] / (exp(r2) - 1.0);
        }
    }
    return;
}   /* B_Planck() */


/***********************************************************
* double T_Planck(double f, double I)
*
* Purpose:
*   Computes Planck brightness temperature as a function of
*   frequency f and radiance I.
*
* Arguments:
*   double f - frequency [GHz]
*   double I - radiance [watt / (cm^2 * GHz * sr)]
*
* Return:
*   Planck brightness temperature [K]
************************************************************/

double T_Planck(double f, double I)
{
    if (I < DBL_MIN)
        return 0.0;
    return H_ON_KB * f / log(1. + TWOH_ON_CSQUARED * f * f * f / I);
}   /* T_Planck() */


/***********************************************************
* double T_Rayleigh_Jeans(double f, double T)
*
* Purpose:
*   Computes the Rayleigh-Jeans brightness temperature at
*   frequency f for a source at physical temperature T
*
* Arguments:
*   double f - frequency [GHz]
*   double T - temperature [K]
*
* Return:
*   Rayleigh-Jeans brightness temperature [K]
************************************************************/

double T_Rayleigh_Jeans(double f, double T)
{
    double x;

    if (T < DBL_MIN)
        return 0.;
    else if (f < DBL_MIN)
        return T;
    else {
        x = H_ON_KB * f / T;
        return T * x / (exp(x) - 1.0);
    }
}   /* T_Rayleigh_Jeans() */


/***********************************************************
* void planck_benchmarks(void)
*
* Purpose:
*   Run timing tests on functions in this file.
************************************************************/

void planck_benchmarks(void)
{
    const gridsize_t BMARK_MIN_NGRID  = 1000;    /* points */
    const gridsize_t BMARK_MAX_NGRID  = 1000000; /* points */
    const gridsize_t BMARK_NGRID_GROW = 10;
    const double BMARK_DF = 0.01;       /* [GHz] */
    const double BMARK_T  = 270.0;       /* [K] */
    const double BMARK_EST_TIME = 0.05;  /* [s] */
    const double BMARK_RUN_TIME = 0.25;  /* [s] */
    double t, tstart;
    double *f = NULL;
    double *B = NULL;
    double *T = NULL;
    unsigned int j;
    unsigned int ncalls;
    gridsize_t i, ngrid;
    /*
     * Start by setting up arrays for f and B
     */
    if (((f = (double *)malloc(BMARK_MAX_NGRID * sizeof(double))) == NULL) ||
        ((B = (double *)malloc(BMARK_MAX_NGRID * sizeof(double))) == NULL) ||
        ((T = (double *)malloc(BMARK_MAX_NGRID * sizeof(double))) == NULL)) {
        free(f);
        free(B);
        free(T);
        fprintf(stderr, "malloc() failed in planck_benchmarks().\n");
        return;
    }
    for (i = 0; i < BMARK_MAX_NGRID; ++i)
        f[i] = (double)i * BMARK_DF;
    /*
     * Print a header line
     */
    printf("%24s  %8s  %8s  %12s\n",
            "function", "ngrid", "ncalls", "t[ns]/point");
    printf("%24s  %8s  %8s  %12s\n",
            "--------", "-----", "------", "-----------");
    /*
     * B_Planck()
     */
    for (ngrid = BMARK_MIN_NGRID;
        ngrid <= BMARK_MAX_NGRID; ngrid *= BMARK_NGRID_GROW) {
        /*
         * estimation loop
         */
        tstart = am_timer(0.0);
        ncalls = 0;
        do {
            B_Planck(B, f, BMARK_DF, ngrid, BMARK_T);
            ++ncalls;
        } while ((t = am_timer(tstart)) < BMARK_EST_TIME);
        ncalls *= (unsigned int)(BMARK_RUN_TIME / t);
        ncalls = ncalls ? ncalls : 1;
        /*
         * measurement loop
         */
        tstart = am_timer(0.0);
        for (j = 0; j < ncalls; ++j)
            B_Planck(B, f, BMARK_DF, ngrid, BMARK_T);
        t = am_timer(tstart);
        printf("%24s  %8ld  %8d  %12.2f\n",
            "B_Planck()",
            (long int)ngrid,
            ncalls,
            1.e9 * t / (ncalls * ngrid));
    }
    /*
     * T_Planck().  Loops and OpenMP directives here are models for
     * similar loops in functions in spectra.c and rt.c, to serve as
     * guide to parallelizing those loops.  For these timings,
     * n_calls corresponds to the number of calls to those
     * functions.
     */
    for (ngrid = BMARK_MIN_NGRID;
        ngrid <= BMARK_MAX_NGRID; ngrid *= BMARK_NGRID_GROW) {
        tstart = am_timer(0.0);
        ncalls = 0;
        do {
            #ifdef _OPENMP
            #pragma omp parallel for schedule(guided, 128) if (ngrid >= 16384)
            #endif
            for(i = 0; i < ngrid; ++i)
                T[i] = T_Planck(f[i], B[i]);
            ++ncalls;
        } while ((t = am_timer(tstart)) < BMARK_EST_TIME);
        ncalls *= (unsigned int)(BMARK_RUN_TIME / t);
        ncalls = ncalls ? ncalls : 1;
        tstart = am_timer(0.0);
        for (j = 0; j < ncalls; ++j) {
            #ifdef _OPENMP
            #pragma omp parallel for schedule(guided, 128) if (ngrid >= 16384)
            #endif
            for(i = 0; i < ngrid; ++i)
                T[i] = T_Planck(f[i], B[i]);
        }
        t = am_timer(tstart);
        printf("%24s  %8ld  %8d  %12.2f\n",
            "T_Planck()",
            (long int)ngrid,
            ncalls,
            1.e9 * t / (ncalls * ngrid));
    }
    free(f);
    free(B);
    free(T);
    return;
}
