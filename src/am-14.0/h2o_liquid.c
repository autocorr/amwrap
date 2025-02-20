/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* h2o_liquid.c               S. Paine rev. 2019 September 26
*
* Liquid water properties and absorption.
************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "am_types.h"
#include "errlog.h"
#include "h2o_liquid.h"
#include "math_const.h"
#include "molecules.h"
#include "phys_const.h"
#include "rayleigh.h"


/*
 * Temperature limits for liquid water.
 *
 * A theoretical homogeneous ice nucleation temperature of
 * approximately
 * -45 C at atmospheric pressure can be inferred by extrapolating
 *  the temperature dependence of properties such as density,
 *  compressibility, etc. as shown in
 *
 *   R. J. Speedy and C. A. Angell 1976, "Isothermal
 *   compressibility of supercooled water and evidence for a
 *   thermodynamic singularity at -45 C."  J. Chem. Phys. 65:851.
 * 
 * This is the lower limit for liquid water in am models;
 * temperatures lower than this will generate an error.
 *
 * In nature, the observed low temperature limit for supercooled
 * liquid clouds is somewhat higher, about -35 C.  See
 *
 *  Yongxiang Hu, et al. 2010, "Occurrence, liquid water content,
 *  and fraction of supercooled water clouds  from combined
 *  CALIOP/IIR/MODIS measurements."  J. Geophys. Res. 115:D00H34.
 * 
 * Models with liquid water below this temperature (but above the
 * error limit) will generate a warning message.
 *
 * The upper temperature limit is the stated upper temperature
 * limit of the Ellison permittivity model.  Temperatures above
 * this will generate an error.
 */

static const double H2O_LIQUID_HIGH_T_ERROR_LIMIT_C  = 100.;
static const double H2O_LIQUID_LOW_T_WARNING_LIMIT_C = -35.;
static const double H2O_LIQUID_LOW_T_ERROR_LIMIT_C   = -45.;

static double H2O_liquid_density(double);
static int    H2O_liquid_permittivity(
        double*, double*, const double*, const gridsize_t, const double);

/***********************************************************
* int lwp_abs_Rayleigh(
*         double *k,
*         const double *f,
*         const gridsize_t ngrid,
*         const double T)
*
* Purpose:
*   Computes the molecular absorption coefficient [cm^2] for
*   liquid water as droplets, in the Rayleigh limit.
*
*   Internally, liquid water path (LWP) is represented as
*   a molecular column density in [cm^-2], but configuration
*   files will typically specify LWP in practical units such
*   as [kg*m^-2] or [mm_pwv].
*
* Arguments:
*   double *k              - absorption coefficient [cm^2]
*   const double *f        - frequency grid [GHz]
*   const gridsize_t ngrid - number of grid points
*   const double T         - temperature [K]
*
* Return:
*   0 if OK, 1 if an error occurred
************************************************************/

int lwp_abs_Rayleigh(
        double *k,
        const double *f,
        const gridsize_t ngrid,
        const double T)
{
    double *epsr = NULL;
    double *epsi = NULL;
    double rho;

    if (T < T_STP + H2O_LIQUID_LOW_T_ERROR_LIMIT_C ||
        T > T_STP + H2O_LIQUID_HIGH_T_ERROR_LIMIT_C) {
        errlog(139, 0);
        return 1;
    }
    if (T < T_STP + H2O_LIQUID_LOW_T_WARNING_LIMIT_C)
        errlog(140, 0);
    /*
     * Allocate temporary arrays and compute the complex permittivity
     * at temperature T.
     */
    if ((epsr = (double*)malloc(ngrid * sizeof(double))) == NULL) {
        errlog(141, 0);
        return 1;
    }
    if ((epsi = (double*)malloc(ngrid * sizeof(double))) == NULL) {
        free(epsr);
        errlog(141, 0);
        return 1;
    }
    if (H2O_liquid_permittivity(epsr, epsi, f, ngrid, T)) {
        free(epsi);
        free(epsr);
        return 1;
    }
    /*
     * Density is converted from [kg / m^3] to [molecules / cm^3]
     */
    rho = H2O_liquid_density(T - T_STP) / (1e6 * MASS_H2O * PCONST_AMU);
    Rayleigh_mol_abscoeff(k, f, ngrid, epsr, epsi, rho);
    free(epsi);
    free(epsr);
    return 0;
}   /* lwp_abs_Rayleigh() */


/***********************************************************
* static double H2O_liquid_density(double t)
*
* Purpose:
*   Computes the density of liquid water [kg / m^3], as a
*   function of temperature at atmospheric pressure, using
*   the rational function approximation in
*
*     G. S. Kell 1975, "Density, Thermal Expansivity, and
*     Compressibility of Liquid Water from 0 to 150 C:
*     Correlations and Tables for Atmospheric Pressure and
*     Saturation Reviewed and Expressed on 1968 Temperature
*     Scale."  Journal of Chemical and Engineering Data
*     20:97.
*
*   Comparison with the measurements of
*
*     D. E. Hare and C. M. Sorensen 1987, "The density of
*     supercooled water.  II. Bulk samples cooled to the
*     homogeneous nucleation limit."  J. Chem. Phys.
*     87:4840.
*
*   shows that the Kell function remains accurate when
*   extended below 0 C, down to the homogeneous ice
*   nucleation temperature.
*
* Arguments:
*   double t - Temperature [C], between -45 C and 150 C.
*     This argument is assumed to have been checked by the
*     calling function.
*
* Return:
*   Density [kg / m^3].
************************************************************/

static double H2O_liquid_density(double t)
{
    static const double a0 =  9.9983952e+02;
    static const double a1 =  1.6945176e+01;
    static const double a2 = -7.9870401e-03;
    static const double a3 = -4.6170461e-05;
    static const double a4 =  1.0556302e-07;
    static const double a5 = -2.8054253e-10;
    static const double b1 =  1.6879850e-02;
    
    return (a0 + t * (a1 + t * (a2 + t * (a3 + t * (a4 + t * a5))))) /
        (1.0 + b1 * t);
}   /* H2O_liquid_density() */


/***********************************************************
* static int H2O_liquid_permittivity(
*         double *epsr,
*         double *epsi,
*         const double *f,
*         const gridsize_t ngrid,
*         const double T)
*
* Purpose:
*   Computes the complex permittivity of water, using a
*   model incorporating three Debye relaxation processes and
*   two resonant absorption processes. The strengths and
*   relaxation time constants of the three Debye processes;
*   the strengths, frequencies, and widths of the two
*   resonances; and their respective temperature dependences
*   are modeled following
*
*     W. J. Ellison 2007, "Permittivity of Pure Water, at
*     Standard Atmospheric Pressure, over the Frequency
*     Range 0 - 25 THz and the Temperature range 0 - 100 C."
*     J. Phys.  Chem. Ref. Data 36:1.
*
*   A more recent model derived using additional field
*   observations of supercooled clouds, but limited to
*   modeling the permittivity at frequencies below 500 GHz
*   using only the two lowest-frequency Debye processes, is
*   the TKC model described in
*
*     D. D. Turner, S. Kneifel, and M. P. Cadeddu 2016, "An
*     Improved Liquid Water Absorption Model at Microwave
*     Frequencies for Supercooled Liquid Water Clouds."
*     Journal of Atmospheric and Oceanic Technology 33:33.
*
*   The TKC model should be an improvement over the Ellison
*   2007 model at supercooled temperatures and below 500
*   GHz.
*
*   This function can be compiled to use one of three
*   possible sets of model coefficients by defining one of
*   the following macros at compile time:
*
*     1. H2O_PERMITTIVITY_ELLISON_2007 
*        Use coefficients from Ellison 2007 only.
*
*     2. H2O_PERMITTIVITY_TKC
*        Use the TKC model coefficients only.  The
*        amplitudes of the third Debye process and two
*        resonances are set to zero.
*
*     3. H2O_PERMITTIVITY_ELLISON_2007_WITH_TKC (default)
*        Use the TKC coefficients for the two
*        lowest-frequency Debye processes and critical
*        temperature tc used in modeling the temperature
*        dependence of the relaxation time constants, and
*        Ellison 2007 for the remaining coefficients.
*   
* Arguments:
*   double *epsr - pointer to array to receive the real part
*                  of the permittivity
*   double *epsi - pointer to array to receive the imaginary
*                  part of the permittivity
*   const double *f        - frequency grid [GHz]
*   const gridsize_t ngrid - number of grid points
*   const double T         - temperature [K]
*
* Return:
*   0 if OK, 1 if an error occurred
************************************************************/

#if !defined H2O_PERMITTIVITY_ELLISON_2007 && \
    !defined H2O_PERMITTIVITY_TKC && \
    !defined H2O_PERMITTIVITY_ELLISON_2007_WITH_TKC
    #define H2O_PERMITTIVITY_ELLISON_2007_WITH_TKC
#endif

static int H2O_liquid_permittivity(
        double *epsr,
        double *epsi,
        const double *f,
        const gridsize_t ngrid,
        const double T)
{
    gridsize_t i;
    const double t = T - T_STP; /* Temperature in C */
    double es;
    double delta1, delta2, delta3;
    double tau1, tau2, tau3;
    double u1, u2, u3, v1, v2, v3, w1, w2, w3;
    double delta4, delta5;
    double tau4, tau5;
    double f0, f1;
    double u4, u5, v4, v5, w4, w5;
    double r1;
#ifdef H2O_PERMITTIVITY_ELLISON_2007
    /*
     * Constants from Ellison 2007 Table 2
     */
    static const double  a1 =  7.923882e+01;
    static const double  b1 =  4.300598e-03;
    static const double  c1 =  1.382264e-13;
    static const double  d1 =  6.527648e+02;
    static const double  a2 =  3.815866e+00;
    static const double  b2 =  1.117295e-02;
    static const double  c2 =  3.510354e-16;
    static const double  d2 =  1.249533e+03;
    static const double  a3 =  1.634967e+00;
    static const double  b3 =  6.841548e-03;
    static const double  c3 =  6.300350e-15;
    static const double  d3 =  4.055169e+02;
    static const double  tc =  1.331383e+02;
    static const double  p0 =  8.379692e-01;
    static const double  p1 = -6.118594e-03;
    static const double  p2 = -1.293680e-05;
    static const double  p3 =  4.235901e+12;
    static const double  p4 = -1.426088e+10;
    static const double  p5 =  2.738157e+08;
    static const double  p6 = -1.246943e+06;
    static const double  p7 =  9.618642e-14;
    static const double  p8 =  1.795786e-16;
    static const double  p9 = -9.310017e-18;
    static const double p10 =  1.655473e-19;
    static const double p11 =  6.165532e-01;
    static const double p12 =  7.238532e-03;
    static const double p13 = -9.523366e-05;
    static const double p14 =  1.598317e+13;
    static const double p15 = -7.441357e+10;
    static const double p16 =  4.974480e+08;
    static const double p17 =  2.882476e-14;
    static const double p18 = -3.142118e-16;
    static const double p19 =  3.528051e-18;
#endif
#ifdef H2O_PERMITTIVITY_TKC
    /*
     * Constants from Turner, Kneifel, and Cadeddu 2016.
     */
    static const double  a1 = 8.111e+01;
    static const double  b1 = 4.434e-03;
    static const double  c1 = 1.302e-13;
    static const double  d1 = 6.627e+02;
    static const double  a2 = 2.025e+00;
    static const double  b2 = 1.073e-02;
    static const double  c2 = 1.012e-14;
    static const double  d2 = 6.089e+02;
    static const double  a3 = 0.0;
    static const double  b3 = 0.0;
    static const double  c3 = 0.0;
    static const double  d3 = 0.0;
    static const double  tc = 1.342e+02;
    static const double  p0 = 0.0;
    static const double  p1 = 0.0;
    static const double  p2 = 0.0;
    static const double  p3 = 0.0;
    static const double  p4 = 0.0;
    static const double  p5 = 0.0;
    static const double  p6 = 0.0;
    static const double  p7 = 0.0;
    static const double  p8 = 0.0;
    static const double  p9 = 0.0;
    static const double p10 = 0.0;
    static const double p11 = 0.0;
    static const double p12 = 0.0;
    static const double p13 = 0.0;
    static const double p14 = 0.0;
    static const double p15 = 0.0;
    static const double p16 = 0.0;
    static const double p17 = 0.0;
    static const double p18 = 0.0;
    static const double p19 = 0.0;
#endif
#ifdef H2O_PERMITTIVITY_ELLISON_2007_WITH_TKC
    /*
     * Ellison 2007 modified per TKC, with the parameters for the
     * two lowest-frequency Debye processes and critical temperature
     * taken from TKC.
     */
    static const double  a1 = 8.111e+01;
    static const double  b1 = 4.434e-03;
    static const double  c1 = 1.302e-13;
    static const double  d1 = 6.627e+02;
    static const double  a2 = 2.025e+00;
    static const double  b2 = 1.073e-02;
    static const double  c2 = 1.012e-14;
    static const double  d2 = 6.089e+02;
    static const double  a3 =  1.634967e+00;
    static const double  b3 =  6.841548e-03;
    static const double  c3 =  6.300350e-15;
    static const double  d3 =  4.055169e+02;
    static const double  tc = 1.342e+02;
    static const double  p0 =  8.379692e-01;
    static const double  p1 = -6.118594e-03;
    static const double  p2 = -1.293680e-05;
    static const double  p3 =  4.235901e+12;
    static const double  p4 = -1.426088e+10;
    static const double  p5 =  2.738157e+08;
    static const double  p6 = -1.246943e+06;
    static const double  p7 =  9.618642e-14;
    static const double  p8 =  1.795786e-16;
    static const double  p9 = -9.310017e-18;
    static const double p10 =  1.655473e-19;
    static const double p11 =  6.165532e-01;
    static const double p12 =  7.238532e-03;
    static const double p13 = -9.523366e-05;
    static const double p14 =  1.598317e+13;
    static const double p15 = -7.441357e+10;
    static const double p16 =  4.974480e+08;
    static const double p17 =  2.882476e-14;
    static const double p18 = -3.142118e-16;
    static const double p19 =  3.528051e-18;
#endif
    /*
     * Static dielectric constant coefficients from Hamelin et al. 1998
     */
    static const double e0 =  8.79144e+01;
    static const double e1 = -4.04399e-01;
    static const double e2 =  9.58726e-04;
    static const double e3 = -1.32802e-06;

    if (t < H2O_LIQUID_LOW_T_ERROR_LIMIT_C ||
        t > H2O_LIQUID_HIGH_T_ERROR_LIMIT_C) {
        return 1;
    }
    es = e0 + t * (e1 + t * (e2 + t * e3));
    delta1 = a1 * exp(-b1 * t);
    delta2 = a2 * exp(-b2 * t);
    delta3 = a3 * exp(-b3 * t);
    r1 = t + tc;
    if (r1 <= 0.)
        return 1;
    r1 = 1.0 / r1;
    tau1   = c1 * exp(d1 * r1);
    tau2   = c2 * exp(d2 * r1);
    tau3   = c3 * exp(d3 * r1);
    delta4 = p0 + t * (p1 + t * p2);
    f0   = p3 + t * (p4 + t * (p5 + t * p6));
    tau4 = p7 + t * (p8 + t * (p9 + t * p10));
    delta5 = p11 + t * (p12 + t * p13);
    f1     = p14 + t * (p15 + t * p16);
    tau5   = p17 + t * (p18 + t * p19);
    u1 = TWOPI * tau1;
    w1 = u1 * delta1;
    u1 *= u1;
    v1 = u1 * delta1;
    u2 = TWOPI * tau2;
    w2 = u2 * delta2;
    u2 *= u2;
    v2 = u2 * delta2;
    u3 = TWOPI * tau3;
    w3 = u3 * delta3;
    u3 *= u3;
    v3 = u3 * delta3;
    u4 = TWOPI * tau4;
    w4 = 0.5 * u4 * delta4;
    v4 = w4 * u4;
    u4 *= u4;
    u5 = TWOPI * tau5;
    w5 = 0.5 * u5 * delta5;
    v5 = w5 * u5;
    u5 *= u5;
    for (i = 0; i < ngrid; ++i) {
        double nu = f[i] * 1.e9; /* [Hz] */
        double nusq = nu * nu;
        double f0p = f0 + nu;
        double f0psq = f0p * f0p;
        double f0m = f0 - nu;
        double f0msq = f0m * f0m;
        double f1p = f1 + nu;
        double f1psq = f1p * f1p;
        double f1m = f1 - nu;
        double f1msq = f1m * f1m;
        epsr[i] = es - nusq * (
            v1 / (1.0 + u1 * nusq) +
            v2 / (1.0 + u2 * nusq) +
            v3 / (1.0 + u3 * nusq)) -
            v4 *
            (nu * f0p / (1.0 + u4 * f0psq) - nu * f0m / (1.0 + u4 * f0msq)) -
            v5 *
            (nu * f1p / (1.0 + u5 * f1psq) - nu * f1m / (1.0 + u5 * f1msq));
        epsi[i] = nu * (
            w1 / (1.0 + u1 * nusq) +
            w2 / (1.0 + u2 * nusq) +
            w3 / (1.0 + u3 * nusq) + 
            w4 * (1.0 / (1.0 + u4 * f0psq) + 1.0 / (1.0 + u4 * f0msq)) +
            w5 * (1.0 / (1.0 + u5 * f1psq) + 1.0 / (1.0 + u5 * f1msq)));
    }
    return 0;
}   /* H2O_liquid_permittivity() */


#ifdef UNIT_TEST

/***********************************************************
* int main(int argc, char **argv)
*
* Purpose:
*   Tests H2O_liquid_density() and H2O_liquid_permittivity()
*   in this file.  Usage is
*
*   a.out T[C] fmin[GHz] fmax[GHz] num_points
*
*   where num_points is the number of logarithmically-
*   spaced frequency grid points.
*
* Example:
*   $ gcc -D UNIT_TEST h2o_liquid.c
*   $ ./a.out -10 0.1 1500 100
************************************************************/


int main(int argc, char **argv)
{
    double T, fmin, fmax, r = 1.0;
    double *f, *epsr, *epsi;
    gridsize_t i;
    int n;

    if (argc < 5) {
        printf("usage: %s T[C] fmin[GHz] fmax[GHz] num_points\n",
            argv[0]);
        return 1;
    }
    T    = atof(argv[1]) + T_STP;
    fmin = atof(argv[2]);
    fmax = atof(argv[3]);
    n    = atoi(argv[4]);
    fmin = fmin < 0.1 ? 0.1 : fmin;
    if (n > 1)
        r = pow(fmax / fmin, 1. / (double)(n - 1));
    f    = (double*)malloc(n * sizeof(double));
    epsr = (double*)malloc(n * sizeof(double));
    epsi = (double*)malloc(n * sizeof(double));
    printf("# rho(%.8g C) = %.8g kg m^-3\n",
        T - T_STP, H2O_liquid_density(T - T_STP));
    f[0] = fmin;
    for (i = 1; i < n; ++i)
        f[i] = f[i-1] * r;
    if (H2O_liquid_permittivity(epsr, epsi, f, n, T))
        return 1;
    printf("#%13s %14s %14s\n", "f[GHz]", "epsr", "epsi");
    for (i = 0; i < n; ++i)
        printf("%14e %14e %14e\n", f[i], epsr[i], epsi[i]);
    free(f);
    free(epsr);
    free(epsi);
    return 0;
}   /* main() */


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
*int Rayleigh_mol_abscoeff(
*    double *k,
*    const double *f,
*    const gridsize_t ngrid,
*    const double *epsr,
*    const double *epsi,
*    const double rho)
*
* Purpose:
*   dummy Rayleigh_mol_abscoeff() function for unit test
************************************************************/

int Rayleigh_mol_abscoeff(
    double *k,
    const double *f,
    const gridsize_t ngrid,
    const double *epsr,
    const double *epsi,
    const double rho)
{
    return 0;
}   /* Rayleigh_mol_abscoeff() */

#endif /* UNIT_TEST */
