/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* cia.c                            S. Paine rev. 2023 May 18
*
* Collision-induced absorption.  The principal references
* for for this file are:
*
*  A. Borysow and L. Frommhold, Collision-induced
*   rototranslational absorption spectra of N2-N2 pairs for
*   temperatures from 50 to 300 K.  Ap. J. 311:1043 (1986
*   Dec 15). (erratum in Ap. J. 320:437 (1987 Sep 1))
*
*  J. Boissoles, C. Boulet, R.H. Tipping, A. Brown, Q. Ma,
*   Theoretical calculation of the translation-rotation
*   collision-induced absorption in N2-N2, O2-O2, and N2-O2
*   pairs.  JQSRT 82:505 (2003)
*
* For useful background see the following, and references
* therein:
*
*  L. Frommhold, _Collision-induced_Absorption_in_Gases_
*   (Cambridge University Press, 1993) (note scale error in
*   Fig. 3.21 - should be 10^-7 cm^-1 amagat^-2)
*
*  R. M. Goody and Y. L. Yung,
*   _Atmospheric_Radiation,_Theoretical Basis_, 2nd ed.
*   (Oxford University Press, 1989) sections 3.4 and 5.2.
*
*  G. Birnbaum and E. R. Cohen, Theory of line shape in
*   pressure- induced absorption.  Can. J. Phys. 54:593
*   (1976).
*
* For N2-N2, this code has been validated from 93 K to 343 K
* with the experimental data in:
*
*  N.W.B. Stone, L.A.A. Read, A. Anderson, I.R. Dagg, and W.
*   Smith, Temperature dependent collision-induced
*   absorption in nitrogen.  Can. J. Phys. 62:338 (1984).
*
*  I.R. Dagg, A. Anderson, S. Yan, W. Smith, and L.A.A. Read,
*   Collision-induced absorption in nitrogen at low
*   temperatures.  Can. J. Phys. 63:625 (1985).
*
*  P. Dore and A. Filabozzi, On the nitrogen-induced
*   far-infrared absorption spectra.  Can. J. Phys. 65:90
*   (1987).
*
* For O2-O2, experimental validation is possible only at 300
* K, using the data of:
*
*  D.R. Bosomworth and H.P. Gush, Collision-induced
*   absorption of compressed gases in the far infrared, Part
*   II.  Can. J. Phys.  43:751 (1965).
************************************************************/

#include <float.h>
#include <math.h>
#ifdef _OPENMP
    #include <omp.h>
#endif

#include "am_sysdep.h"
#include "am_types.h"
#include "errlog.h"
#include "math_const.h"
#include "phys_const.h"
#include "specfunc.h"

/*
 * EBC_pars_t - Structure for CIA lineshape parameters.
 *
 * The rototranslational spectra are modeled using the extended
 * Birnbaum-Cohen (EBC) lineshape, as defined in Borysow and
 * Frommhold, Eq. 4.  The line shape depends on up to six
 * temperature-dependent parameters: S, tau1, tau2, eps, tau1p,
 * and tau2p.  For each parameter X, the temperature dependence
 * is modeled as
 *
 *  X(T) = X_A * exp(X_B * log(T) + X_C * (log(T))^2).
 *
 * The following structure holds the A, B, and C coefficients for
 * each line shape parameter, and saved intermediate results for
 * efficient computation of the line shape function at a given
 * temperature.  Units are time in ns, frequency in GHz, and S
 * has the dimensions [cm^5] of a binary absorption coefficient.
 * (S[erg cm^6] as defined in Borysow and Frommhold is multiplied
 * by (4 pi^3) / (3 hbar c) to obtain S[cm^5]).  The coefficients
 * here are derived from Borysow and Frommhold, Table 5, with
 * adjustments as noted below.
 */
typedef struct EBC_pars_t {
    double S, S_A, S_B, S_C;
    double tau1, tau1_A, tau1_B, tau1_C;
    double tau2, tau2_A, tau2_B, tau2_C;
    double eps, eps_A, eps_B, eps_C;
    double tau1p, tau1p_A, tau1p_B, tau1p_C;
    double tau2p, tau2p_A, tau2p_B, tau2p_C;
    double tau0;    /* hbar / (2kT)                     */
    double r0;      /* S / (pi * (1 + eps))             */
    double r1;      /* 2 * pi * tau0                    */
    double r2;      /* (2 * pi * tau1)^2                */
    double r3;      /* tau2 / tau 1                     */
    double r4;      /* sqrt(tau0^2 + tau2^2) / tau1     */
    double r5;      /* tau1p / tau2p                    */
    double r6;      /* sqrt(tau0^2 + tau1p^2) / tau2p   */
    double r7;      /* (2 * pi * tau2p)^2               */
} EBCpars_t;

static double clebsch_squared(const int, const int, const int);
static double Erot(const int, const double, const double);
static void   initialize_EBCpars(EBCpars_t*, const double);
static double EBCfunc(const EBCpars_t*, const double);

/*
 * Rotational constants [GHz].
 */

/*
 * N2, from J. Bendtsen, J. Raman Spectroscopy 2:133 (1974)
 */
static const double N2_B_ROT = 59.6459;
static const double N2_D_ROT = 1.727e-4;
/*
 * O2, from G. Yu. Golubiatnikov and A. F. Krupnov, J. Mol. Spec.
 * 217:282 (2003)
 */
static const double O2_B_ROT = 43.1004438;
static const double O2_D_ROT = 1.45115e-4;

/*
 * Maximum rotational quantum number for the initial state in
 * the line-by-line sums.
 */
enum {
    N2_JMAX = 30,
    O2_JMAX = 35
};

/*
 * Coefficients for the series approximation of the rotational
 * partition function,
 *  Q = U0 * x + U1 + U2 / x, where  x = kT / hB
 */
/* N2, rel. error < 4e-5 for 120 K < T < 320 K, 3e-4 at 50 K */
static const double N2_U0 = 4.50454;
static const double N2_U1 = 1.27450;
static const double N2_U2 = 3.40814;
/* O2, rel. error < 6e-5 for 120 K < T < 320 K, 5e-4 at 50 K */
static const double O2_U0 = 0.50082;
static const double O2_U1 = 0.11041;
static const double O2_U2 = 1.1058;

/*
 * Extended Birnbaum-Cohen (EBC) lineshape data.
 */
/*
 * N2-N2, temperature ranges fitted in Borysow and Frommhold, Ap.
 * J. 311:1043 Table 5.  The temperature range of the fits there
 * is 50 K - 300 K, but the output of am has been validated
 * against experimental data up to 350 K.
 */
static const double N2_TMIN = 50.;
static const double N2_TMID = 140.;
static const double N2_TMAX = 350.;
/*
 * N2-N2, quadrupolar induction {3220 + 3202}, 50 K - 350 K From
 * Borysow and Frommhold, Table 5, as corrected in Ap. J. 320:437
 * (1987).  The parameters for tau1 and tau2 are unchanged here,
 * but the parameters for S have been adjusted to obtain a better
 * fit to the experimental data referenced at the top of this
 * file.  This adjustment is needed because am includes only the
 * free-free and bound- free components of the translational
 * spectral function; the bound-bound component included in the
 * analysis of Borysow and Frommhold is omitted.  This is a
 * computational simplification, and also avoids the need to
 * postulate a broadening mechanism for the sharp bound-bound
 * spectral lines.  With the parameters below, the increase in
 * the free-free plus bound-free quadrupole component is +7.2% at
 * 300K, +%23 at 150 K, and +48% at 100 K.  The total spectra
 * agree, within experimental errors, with the spectra cited
 * above, which cover a temperature range from 90 K to 350 K.
 */
static const EBCpars_t N2_EBC3220_INIT = {
/*  original Borysow-Frommhold parameters       */
/*  0.0, 3.88672e-42, -0.99569,  0.09464,       */
/*  Refitted per above without b-b component    */
    0.0, 1.31362e-38, -3.79195,  0.33728,
    0.0, 0.12962e-02, -0.13048, -0.03128,
    0.0, 0.37969e-04,  1.03681, -0.14336,
    0.0, 0.0,          0.0,      0.0,
    0.0, 0.0,          0.0,      0.0,
    0.0, 0.0,          0.0,      0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
};

/*
 * N2-N2, hexadecapolar induction {5440 + 5404}, 50 K - 140 K
 */
static const EBCpars_t N2_EBC5440_LT_INIT = {
    0.0, 2.36588e-43, -1.69153,  0.18605,
    0.0, 0.66017e-06,  2.59982, -0.31831,
    0.0, 0.12481e-02, -0.57028,  0.05983,
    0.0, 0.3,          0.0,      0.0,
    0.0, 0.52681e-03, -0.24719,  0.00519,
    0.0, 0.27518e+25,-25.38969,  2.46542,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
};

/*
 * N2-N2, hexadecapolar induction {5440 + 5404}, 140 K - 350 K
 */
static const EBCpars_t N2_EBC5440_HT_INIT = {
    0.0, 1.07920e-43, -1.25562,  0.12981,
    0.0, 0.36611e-05,  1.47688, -0.16537,
    0.0, 0.61264e+00, -2.25011,  0.15289,
    0.0, 0.3,          0.0,      0.0,
    0.0, 0.79820e+00, -2.76152,  0.21847,
    0.0, 0.52868e-12,  7.66253, -0.77527,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
};

/*
 * O2-O2 spectral functions are just scaled versions of the N2-N2
 * functions, following Boissoles et al., referenced above.
 * Given the lack of experimental data at temperatures other than
 * 300 K, the same temperature dependence as for N2-N2 is
 * assumed.
 *
 * For am-12.3+, the initial band strength estimates from
 * Boissoles et al., and the coefficient governing the admixture
 * of BC and K0 lineshapes in the EBC lineshape, were adjusted
 * slightly to minimize the squared residuals relative to the 300
 * K spectrum of Bosomworth and Gush (1965).  The values before
 * this adjustment was made are commented out below.  This
 * adjustment reduces the standard deviation of these residuals
 * by a factor 0.75.
 */
static const double O2_TMIN = 50.;
static const double O2_TMID = 140.;
static const double O2_TMAX = 350.;
/*
 * O2-O2, quadrupolar induction {3220 + 3202}, 50 K - 350 K
 */
static const EBCpars_t O2_EBC3220_INIT = {
/*  0.0, 9.75000e-40, -3.79195,  0.33728, */
    0.0, 1.09900e-39, -3.79195,  0.33728,
    0.0, 0.12962e-02, -0.13048, -0.03128,
    0.0, 0.37969e-04,  1.03681, -0.14336,
    0.0, 0.0,          0.0,      0.0,
    0.0, 0.0,          0.0,      0.0,
    0.0, 0.0,          0.0,      0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
};

/*
 * O2-O2, hexadecapolar induction {5440 + 5404}, 50 K - 140 K
 */
static const EBCpars_t O2_EBC5440_LT_INIT = {
/*  0.0, 5.81000e-43, -1.69153,  0.18605,  */
    0.0, 5.38800e-43, -1.69153,  0.18605,
    0.0, 0.66017e-06,  2.59982, -0.31831,
    0.0, 0.12481e-02, -0.57028,  0.05983,
/*  0.0, 0.3,          0.0,      0.0,   */
    0.0, 0.520,        0.0,      0.0,
    0.0, 0.52681e-03, -0.24719,  0.00519,
    0.0, 0.27518e+25,-25.38969,  2.46542,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
};

/*
 * O2-O2, hexadecapolar induction {5440 + 5404}, 140 K - 350 K
 */
static const EBCpars_t O2_EBC5440_HT_INIT = {
/*  0.0, 2.65000e-43, -1.25562,  0.12981,   */
    0.0, 2.45800e-43, -1.25562,  0.12981,
    0.0, 0.36611e-05,  1.47688, -0.16537,
    0.0, 0.61264e+00, -2.25011,  0.15289,
/*  0.0, 0.3,          0.0,      0.0,   */
    0.0, 0.520,        0.0,      0.0,
    0.0, 0.79820e+00, -2.76152,  0.21847,
    0.0, 0.52868e-12,  7.66253, -0.77527,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
};


/***********************************************************
* int N2N2_cia(
*   double *k,
*   const double *f,
*   const gridsize_t ngrid,
*   const double T)
*
* Purpose:
*   Computes the N2-N2 collision-induced absorption
*   coefficient, following the treatment of Borysow and
*   Frommhold referenced above.  Here, only the dominant
*   quadrupole and hexadecapole induced rototranslational
*   spectra are included.  Bound-bound and double
*   transitions are neglected.
*
*   k [cm^5], is a binary absorption coefficient, such that
*   if n [cm^-3] is a molecular number density, then
*   av [cm-1] = k * n^2 is a volume absorption coefficient.
*
* Arguments:
*   double *k               - binary spectral absorption
*                             coefficient [cm^5]
*   const double *f         - frequency grid [GHz]
*   const gridsize_t ngrid  - number of grid points
*   const double T          - temperature [K]
*
* Return:
*   0 if OK, 1 if an error occurred
************************************************************/

int N2N2_cia(
        double *k,
        const double *f,
        const gridsize_t ngrid,
        const double T)
{
    double pop[N2_JMAX + 1];    /* rotational level populations     */
    double Qrot;                /* rotational partition function    */
    double r0;
    EBCpars_t EBCpars_3220, EBCpars_5440;
    gridsize_t offset, blocksize;
    int j;

    /*
     * Log a warning if T is outside the fitted range for the EBC
     * line shape parameters.
     */
    if (T < N2_TMIN || T > N2_TMAX) {
        errlog(32, 0);
    }
    /*
     * Calculate the rotational level populations.
     */
    r0 = T / (H_ON_KB * N2_B_ROT);
    Qrot = N2_U0 * r0 + N2_U1 + N2_U2 / r0;
    r0 = H_ON_KB / T;
    for (j = 0; j <= N2_JMAX; ++j) {
        pop[j] = (j % 2 == 0 ? 6.0 : 3.0);
        pop[j] *= 2.0 * (double)j + 1.0;
        pop[j] *= am_exp(-Erot(j, N2_B_ROT, N2_D_ROT) * r0) / Qrot;
    }
    /*
     * Initialize the lineshape parameters.
     */
    EBCpars_3220 = N2_EBC3220_INIT;
    initialize_EBCpars(&EBCpars_3220, T);
    if (T > N2_TMID)
        EBCpars_5440 = N2_EBC5440_HT_INIT;
    else
        EBCpars_5440 = N2_EBC5440_LT_INIT;
    initialize_EBCpars(&EBCpars_5440, T);
    /*
     * Break the computation over the frequency grid into blocks
     * to fit in L1 cache.  Under OpenMP, make at least one block
     * for each thread, and do them in parallel.
     */
    blocksize = ((L1_CACHE_BYTES) == 0) ? ngrid :
        (gridsize_t)((L1_CACHE_BYTES) / ((L1_CACHE_WAYS) * sizeof(double)));
    #ifdef _OPENMP
    {
        gridsize_t thd_blocksize;
        thd_blocksize = ngrid / omp_get_max_threads();
        blocksize = (blocksize < thd_blocksize) ? blocksize : thd_blocksize;
    }
    #endif
    blocksize = blocksize < 1 ? 1 : blocksize;
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 1)
    #endif
    for (offset = 0; offset < ngrid; offset += blocksize) {
        double r1;
        double *kb, *fb;
        gridsize_t i, nblock;
        int j_1, j_1p;  /* angular momenta  */
        nblock = ngrid - offset;
        if (nblock > blocksize)
            nblock = blocksize;
        kb = k + offset;
        fb = (double*)f + offset; 
        /*
         * Initialize the absorption coefficient array.
         */
        for (i = 0; i < nblock; ++i)
            kb[i] = 0.0;
        /*
         * Line-by-line sum for quadrupole-induced single
         * transitions, Note that for single transitions, the sum
         * over the states of the second molecule is unity, so it
         * is suppressed here.
         */
        for (j_1 = 0; j_1 <= N2_JMAX; ++j_1) {
            double f1;
            f1 = Erot(j_1, N2_B_ROT, N2_D_ROT);
            for (j_1p = j_1 - 2; j_1p <= j_1 + 2; j_1p += 2) {
                double f0, r2;
                if (j_1p < 0)
                    continue;
                r2 = clebsch_squared(j_1, 2, j_1p);
                r2 *= pop[j_1];
                if (r2 < DBL_EPSILON)
                    continue;
                f0 = Erot(j_1p, N2_B_ROT, N2_D_ROT) - f1;
                for (i = 0; i < nblock; ++i)
                    kb[i] += r2 * EBCfunc(&EBCpars_3220, fb[i] - f0);
            }
        }
        /*
         * Line-by-line sum for hexadecapole-induced single
         * transitions.
         */
        for (j_1 = 0; j_1 < N2_JMAX; ++j_1) {
            double f1;
            f1 = Erot(j_1, N2_B_ROT, N2_D_ROT);
            for (j_1p = j_1 - 4; j_1p <= j_1 + 4; j_1p += 2) {
                double f0, r2;
                if (j_1p < 0)
                    continue;
                r2 = clebsch_squared(j_1, 4, j_1p);
                r2 *= pop[j_1];
                if (r2 < DBL_EPSILON)
                    continue;
                f0 = Erot(j_1p, N2_B_ROT, N2_D_ROT) - f1;
                for (i = 0; i < nblock; ++i)
                    kb[i] += r2 * EBCfunc(&EBCpars_5440, fb[i] - f0);
            }
        }
        /*
         * Multiply k through by factors outside the summation.
         */
        r1 = H_ON_KB / T;
        for (i = 0; i < nblock; ++i)
            kb[i] *= fb[i] * (1.0 - am_exp(-fb[i] * r1));
    }
    return 0;
}   /* N2N2_cia() */


/***********************************************************
* int O2O2_cia(
*   double *k,
*   const double *f,
*   const gridsize_t ngrid,
*   const double T)
*
* Purpose:
*   Computes the O2-O2 collision-induced absorption
*   coefficient.  Following Boissoles, et al., referenced
*   above, the Borysow- Frommhold lineshape used for the
*   N2-N2 CIA is used, rescaled to reflect the differences
*   in the multipole moments and polarizability of O2 versus
*   N2.
*
*   k [cm^5], is a binary absorption coefficient, such that
*   if n [cm^-3] is a molecular number density, then
*   av [cm-1] = k * n^2 is a volume absorption coefficient.
*
* Arguments:
*   double *k               - binary spectral absorption
*                             coefficient [cm^5]
*   const double *f         - frequency grid [GHz]
*   const gridsize_t ngrid  - number of grid points
*   const double T          - temperature [K]
*
* Return:
*   0 if OK, 1 if an error occurred
************************************************************/

int O2O2_cia(
        double *k,
        const double *f,
        const gridsize_t ngrid,
        const double T)
{
    double pop[O2_JMAX + 1];    /* rotational level populations     */
    double Qrot;                /* rotational partition function    */
    double r0;
    EBCpars_t EBCpars_3220, EBCpars_5440;
    gridsize_t offset, blocksize;
    int j;

    /*
     * Log a warning if T is outside the fitted range for the EBC
     * line shape parameters.
     */
    if (T < O2_TMIN || T > O2_TMAX) {
        errlog(32, 0);
    }
    /*
     * Calculate the rotational level populations.
     */
    r0 = T / (H_ON_KB * O2_B_ROT);
    Qrot = O2_U0 * r0 + O2_U1 + O2_U2 / r0;
    r0 = H_ON_KB / T;
    for (j = 0; j <= O2_JMAX; ++j) {
        pop[j] = (j % 2 == 0 ? 0.0 : 1.0);
        pop[j] *= 2.0 * (double)j + 1.0;
        pop[j] *= am_exp(-Erot(j, O2_B_ROT, O2_D_ROT) * r0) / Qrot;
    }
    /*
     * Initialize the lineshape parameters.
     */
    EBCpars_3220 = O2_EBC3220_INIT;
    initialize_EBCpars(&EBCpars_3220, T);
    if (T > O2_TMID)
        EBCpars_5440 = O2_EBC5440_HT_INIT;
    else
        EBCpars_5440 = O2_EBC5440_LT_INIT;
    initialize_EBCpars(&EBCpars_5440, T);
    /*
     * Break the computation over the frequency grid into blocks
     * to fit in L1 cache.  Under OpenMP, make at least one block
     * for each thread, and do them in parallel.
     */
    blocksize = ((L1_CACHE_BYTES) == 0) ?
        ngrid : (gridsize_t)((L1_CACHE_BYTES) / (3 * sizeof(double)));
    #ifdef _OPENMP
    {
        gridsize_t thd_blocksize;
        thd_blocksize = ngrid / omp_get_max_threads();
        blocksize = (blocksize < thd_blocksize) ? blocksize : thd_blocksize;
    }
    #endif
    blocksize = blocksize < 1 ? 1 : blocksize;
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 1)
    #endif
    for (offset = 0; offset < ngrid; offset += blocksize) {
        double r1;
        double *kb, *fb;
        gridsize_t i, nblock;
        int j_1, j_1p;  /* angular momenta  */
        nblock = ngrid - offset;
        if (nblock > blocksize)
            nblock = blocksize;
        kb = k + offset;
        fb = (double*)f + offset; 
        /*
         * Initialize the absorption coefficient array.
         */
        for (i = 0; i < nblock; ++i)
            kb[i] = 0.0;
        /*
         * Line-by-line sum for quadrupole-induced single
         * transitions, Note that for single transitions, the sum
         * over the states of the second molecule is unity, so it
         * is suppressed here.
         */
        for (j_1 = 0; j_1 <= O2_JMAX; ++j_1) {
            double f1;
            f1 = Erot(j_1, O2_B_ROT, O2_D_ROT);
            for (j_1p = j_1 - 2; j_1p <= j_1 + 2; j_1p += 2) {
                double f0, r2;
                if (j_1p < 0)
                    continue;
                r2 = clebsch_squared(j_1, 2, j_1p);
                r2 *= pop[j_1];
                if (r2 < DBL_EPSILON)
                    continue;
                f0 = Erot(j_1p, O2_B_ROT, O2_D_ROT) - f1;
                for (i = 0; i < nblock; ++i)
                    kb[i] += r2 * EBCfunc(&EBCpars_3220, fb[i] - f0);
            }
        }
        /*
         * Line-by-line sum for hexadecapole-induced single transitions.
         */
        for (j_1 = 0; j_1 < O2_JMAX; ++j_1) {
            double f1;
            f1 = Erot(j_1, O2_B_ROT, O2_D_ROT);
            for (j_1p = j_1 - 4; j_1p <= j_1 + 4; j_1p += 2) {
                double f0, r2;
                if (j_1p < 0)
                    continue;
                r2 = clebsch_squared(j_1, 4, j_1p);
                r2 *= pop[j_1];
                if (r2 < DBL_EPSILON)
                    continue;
                f0 = Erot(j_1p, O2_B_ROT, O2_D_ROT) - f1;
                for (i = 0; i < nblock; ++i)
                    kb[i] += r2 * EBCfunc(&EBCpars_5440, fb[i] - f0);
            }
        }
        /*
         * Multiply k through by factors outside the summation.
         */
        r1 = H_ON_KB / T;
        for (i = 0; i < nblock; ++i)
            kb[i] *= fb[i] * (1.0 - am_exp(-fb[i] * r1));
    }
    return 0;
}   /* O2O2_cia() */


/***********************************************************
* static double clebsch_squared(
*         const int j, const int lambda, const int jp)
*
* Purpose:
*   Computes a squared Clebsch-Gordan coefficient of the
*   form C(j lambda jp; 00)^2 needed for computing
*   collision-induced dipole moments.  Here, lambda is the
*   multipole order of the transition moment of the active
*   molecule, making a rotational transition j -> jp.  See
*   Birnbaum and Cohen, Can. J. Phys. 54:593, Table 1.  Also
*   Landau and Lifshitz, _Quantum_Mechanics_, sec. 106.
*
* Arguments:
*   const int j         - initial state rotational quantum
*                         number
*   const int lambda    - transition multipole order, 2 for
*                         quadrupole, 4 for hexadecapole.
*   const int jp        - final state rotational quantum
*                         number
*
* Return:
*   Squared Clebsch-Gordan coefficient C(j lambda jp; 00)^2
************************************************************/
static double clebsch_squared(
        const int j, const int lambda, const int jp)
{
    double csq, j_1, j_2;

    j_1 = (double)j;
    j_2 = j_1 * 2.;
    if (j < 0 || jp < 0) {
        errlog(80, 0);
        return 0.;
    }
    if (lambda == 2) {
        switch (jp - j) {
        case -2:
            csq = 3. * (j_1 - 1.) * j_1;
            csq /= 2. * (j_2 - 1.) * (j_2 + 1.);
            break;
        case 0:
            csq = j_1 * (j_1 + 1.);
            csq /= (j_2 - 1.) * (j_2 + 3.);
            break;
        case 2:
            csq = 3. * (j_1 + 1.) * (j_1 + 2.);
            csq /= 2. * (j_2 + 1.) * (j_2 + 3.);
            break;
        default:
            csq = 0.;
            break;
        }
    } else if (lambda == 4) {
        switch (jp - j) {
        case -4:
            csq = 35. * (j_1 - 3.) * (j_1 - 2.) * (j_1 - 1.) * j_1;
            csq /= 8. * (j_2 - 5.) * (j_2 - 3.) * (j_2 - 1.) * (j_2 + 1);
            break;
        case -2:
            csq = 5. * (j_1 - 2.) * (j_1 - 1.) * j_1 * (j_1 + 1.);
            csq /= 2. * (j_2 - 5.) * (j_2 - 1.) * (j_2 + 1.) * (j_2 + 3.);
            break;
        case 0:
            csq = 9. * (j_1 - 1.) * j_1 * (j_1 + 1.) * (j_1 + 2.);
            csq /= 4. * (j_2 - 3.) * (j_2 - 1.) * (j_2 + 3.) * (j_2 + 5.);
            break;
        case 2:
            csq = 5. * j_1 * (j_1 + 1.) * (j_1 + 2.) * (j_1 + 3.);
            csq /= 2. * (j_2 - 1.) * (j_2 + 1.) * (j_2 + 3.) * (j_2 + 7.);
            break;
        case 4:
            csq = 35. * (j_1 + 1.) * (j_1 + 2.) * (j_1 + 3.) * (j_1 + 4);
            csq /= 8. * (j_2 + 1.) * (j_2 + 3.) * (j_2 + 5.) * (j_2 + 7);
            break;
        default:
            csq = 0.;
            break;
        }
    } else {
        errlog(79, lambda);
        return 0.;
    }
    return csq;
}   /* clebsch_squared() */


/***********************************************************
* static double Erot(
*         const int j, const double B, const double D)
*
* Purpose:
*   Calculates rotational energy levels of a diatomic
*   molecule:
*
*   E(j) = B * j(j+1) - D * [j(j+1)]^2
*
* Arguments:
*   int j    - rotational quantum number
*   double B - rotational constant
*   double D - rotational constant
*
* Return:
*   rotational energy, in the same units as B, D
************************************************************/

static double Erot(
        const int j, const double B, const double D)
{
    double x;

    x  = (double)j;
    x *= x + 1.0;
    return x * (B - D * x);
}   /* Erot() */


/***********************************************************
* static void initialize_EBCpars(
*         EBCpars_t *EBCpars, const double T)
*
* Purpose:
*   Initializes an EBCpars_t structure for computation of
*   the extended Birnbaum-Cohen lineshape at temperature T.
*
* Arguments:
*   EBCpars_t *EBCpars  - pointer to line shape parameter
*                         structure
*   const double T      - Temperature [K]
************************************************************/

static void initialize_EBCpars(
        EBCpars_t *EBCpars, const double T)
{
    double x, x2;

    x = log(T);
    x2 = x * x;
    EBCpars->S = EBCpars->S_A * am_exp(EBCpars->S_B * x + EBCpars->S_C * x2);
    EBCpars->tau1 =
        EBCpars->tau1_A * am_exp(EBCpars->tau1_B * x + EBCpars->tau1_C * x2);
    EBCpars->tau2 =
        EBCpars->tau2_A * am_exp(EBCpars->tau2_B * x + EBCpars->tau2_C * x2);
    EBCpars->eps =
        EBCpars->eps_A * am_exp(EBCpars->eps_B * x + EBCpars->eps_C * x2);
    EBCpars->tau0 = H_ON_KB / (4.0 * PI * T);
    EBCpars->r0  = EBCpars->S / (PI * (1.0 + EBCpars->eps));
    EBCpars->r1  = TWOPI * EBCpars->tau0;
    EBCpars->r2  = TWOPI * EBCpars->tau1;
    EBCpars->r2 *= EBCpars->r2;
    EBCpars->r3  = EBCpars->tau2 / EBCpars->tau1;
    EBCpars->r4 =
        sqrt(EBCpars->tau0 * EBCpars->tau0 + EBCpars->tau2 * EBCpars->tau2);
    EBCpars->r4 /= EBCpars->tau1;
    if (EBCpars->eps < DBL_EPSILON)
        return;
    EBCpars->tau1p =
        EBCpars->tau1p_A *
        am_exp(EBCpars->tau1p_B * x + EBCpars->tau1p_C * x2);
    EBCpars->tau2p =
        EBCpars->tau2p_A *
        am_exp(EBCpars->tau2p_B * x + EBCpars->tau2p_C * x2);
    EBCpars->r5 = EBCpars->tau1p / EBCpars->tau2p;
    EBCpars->r6 =
        sqrt(EBCpars->tau0 * EBCpars->tau0 + EBCpars->tau1p * EBCpars->tau1p);
    EBCpars->r6 /= EBCpars->tau2p;
    EBCpars->r7  = TWOPI * EBCpars->tau2p;
    EBCpars->r7 *= EBCpars->r7;
    return;
}   /* initialize_EBCpars() */


/***********************************************************
* static double EBCfunc(
*         const EBCpars_t *EBCpars, const double f)
*
* Purpose:
*   Computes the EBC lineshape function at frequency f,
*   given a table of parameters and precomputed intermediate
*   results in the structure EBCpars.
*
* Arguments:
*   const EBCpars_t *EBCpars    - pointer to line shape
*                                 parameter structure
*   const double f              - frequency [GHz]
*
* Return:
*   EBC lineshape function at frequency f
************************************************************/

static double EBCfunc(
        const EBCpars_t *EBCpars, const double f)
{
    double f2, g1, g2, z;

    f2  = f * f;
    g1  = 1.0 + f2 * EBCpars->r2;
    z   = EBCpars->r4 * sqrt(g1);
    g1  = EBCpars->tau1 * am_xK1(z) / g1;
    g1 *= am_exp(EBCpars->r3 + EBCpars->r1 * f);
    if (EBCpars->eps >= DBL_EPSILON) {
        z   = EBCpars->r6 * sqrt(1.0 + f2 * EBCpars->r7);
        g2  = EBCpars->tau1p * am_K0(z);
        g2 *= am_exp(EBCpars->r5 + EBCpars->r1 * f);
        g1 += EBCpars->eps * g2;
    }
    return EBCpars->r0 * g1;
}   /* EBCfunc() */
