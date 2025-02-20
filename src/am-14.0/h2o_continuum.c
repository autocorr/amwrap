/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* h2o_continuum.c               S. Paine rev. 2023 August 30
*
* H2O continuum absorption, from MT_CKD.  MT_CKD is
* maintained by the Radiative Transfer Group at AER, Inc.
* For MT_CKD copyright, reference, and version information,
* see the file mt_ckd.c.
************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "am_types.h"
#include "errlog.h"
#include "h2o_continuum.h"
#include "mt_ckd.h"
#include "phys_const.h"
#include "units.h"

int interpolate_continuum_table(
    double*,
    const double*,
    const gridsize_t,
    const double,
    const double,
    const double,
    const double*,
    const double*,
    const double*,
    const int);

/***********************************************************
* int H2O_air_continuum(
*   double *kb,
*   const double *f,
*   const gridsize_t ngrid,
*   const double T)
*
* Purpose:
*   Computes the H2O air-induced continuum absorption
*   coefficient.
*
* Arguments:
*   double *kb             - binary absorption coeff. [cm^5]
*   const double *f        - frequency grid [GHz]
*   const gridsize_t ngrid - number of grid points
*   const double T         - temperature [K]
*
* Return:
*   0 if OK, 1 if an error occurred
************************************************************/

int H2O_air_continuum(
    double *kb,
    const double *f,
    const gridsize_t ngrid,
    const double T)
{
    return interpolate_continuum_table(
            kb,
            f,
            ngrid,
            T,
            mt_ckd_ref_press,
            mt_ckd_ref_temp,
            mt_ckd_for_absco_ref,
            NULL,
            mt_ckd_wavenumbers,
            mt_ckd_ntab);
}   /* H2O_air_continuum() */


/***********************************************************
* int H2O_self_continuum(
*   double *kb,
*   const double *f,
*   const gridsize_t ngrid,
*   const double T)
*
* Purpose:
*   Computes the H2O self-induced continuum absorption
*   coefficient.
*
* Arguments:
*   double *kb             - binary absorption coeff. [cm^5]
*   const double *f        - frequency grid [GHz]
*   const gridsize_t ngrid - number of grid points
*   const double T         - temperature [K]
*
* Return:
*   0 if OK, 1 if an error occurred
************************************************************/

int H2O_self_continuum(
    double *kb,
    const double *f,
    const gridsize_t ngrid,
    const double T)
{
    return interpolate_continuum_table(
            kb,
            f,
            ngrid,
            T,
            mt_ckd_ref_press,
            mt_ckd_ref_temp,
            mt_ckd_self_absco_ref,
            mt_ckd_self_texp,
            mt_ckd_wavenumbers,
            mt_ckd_ntab);
}   /* H2O_self_continuum() */


/***********************************************************
* int interpolate_continuum_table(
*   double *kb,
*   const double *f,
*   const gridsize_t ngrid,
*   const double T,
*   const double ref_press,
*   const double ref_temp,
*   const double *absco_ref,
*   const double *texp,
*   const double *wavenumbers,
*   const int ntab)
*
* Purpose:
*   Computes a binary spectral absorption coefficient
*   kb[] [cm^5] from a tabulated MT_CKD continuum 
*   spectral absorption coefficient [1/(cm^-1 molec/cm^2)],
*   and interpolates the result in frequency onto the
*   model frequency grid f[].
*
* Arguments:
*   double *kb              - binary absorption coefficient
*                             [cm^5]
*   const double *f         - model frequency grid [GHz]
*   const gridsize_t ngrid  - number of model grid points
*   const double T          - temperature [K]
*   const double ref_press  - reference pressure [mbar]
*   const double ref_temp   - reference temperature [K]
*   const double *absco_ref - continuum absorption coeff.
*                             [1/(cm^-1 molec/cm^2)] at
*                             reference press. and temp.
*   const double *texp      - temperature dependence
*                             exponent for continuum
*                             (NULL if none)
*   const double *wavenumbers
*                           - tabulated wavenumbers [cm^-1]
*   const int ntab          - number of table entries
*
* Return:
*   0 if OK, 1 if an error occurred
************************************************************/

int interpolate_continuum_table(
    double *kb,
    const double *f,
    const gridsize_t ngrid,
    const double T,
    const double ref_press,
    const double ref_temp,
    const double *absco_ref,
    const double *texp,
    const double *wavenumbers,
    const int ntab)
{
    gridsize_t i;
    int j, jmin, jmax;
    double df;
    double r, rho0;
    double *kb_tab = NULL, *m_tab = NULL;

    /*
     * Compute the table range which will be needed.  Fail if
     * the model frequency grid runs beyond the table range.
     */
    df   = (wavenumbers[1] - wavenumbers[0]) * C_CM_GHZ;
    jmin = (int)floor(f[0] / df);
    jmax = (int)ceil(f[ngrid-1] / df);
    if (jmax > ntab - 2) {
        errlog(102, 0);
        return 1;
    }
    /*
     * Allocate temporary space for interpolation tables for the
     * binary spectral absorption coefficient adjusted to
     * temperature T.
     */
    if (    (kb_tab = (double*)malloc(ntab * sizeof(double))) == NULL ||
            (m_tab  = (double*)malloc(ntab * sizeof(double))) == NULL) {
        free(kb_tab);
        free(m_tab);
        errlog(201, 0);
        return 1;
    }
    /*
     * Multiplying the tabulated continuum coefficients in 
     * [1 / (cm^-1 molec/cm^2)] by the radiation term  
     *
     *   nu * tanh((h c nu) / (2 k T)),
     * 
     * with nu in wavenumbers [cm^-1], yields a molecular
     * absorption coefficient [cm^2] at a reference density rho0
     * that corresponds to the reference pressure and temperature
     * for MT_CKD.  Dividing by rho0 yields a binary absorption
     * coefficient [cm^5] as used in am. 
     */
    r    = H_ON_KB * C_CM_GHZ / (2.0 * T);
    rho0 = N_STP * (ref_press / P_STP) * (T_STP / ref_temp);
    for (j = jmin > 0 ? jmin - 1 : 0; j <= jmax + 1; ++j) {
        kb_tab[j] =
            wavenumbers[j] * tanh(r * wavenumbers[j]) * absco_ref[j] / rho0; 
            
    }
    /*
     * If the spectral temperature dependence exponent texp is
     * non-NULL, adjust the absorption coefficient by a factor
     *
     *   (ref_temp / T)^texp.
     */
    if (texp != NULL) {
        r = ref_temp / T;
        for (j = jmin > 0 ? jmin - 1 : 0; j <= jmax + 1; ++j)
            kb_tab[j] *= pow(r, texp[j]);
    }
    /*
     * Compute slopes for Catmull-Rom cubic Hermite spline
     * interpolation of kb_tab[].
     */
    if (jmin == 0)
        m_tab[jmin] = 0.0;
    else
        m_tab[jmin] = 0.5 * (kb_tab[jmin+1] - kb_tab[jmin-1]);
    for (j = jmin + 1; j <= jmax; ++j)
        m_tab[j] = 0.5 * (kb_tab[j+1] - kb_tab[j-1]);
    /*
     * Interpolate from tabulated frequencies onto the model
     * frequency grid.
     */
    for (i = 0; i < ngrid; ++i) {
        double p, jd;
        double a0, a1, b0, b1;
        p = modf(f[i] / df, &jd);
        j = (int)jd;
        a0 = p - 1.0;
        a1 = b1 = p * p;
        a1 *= (3.0 - 2.0 * p);
        b1 *= a0;
        a0 *= a0;
        b0 = p * a0;
        a0 *= (1.0 + 2.0 * p);
        kb[i] = a0 * kb_tab[j]   + b0 * m_tab[j]
              + a1 * kb_tab[j+1] + b1 * m_tab[j+1];
    }
    free(kb_tab);
    free(m_tab);
    return 0;
}   /* interpolate_continuum_table() */
