/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* h2o_psat.c                       S. Paine rev. 2021 June 9
*
* Saturation vapor pressure for liquid water and ice.
************************************************************/

#include <math.h>
#include <stdio.h>

#include "am_types.h"
#include "column.h"
#include "errlog.h"
#include "h2o_psat.h"
#include "phys_const.h"


/***********************************************************
* double H2O_liquid_Psat(double T, int col_typenum)
*
* Purpose:
*   Computes H2O saturation vapor pressure over a flat
*   liquid surface at temperature T.  This is the customary
*   meteorological definition for reporting RH, even at
*   sub-freezing temperature, and for vapor in equilibrium
*   with small droplets.  The formulation used in this
*   function is  Eq. 10 of
*
*     D. M. Murphy and T. Koop 2005, "Review of the vapour
*     pressures of ice and supercooled water for atmospheric
*     applications."  Q. J. R. Meteorol. Soc. 131:1539.
*
*   which has a stated range of applicability of 123 K < T <
*   332 K.
*
*   Standard isotopic composition is assumed, and
*   differences in equilibrium isotopic fraction between
*   phases are ignored.  This is an approximation.  In
*   reality, in equilibrium, the condensed phase is enriched
*   in heavier species relative to the vapor phase by a
*   factor which increases with decreasing temperature.  At
*   the spontaneous nucleation point near -40 C, this factor
*   is approximately 1.02 for H2_18O / H2_16O, and
*   approximately 1.4 for HD_16O / H2_16O.  See
*
*     M. Bolot, B. Legras, and E. J. Moyer 2013, "Modelling
*     and interpreting the isotopic composition of water
*     vapour in convective updrafts."  Atmos. Chem. Phys.
*     13:7903.
*
* Arguments:
*   T           - temperature [K]
*   col_typenum - index number of column type
*
* Return:
*   Psat [mbar]
************************************************************/

double H2O_liquid_Psat(double T, int col_typenum)
{
    double r1, r2, frac;
    double ln_Psat, ln_T;

    /*
     * Log a warning if T is outside the applicable range of
     * the Murphy-Koop formula.
     */
    if (T < 123. || T > 332.)
        errlog(116, 0);
    /*
     * The following formula gives log(Psat), with Psat in Pa.
     */
    ln_T = log(T);
    r1 = 54.842763 - 6763.22 / T - 4.210   * ln_T + 0.000367 * T;
    r2 = 53.878    - 1331.22 / T - 9.44523 * ln_T + 0.014025 * T;
    ln_Psat = r1 + r2 * tanh(0.0415 * (T - 218.8));
    switch (col_typenum) {
    case COL_TYPE_H2O:
    case COL_TYPE_H2O_OPTICAL_REFRACTIVITY:
    case COL_TYPE_H2O_CONTINUUM:
    case COL_TYPE_H2O_AIR_CONTINUUM:
    case COL_TYPE_H2O_SELF_CONTINUUM:
    case COL_TYPE_H2O_LINES:
        frac = h2o_abundance_tab[0];
        break;
    case COL_TYPE_H2_16O:
        frac = h2o_abundance_tab[1];
        break;
    case COL_TYPE_H2_18O:
        frac = h2o_abundance_tab[2];
        break;
    case COL_TYPE_H2_17O:
        frac = h2o_abundance_tab[3];
        break;
    case COL_TYPE_HD_16O:
        frac = h2o_abundance_tab[4];
        break;
    case COL_TYPE_HD_18O:
        frac = h2o_abundance_tab[5];
        break;
    case COL_TYPE_HD_17O:
        frac = h2o_abundance_tab[6];
        break;
    default:
        errlog(117, col_typenum);
        frac = 1.0;
        break;
    }
    /*
     * Convert Psat to mbar.
     */
    return frac * 0.01 * exp(ln_Psat);
}   /* H2O_liquid_Psat() */


/***********************************************************
* double H2O_ice_Psat(double T, int col_typenum)
*
* Purpose:
*   Computes H2O saturation vapor pressure over ice at
*   temperature T.  The formulation used in this function is
*   Eq. 7 of
*
*     D. M. Murphy and T. Koop 2005, "Review of the vapour
*     pressures of ice and supercooled water for atmospheric
*     applications."  Q. J. R. Meteorol. Soc. 131:1539.
*
*   This formula applies to hexagonal ice Ih, and is stated
*   to be valid from 111 K to the triple point at 273.16 K.
*
*   Standard isotopic composition is assumed, and
*   differences in equilibrium isotopic fraction between
*   phases are ignored.  This is an approximation.  In
*   reality, in equilibrium, the condensed phase is enriched
*   in heavier species relative to the vapor phase by a
*   factor which increases with decreasing temperature.  At
*   -80 C, this factor is approximately 1.03 for H2_18O /
*   H2_16O, and approximately 1.4 for HD_16O / H2_16O.  See
*
*     M. Bolot, B. Legras, and E. J. Moyer 2013, "Modelling
*     and interpreting the isotopic composition of water
*     vapour in convective updrafts."  Atmos. Chem. Phys.
*     13:7903.
*
* Arguments:
*   T           - temperature [K]
*   col_typenum - index number of column type
*
* Return:
*   Psat [mbar]
************************************************************/

double H2O_ice_Psat(double T, int col_typenum)
{
    double frac;
    double ln_Psat, ln_T;

    /*
     * Log a warning if T is outside the applicable range of
     * the Murphy-Koop formula.
     */
    if (T < 111. || T > 273.16)
        errlog(126, 0);
    /*
     * The following formula gives log(Psat), with Psat in Pa.
     */
    ln_T = log(T);
    ln_Psat = 9.550426 - 5723.265 / T + 3.53068 * ln_T - 0.00728332 * T;
    switch (col_typenum) {
    case COL_TYPE_H2O:
    case COL_TYPE_H2O_OPTICAL_REFRACTIVITY:
    case COL_TYPE_H2O_CONTINUUM:
    case COL_TYPE_H2O_AIR_CONTINUUM:
    case COL_TYPE_H2O_SELF_CONTINUUM:
    case COL_TYPE_H2O_LINES:
        frac = h2o_abundance_tab[0];
        break;
    case COL_TYPE_H2_16O:
        frac = h2o_abundance_tab[1];
        break;
    case COL_TYPE_H2_18O:
        frac = h2o_abundance_tab[2];
        break;
    case COL_TYPE_H2_17O:
        frac = h2o_abundance_tab[3];
        break;
    case COL_TYPE_HD_16O:
        frac = h2o_abundance_tab[4];
        break;
    case COL_TYPE_HD_18O:
        frac = h2o_abundance_tab[5];
        break;
    case COL_TYPE_HD_17O:
        frac = h2o_abundance_tab[6];
        break;
    default:
        errlog(127, col_typenum);
        frac = 1.0;
        break;
    }
    /*
     * Convert Psat to mbar.
     */
    return frac * 0.01 * exp(ln_Psat);
}   /* H2O_ice_Psat() */
