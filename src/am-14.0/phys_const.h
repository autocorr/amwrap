/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* phys_const.h                    S. Paine rev. 2024 July 15
*
* Physical constants
************************************************************/

#ifndef AM_PHYS_CONST_H
#define AM_PHYS_CONST_H

/*
 * Physical Constants (CODATA 2014) in SI units.  See
 *
 *  P. J. Mohr, D. B. Newell, and B. N. Taylor 2016, "CODATA
 *  recommended values of the fundamental physical constants:
 *  2014."  Rev. Mod. Phys. 88:035009.
 *
 * Apart from these, the other physical constants defined
 * in this file are in am's internal units.
 */
#define PCONST_C    299792458.       /* speed of light          [m s^-1] */
#define PCONST_H    6.626070040e-34  /* Planck constant            [J s] */
#define PCONST_KB   1.38064852e-23   /* Boltzmann constant      [J K^-1] */
#define PCONST_AMU  1.660539040e-27  /* Unified atomic mass unit    [kg] */

/*
 * Standard temperature and pressure
 */
#define T_STP  273.15   /* Standard temperature     [K] */
#define P_STP  1013.25  /* Standard pressure     [mbar] */

/*
 * Earth
 */
#define G_STD   980.665 /* Standard acceleration of gravity [cm s^-2] */
#define R_EARTH 6.37e8  /* Mean radius [cm]                          */

/*
 * Latent heats of sublimation and vaporization for water at 0 C.
 * From Eq.(5) and Eq.(9) of:
 *
 *     D. M. Murphy and T. Koop 2005, "Review of the vapour
 *     pressures of ice and supercooled water for atmospheric
 *     applications."  Q. J. R. Meteorol. Soc. 131:1539.
 *
 * converted from [J mol^-1] to [J kg^-1]
 */
#define H2O_LATENT_HEAT_SUB_0C 2.834e6 /* water ice at 0 C    [J kg^-1] */ 
#define H2O_LATENT_HEAT_VAP_0C 2.500e6 /* liquid water at 0 C [J kg^-1] */ 

/*
 * Radio refractivity constants.  See
 *
 *   M. Bevis, et al. "GPS Meteorology: Mapping Zenith Wet Delays
 *   onto Precipitable Water," J. Appl. Meteo. 33:379 (1994).
 *
 * These are updated coefficients for the Smith-Weintraub
 * equation:
 *
 *   n - 1 = k1 * (Pdry / T) + k2 * (Ph2o / T) + k2 * (Ph2o / T^2)
 *
 * where Pdry, Ph2o are in mbar, and T is in K.  See
 *
 *   E. K. Smith and S. Weintraub 1953, "The Constants in the
 *   Equation of Atmospheric Refractive Index at Radio
 *   Frequencies," Proc. IRE 41:1035 (1953).
 */
#define RADIO_REFRACTIVITY_K1 77.6e-6   /* [K mbar^-1]   */
#define RADIO_REFRACTIVITY_K2 70.4e-6   /* [K mbar^-1]   */
#define RADIO_REFRACTIVITY_K3 3.739e-1  /* [K^2 mbar^-1] */

/*
 * Optical refractivity is approximated by the first two
 * Smith-Weintraub terms, with an adjustment to the dry term.
 * This adjustment mainly accounts for the refractivity added by
 * the O2 spin-rotation band.  See
 *
 *   G. D. Thayer, "An improved equation for the radio refractive
 *   index of air," Radio Science 9:803 (Oct 1974).
 * 
 * and
 *
 *   R. J. Hill, R. S. Lawrence, and J. T. Priestley,
 *   "Theoretical and calculational aspects of the radio
 *   refractive index of water vapor," Radio Science 17:1251
 *   (Sep-Oct 1982).
 *
 * Thayer notes that for dry air at STP
 *
 *   10^6 * (n_radio - 1) = Na = 288.04(5)
 *
 * and
 *
 *   10^6 * (n_optical - 1) = Na' = 287.82(5)
 *
 * The difference is Na - Na' = 0.22(7) This is consistent with
 * the dry refractivity difference between 0 Hz and 15 THz
 * computed by Hilbert transform using am (0.231), which is used
 * for the adjustment here.
 */
#define OPTICAL_REFRACTIVITY_K1 (RADIO_REFRACTIVITY_K1 * (1.0 - 0.231/288.04))
#define OPTICAL_REFRACTIVITY_K2 RADIO_REFRACTIVITY_K2

/*
 * Optical refractivity coefficients [mbar K^-1] are converted to
 * volume refractivity in [cm^3] by multiplying by Boltzmann's
 * constant in suitable units, noting that
 *
 *      kB [mbar cm^3 K^-1] = 1e4 kB [Pa m^3 K^-1] = 1e4 kB [J K^-1]
 */
#define DRY_OPTICAL_VOLUME_REFRACTIVITY \
    (1e4 * PCONST_KB * OPTICAL_REFRACTIVITY_K1)
#define H2O_OPTICAL_VOLUME_REFRACTIVITY \
    (1e4 * PCONST_KB * OPTICAL_REFRACTIVITY_K2)

/*
 * Below are combinations of constants and internal unit conversions.
 */

/*
 * Loschmidt constant [cm^-3]
 */
#define N_STP (1.0e-4 * P_STP / (PCONST_KB * T_STP))

/*
 * (h / k)  [K GHz^-1]
 */
#define H_ON_KB (1.0e9 * PCONST_H / PCONST_KB)

/*
 * (2h / c^2)  [watt cm^-2 GHz^-4]
 */
#define TWOH_ON_CSQUARED (2.0e32 * PCONST_H / (PCONST_C * PCONST_C))

/*
 * (2k / c^2)  [watt K^-1 cm^-2 GHz^-3]
 */
#define TWOKB_ON_CSQUARED (2.0e23 * PCONST_KB / (PCONST_C * PCONST_C))

/*
 * (2k / (1 amu * c^2))  [K^-1]
 */
#define TWOKB_ON_AMU_CSQUARED \
                         (2.0 * PCONST_KB / (PCONST_AMU * PCONST_C * PCONST_C))

/*
 * c [cm GHz]
 */
#define C_CM_GHZ        (1e-7 * PCONST_C)

/*
 * c / (4 pi)  [cm GHz]
 */
#define C_ON_FOURPI      (7.957747154594767e-9 * PCONST_C)

/*
 * (4 pi) / c  [(cm GHz)^-1]
 */
#define FOURPI_ON_C      (1.256637061435917e8 / PCONST_C)

/*
 * Some expressions do not evaluate correctly in the low
 * frequency limit f->0.  F_EPSILON is the threshold for special
 * handling in several such instances.
 */
#ifndef F_EPSILON
    #define F_EPSILON   1.0e-6 /* GHz */
#endif

#endif /* AM_PHYS_CONST_H */
