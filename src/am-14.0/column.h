/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* column.h                       S. Paine rev. 2024 April 20
*
* Constants and declarations for column.c
************************************************************/

#ifndef AM_COLUMN_H
#define AM_COLUMN_H

#include "am_types.h"

/*
 * Column type properties are listed in the table col_type[],
 * declared below and defined in column.c.
 *
 * common_unit is the unit index number for a customary unit,
 * used when reporting total column densities in the function
 * write_model_config_data().
 *
 * default_vmr is the default volume mixing ratio for the column
 * type For most species other than the constituents of dry air,
 * default_vmr = 0.0
 *
 * Mr is the relative molar mass for molecular column types.  For
 * blends of isotopologues, it is the abundance-weighted average.
 * For non-molecular column types, Mr == 0.
 *
 * k_type[] is an array of type numbers for the absorption
 * coefficients needed to compute the opacity of this column.
 * The maximum number of absorption coefficients which can be
 * associated with a column type is MAX_K_TYPES_PER_COL. 
 *
 * pmr[] is an array of partial mixing ratios (relative to total
 * mixing ratio) to support mixed species like dry_air.
 */

enum {
    MAX_K_TYPES_PER_COL = 10
};

extern struct col_type_tabentry {
    const char *name;   /* text name of column type                       */
    int common_unit;    /* unit number of commonly-used unit              */
    double default_vmr; /* default volume mixing ratio in dry air         */
    double Mr;  /* relative molar mass for molecular species, else 0      */
    int flags;  /* miscellaneous flags, see enum below                    */
    int zadep;  /* dependence of opacity on zenith angle, see enum below  */
    int n_abscoeffs;    /* number of absorption coeffs for this col type  */
    int k_type[MAX_K_TYPES_PER_COL]; /* list of absorption coeff typenums */
    double pmr[MAX_K_TYPES_PER_COL]; /* partial mixing ratios             */
} col_type[];

enum {
    COL_TYPE_NONE,
    COL_TYPE_ATTEN,
    COL_TYPE_ATTEN_AIRMASS,
    COL_TYPE_ATTEN_SEC_ZA,
    COL_TYPE_ATTEN_SIN_SQUARED_ZA,
    COL_TYPE_ATTEN_SIN_SQUARED_2ZA,
    COL_TYPE_CH4,
        COL_TYPE_12CH4,
        COL_TYPE_13CH4,
        COL_TYPE_12CH3D,
    COL_TYPE_CH3CN,
        COL_TYPE_12CH3_12C14N,
    COL_TYPE_CH3OH,
        COL_TYPE_12CH3_16OH,
    COL_TYPE_CO,
        COL_TYPE_12C_16O,
        COL_TYPE_13C_16O,
        COL_TYPE_12C_18O,
        COL_TYPE_12C_17O,
        COL_TYPE_13C_18O,
        COL_TYPE_13C_17O,
    COL_TYPE_CO2,
        COL_TYPE_12C_16O2,
        COL_TYPE_13C_16O2,
        COL_TYPE_16O_12C_18O,
        COL_TYPE_16O_12C_17O,
        COL_TYPE_16O_13C_18O,
        COL_TYPE_16O_13C_17O,
        COL_TYPE_12C_18O2,
    COL_TYPE_ClO,
        COL_TYPE_35Cl_16O,
        COL_TYPE_37Cl_16O,
    COL_TYPE_DRY_AIR,
    COL_TYPE_DRY_AIR_STD,
    COL_TYPE_DRY_AIR_OPTICAL_REFRACTIVITY,
    COL_TYPE_HBr,
        COL_TYPE_H_79Br,
        COL_TYPE_H_81Br,
    COL_TYPE_HCN,
        COL_TYPE_H_12C_14N,
        COL_TYPE_H_13C_14N,
        COL_TYPE_H_12C_15N,
    COL_TYPE_H2CO,
        COL_TYPE_H2_12C_16O,
        COL_TYPE_H2_13C_16O,
        COL_TYPE_H2_12C_18O,
    COL_TYPE_HCl,
        COL_TYPE_H_35Cl,
        COL_TYPE_H_37Cl,
    COL_TYPE_HF,
        COL_TYPE_H_19F,
    COL_TYPE_HNO3,
        COL_TYPE_H_14N_16O3,
    COL_TYPE_H2O,
    COL_TYPE_H2_16O_LINES_PLUS_CONTINUUM,
    COL_TYPE_H2O_OPTICAL_REFRACTIVITY,
    COL_TYPE_H2O_CONTINUUM,
    COL_TYPE_H2O_AIR_CONTINUUM,
    COL_TYPE_H2O_SELF_CONTINUUM,
    COL_TYPE_H2O_LINES,
        COL_TYPE_H2_16O,
        COL_TYPE_H2_18O,
        COL_TYPE_H2_17O,
        COL_TYPE_HD_16O,
        COL_TYPE_HD_18O,
        COL_TYPE_HD_17O,
    COL_TYPE_H2O2,
        COL_TYPE_H2_16O2,
    COL_TYPE_HO2,
        COL_TYPE_H_16O2,
    COL_TYPE_HOCl,
        COL_TYPE_H_16O_35Cl,
        COL_TYPE_H_16O_37Cl,
    COL_TYPE_H2S,
        COL_TYPE_H2_32S,
        COL_TYPE_H2_34S,
        COL_TYPE_H2_33S,
    COL_TYPE_IWP_ABS_RAYLEIGH,
    COL_TYPE_LWP_ABS_RAYLEIGH,
    COL_TYPE_NH3,
        COL_TYPE_14NH3,
        COL_TYPE_15NH3,
    COL_TYPE_N2N2,
    COL_TYPE_N2AIR,
    COL_TYPE_N2O,
        COL_TYPE_14N2_16O,
        COL_TYPE_14N_15N_16O,
        COL_TYPE_15N_14N_16O,
        COL_TYPE_14N2_18O,
        COL_TYPE_14N2_17O,
    COL_TYPE_NO,
        COL_TYPE_14N_16O,
    COL_TYPE_NO2,
        COL_TYPE_14N_16O2,
    COL_TYPE_O,
        COL_TYPE_16O,
    COL_TYPE_O2,
        COL_TYPE_16O2,
        COL_TYPE_16O_18O,
        COL_TYPE_16O_17O,
    COL_TYPE_O2_COUPLED,
        COL_TYPE_16O2_COUPLED,
    COL_TYPE_O2_UNCOUPLED,
        COL_TYPE_16O2_UNCOUPLED,
        COL_TYPE_16O_18O_UNCOUPLED,
        COL_TYPE_16O_17O_UNCOUPLED,
    COL_TYPE_O2O2,
    COL_TYPE_O2AIR,
    COL_TYPE_O3,
        COL_TYPE_16O3,
        COL_TYPE_16O2_18O,
        COL_TYPE_16O_18O_16O,
        COL_TYPE_16O2_17O,
        COL_TYPE_16O_17O_16O,
    COL_TYPE_OCS,
        COL_TYPE_16O_12C_32S,
        COL_TYPE_16O_12C_34S,
        COL_TYPE_16O_13C_32S,
        COL_TYPE_16O_12C_33S,
        COL_TYPE_18O_12C_32S,
    COL_TYPE_OH,
        COL_TYPE_16O_H,
        COL_TYPE_18O_H,
        COL_TYPE_16O_D,
    COL_TYPE_SO2,
        COL_TYPE_32S_16O2,
        COL_TYPE_34S_16O2,
    COL_TYPE_ONELINE,
    COL_TYPE_END_OF_TABLE
};

/*
 * Column flag bits
 *
 * COL_PARAMETRIC is a constant flag indicating that this is a
 * parametric column type, such as atten, for which N is
 * interpreted as a parameter rather than as a column density.
 *
 * COL_OPTICAL_DELAY is a constant flag for column types which
 * include a constant excess delay, proportional to N.
 *
 * COL_NO_ARRAYS is a flag to skip array allocation on this
 * column, as for delay-only column types.
 *
 * COL_RH_ALLOWED is a flag that indicates that a relative
 * humidity specification is allowed for this column.
 */
enum {
    COL_PARAMETRIC    = 0x1,
    COL_OPTICAL_DELAY = 0x2,
    COL_NO_ARRAYS     = 0x4,
    COL_RH_ALLOWED    = 0x8
};

/*
 * Zenith angle dependence of line-of-sight opacity vs. zenith
 * opacity for a given column type. The common case is
 * ZADEP_AIRMASS, which is synonymous with ZADEP_SEC_ZA for
 * plane-parallel geometry.
 */
enum {
    ZADEP_NONE,           /* tau(los) = tau(zenith)                  */
    ZADEP_AIRMASS,        /* tau(los) = tau(zenith) *                */
                          /*            geometric_path / zenith_path */
    ZADEP_SEC_ZA,         /* tau(los) = tau(zenith) * sec(za)        */
    ZADEP_SIN_SQUARED_ZA, /* tau(los) = tau(zenith) * sin^2(za)      */
    ZADEP_SIN_SQUARED_2ZA /* tau(los) = tau(zenith) * sin^2(2za)     */
};

/*
 * Bits in column_t.N_mode, which controls column density
 * computation and presentation in output config data.
 */
enum {
    N_EXPLICIT     = 0x1, /* Column density specified explicitly        */
    N_BY_DISTANCE  = 0x2, /* by vmr, using layer height                 */
    N_HYDROSTATIC  = 0x4, /* by vmr, using hydrostatic column density   */
    N_INTERPOLATED = 0x8  /* N was split across interpolated layers     */
};

enum {
    N_BY_VMR = (N_BY_DISTANCE | N_HYDROSTATIC)
};

int    column_opacity(model_t*, model_t*, int, int, int*);
double column_optical_delay(model_t*, int, int);

#endif /* AM_COLUMN_H */
