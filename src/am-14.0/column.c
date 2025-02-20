/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* column.c                       S. Paine rev. 2024 April 20
*
* Definitions of column types, and computation of column
* opacities.
************************************************************/

#include <float.h>
#include <math.h>
#include <stdlib.h>

#include "abscoeff.h"
#include "am_types.h"
#include "cia.h"
#include "column.h"
#include "errlog.h"
#include "layer.h"
#include "molecules.h"
#include "phys_const.h"
#include "spectra.h"
#include "units.h"


/*
 * Column types
 *
 * Each column can be associated with up to MAX_K_TYPES_PER_COL
 * (defined in column.h) absorption coefficients.  The array k_type[]
 * in the col_type_tabentry structure for each column type contains
 * the associated absorption coefficient type numbers.  Uninitialized
 * elements of ktype[] default to 0, which is equal to K_TYPE_NONE.
 *
 * The array pmr[] is an array of partial mixing ratios corresponding
 * to each of the absorption coefficients listed in k_type[].
 * This supports composite column types such as dry_air; for such
 * column types, the overall mixing ratio multiplied by the partial
 * mixing ratio for each constituent is the actual mixing ratio
 * for that constituent.
 */
struct col_type_tabentry col_type[] = {
    {
        "none",        /* name          */
        UNIT_NONE,     /* common_unit   */
        0.0,           /* default_vmr   */
        0.0,           /* Mr            */
        0,             /* flags         */
        ZADEP_NONE,    /* zadep         */
        0,             /* n_abscoeffs   */
        {K_TYPE_NONE}, /* k_type[]      */
        {1.0}          /* pmr[]         */
    },
    {
        "atten",
        UNIT_NEPER,
        0.0,
        0.0,
        COL_PARAMETRIC,
        ZADEP_NONE,
        0,
        {K_TYPE_NONE},
        {1.0}
    },
    {
        "atten_airmass",
        UNIT_NEPER,
        0.0,
        0.0,
        COL_PARAMETRIC,
        ZADEP_AIRMASS,
        0,
        {K_TYPE_NONE},
        {1.0}
    },
    {
        "atten_sec(za)",
        UNIT_NEPER,
        0.0,
        0.0,
        COL_PARAMETRIC,
        ZADEP_SEC_ZA,
        0,
        {K_TYPE_NONE},
        {1.0}
    },
    {
        "atten_sin^2(za)",
        UNIT_NEPER,
        0.0,
        0.0,
        COL_PARAMETRIC,
        ZADEP_SIN_SQUARED_ZA,
        0,
        {K_TYPE_NONE},
        {1.0}
    },
    {
        "atten_sin^2(2za)",
        UNIT_NEPER,
        0.0,
        0.0,
        COL_PARAMETRIC,
        ZADEP_SIN_SQUARED_2ZA,
        0,
        {K_TYPE_NONE},
        {1.0}
    },
    {
        "ch4",
        AM_UNIT_COL_DENSITY,
        FRAC_CH4,
        MASS_CH4,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_CH4},
        {1.0}
    },
    {
        "12ch4",
        AM_UNIT_COL_DENSITY,
        FRAC_12CH4,
        MASS_12CH4,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_12CH4},
        {1.0}
    },
    {
        "13ch4",
        AM_UNIT_COL_DENSITY,
        FRAC_13CH4,
        MASS_13CH4,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_13CH4},
        {1.0}
    },
    {
        "12ch3d",
        AM_UNIT_COL_DENSITY,
        FRAC_12CH3D,
        MASS_12CH3D,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_12CH3D},
        {1.0}
    },
    {
        "ch3cn",
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_CH3CN,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_CH3CN},
        {1.0}
    },
    {
        "12ch3_12c14n",
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_12CH3_12C14N,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_12CH3_12C14N},
        {1.0}
    },
    {
        "ch3oh",
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_CH3OH,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_CH3OH},
        {1.0}
    },
    {
        "12ch3_16oh",
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_12CH3_16OH,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_12CH3_16OH},
        {1.0}
    },
    {
        "co",              
        AM_UNIT_COL_DENSITY,
        FRAC_CO,
        MASS_CO,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_CO},
        {1.0}
    },
    {
        "12c_16o",              
        AM_UNIT_COL_DENSITY,
        FRAC_12C_16O,
        MASS_12C_16O,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_12C_16O},
        {1.0}
    },
    {
        "13c_16o",              
        AM_UNIT_COL_DENSITY,
        FRAC_13C_16O,
        MASS_13C_16O,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_13C_16O},
        {1.0}
    },
    {
        "12c_18o",              
        AM_UNIT_COL_DENSITY,
        FRAC_12C_18O,
        MASS_12C_18O,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_12C_18O},
        {1.0}
    },
    {
        "12c_17o",              
        AM_UNIT_COL_DENSITY,
        FRAC_12C_17O,
        MASS_12C_17O,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_12C_17O},
        {1.0}
    },
    {
        "13c_18o",              
        AM_UNIT_COL_DENSITY,
        FRAC_13C_18O,
        MASS_13C_18O,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_13C_18O},
        {1.0}
    },
    {
        "13c_17o",              
        AM_UNIT_COL_DENSITY,
        FRAC_13C_17O,
        MASS_13C_17O,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_13C_17O},
        {1.0}
    },
    {
        "co2",              
        AM_UNIT_COL_DENSITY,
        FRAC_CO2,
        MASS_CO2,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_CO2},
        {1.0}
    },
    {
        "12c_16o2",              
        AM_UNIT_COL_DENSITY,
        FRAC_12C_16O2,
        MASS_12C_16O2,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_12C_16O2},
        {1.0}
    },
    {
        "13c_16o2",              
        AM_UNIT_COL_DENSITY,
        FRAC_13C_16O2,
        MASS_13C_16O2,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_13C_16O2},
        {1.0}
    },
    {
        "16o_12c_18o",              
        AM_UNIT_COL_DENSITY,
        FRAC_16O_12C_18O,
        MASS_16O_12C_18O,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_16O_12C_18O},
        {1.0}
    },
    {
        "16o_12c_17o",              
        AM_UNIT_COL_DENSITY,
        FRAC_16O_12C_17O,
        MASS_16O_12C_17O,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_16O_12C_17O},
        {1.0}
    },
    {
        "16o_13c_18o",              
        AM_UNIT_COL_DENSITY,
        FRAC_16O_13C_18O,
        MASS_16O_13C_18O,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_16O_13C_18O},
        {1.0}
    },
    {
        "16o_13c_17o",              
        AM_UNIT_COL_DENSITY,
        FRAC_16O_13C_17O,
        MASS_16O_13C_17O,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_16O_13C_17O},
        {1.0}
    },
    {
        "12c_18o2",              
        AM_UNIT_COL_DENSITY,
        FRAC_12C_18O2,
        MASS_12C_18O2,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_12C_18O2},
        {1.0}
    },
    {
        "clo",              
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_ClO,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_ClO},
        {1.0}
    },
    {
        "35cl_16o",              
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_35Cl_16O,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_35Cl_16O},
        {1.0}
    },
    {
        "37cl_16o",              
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_37Cl_16O,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_37Cl_16O},
        {1.0}
    },
    {
        "dry_air",
        AM_UNIT_COL_DENSITY,
        1.0, /* default mixing ratio for dry air */
        MASS_DRY_AIR,
        COL_OPTICAL_DELAY,
        ZADEP_AIRMASS,
        8,
        {K_TYPE_CH4, K_TYPE_CO, K_TYPE_CO2, K_TYPE_N2O,
            K_TYPE_O2_COUPLED, K_TYPE_O2_UNCOUPLED, K_TYPE_N2N2, K_TYPE_O2O2},
        {FRAC_CH4, FRAC_CO, FRAC_CO2, FRAC_N2O,
            FRAC_O2, FRAC_O2, FRAC_N2, FRAC_O2}
    },
    {
        "dry_air_std",
        AM_UNIT_COL_DENSITY,
        1.0, /* default mixing ratio for dry air */
        MASS_DRY_AIR,
        COL_OPTICAL_DELAY,
        ZADEP_AIRMASS,
        5,
        {K_TYPE_CO2,
            K_TYPE_O2_COUPLED, K_TYPE_O2_UNCOUPLED, K_TYPE_N2N2, K_TYPE_O2O2},
        {FRAC_CO2,
            FRAC_O2, FRAC_O2, FRAC_N2, FRAC_O2}
    },
    {
        "dry_air_optical_refractivity",
        AM_UNIT_COL_DENSITY,
        1.0,
        MASS_DRY_AIR,
        COL_OPTICAL_DELAY | COL_NO_ARRAYS,
        ZADEP_AIRMASS,
        0,
        {K_TYPE_NONE},
        {1.0}
    },
    {
        "hbr",              
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_HBr,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_HBr},
        {1.0}
    },
    {
        "h_79br",              
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_H_79Br,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_H_79Br},
        {1.0}
    },
    {
        "h_81br",              
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_H_81Br,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_H_81Br},
        {1.0}
    },
    {
        "hcn",              
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_HCN,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_HCN},
        {1.0}
    },
    {
        "h_12c_14n",              
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_H_12C_14N,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_H_12C_14N},
        {1.0}
    },
    {
        "h_13c_14n",              
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_H_13C_14N,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_H_13C_14N},
        {1.0}
    },
    {
        "h_12c_15n",              
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_H_12C_15N,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_H_12C_15N},
        {1.0}
    },
    {
        "h2co",              
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_H2CO,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_H2CO},
        {1.0}
    },
    {
        "h2_12c_16o",              
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_H2_12C_16O,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_H2_12C_16O},
        {1.0}
    },
    {
        "h2_13c_16o",              
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_H2_13C_16O,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_H2_13C_16O},
        {1.0}
    },
    {
        "h2_12c_18o",              
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_H2_12C_18O,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_H2_12C_18O},
        {1.0}
    },
    {
        "hcl",              
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_HCl,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_HCl},
        {1.0}
    },
    {
        "h_35cl",              
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_H_35Cl,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_H_35Cl},
        {1.0}
    },
    {
        "h_37cl",              
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_H_37Cl,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_H_37Cl},
        {1.0}
    },
    {
        "hf",              
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_HF,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_HF},
        {1.0}
    },
    {
        "h_19f",              
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_H_19F,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_H_19F},
        {1.0}
    },
    {
        "hno3",              
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_HNO3,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_HNO3},
        {1.0}
    },
    {
        "h_14n_16o3",              
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_H_14N_16O3,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_H_14N_16O3},
        {1.0}
    },
    {
        "h2o",
        UNIT_UM_PWV,
        0.0,
        MASS_H2O,
        COL_OPTICAL_DELAY | COL_RH_ALLOWED,
        ZADEP_AIRMASS,
        3,
        {K_TYPE_H2O_LINES, K_TYPE_H2O_AIR_CONTINUUM, K_TYPE_H2O_SELF_CONTINUUM},
        {1.0, 1.0, 1.0}
    },
    {
        "h2_16o_lines_plus_continuum",
        UNIT_UM_PWV,
        0.0,
        MASS_H2_16O,
        COL_OPTICAL_DELAY | COL_RH_ALLOWED,
        ZADEP_AIRMASS,
        3,
        {K_TYPE_H2_16O, K_TYPE_H2O_AIR_CONTINUUM, K_TYPE_H2O_SELF_CONTINUUM},
        {1.0, 1.0, 1.0}
    },
    {
        "h2o_optical_refractivity",
        UNIT_UM_PWV,
        0.0,
        MASS_H2O,
        COL_OPTICAL_DELAY | COL_NO_ARRAYS | COL_RH_ALLOWED,
        ZADEP_AIRMASS,
        0,
        {K_TYPE_NONE},
        {1.0}
    },
    {
        "h2o_continuum",
        UNIT_UM_PWV,
        0.0,
        MASS_H2O,
        COL_RH_ALLOWED,
        ZADEP_AIRMASS,
        2,
        {K_TYPE_H2O_AIR_CONTINUUM, K_TYPE_H2O_SELF_CONTINUUM},
        {1.0, 1.0}
    },
    {
        "h2o_air_continuum",
        UNIT_UM_PWV,
        0.0,
        MASS_H2O,
        COL_RH_ALLOWED,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_H2O_AIR_CONTINUUM},
        {1.0}
    },
    {
        "h2o_self_continuum",
        UNIT_UM_PWV,
        0.0,
        MASS_H2O,
        COL_RH_ALLOWED,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_H2O_SELF_CONTINUUM},
        {1.0}
    },
    {
        "h2o_lines",
        UNIT_UM_PWV,
        0.0,
        MASS_H2O,
        COL_RH_ALLOWED,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_H2O_LINES},
        {1.0}
    },
    {
        "h2_16o",
        UNIT_UM_PWV,
        0.0,
        MASS_H2_16O,
        COL_RH_ALLOWED,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_H2_16O},
        {1.0}
    },
    {
        "h2_18o",
        UNIT_UM_PWV,
        0.0,
        MASS_H2_18O,
        COL_RH_ALLOWED,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_H2_18O},
        {1.0}
    },
    {
        "h2_17o",
        UNIT_UM_PWV,
        0.0,
        MASS_H2_17O,
        COL_RH_ALLOWED,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_H2_17O},
        {1.0}
    },
    {
        "hd_16o",
        UNIT_UM_PWV,
        0.0,
        MASS_HD_16O,
        COL_RH_ALLOWED,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_HD_16O},
        {1.0}
    },
    {
        "hd_18o",
        UNIT_UM_PWV,
        0.0,
        MASS_HD_18O,
        COL_RH_ALLOWED,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_HD_18O},
        {1.0}
    },
    {
        "hd_17o",
        UNIT_UM_PWV,
        0.0,
        MASS_HD_17O,
        COL_RH_ALLOWED,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_HD_17O},
        {1.0}
    },
    {
        "h2o2",
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_H2O2,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_H2O2},
        {1.0}
    },
    {
        "h2_16o2",
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_H2_16O2,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_H2_16O2},
        {1.0}
    },
    {
        "ho2",
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_HO2,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_HO2},
        {1.0}
    },
    {
        "h_16o2",
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_H_16O2,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_H_16O2},
        {1.0}
    },
    {
        "hocl",
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_HOCl,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_HOCl},
        {1.0}
    },
    {
        "h_16o_35cl",
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_H_16O_35Cl,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_H_16O_35Cl},
        {1.0}
    },
    {
        "h_16o_37cl",
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_H_16O_37Cl,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_H_16O_37Cl},
        {1.0}
    },
    {
        "h2s",
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_H2S,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_H2S},
        {1.0}
    },
    {
        "h2_32s",
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_H2_32S,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_H2_32S},
        {1.0}
    },
    {
        "h2_34s",
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_H2_34S,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_H2_34S},
        {1.0}
    },
    {
        "h2_33s",
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_H2_33S,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_H2_33S},
        {1.0}
    },
    {
        "iwp_abs_Rayleigh",
        UNIT_KG_ON_MSQUARED,
        0.0,
        0.0,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_IWP_ABS_RAYLEIGH},
        {1.0}
    },
    {
        "lwp_abs_Rayleigh",
        UNIT_KG_ON_MSQUARED,
        0.0,
        0.0,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_LWP_ABS_RAYLEIGH},
        {1.0}
    },
    {
        "nh3",
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_NH3,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_NH3},
        {1.0}
    },
    {
        "14nh3",
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_14NH3,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_14NH3},
        {1.0}
    },
    {
        "15nh3",
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_15NH3,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_15NH3},
        {1.0}
    },
    {
        "n2n2",
        AM_UNIT_COL_DENSITY,
        FRAC_N2,
        MASS_N2,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_N2N2},
        {1.0}
    },
    {
        "n2air", 
        AM_UNIT_COL_DENSITY,
        FRAC_N2,
        MASS_N2,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_N2N2},
        {1.0}
    },
    {
        "n2o",
        AM_UNIT_COL_DENSITY,
        FRAC_N2O,
        MASS_N2O,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_N2O},
        {1.0}
    },
    {
        "14n2_16o",
        AM_UNIT_COL_DENSITY,
        FRAC_14N2_16O,
        MASS_14N2_16O,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_14N2_16O},
        {1.0}
    },
    {
        "14n_15n_16o",
        AM_UNIT_COL_DENSITY,
        FRAC_14N_15N_16O,
        MASS_14N_15N_16O,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_14N_15N_16O},
        {1.0}
    },
    {
        "15n_14n_16o",
        AM_UNIT_COL_DENSITY,
        FRAC_15N_14N_16O,
        MASS_15N_14N_16O,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_15N_14N_16O},
        {1.0}
    },
    {
        "14n2_18o",
        AM_UNIT_COL_DENSITY,
        FRAC_14N2_18O,
        MASS_14N2_18O,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_14N2_18O},
        {1.0}
    },
    {
        "14n2_17o",
        AM_UNIT_COL_DENSITY,
        FRAC_14N2_17O,
        MASS_14N2_17O,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_14N2_17O},
        {1.0}
    },
    {
        "no",              
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_NO,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_NO},
        {1.0}
    },
    {
        "14n_16o",              
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_14N_16O,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_14N_16O},
        {1.0}
    },
    {
        "no2",              
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_NO2,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_NO2},
        {1.0}
    },
    {
        "14n_16o2",              
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_14N_16O2,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_14N_16O2},
        {1.0}
    },
    {
        "o",
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_O,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_O},
        {1.0}
    },
    {
        "16o",
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_16O,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_16O},
        {1.0}
    },
    {
        "o2",
        AM_UNIT_COL_DENSITY,
        FRAC_O2,
        MASS_O2,
        0,
        ZADEP_AIRMASS,
        2,
        {K_TYPE_O2_COUPLED, K_TYPE_O2_UNCOUPLED},
        {1.0, 1.0}
    },
    {
        "16o2",
        AM_UNIT_COL_DENSITY,
        FRAC_16O2,
        MASS_16O2,
        0,
        ZADEP_AIRMASS,
        2,
        {K_TYPE_16O2_COUPLED, K_TYPE_16O2_UNCOUPLED},
        {1.0, 1.0}
    },
    {
        "16o_18o",
        AM_UNIT_COL_DENSITY,
        FRAC_16O_18O,
        MASS_16O_18O,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_16O_18O_UNCOUPLED},
        {1.0}
    },
    {
        "16o_17o",
        AM_UNIT_COL_DENSITY,
        FRAC_16O_17O,
        MASS_16O_17O,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_16O_17O_UNCOUPLED},
        {1.0}
    },
    {
        "o2_coupled",
        AM_UNIT_COL_DENSITY,
        FRAC_O2,
        MASS_O2,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_O2_COUPLED},
        {1.0}
    },
    {
        "16o2_coupled",
        AM_UNIT_COL_DENSITY,
        FRAC_16O2,
        MASS_16O2,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_16O2_COUPLED},
        {1.0}
    },
    {
        "o2_uncoupled",
        AM_UNIT_COL_DENSITY,
        FRAC_O2,
        MASS_O2,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_O2_UNCOUPLED},
        {1.0}
    },
    {
        "16o2_uncoupled",
        AM_UNIT_COL_DENSITY,
        FRAC_16O2,
        MASS_16O2,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_16O2_UNCOUPLED},
        {1.0}
    },
    {
        "16o_18o_uncoupled",
        AM_UNIT_COL_DENSITY,
        FRAC_16O_18O,
        MASS_16O_18O,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_16O_18O_UNCOUPLED},
        {1.0}
    },
    {
        "16o_17o_uncoupled",
        AM_UNIT_COL_DENSITY,
        FRAC_16O_17O,
        MASS_16O_17O,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_16O_17O_UNCOUPLED},
        {1.0}
    },
    {
        "o2o2",
        AM_UNIT_COL_DENSITY,
        FRAC_O2,
        MASS_O2,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_O2O2},
        {1.0}
    },
    {
        "o2air",
        AM_UNIT_COL_DENSITY,
        FRAC_O2,
        MASS_O2,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_O2O2},
        {1.0}
    },
    {
        "o3",
        UNIT_DU,
        0.0,
        MASS_O3,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_O3},
        {1.0}
    },
    {
        "16o3",
        UNIT_DU,
        0.0,
        MASS_16O3,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_16O3},
        {1.0}
    },
    {
        "16o2_18o",
        UNIT_DU,
        0.0,
        MASS_16O2_18O,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_16O2_18O},
        {1.0}
    },
    {
        "16o_18o_16o",
        UNIT_DU,
        0.0,
        MASS_16O_18O_16O,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_16O_18O_16O},
        {1.0}
    },
    {
        "16o2_17o",
        UNIT_DU,
        0.0,
        MASS_16O2_17O,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_16O2_17O},
        {1.0}
    },
    {
        "16o_17o_16o",
        UNIT_DU,
        0.0,
        MASS_16O_17O_16O,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_16O_17O_16O},
        {1.0}
    },
    {
        "ocs",              
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_OCS,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_OCS},
        {1.0}
    },
    {
        "16o_12c_32s",              
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_16O_12C_32S,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_16O_12C_32S},
        {1.0}
    },
    {
        "16o_12c_34s",              
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_16O_12C_34S,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_16O_12C_34S},
        {1.0}
    },
    {
        "16o_13c_32s",              
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_16O_13C_32S,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_16O_13C_32S},
        {1.0}
    },
    {
        "16o_12c_33s",              
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_16O_12C_33S,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_16O_12C_33S},
        {1.0}
    },
    {
        "18o_12c_32s",              
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_18O_12C_32S,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_18O_12C_32S},
        {1.0}
    },
    {
        "oh",
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_OH,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_OH},
        {1.0}
    },
    {
        "16o_h",
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_16O_H,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_16O_H},
        {1.0}
    },
    {
        "18o_h",
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_18O_H,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_18O_H},
        {1.0}
    },
    {
        "16o_d",
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_16O_D,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_16O_D},
        {1.0}
    },
    {
        "so2",
        UNIT_DU,
        0.0,
        MASS_SO2,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_SO2},
        {1.0}
    },
    {
        "32s_16o2",
        UNIT_DU,
        0.0,
        MASS_32S_16O2,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_32S_16O2},
        {1.0}
    },
    {
        "34s_16o2",
        UNIT_DU,
        0.0,
        MASS_34S_16O2,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_34S_16O2},
        {1.0}
    },
    {
        "oneline",
        AM_UNIT_COL_DENSITY,
        0.0,
        MASS_ONELINE,
        0,
        ZADEP_AIRMASS,
        1,
        {K_TYPE_ONELINE},
        {1.0}
    },
    {
        "END",
        UNIT_NONE,
        0.0,
        0.0,
        0,
        ZADEP_AIRMASS,
        0,
        {K_TYPE_NONE},
        {1.0}
    }
};


/***********************************************************
* int column_opacity(
*   model_t *model,
*   model_t *lmodel,
*   int lnum,
*   int cnum,
*   int *updateflag)
*
* Purpose:
*   Computes the opacity for column cnum in layer lnum.
*
*   If lmodel != NULL, then this is an update of a prior computation,
*   and lmodel holds the values of the model parameters at the time
*   of the prior computation.  The prior model parameters are used
*   to determine whether or not expensive prior computations need
*   to be redone.
*
*   The flag updateflag tells the calling function whether or not
*   this function needed to update the column opacity.  If the
*   opacity was updated (or if lmodel == NULL), updateflag is set
*   to 1.  Otherwise, updateflag is left in its initial state.
*
* Arguments:
*   model_t *model - pointer to model structure
*   model_t *lmodel - pointer to a model data structure containing
*       scalar data from a prior computation.
*   int lnum - layer number
*   int cnum - column number
*   char *updateflag - gets set to 1 if any quantity gets (re)computed,
*       otherwise retains initial state.
*
* Return:
*   0 if OK, 1 if an error occurred
************************************************************/

int column_opacity(
    model_t *model,
    model_t *lmodel,
    int lnum,
    int cnum,
    int *updateflag
    )
{
    int knum;
    int col_var_bitmask;
    int update_this_col;
    gridsize_t i;
    double n, Nn;
    layer_t *layer;
    layer_t *llayer;
    column_t *column;
    column_t *lcolumn;

    layer = model->layer[lnum];
    column = layer->column[cnum];
    if (lmodel != NULL) {
        llayer = lmodel->layer[lnum];
        lcolumn = llayer->column[cnum];
    } else {
        llayer = NULL;
        lcolumn = NULL;
    }
    /*
     * Check for negative variables and quit if they are encountered, since
     * this will cause numerical failures.  This could happen during a fit if
     * the user has disabled log mapping of fit variables, or could arise
     * from a coding error.
     */
    if (layer->P <= 0.0) {
        errlog(0, lnum);
        return 1;
    }
    if (layer->T <= 0.0) {
        errlog(1, lnum);
        return 1;
    }
    if (column->N_scaled < 0.0) {
        errlog(2, lnum);
        return 1;
    }
    /*
     * Mixing ratios cannot be negative, but they are allowed to exceed 1.0,
     * since this could be a transient condition during a fit.
     */
    if (column->vmr_scaled < 0.0) {
        errlog(3, lnum);
        return 1;
    }
    /*
     * Set bitmask indicating changed variables.  If any of these variables
     * have changed, the column opacity will need to be be recomputed.
     */
    update_this_col = 0;
    if (lmodel == NULL) {
        col_var_bitmask = DEP_ON_ALL;
    } else {
        col_var_bitmask = 0x0;
        if (fabs(layer->P - llayer->P) > (DBL_EPSILON * layer->P))
            col_var_bitmask |= DEP_ON_P;
        if (fabs(layer->T - llayer->T) > (DBL_EPSILON * layer->T))
            col_var_bitmask |= DEP_ON_T;
        if ((fabs(column->N_scaled - lcolumn->N_scaled)
            > (DBL_EPSILON * column->N_scaled)))
            col_var_bitmask |= DEP_ON_N;
        if (fabs(column->vmr_scaled - lcolumn->vmr_scaled)
            > (DBL_EPSILON * column->vmr_scaled))
            col_var_bitmask |= DEP_ON_VMR_SCALED;
    }
    if (col_var_bitmask)
        update_this_col = 1;
    /*
     * Additionally, any absorption coefficients that depend on changed
     * column variables need to be recomputed.
     */
    for (knum = 0; knum < column->n_abscoeffs; ++knum) {
        int dep_bitmask;
        int dep_flags = k_type[column->abscoeff[knum]->k_typenum].dep_flags;
        if (lmodel == NULL) {
            dep_bitmask = DEP_ON_ALL;
        } else {
            dep_bitmask = col_var_bitmask;
            if (fabs(column->abscoeff[knum]->vmr_selfbroad -
                lcolumn->abscoeff[knum]->vmr_selfbroad) >
                (DBL_EPSILON * column->abscoeff[knum]->vmr_selfbroad)) {
                dep_bitmask |= DEP_ON_VMR_SELFBROAD;
            }
        }
        if (dep_flags & dep_bitmask) {
            update_this_col = 1;
            if (get_absorption_coefficient(model, lmodel, lnum, cnum, knum))
                return 1;
        }
    }
    /*
     * The column opacity also needs to be recomputed if the model subgrid
     * ranges have changed, for example if the LO frequency has changed, or
     * if the output frequency range has changed.
     */
    if (compare_spectral_subgrid_ranges(model, lmodel))
        update_this_col = 1;
    /*
     * If a recomputation has been triggered by any of the criteria above,
     * set the update flag and proceed with the computation.  Otherwise,
     * leave the update flag in its original state (which may have been
     * set by the caller) and return.
     */
    if (update_this_col)
        *updateflag = 1;
    else
        return 0;   /* no changes on this column, return */
    /*
     * Sum up the unresolved line counts from all the line-by-line absorption
     * coefficient computations that go into this column's opacity
     * computation.
     */
    column->unres_lines = 0;
    for (knum = 0; knum < column->n_abscoeffs; ++knum)
        column->unres_lines += column->abscoeff[knum]->unres_lines;
    /*
     * Compute the column opacity.
     */
    switch (column->col_typenum) {
    /*
     * Parametric column types with no absorption coefficient:
     */
    case COL_TYPE_ATTEN:
    case COL_TYPE_ATTEN_AIRMASS:
    case COL_TYPE_ATTEN_SEC_ZA:
    case COL_TYPE_ATTEN_SIN_SQUARED_ZA:
    case COL_TYPE_ATTEN_SIN_SQUARED_2ZA:
        for (i = model->isub[0].min; i <= model->isub[0].max; ++i)
            column->ztau[i] = column->N_scaled;
        for (i = model->isub[1].min; i <= model->isub[1].max; ++i)
            column->ztau[i] = column->N_scaled;
        break;
    /*
     * Simple line-by-line column types depending on a single
     * molecular absorption coefficient.  Note that 16o_18o and
     * 16o_17o are in this case, but 16o2 is not.
     */
    case COL_TYPE_CH4:
    case COL_TYPE_12CH4:
    case COL_TYPE_13CH4:
    case COL_TYPE_12CH3D:
    case COL_TYPE_CH3CN:
    case COL_TYPE_12CH3_12C14N:
    case COL_TYPE_CH3OH:
    case COL_TYPE_12CH3_16OH:
    case COL_TYPE_CO:
    case COL_TYPE_12C_16O:
    case COL_TYPE_13C_16O:
    case COL_TYPE_12C_18O:
    case COL_TYPE_12C_17O:
    case COL_TYPE_13C_18O:
    case COL_TYPE_13C_17O:
    case COL_TYPE_CO2:
    case COL_TYPE_12C_16O2:
    case COL_TYPE_13C_16O2:
    case COL_TYPE_16O_12C_18O:
    case COL_TYPE_16O_12C_17O:
    case COL_TYPE_16O_13C_18O:
    case COL_TYPE_16O_13C_17O:
    case COL_TYPE_12C_18O2:
    case COL_TYPE_ClO:
    case COL_TYPE_35Cl_16O:
    case COL_TYPE_37Cl_16O:
    case COL_TYPE_HBr:
    case COL_TYPE_H_79Br:
    case COL_TYPE_H_81Br:
    case COL_TYPE_HCN:
    case COL_TYPE_H_12C_14N:
    case COL_TYPE_H_13C_14N:
    case COL_TYPE_H_12C_15N:
    case COL_TYPE_H2CO:
    case COL_TYPE_H2_12C_16O:
    case COL_TYPE_H2_13C_16O:
    case COL_TYPE_H2_12C_18O:
    case COL_TYPE_HCl:
    case COL_TYPE_H_35Cl:
    case COL_TYPE_H_37Cl:
    case COL_TYPE_HF:
    case COL_TYPE_H_19F:
    case COL_TYPE_HNO3:
    case COL_TYPE_H_14N_16O3:
    case COL_TYPE_H2O_LINES:
    case COL_TYPE_H2_16O:
    case COL_TYPE_H2_18O:
    case COL_TYPE_H2_17O:
    case COL_TYPE_HD_16O:
    case COL_TYPE_HD_18O:
    case COL_TYPE_HD_17O:
    case COL_TYPE_H2O2:
    case COL_TYPE_H2_16O2:
    case COL_TYPE_HO2:
    case COL_TYPE_H_16O2:
    case COL_TYPE_HOCl:
    case COL_TYPE_H_16O_35Cl:
    case COL_TYPE_H_16O_37Cl:
    case COL_TYPE_H2S:
    case COL_TYPE_H2_32S:
    case COL_TYPE_H2_34S:
    case COL_TYPE_H2_33S:
    case COL_TYPE_NH3:
    case COL_TYPE_14NH3:
    case COL_TYPE_15NH3:
    case COL_TYPE_N2O:
    case COL_TYPE_14N2_16O:
    case COL_TYPE_14N_15N_16O:
    case COL_TYPE_15N_14N_16O:
    case COL_TYPE_14N2_18O:
    case COL_TYPE_14N2_17O:
    case COL_TYPE_NO:
    case COL_TYPE_14N_16O:
    case COL_TYPE_NO2:
    case COL_TYPE_14N_16O2:
    case COL_TYPE_O:
    case COL_TYPE_16O:
    case COL_TYPE_16O_18O:
    case COL_TYPE_16O_17O:
    case COL_TYPE_O2_COUPLED:
    case COL_TYPE_16O2_COUPLED:
    case COL_TYPE_O2_UNCOUPLED:
    case COL_TYPE_16O2_UNCOUPLED:
    case COL_TYPE_16O_18O_UNCOUPLED:
    case COL_TYPE_16O_17O_UNCOUPLED:
    case COL_TYPE_O3:
    case COL_TYPE_16O3:
    case COL_TYPE_16O2_18O:
    case COL_TYPE_16O_18O_16O:
    case COL_TYPE_16O2_17O:
    case COL_TYPE_16O_17O_16O:
    case COL_TYPE_OCS:
    case COL_TYPE_16O_12C_32S:
    case COL_TYPE_16O_12C_34S:
    case COL_TYPE_16O_13C_32S:
    case COL_TYPE_16O_12C_33S:
    case COL_TYPE_18O_12C_32S:
    case COL_TYPE_OH:
    case COL_TYPE_16O_H:
    case COL_TYPE_18O_H:
    case COL_TYPE_16O_D:
    case COL_TYPE_SO2:
    case COL_TYPE_32S_16O2:
    case COL_TYPE_34S_16O2:
    case COL_TYPE_ONELINE:
        for (i = model->isub[0].min; i <= model->isub[0].max; ++i)
            column->ztau[i] = column->N_scaled * column->abscoeff[0]->k[i];
        for (i = model->isub[1].min; i <= model->isub[1].max; ++i)
            column->ztau[i] = column->N_scaled * column->abscoeff[0]->k[i];
        break;
    /*
     * Type h2o includes line-by-line, air-induced continuum, and
     * self-induced continuum.  Type h2_16o_lines_plus_continuum is
     * similar, but only includes lines for h2_16o.
     */
    case COL_TYPE_H2O:
    case COL_TYPE_H2_16O_LINES_PLUS_CONTINUUM:
        /* line-by-line */
        for (i = model->isub[0].min; i <= model->isub[0].max; ++i)
            column->ztau[i] = column->N_scaled * column->abscoeff[0]->k[i];
        for (i = model->isub[1].min; i <= model->isub[1].max; ++i)
            column->ztau[i] = column->N_scaled * column->abscoeff[0]->k[i];
        /* For continua, warn if h2o vmr has defaulted to 0.0 */
        if ((column->N_scaled > 0.0) && (column->vmr_stat & VMR_ZERO_DEFAULT)) {
            errlog(57, 0);
            column->vmr_stat |= VMR_ERROR;
        }
        /* air-induced continuum */
        n = N_STP * (layer->P / P_STP) * (T_STP / layer->T);
        Nn = n * column->N_scaled * (1.0 - column->vmr_scaled);
        for (i = model->isub[0].min; i <= model->isub[0].max; ++i)
            column->ztau[i] += Nn * column->abscoeff[1]->k[i];
        for (i = model->isub[1].min; i <= model->isub[1].max; ++i)
            column->ztau[i] += Nn * column->abscoeff[1]->k[i];
        /* self-induced continuum */
        Nn = n * column->N_scaled * column->vmr_scaled;
        if (Nn != 0.0) {
            for (i = model->isub[0].min; i <= model->isub[0].max; ++i)
                column->ztau[i] += Nn * column->abscoeff[2]->k[i];
            for (i = model->isub[1].min; i <= model->isub[1].max; ++i)
                column->ztau[i] += Nn * column->abscoeff[2]->k[i];
        }
        break;
    /*
     * Type h2o_continuum is air- and self-induced h2o continuum, without
     * line-by-line absorption.
     */
    case COL_TYPE_H2O_CONTINUUM:
        /* For continua, warn if h2o vmr has defaulted to 0.0 */
        if ((column->N_scaled > 0.0) && (column->vmr_stat & VMR_ZERO_DEFAULT)) {
            errlog(57, 0);
            column->vmr_stat |= VMR_ERROR;
        }
        /* air-induced continuum*/
        n = N_STP * (layer->P / P_STP) * (T_STP / layer->T);
        Nn = n * column->N_scaled * (1.0 - column->vmr_scaled);
        for (i = model->isub[0].min; i <= model->isub[0].max; ++i)
            column->ztau[i] = Nn * column->abscoeff[0]->k[i];
        for (i = model->isub[1].min; i <= model->isub[1].max; ++i)
            column->ztau[i] = Nn * column->abscoeff[0]->k[i];
        /* self-induced continuum */
        Nn = n * column->N_scaled * column->vmr_scaled;
        if (Nn != 0.0) {
            for (i = model->isub[0].min; i <= model->isub[0].max; ++i)
                column->ztau[i] += Nn * column->abscoeff[1]->k[i];
            for (i = model->isub[1].min; i <= model->isub[1].max; ++i)
                column->ztau[i] += Nn * column->abscoeff[1]->k[i];
        }
        break;
    /*
     * Air-induced h2o continuum only
     */
    case COL_TYPE_H2O_AIR_CONTINUUM:
        /* Warn if h2o vmr has defaulted to 0.0 */
        if ((column->N_scaled > 0.0) && (column->vmr_stat & VMR_ZERO_DEFAULT)) {
            errlog(57, 0);
            column->vmr_stat |= VMR_ERROR;
        }
        n = N_STP * (layer->P / P_STP) * (T_STP / layer->T);
        Nn = n * column->N_scaled * (1.0 - column->vmr_scaled);
        for (i = model->isub[0].min; i <= model->isub[0].max; ++i)
            column->ztau[i] = Nn * column->abscoeff[0]->k[i];
        for (i = model->isub[1].min; i <= model->isub[1].max; ++i)
            column->ztau[i] = Nn * column->abscoeff[0]->k[i];
        break;
    /*
     * Self-induced h2o continuum only
     */
    case COL_TYPE_H2O_SELF_CONTINUUM:
        /* Warn if h2o vmr has defaulted to 0.0 */
        if ((column->N_scaled > 0.0) && (column->vmr_stat & VMR_ZERO_DEFAULT)) {
            errlog(57, 0);
            column->vmr_stat |= VMR_ERROR;
        }
        n = N_STP * (layer->P / P_STP) * (T_STP / layer->T);
        Nn = n * column->N_scaled * column->vmr_scaled;
        if (Nn != 0.0) {
            for (i = model->isub[0].min; i <= model->isub[0].max; ++i)
                column->ztau[i] = Nn * column->abscoeff[0]->k[i];
            for (i = model->isub[1].min; i <= model->isub[1].max; ++i)
                column->ztau[i] = Nn * column->abscoeff[0]->k[i];
        }
        break;
    /*
     * Types lwp_abs_Rayleigh and lwp_abs_Rayleigh are water ice
     * particle and liquid water droplet absorption in the Rayleigh
     * limit.  The absorption coefficient is stored internally as a
     * molecular absorption coefficient in [cm^2], and the column
     * density as a molecular column density in [molecules cm^-2].
     * Consequently, the computation here looks just like that for
     * gaseous species.
     */
    case COL_TYPE_IWP_ABS_RAYLEIGH:
    case COL_TYPE_LWP_ABS_RAYLEIGH:
        for (i = model->isub[0].min; i <= model->isub[0].max; ++i)
            column->ztau[i] = column->N_scaled * column->abscoeff[0]->k[i];
        for (i = model->isub[1].min; i <= model->isub[1].max; ++i)
            column->ztau[i] = column->N_scaled * column->abscoeff[0]->k[i];
        break;
    case COL_TYPE_N2N2:
    case COL_TYPE_O2O2:
        /*
         * Use layer P, T, and vmr to compute the number density n [cm^-3].
         * The opacity is then tau = N [cm^-2] * n [cm^-3] * k [cm^5].
         */
        n = N_STP * (layer->P / P_STP) * (T_STP / layer->T);
        Nn = n * column->vmr_scaled * column->N_scaled;
        for (i = model->isub[0].min; i <= model->isub[0].max; ++i)
            column->ztau[i] = Nn * column->abscoeff[0]->k[i];
        for (i = model->isub[1].min; i <= model->isub[1].max; ++i)
            column->ztau[i] = Nn * column->abscoeff[0]->k[i];
        break;
    case COL_TYPE_N2AIR:
        /*
         * N2-air is combined N2-N2 and N2-O2 CIA, where N2 makes the
         * rotational transition.  The assumed O2 mixing ratio is the
         * o2 dry air default.
         */
        n = N_STP * (layer->P / P_STP) * (T_STP / layer->T);
        n *= column->vmr_scaled + S_N2O2_ON_S_N2N2 * FRAC_O2;
        Nn = n * column->N_scaled;
        for (i = model->isub[0].min; i <= model->isub[0].max; ++i)
            column->ztau[i] = Nn * column->abscoeff[0]->k[i];
        for (i = model->isub[1].min; i <= model->isub[1].max; ++i)
            column->ztau[i] = Nn * column->abscoeff[0]->k[i];
        break;
    case COL_TYPE_O2AIR:
        /*
         * O2-air is combined O2-O2 and O2-N2 CIA, where O2 makes the
         * rotational transition.  The assumed N2 mixing ratio is the
         * n2 dry air default.
         */
        n = N_STP * (layer->P / P_STP) * (T_STP / layer->T);
        n *= column->vmr_scaled + S_O2N2_ON_S_O2O2 * FRAC_N2;
        Nn = n * column->N_scaled;
        for (i = model->isub[0].min; i <= model->isub[0].max; ++i)
            column->ztau[i] = Nn * column->abscoeff[0]->k[i];
        for (i = model->isub[1].min; i <= model->isub[1].max; ++i)
            column->ztau[i] = Nn * column->abscoeff[0]->k[i];
        break;
    /*
     * Type o2 combines o2_coupled and o2_uncoupled line-by-line
     * absorption.  Type 16o2 combines 16o2_coupled and 16o2_uncoupled
     * line-by-line absorption.
     */
    case COL_TYPE_O2:
    case COL_TYPE_16O2:
        for (i = model->isub[0].min; i <= model->isub[0].max; ++i)
            column->ztau[i] = column->N_scaled *
                  (column->abscoeff[0]->k[i] + column->abscoeff[1]->k[i]);
        for (i = model->isub[1].min; i <= model->isub[1].max; ++i)
            column->ztau[i] = column->N_scaled *
                  (column->abscoeff[0]->k[i] + column->abscoeff[1]->k[i]);
        break;
    case COL_TYPE_DRY_AIR:
        /*
         * Dry air consisting of N2, O2, (Ar), CO2, and nominal
         * concentrations of CH4, N2O, and CO.  See molecules.h for
         * more information.  The absorption coefficients are indexed
         * in the order given in the table above, which is
         * {K_TYPE_CH4, K_TYPE_CO, K_TYPE_CO2, K_TYPE_N2O,
         * K_TYPE_O2_COUPLED, K_TYPE_O2_UNCOUPLED, K_TYPE_N2N2,
         * K_TYPE_O2O2}
         */
        {
            double Nn_n2n2, Nn_o2o2;
            n = N_STP * (layer->P / P_STP) * (T_STP / layer->T);
            Nn_n2n2 = FRAC_N2 * n
                * (FRAC_N2 * column->vmr_scaled + S_N2O2_ON_S_N2N2 * FRAC_O2);
            Nn_o2o2 = FRAC_O2 * n
                * (FRAC_O2 * column->vmr_scaled + S_O2N2_ON_S_O2O2 * FRAC_N2);
            for (i = model->isub[0].min; i <= model->isub[0].max; ++i)
                column->ztau[i] = column->N_scaled * (
                      FRAC_CH4 * column->abscoeff[0]->k[i]
                    +  FRAC_CO * column->abscoeff[1]->k[i]
                    + FRAC_CO2 * column->abscoeff[2]->k[i]
                    + FRAC_N2O * column->abscoeff[3]->k[i]
                    );
            for (i = model->isub[1].min; i <= model->isub[1].max; ++i)
                column->ztau[i] = column->N_scaled * (
                      FRAC_CH4 * column->abscoeff[0]->k[i]
                    +  FRAC_CO * column->abscoeff[1]->k[i]
                    + FRAC_CO2 * column->abscoeff[2]->k[i]
                    + FRAC_N2O * column->abscoeff[3]->k[i]
                    );
            for (i = model->isub[0].min; i <= model->isub[0].max; ++i)
                column->ztau[i] += column->N_scaled * FRAC_O2 * (
                      column->abscoeff[4]->k[i]
                    + column->abscoeff[5]->k[i]
                    );
            for (i = model->isub[1].min; i <= model->isub[1].max; ++i)
                column->ztau[i] += column->N_scaled * FRAC_O2 * (
                      column->abscoeff[4]->k[i]
                    + column->abscoeff[5]->k[i]
                    );
            for (i = model->isub[0].min; i <= model->isub[0].max; ++i)
                column->ztau[i] += column->N_scaled * (
                      Nn_n2n2 * column->abscoeff[6]->k[i]
                    + Nn_o2o2 * column->abscoeff[7]->k[i]
                    );
            for (i = model->isub[1].min; i <= model->isub[1].max; ++i)
                column->ztau[i] += column->N_scaled * (
                      Nn_n2n2 * column->abscoeff[6]->k[i]
                    + Nn_o2o2 * column->abscoeff[7]->k[i]
                    );
        }
        break;
    case COL_TYPE_DRY_AIR_STD:
        /*
         * "Standard" dry air, consisting of N2, O2, (Ar), and CO2
         * only.  See molecules.h for more information.  The
         * absorption coefficients are indexed in the order given in
         * the table above, which is {K_TYPE_CO2, K_TYPE_O2_COUPLED,
         * K_TYPE_O2_UNCOUPLED, K_TYPE_N2N2, K_TYPE_O2O2}
         */
        {
            double Nn_n2n2, Nn_o2o2;
            n = N_STP * (layer->P / P_STP) * (T_STP / layer->T);
            Nn_n2n2 = FRAC_N2 * n
                * (FRAC_N2 * column->vmr_scaled + S_N2O2_ON_S_N2N2 * FRAC_O2);
            Nn_o2o2 = FRAC_O2 * n
                * (FRAC_O2 * column->vmr_scaled + S_O2N2_ON_S_O2O2 * FRAC_N2);
            for (i = model->isub[0].min; i <= model->isub[0].max; ++i)
                column->ztau[i] = column->N_scaled * (
                    FRAC_CO2 * column->abscoeff[0]->k[i]
                    );
            for (i = model->isub[1].min; i <= model->isub[1].max; ++i)
                column->ztau[i] = column->N_scaled * (
                    FRAC_CO2 * column->abscoeff[0]->k[i]
                    );
            for (i = model->isub[0].min; i <= model->isub[0].max; ++i)
                column->ztau[i] += column->N_scaled * FRAC_O2 * (
                      column->abscoeff[1]->k[i]
                    + column->abscoeff[2]->k[i]
                    );
            for (i = model->isub[1].min; i <= model->isub[1].max; ++i)
                column->ztau[i] += column->N_scaled * FRAC_O2 * (
                      column->abscoeff[1]->k[i]
                    + column->abscoeff[2]->k[i]
                    );
            for (i = model->isub[0].min; i <= model->isub[0].max; ++i)
                column->ztau[i] += column->N_scaled * (
                      Nn_n2n2 * column->abscoeff[3]->k[i]
                    + Nn_o2o2 * column->abscoeff[4]->k[i]
                    );
            for (i = model->isub[1].min; i <= model->isub[1].max; ++i)
                column->ztau[i] += column->N_scaled * (
                      Nn_n2n2 * column->abscoeff[3]->k[i]
                    + Nn_o2o2 * column->abscoeff[4]->k[i]
                    );
        }
        break;
    case COL_TYPE_DRY_AIR_OPTICAL_REFRACTIVITY:
    case COL_TYPE_H2O_OPTICAL_REFRACTIVITY:
        break;
    default:
        errlog(49, column->col_typenum);
        return 1;
    }
    return 0;
}   /* column_opacity() */


/***********************************************************
* double column_optical_delay(
*   model_t *model,
*   int lnum,
*   int cnum)
*
* Purpose:
*   Computes and returns the optical delay for column cnum
*   in layer lnum.
*
* Arguments:
*   model_t *model - pointer to model structure
*   int lnum - layer number
*   int cnum - column number
*
* Return:
*   optical delay, as an excess path in cm (as a double)
************************************************************/

double column_optical_delay(
    model_t *model,
    int lnum,
    int cnum)
{
    layer_t *layer;
    column_t *column;

    layer = model->layer[lnum];
    column = layer->column[cnum];
    if (col_type[column->col_typenum].flags & COL_OPTICAL_DELAY) {
        switch (column->col_typenum) {
        case COL_TYPE_DRY_AIR:
        case COL_TYPE_DRY_AIR_OPTICAL_REFRACTIVITY:
            return column->N_scaled * layer->m *
                DRY_OPTICAL_VOLUME_REFRACTIVITY;
        case COL_TYPE_H2O:
        case COL_TYPE_H2_16O_LINES_PLUS_CONTINUUM:
        case COL_TYPE_H2O_OPTICAL_REFRACTIVITY:
            return column->N_scaled * layer->m *
                H2O_OPTICAL_VOLUME_REFRACTIVITY;
        default:
            errlog(103, column->col_typenum);
            return 1;
        }
    }
    return 0.0;
}   /* column_optical_delay() */
