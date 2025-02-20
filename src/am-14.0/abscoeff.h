/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* abscoeff.h                     S. Paine rev. 2024 April 20
*
* Constants and declarations for abscoeff.c
************************************************************/

#ifndef AM_ABSCOEFF_H
#define AM_ABSCOEFF_H

#include "am_types.h"

/*
 * Absorption coefficient type properties are defined in the
 * following table, defined in abscoeff.c.
 */
extern struct k_type_tabentry {
    const char *name;  /* text name of absorption coefficient type      */
    const int unitnum; /* native unit number for absorption coefficient */
    const int iso;     /* isotopologue for line-by-line types, 0 = all  */
    const int default_lineshape; /* for line-by-line abscoeff_types     */
    const int comp_flags; /* flags relating to computation              */
    const int dep_flags;  /* dependency flags                           */
} k_type[];

/*
 * Additions to this list between major program version number changes
 * must be made at the end, since these values are part of the cache
 * file header data.
 */
enum {
    K_TYPE_NONE,
    K_TYPE_CH4,
        K_TYPE_12CH4,
        K_TYPE_13CH4,
        K_TYPE_12CH3D,
    K_TYPE_CH3CN,
        K_TYPE_12CH3_12C14N,
    K_TYPE_CH3OH,
        K_TYPE_12CH3_16OH,
    K_TYPE_CO,
        K_TYPE_12C_16O,
        K_TYPE_13C_16O,
        K_TYPE_12C_18O,
        K_TYPE_12C_17O,
        K_TYPE_13C_18O,
        K_TYPE_13C_17O,
    K_TYPE_CO2,
        K_TYPE_12C_16O2,
        K_TYPE_13C_16O2,
        K_TYPE_16O_12C_18O,
        K_TYPE_16O_12C_17O,
        K_TYPE_16O_13C_18O,
        K_TYPE_16O_13C_17O,
        K_TYPE_12C_18O2,
    K_TYPE_ClO,
        K_TYPE_35Cl_16O,
        K_TYPE_37Cl_16O,
    K_TYPE_HBr,
        K_TYPE_H_79Br,
        K_TYPE_H_81Br,
    K_TYPE_HCN,
        K_TYPE_H_12C_14N,
        K_TYPE_H_13C_14N,
        K_TYPE_H_12C_15N,
    K_TYPE_H2CO,
        K_TYPE_H2_12C_16O,
        K_TYPE_H2_13C_16O,
        K_TYPE_H2_12C_18O,
    K_TYPE_HCl,
        K_TYPE_H_35Cl,
        K_TYPE_H_37Cl,
    K_TYPE_HF,
        K_TYPE_H_19F,
    K_TYPE_HNO3,
        K_TYPE_H_14N_16O3,
    K_TYPE_H2O_LINES,
        K_TYPE_H2_16O,
        K_TYPE_H2_18O,
        K_TYPE_H2_17O,
        K_TYPE_HD_16O,
        K_TYPE_HD_18O,
        K_TYPE_HD_17O,
    K_TYPE_H2O_AIR_CONTINUUM,
    K_TYPE_H2O_SELF_CONTINUUM,
    K_TYPE_H2O2,
        K_TYPE_H2_16O2,
    K_TYPE_HO2,
        K_TYPE_H_16O2,
    K_TYPE_HOCl,
        K_TYPE_H_16O_35Cl,
        K_TYPE_H_16O_37Cl,
    K_TYPE_H2S,
        K_TYPE_H2_32S,
        K_TYPE_H2_34S,
        K_TYPE_H2_33S,
    K_TYPE_IWP_ABS_RAYLEIGH,
    K_TYPE_LWP_ABS_RAYLEIGH,
    K_TYPE_NH3,
        K_TYPE_14NH3,
        K_TYPE_15NH3,
    K_TYPE_N2N2,
    K_TYPE_N2O,
        K_TYPE_14N2_16O,
        K_TYPE_14N_15N_16O,
        K_TYPE_15N_14N_16O,
        K_TYPE_14N2_18O,
        K_TYPE_14N2_17O,
    K_TYPE_NO,
        K_TYPE_14N_16O,
    K_TYPE_NO2,
        K_TYPE_14N_16O2,
    K_TYPE_O,
        K_TYPE_16O,
    K_TYPE_O2_COUPLED,
        K_TYPE_16O2_COUPLED,
    K_TYPE_O2_UNCOUPLED,
        K_TYPE_16O2_UNCOUPLED,
        K_TYPE_16O_18O_UNCOUPLED,
        K_TYPE_16O_17O_UNCOUPLED,
    K_TYPE_O2O2,
    K_TYPE_O3,
        K_TYPE_16O3,
        K_TYPE_16O2_18O,
        K_TYPE_16O_18O_16O,
        K_TYPE_16O2_17O,
        K_TYPE_16O_17O_16O,
    K_TYPE_OCS,
        K_TYPE_16O_12C_32S,
        K_TYPE_16O_12C_34S,
        K_TYPE_16O_13C_32S,
        K_TYPE_16O_12C_33S,
        K_TYPE_18O_12C_32S,
    K_TYPE_OH,
        K_TYPE_16O_H,
        K_TYPE_18O_H,
        K_TYPE_16O_D,
    K_TYPE_SO2,
        K_TYPE_32S_16O2,
        K_TYPE_34S_16O2,
    K_TYPE_ONELINE,
    K_TYPE_END_OF_TABLE
};

/*
 * Computation flags for absorption coefficients.
 */
enum {
    CACHEABLE = 0x1
};

/*
 * Dependency flag bits for absorption coefficients and opacities.
 */
enum {
    DEP_ON_P = 0x1,
    DEP_ON_T = 0x2,
    DEP_ON_VMR_SELFBROAD = 0x4,
    DEP_ON_VMR_SCALED = 0x8,
    DEP_ON_N = 0x10,
    DEP_ON_LSHAPE = 0x20,
    DEP_ON_TOL = 0x40,
    DEP_ON_ALL = ~0x0
};

enum {
    LINESUM_DEP_FLAGS = (DEP_ON_P | DEP_ON_T | DEP_ON_VMR_SELFBROAD |
    DEP_ON_LSHAPE | DEP_ON_TOL)
};


int get_absorption_coefficient(model_t*, model_t*, int, int, int);

#endif /* AM_ABSCOEFF_H */
