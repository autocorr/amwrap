/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* units.h                           S. Paine rev. 2024 May 9
*
* Constants and declarations for units.c
************************************************************/

#ifndef AM_UNITS_H
#define AM_UNITS_H

#include <stdio.h>

/*
 * Units and unit conversion data are contained in an array of the
 * following structures, defined in units.c
 *
 * The group field controls which unit conversions are allowed.  This
 * does not necessarily mean the dimensions are strictly the same.
 * Conventional conversions via the speed of light are allowed between
 * frequency and wavenumber, and between distance and time delay.
 * Customary units for H2O and O3 column density convert to their
 * molecular column density equivalents.  The purpose of unit groups
 * is to catch definite wrong unit errors instead of silently
 * performing an inappropriate conversion.
 */
extern struct unit_tabentry {
    const char    *name; /* text name of unit */
    const int     group; /* unit conversions are only allowed within group */
    const double offset; /* am_native_unit = factor * (user_unit + offset) */
    const double factor;
} unit_tab[];

/*
 * unit_tab[] is indexed using the following enums.
 */
enum {
    /*
     * dimensionless quantity
     */
    UNIT_NONE,
    /*
     * unit assigned at runtime
     */
    UNIT_AUTO,
    /*
     * temperature
     */
    UNIT_KELVIN,
    UNIT_CELSIUS,
    /*
     * column density
     */
    UNIT_ONE_ON_CMSQUARED,
    UNIT_ONE_ON_CMSQUARED_ALT,
    UNIT_ONE_ON_MSQUARED,
    UNIT_ONE_ON_MSQUARED_ALT,
    UNIT_CM_PWV,
    UNIT_MM_PWV,
    UNIT_G_ON_CMSQUARED,
    UNIT_G_ON_CMSQUARED_ALT,
    UNIT_G_ON_MSQUARED,
    UNIT_G_ON_MSQUARED_ALT,
    UNIT_KG_ON_MSQUARED,
    UNIT_KG_ON_MSQUARED_ALT,
    UNIT_UM_PWV,
    UNIT_DOBSON,
    UNIT_DU,
    /*
     * pressure
     */
    UNIT_BAR,
    UNIT_MBAR,
    UNIT_HPA,
    UNIT_PA,
    UNIT_ATM,
    UNIT_TORR,
    /*
     * frequency
     */
    UNIT_GHZ,
    UNIT_HZ,
    UNIT_KHZ,
    UNIT_MHZ,
    UNIT_THZ,
    UNIT_WAVENUM,
    UNIT_WAVENUM_ALT,
    /*
     * angle
     */
    UNIT_RAD,
    UNIT_RADIAN,
    UNIT_DEG,
    UNIT_DEGREE,
    UNIT_ARCMIN,
    UNIT_ARCSEC,
    /*
     * spectral radiance
     */
    UNIT_WATT_ON_CMSQUARED_GHZ_SR,
    UNIT_WATT_ON_CMSQUARED_GHZ_SR_ALT,
    UNIT_MILLIWATT_ON_MSQUARED_WAVENUM_SR,
    UNIT_MILLIWATT_ON_MSQUARED_WAVENUM_SR_ALT,
    UNIT_RU,
    UNIT_WATT_ON_MSQUARED_HERTZ_SR,
    UNIT_WATT_ON_MSQUARED_HERTZ_SR_ALT,
    UNIT_JY_ON_ARCSECSQUARED,
    UNIT_JY_ON_ARCSECSQUARED_ALT,
    UNIT_JY_ON_STERADIAN,
    UNIT_JY_ON_STERADIAN_ALT,
    /*
     * geometric delay and distance
     */
    UNIT_CM,
    UNIT_UM,
    UNIT_MM,
    UNIT_M,
    UNIT_KM,
    /*
     * time delay equivalent to geometric delay
     */
    UNIT_PS,
    UNIT_FS,
    UNIT_S,
    /*
     * gravitational acceleration
     */
    UNIT_CM_ON_SSQUARED,
    UNIT_CM_ON_SSQUARED_ALT,
    UNIT_GAL,
    UNIT_M_ON_SSQUARED,
    UNIT_M_ON_SSQUARED_ALT,
    /*
     * gravity gradient
     */
    UNIT_ONE_ON_SSQUARED,
    UNIT_ONE_ON_SSQUARED_ALT,
    UNIT_GAL_ON_CM,
    UNIT_GAL_ON_CM_ALT,
    UNIT_MILLIGAL_ON_M,
    UNIT_MILLIGAL_ON_M_ALT,
    /*
     * optical depth (opacity)
     */
    UNIT_NEPER,
    UNIT_DB,
    /*
     * molecular absorption coefficient
     */
    UNIT_CM2,
    /*
     * binary absorption coefficient
     */
    UNIT_CM5,
    UNIT_END_OF_TABLE
};

enum { /* names for groups of interconvertable units */
    UGROUP_NONE,
    UGROUP_AUTO,
    UGROUP_TEMPERATURE,
    UGROUP_COL_DENSITY,
    UGROUP_PRESSURE,
    UGROUP_FREQUENCY,
    UGROUP_ANGLE,
    UGROUP_RADIANCE,
    UGROUP_DELAY_DIST,
    UGROUP_ACCEL,
    UGROUP_ACCEL_GRADIENT,
    UGROUP_OPACITY,
    UGROUP_MOL_ABSCOEFF,
    UGROUP_BIN_ABSCOEFF
};

enum { /* names for am native units */
    AM_UNIT_NONE            = UNIT_NONE,
    AM_UNIT_TEMPERATURE     = UNIT_KELVIN,
    AM_UNIT_COL_DENSITY     = UNIT_ONE_ON_CMSQUARED,
    AM_UNIT_PRESSURE        = UNIT_MBAR,
    AM_UNIT_FREQUENCY       = UNIT_GHZ,
    AM_UNIT_ANGLE           = UNIT_RAD,
    AM_UNIT_RADIANCE        = UNIT_WATT_ON_CMSQUARED_GHZ_SR,
    AM_UNIT_DELAY           = UNIT_CM,
    AM_UNIT_ACCEL           = UNIT_CM_ON_SSQUARED,
    AM_UNIT_ACCEL_GRADIENT  = UNIT_ONE_ON_SSQUARED,
    AM_UNIT_OPACITY         = UNIT_NEPER,
    AM_UNIT_DISTANCE        = UNIT_CM,
    AM_UNIT_MOL_ABSCOEFF    = UNIT_CM2,
    AM_UNIT_BIN_ABSCOEFF    = UNIT_CM5
};

/*
 * Functions in units.c
 */

int convert_to_am_unit(double*, const int, const int);
int is_unit(const char*);
int get_unitnum(const char*);
int print_differential_with_unit(FILE*, const char*, double, int);
int print_with_unit(FILE*, const char*, double, int);
int snprint_with_unit(char*, size_t, const char*, double, int);

#endif /* AM_UNITS_H */
