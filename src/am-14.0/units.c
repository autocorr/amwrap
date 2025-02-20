/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* units.c                           S. Paine rev. 2024 May 9
*
* Units and unit conversions.
************************************************************/

#include <string.h>

#include "errlog.h"
#include "units.h"

/*
 * Quantities are converted from user units to native units as
 * follows:
 *
 *  x [native_units] = factor * (x' [user_units] + offset)
 *
 * where factor and offset for each unit are found in the table below.
 *
 * The integer tag after each unit name defines a unit group within
 * which unit conversion is allowed.
 *
 * Changes to this table require corresponding changes to the enum
 * table for units in units.h.
 */

struct unit_tabentry unit_tab[] = {
    /*
     * dimensionless quantity
     */
    {"none",                 UGROUP_NONE,           0.0, 1.0},
    /*
     * unit determined at runtime
     */
    {"auto",                 UGROUP_AUTO,           0.0, 1.0},
    /*
     * temperature
     */
    {"K",                    UGROUP_TEMPERATURE,    0.0, 1.0},
    {"C",                    UGROUP_TEMPERATURE, 273.15, 1.0},
    /*
     * column density
     */
    {"cm^-2",                UGROUP_COL_DENSITY,    0.0, 1.0},
    {"cm-2",                 UGROUP_COL_DENSITY,    0.0, 1.0},
    {"m^-2",                 UGROUP_COL_DENSITY,    0.0, 1.0e-4},
    {"m-2",                  UGROUP_COL_DENSITY,    0.0, 1.0e-4},
    {"cm_pwv",               UGROUP_COL_DENSITY,    0.0, 3.3427e22},
    {"mm_pwv",               UGROUP_COL_DENSITY,    0.0, 3.3427e21},
    {"g*cm^-2",              UGROUP_COL_DENSITY,    0.0, 3.3427e22},
    {"g*cm-2",               UGROUP_COL_DENSITY,    0.0, 3.3427e22},
    {"g*m^-2",               UGROUP_COL_DENSITY,    0.0, 3.3427e18},
    {"g*m-2",                UGROUP_COL_DENSITY,    0.0, 3.3427e18},
    {"kg*m^-2",              UGROUP_COL_DENSITY,    0.0, 3.3427e21},
    {"kg*m-2",               UGROUP_COL_DENSITY,    0.0, 3.3427e21},
    {"um_pwv",               UGROUP_COL_DENSITY,    0.0, 3.3427e18},
    {"dobson",               UGROUP_COL_DENSITY,    0.0, 2.6868e16},
    {"DU",                   UGROUP_COL_DENSITY,    0.0, 2.6868e16},
    /*
     * pressure
     */
    {"bar",                     UGROUP_PRESSURE,    0.0, 1.0e3},
    {"mbar",                    UGROUP_PRESSURE,    0.0, 1.0},
    {"hPa",                     UGROUP_PRESSURE,    0.0, 1.0},
    {"Pa",                      UGROUP_PRESSURE,    0.0, 1.0e-2},
    {"atm",                     UGROUP_PRESSURE,    0.0, 1013.25},
    {"Torr",                    UGROUP_PRESSURE,    0.0, 1.333223684210526},
    /*
     * frequency
     */
    {"GHz",                    UGROUP_FREQUENCY,    0.0, 1.0},
    {"Hz",                     UGROUP_FREQUENCY,    0.0, 1.0e-9},
    {"kHz",                    UGROUP_FREQUENCY,    0.0, 1.0e-6},
    {"MHz",                    UGROUP_FREQUENCY,    0.0, 1.0e-3},
    {"THz",                    UGROUP_FREQUENCY,    0.0, 1.0e3},
    {"cm^-1",                  UGROUP_FREQUENCY,    0.0, 29.9792458},
    {"cm-1",                   UGROUP_FREQUENCY,    0.0, 29.9792458},
    /*
     * angle
     */
    {"rad",                        UGROUP_ANGLE,    0.0, 1.0},
    {"radian",                     UGROUP_ANGLE,    0.0, 1.0},
    {"deg",                        UGROUP_ANGLE,    0.0, 1.7453292519943295e-2},
    {"degree",                     UGROUP_ANGLE,    0.0, 1.7453292519943295e-2},
    {"arcmin",                     UGROUP_ANGLE,    0.0, 2.908882086657216e-4},
    {"arcsec",                     UGROUP_ANGLE,    0.0, 4.848136811095360e-6},
    /*
     * spectral radiance
     */
    {"watt*cm^-2*GHz^-1*sr^-1", UGROUP_RADIANCE,    0.0, 1.0},
    {"watt*cm-2*GHz-1*sr-1",    UGROUP_RADIANCE,    0.0, 1.0},
    {"mW*m^-2*(cm^-1)^-1*sr^-1",
                                UGROUP_RADIANCE,    0.0, 3.33564095198152e-9},
    {"mW*m-2*(cm-1)-1*sr-1",    UGROUP_RADIANCE,    0.0, 3.33564095198152e-9},
    {"RU",                      UGROUP_RADIANCE,    0.0, 3.33564095198152e-9},
    {"watt*m^-2*Hz^-1*sr^-1",   UGROUP_RADIANCE,    0.0, 1.0e5},
    {"watt*m-2*Hz-1*sr-1",      UGROUP_RADIANCE,    0.0, 1.0e5},
    {"Jy*arcsec^-2",            UGROUP_RADIANCE,    0.0, 4.254517029e-11},
    {"Jy*arcsec-2",             UGROUP_RADIANCE,    0.0, 4.254517029e-11},
    {"Jy*sr^-1",                UGROUP_RADIANCE,    0.0, 1.0e-21},
    {"Jy*sr-1",                 UGROUP_RADIANCE,    0.0, 1.0e-21},
    /*
     * geometric delay and distance
     */
    {"cm",                    UGROUP_DELAY_DIST,    0.0, 1.0},
    {"um",                    UGROUP_DELAY_DIST,    0.0, 1.0e-4},
    {"mm",                    UGROUP_DELAY_DIST,    0.0, 0.1},
    {"m",                     UGROUP_DELAY_DIST,    0.0, 1.0e2},
    {"km",                    UGROUP_DELAY_DIST,    0.0, 1.0e5},
    /*
     * time delay equivalent to geometric delay (in same unit group as
     * geometric delay and distance, above).
     */
    {"ps",                    UGROUP_DELAY_DIST,    0.0, 2.99792458e-2},
    {"fs",                    UGROUP_DELAY_DIST,    0.0, 2.99792458e-5},
    {"s",                     UGROUP_DELAY_DIST,    0.0, 2.99792458e10},
    /*
     * gravitational acceleration
     */
    {"cm*s^-2",                    UGROUP_ACCEL,    0.0, 1.0},
    {"cm*s-2",                     UGROUP_ACCEL,    0.0, 1.0},
    {"Gal",                        UGROUP_ACCEL,    0.0, 1.0},
    {"m*s^-2",                     UGROUP_ACCEL,    0.0, 1.0e2},
    {"m*s-2",                      UGROUP_ACCEL,    0.0, 1.0e2},
    /*
     * gravity gradient
     */
    {"s^-2",              UGROUP_ACCEL_GRADIENT,    0.0, 1.0},
    {"s-2",               UGROUP_ACCEL_GRADIENT,    0.0, 1.0},
    {"Gal*cm^-1",         UGROUP_ACCEL_GRADIENT,    0.0, 1.0},
    {"Gal*cm-1",          UGROUP_ACCEL_GRADIENT,    0.0, 1.0},
    {"mGal*m^-1",         UGROUP_ACCEL_GRADIENT,    0.0, 1.0e-5},
    {"mGal*m-1",          UGROUP_ACCEL_GRADIENT,    0.0, 1.0e-5},
    /*
     * optical depth (opacity)
     */
    {"neper",                    UGROUP_OPACITY,    0.0, 1.0},
    {"dB",                       UGROUP_OPACITY,    0.0, 2.30258509299405e-1},
    /*
     * molecular absorption coefficient
     */
    {"cm^2",                UGROUP_MOL_ABSCOEFF,    0.0, 1.0},
    /*
     * binary absorption coefficient
     */
    {"cm^5",                UGROUP_BIN_ABSCOEFF,    0.0, 1.0},
    /*
     * end of table
     */
    {"END",                         UGROUP_NONE,    0.0, 1.0},
};


/***********************************************************
* int convert_to_am_unit(
*         double *q, const int unitnum, const int * diff_flag)
*
* Purpose:
*   Converts a from user units to am native units, modifying
*   the stored value of the quantity q.  If diff_flag is
*   nonzero, this is a differential conversion, and any
*   zero offset, as for degrees C, is not applied.
*
* Arguments:
*   double    q      a  - pointer to converted variable
*   const int unitnuma  - user unit number converted from
*   const int diff_flag - flag for differential conversion
*
* Return:
*   0 if OK
*   1 on error
************************************************************/

int convert_to_am_unit(
        double *q, const int unitnum, const int diff_flag)
{
    if (unitnum < 0 || unitnum >= UNIT_END_OF_TABLE) {
        errlog(12, 0);
        return 1;
    }
    if (!diff_flag)
        *q += unit_tab[unitnum].offset;
    *q *= unit_tab[unitnum].factor;
    return 0;
}   /* convert_to_am_unit() */


/***********************************************************
* int is_unit(const char *s)
*
* Purpose:
*   Checks whether the string *s is a unit name.
*
* Arguments:
*   const char *s
*
* Return:
*   1 if *s is the name of a unit
*   0 otherwise
************************************************************/

int is_unit(const char *s)
{
    return (get_unitnum(s) < UNIT_END_OF_TABLE);
}   /* is_unit() */


/***********************************************************
* int get_unitnum(const char *unitname)
*
* Purpose:
*   Finds a unit by name in the unit conversion table, and
*   returns the index number.  If the unit is not found,
*   UNIT_END_OF_TABLE is returned.
*
* Arguments:
*   const char *unitname - name of unit
*
* Return:
*   index number of unit, or UNIT_END_OF_TABLE if unit name
*   not found
************************************************************/

int get_unitnum(const char *unitname)
{
    int i;
    for (i = 0; i < UNIT_END_OF_TABLE; ++i) {
        if (strcmp(unit_tab[i].name, unitname) == 0) {
            return i;
        }
    }
    return i;
}   /* get_unitnum() */


/***********************************************************
* int print_differential_with_unit(
*   FILE *stream, const char *fmt, double q, int unitnum)
*
* Purpose:
*   Same as print_with_unit(), above, but does not include any
*   offset in the unit conversion.  In particular, degrees C
*   should not be converted with an offset when a differential
*   value is being printed.
*
* Arguments:
*   FILE *stream - stream to be printed to
*   const char *fmt - format string, including %s for the unit
*   double q - quantity to be printed
*   int unitnum - table index number of unit
*
* Return:
*   number of characters printed, or negative on error.
************************************************************/

int print_differential_with_unit(
    FILE *stream, const char *fmt, double q, int unitnum)
{
    return fprintf(
            stream,
            fmt,
            q / unit_tab[unitnum].factor,
            unit_tab[unitnum].name);
}   /* print_differential_with_unit() */


/***********************************************************
* int print_with_unit(FILE *stream, const char *fmt, double q, int unitnum)
*
* Purpose:
*   Prints a dimensioned quantity q to a specified stream, using a
*   specified format and specified unit.  The quantity q, assumed
*   to be stored in am native units, is converted to the specified
*   unit for printing, but the stored value is not modified.
*
* Arguments:
*   FILE *stream - stream to be printed to
*   const char *fmt - format string, including %s for the unit
*   double q - quantity to be printed
*   int unitnum - table index number of unit
*
* Return:
*   number of characters printed, or negative on error.
************************************************************/

int print_with_unit(FILE *stream, const char *fmt, double q, int unitnum)
{
    return fprintf(
            stream,
            fmt,
            (q / unit_tab[unitnum].factor) - unit_tab[unitnum].offset,
            unit_tab[unitnum].name);
}   /* print_with_unit() */


/***********************************************************
* int snprint_with_unit(
*         char *s,
*         size_t n,
*         const char *fmt,
*         double q, int unitnum)
*
* Purpose:
*   Prints a dimensioned quantity q to a string s, using a
*   specified format and unit.  The quantity q, assumed to
*   be stored in am native units, is converted to the
*   specified unit for printing, but the stored value is not
*   modified.  At most n-1 characters will be printed to s.
*
* Arguments:
*   char *s         - string to be printed to
*   size_t n        - 1 + maximum number of characters printed
*   const char *fmt - format string, including %s for the unit
*   double q        - quantity to be printed
*   int unitnum     - table index number of unit
*
* Return:
*   If fewer than n characters were printed, the return is
*   the number of characters printed to s.  If the return
*   value is n or greater, it indicates the number of
*   characters that were discarded.
************************************************************/

int snprint_with_unit(
        char *s,
        size_t n,
        const char *fmt,
        double q,
        int unitnum)
{
    return snprintf(
            s,
            n,
            fmt,
            (q / unit_tab[unitnum].factor) - unit_tab[unitnum].offset,
            unit_tab[unitnum].name);
}   /* snprint_with_unit() */
