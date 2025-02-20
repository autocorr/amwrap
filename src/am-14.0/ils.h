/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* ils.h                          S. Paine rev. 2020 August 5
*
* Constants and declarations for ils.c
************************************************************/

#ifndef AM_ILS_H
#define AM_ILS_H

#include "am_types.h"

/*
 * Instrumental line shape properties are defined in the following
 * table, defined in ils.c.
 */
extern struct ils_type_tabentry {
    const char *name; /* text name of ils type */
    double     fwhm;  /* fwhm, in normal units */
} ils_type[];

enum {
    ILS_NONE,
    ILS_SINC,
    ILS_SINC_SQUARED,
    ILS_BARTLETT,
    ILS_HANN,
    ILS_HAMMING,
    ILS_BLACKMAN,
    ILS_NORTON_BEER_I1,
    ILS_NORTON_BEER_I2,
    ILS_NORTON_BEER_I3,
    ILS_GAUSSIAN,
    ILS_RECTANGLE,
    ILS_END_OF_TABLE
};

/*
 * ILS mode control bits.
 */
enum {
    ILSMODE_NORMAL = 0x1,
    ILSMODE_DSB    = 0x2,
    ILSMODE_USB    = 0x4,
    ILSMODE_LSB    = 0x8
};

int apply_instrumental_line_shape(model_t*, model_t*);
int initialize_ils(model_t*);

#endif /* AM_ILS_H */
