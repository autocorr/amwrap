/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* linesum.h                       S. Paine rev. 2020 July 20
*
* Constants and declarations for linesum.c
************************************************************/

#ifndef AM_LINESUM_H
#define AM_LINESUM_H

#include "am_types.h"

/*
 * MAX_ISO is the maximum number of isotopologues for any molecule.
 */
#define MAX_ISO 16

/*
 * Line shapes are listed in the following table, defined in linesum.c,
 * and indexed using the enums below.
 */
extern struct lineshape_type_tabentry {
    const char *name; /* text name of lineshape type */
    double f_cutoff;  /* cutoff frequency [GHz] for truncated lineshapes */
} lineshape_type[];

enum {
    LINESHAPE_NONE,
    LINESHAPE_DOPPLER,
    LINESHAPE_FULL_LORENTZ,
    LINESHAPE_GROSS,
    LINESHAPE_LORENTZ,
    LINESHAPE_VVW,
    LINESHAPE_VVW_COUPLED,
    LINESHAPE_VVH,
    LINESHAPE_VVH_750,
    LINESHAPE_VOIGT_KIELKOPF,
    LINESHAPE_END_OF_TABLE
};

int linesum(
    double*,
    const double*,
    const double*,
    const double,
    const gridsize_t,
    const cat_entry_t*,
    const line_coupling_table_entry_t*,
    const int,
    const int,
    const int,
    const double,
    const double,
    const double,
    const double,
    const double*,
    const double*,
    const double*,
    const double,
    int*);
int Qratio(
    const double,
    const double,
    const double*,
    const int,
    const int,
    double *);

#endif /* AM_LINESUM_H */
