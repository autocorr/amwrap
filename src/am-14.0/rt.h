/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* rt.h                          S. Paine rev. 2019 January 9
*
* Declarations for rt.c
************************************************************/

#ifndef AM_RT_H
#define AM_RT_H

#include "am_types.h"

/*
 * Directions for radiative transfer.
 */
enum {
    RT_DOWN,
    RT_UP
};

void compute_spectral_radiance(model_t*, model_t*);
int  compute_opacity_spectrum(model_t*, model_t*);

#endif /* AM_RT_H */
