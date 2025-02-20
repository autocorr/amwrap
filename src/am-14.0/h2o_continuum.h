/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* h2o_continuum.h               S. Paine rev. 2019 January 9
*
* Declarations for h2o_continuum.c
************************************************************/

#ifndef AM_H2O_CONTINUUM_H
#define AM_H2O_CONTINUUM_H

#include "am_types.h"

int H2O_air_continuum(
    double *kb,
    const double *f,
    const gridsize_t ngrid,
    const double T);

int H2O_self_continuum(
    double *kb,
    const double *f,
    const gridsize_t ngrid,
    const double T);

#endif /* AM_H2O_CONTINUUM_H */
