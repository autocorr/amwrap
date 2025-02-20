/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* jacobian.h                     S. Paine rev. 2021 April 15
*
* Declarations for jacobian.c
************************************************************/

#ifndef AM_JACOBIAN_H
#define AM_JACOBIAN_H

#include "am_types.h"

int compute_jacobians(model_t*, model_t*, simplex_t*);

#endif /* AM_JACOBIAN_H */
