/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* simplex.h                     S. Paine rev. 2020 October 2
*
* Declarations for simplex.c
************************************************************/

#ifndef AM_SIMPLEX_H
#define AM_SIMPLEX_H

#include "am_types.h"

int add_var_to_simplex(
        simplex_t*,
        const char *,
        double*,
        double,
        double,
        int,
        int);
int          create_null_simplex(simplex_t*);
unsigned int get_simplex_variable_index(simplex_t*, double*);
int          isvar(simplex_t*, double*);
void         free_simplex_entities(simplex_t*);
void         reset_simplex_vertices(simplex_t*);
double       simplex_scaled_diameter(simplex_t*);
double       simplex_scaled_distance(simplex_t*, double*, double*);
double       simplex_variable_range(simplex_t*, unsigned int);

#endif /* AM_SIMPLEX_H */
