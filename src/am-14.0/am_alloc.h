/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* am_alloc.h                      S. Paine rev. 2021 March 3
*
* Declarations for am_alloc.c
************************************************************/

#ifndef AM_AM_ALLOC_H
#define AM_AM_ALLOC_H

#include "am_types.h"

int  add_column(layer_t*, int);
int  add_layer(model_t*, int);
int  alloc_layer_arrays(model_t*, int);
int  alloc_model_arrays(model_t*);
int  alloc_jacobians(model_t*, simplex_t*);
int  clear_layer_kcache_entries(model_t*, int);
int  copy_layer_allocations(model_t*, layer_t*, layer_t*);
int  copy_layer_dimensions(layer_t*, layer_t*);
int  copy_model_dimensions(model_t*, model_t*);
int  delete_layer(model_t*, int);
int  init_kcache(abscoeff_t*, int);
void free_fit_data_entities(fit_data_t*);
void free_model_entities(model_t*);
void free_jacobians(simplex_t*);
int  grow_fit_data_arrays(fit_data_t*);

#endif /* AM_AM_ALLOC_H */
