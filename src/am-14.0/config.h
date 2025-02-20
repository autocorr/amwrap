/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* config.h                      S. Paine rev. 2019 January 9
*
* Declarations and enums for config.c
************************************************************/

#ifndef AM_CONFIG_H
#define AM_CONFIG_H

#include "am_types.h"

int parse_config_file(int, char**, model_t*, fit_data_t*, simplex_t*);
int set_config_parameter(char*, char*, int, model_t*, fit_data_t*, simplex_t*);

#endif /* AM_CONFIG_H */
