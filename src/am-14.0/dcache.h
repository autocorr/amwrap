/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* dcache.h                      S. Paine rev. 2019 January 9
*
* Constants and declarations for dcache.c 
************************************************************/

#ifndef AM_DCACHE_H
#define AM_DCACHE_H

#include <stdio.h>

#include "am_types.h"

/*
 * Constants used in dcache_log()
 */
enum {
    DCACHE_HIT,
    DCACHE_MISS,
    DCACHE_DISCARD,
    DCACHE_CLEAR_LOG
};


void dcache_log(int);
int dcache_lookup_absorption_coefficient(
    model_t*,
    int,
    int,
    int,
    double*,
    double);
int dcache_save_absorption_coefficient(
    model_t*,
    int,
    int,
    int,
    double*,
    double);
void report_dcache_env_info(FILE *);
void report_dcache_log_data(FILE *);

#endif /* AM_DCACHE_H */
