/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* kcache.h                      S. Paine rev. 2019 January 9
*
* Constants and declarations for kcache.c
************************************************************/

#ifndef AM_KCACHE_H
#define AM_KCACHE_H

#include <stdio.h>

#include "am_types.h"

/*
 * Constants used in kcache_log()
 */
enum {
    KCACHE_HIT,
    KCACHE_MISS,
    KCACHE_DISCARD,
    KCACHE_CLEAR_LOG
};

double *kcache_alloc(model_t*);
void kcache_free(double *);
void kcache_free_all(void);
void kcache_lock_entry(double *);
void kcache_log(int);
void kcache_unlock_entry(double *);
void report_kcache_env_info(FILE *);
void report_kcache_log_data(FILE *);

#endif /* AM_KCACHE_H */
