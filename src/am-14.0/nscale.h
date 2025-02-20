/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* nscale.h                      S. Paine rev. 2019 January 9
*
* Declarations for nscale.c
************************************************************/

#ifndef AM_NSCALE_H
#define AM_NSCALE_H

#include "am_types.h"

typedef struct Nscale_list_t Nscale_list_t;

struct Nscale_list_t {
    Nscale_list_t *next;
    int col_typenum;
    int tagnum;
    double Nscale;
};

Nscale_list_t *create_Nscale_list_entry(const int, const int, const double);
Nscale_list_t *find_Nscale_list_entry(const int, const int);
double lookup_Nscale(const int, const int);
void free_Nscale_list(void);
Nscale_list_t *Nscale_list_head(void);

#endif /* AM_NSCALE_H */
