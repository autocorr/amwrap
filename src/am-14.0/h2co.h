/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* h2co.h                        S. Paine rev. 2024 August 06
*
* External declarations for catalog data defined in h2co.c.
************************************************************/

#ifndef AM_H2CO_H
#define AM_H2CO_H

#include "am_types.h"

extern const double h2co_abundance_tab[];
extern const double h2co_mass_tab[];

extern const double h2co_Tref;
extern const cat_entry_t h2co_cat[];
extern const int h2co_num_lines;

extern const double h2co_Qtab[];
extern const int h2co_Qtab_cols;
extern const int h2co_Qtab_rows;

#endif /* AM_H2CO_H */
