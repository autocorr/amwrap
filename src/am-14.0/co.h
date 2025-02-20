/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* co.h                          S. Paine rev. 2024 August 06
*
* External declarations for catalog data defined in co.c.
************************************************************/

#ifndef AM_CO_H
#define AM_CO_H

#include "am_types.h"

extern const double co_abundance_tab[];
extern const double co_mass_tab[];

extern const double co_Tref;
extern const cat_entry_t co_cat[];
extern const int co_num_lines;

extern const double co_Qtab[];
extern const int co_Qtab_cols;
extern const int co_Qtab_rows;

#endif /* AM_CO_H */
