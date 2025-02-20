/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* no.h                          S. Paine rev. 2024 August 06
*
* External declarations for catalog data defined in no.c.
************************************************************/

#ifndef AM_NO_H
#define AM_NO_H

#include "am_types.h"

extern const double no_abundance_tab[];
extern const double no_mass_tab[];

extern const double no_Tref;
extern const cat_entry_t no_cat[];
extern const int no_num_lines;

extern const double no_Qtab[];
extern const int no_Qtab_cols;
extern const int no_Qtab_rows;

#endif /* AM_NO_H */
