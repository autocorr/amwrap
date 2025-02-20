/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* clo.h                         S. Paine rev. 2024 August 06
*
* External declarations for catalog data defined in clo.c.
************************************************************/

#ifndef AM_ClO_H
#define AM_ClO_H

#include "am_types.h"

extern const double clo_abundance_tab[];
extern const double clo_mass_tab[];

extern const double clo_Tref;
extern const cat_entry_t clo_cat[];
extern const int clo_num_lines;

extern const double clo_Qtab[];
extern const int clo_Qtab_cols;
extern const int clo_Qtab_rows;

#endif /* AM_ClO_H */
