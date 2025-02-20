/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* so2.h                         S. Paine rev. 2024 August 06
*
* External declarations for catalog data defined in so2.c.
************************************************************/

#ifndef AM_SO2_H
#define AM_SO2_H

#include "am_types.h"

extern const double so2_abundance_tab[];
extern const double so2_mass_tab[];

extern const double so2_Tref;
extern const cat_entry_t so2_cat[];
extern const int so2_num_lines;

extern const double so2_Qtab[];
extern const int so2_Qtab_cols;
extern const int so2_Qtab_rows;

#endif /* AM_SO2_H */
