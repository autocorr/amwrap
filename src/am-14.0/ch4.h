/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* ch4.h                         S. Paine rev. 2024 August 06
*
* External declarations for catalog data defined in ch4.c.
************************************************************/

#ifndef AM_CH4_H
#define AM_CH4_H

#include "am_types.h"

extern const double ch4_abundance_tab[];
extern const double ch4_mass_tab[];

extern const double ch4_Tref;
extern const cat_entry_t ch4_cat[];
extern const int ch4_num_lines;

extern const double ch4_Qtab[];
extern const int ch4_Qtab_cols;
extern const int ch4_Qtab_rows;

#endif /* AM_CH4_H */
