/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* no2.h                         S. Paine rev. 2024 August 06
*
* External declarations for catalog data defined in no2.c.
************************************************************/

#ifndef AM_NO2_H
#define AM_NO2_H

#include "am_types.h"

extern const double no2_abundance_tab[];
extern const double no2_mass_tab[];

extern const double no2_Tref;
extern const cat_entry_t no2_cat[];
extern const int no2_num_lines;

extern const double no2_Qtab[];
extern const int no2_Qtab_cols;
extern const int no2_Qtab_rows;

#endif /* AM_NO2_H */
