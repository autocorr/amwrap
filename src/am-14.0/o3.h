/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* o3.h                          S. Paine rev. 2024 August 06
*
* External declarations for catalog data defined in o3.c.
************************************************************/

#ifndef AM_O3_H
#define AM_O3_H

#include "am_types.h"

extern const double o3_abundance_tab[];
extern const double o3_mass_tab[];

extern const double o3_Tref;
extern const cat_entry_t o3_cat[];
extern const int o3_num_lines;

extern const double o3_Qtab[];
extern const int o3_Qtab_cols;
extern const int o3_Qtab_rows;

#endif /* AM_O3_H */
