/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* o.h                           S. Paine rev. 2024 August 06
*
* External declarations for catalog data defined in o.c.
************************************************************/

#ifndef AM_O_H
#define AM_O_H

#include "am_types.h"

extern const double o_abundance_tab[];
extern const double o_mass_tab[];

extern const double o_Tref;
extern const cat_entry_t o_cat[];
extern const int o_num_lines;

extern const double o_Qtab[];
extern const int o_Qtab_cols;
extern const int o_Qtab_rows;

#endif /* AM_O_H */
