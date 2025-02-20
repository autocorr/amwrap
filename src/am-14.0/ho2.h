/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* ho2.h                         S. Paine rev. 2024 August 06
*
* External declarations for catalog data defined in ho2.c.
************************************************************/

#ifndef AM_HO2_H
#define AM_HO2_H

#include "am_types.h"

extern const double ho2_abundance_tab[];
extern const double ho2_mass_tab[];

extern const double ho2_Tref;
extern const cat_entry_t ho2_cat[];
extern const int ho2_num_lines;

extern const double ho2_Qtab[];
extern const int ho2_Qtab_cols;
extern const int ho2_Qtab_rows;

#endif /* AM_HO2_H */
