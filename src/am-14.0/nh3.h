/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* nh3.h                         S. Paine rev. 2024 August 06
*
* External declarations for catalog data defined in nh3.c.
************************************************************/

#ifndef AM_NH3_H
#define AM_NH3_H

#include "am_types.h"

extern const double nh3_abundance_tab[];
extern const double nh3_mass_tab[];

extern const double nh3_Tref;
extern const cat_entry_t nh3_cat[];
extern const int nh3_num_lines;

extern const double nh3_Qtab[];
extern const int nh3_Qtab_cols;
extern const int nh3_Qtab_rows;

#endif /* AM_NH3_H */
