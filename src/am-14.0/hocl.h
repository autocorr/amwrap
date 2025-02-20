/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* hocl.h                        S. Paine rev. 2024 August 06
*
* External declarations for catalog data defined in hocl.c.
************************************************************/

#ifndef AM_HOCl_H
#define AM_HOCl_H

#include "am_types.h"

extern const double hocl_abundance_tab[];
extern const double hocl_mass_tab[];

extern const double hocl_Tref;
extern const cat_entry_t hocl_cat[];
extern const int hocl_num_lines;

extern const double hocl_Qtab[];
extern const int hocl_Qtab_cols;
extern const int hocl_Qtab_rows;

#endif /* AM_HOCl_H */
