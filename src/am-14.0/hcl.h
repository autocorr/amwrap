/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* hcl.h                         S. Paine rev. 2024 August 06
*
* External declarations for catalog data defined in hcl.c.
************************************************************/

#ifndef AM_HCl_H
#define AM_HCl_H

#include "am_types.h"

extern const double hcl_abundance_tab[];
extern const double hcl_mass_tab[];

extern const double hcl_Tref;
extern const cat_entry_t hcl_cat[];
extern const int hcl_num_lines;

extern const double hcl_Qtab[];
extern const int hcl_Qtab_cols;
extern const int hcl_Qtab_rows;

#endif /* AM_HCl_H */
