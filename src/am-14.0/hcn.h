/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* hcn.h                         S. Paine rev. 2024 August 06
*
* External declarations for catalog data defined in hcn.c.
************************************************************/

#ifndef AM_HCN_H
#define AM_HCN_H

#include "am_types.h"

extern const double hcn_abundance_tab[];
extern const double hcn_mass_tab[];

extern const double hcn_Tref;
extern const cat_entry_t hcn_cat[];
extern const int hcn_num_lines;

extern const double hcn_Qtab[];
extern const int hcn_Qtab_cols;
extern const int hcn_Qtab_rows;

#endif /* AM_HCN_H */
