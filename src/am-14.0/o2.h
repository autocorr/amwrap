/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* o2.h                          S. Paine rev. 2024 August 06
*
* External declarations for catalog data defined in o2.c.
************************************************************/

#ifndef AM_O2_H
#define AM_O2_H

#include "am_types.h"

extern const double o2_abundance_tab[];
extern const double o2_mass_tab[];

extern const double o2_Tref;
extern const cat_entry_t o2_coupled_cat[];
extern const int o2_num_coupled_lines;
extern const line_coupling_table_entry_t o2_line_coupling_coeffs[];
extern const cat_entry_t o2_uncoupled_cat[];
extern const int o2_num_uncoupled_lines;

extern const double o2_Qtab[];
extern const int o2_Qtab_cols;
extern const int o2_Qtab_rows;

#endif /* AM_O2_H */
