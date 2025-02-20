/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* n2o.h                         S. Paine rev. 2024 August 06
*
* External declarations for catalog data defined in n2o.c.
************************************************************/

#ifndef AM_N2O_H
#define AM_N2O_H

#include "am_types.h"

extern const double n2o_abundance_tab[];
extern const double n2o_mass_tab[];

extern const double n2o_Tref;
extern const cat_entry_t n2o_cat[];
extern const int n2o_num_lines;

extern const double n2o_Qtab[];
extern const int n2o_Qtab_cols;
extern const int n2o_Qtab_rows;

#endif /* AM_N2O_H */
