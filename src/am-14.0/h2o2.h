/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* h2o2.h                        S. Paine rev. 2024 August 06
*
* External declarations for catalog data defined in h2o2.c.
************************************************************/

#ifndef AM_H2O2_H
#define AM_H2O2_H

#include "am_types.h"

extern const double h2o2_abundance_tab[];
extern const double h2o2_mass_tab[];

extern const double h2o2_Tref;
extern const cat_entry_t h2o2_cat[];
extern const int h2o2_num_lines;

extern const double h2o2_Qtab[];
extern const int h2o2_Qtab_cols;
extern const int h2o2_Qtab_rows;

#endif /* AM_H2O2_H */
