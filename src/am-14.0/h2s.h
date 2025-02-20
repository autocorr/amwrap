/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* h2s.h                         S. Paine rev. 2024 August 06
*
* External declarations for catalog data defined in h2s.c.
************************************************************/

#ifndef AM_H2S_H
#define AM_H2S_H

#include "am_types.h"

extern const double h2s_abundance_tab[];
extern const double h2s_mass_tab[];

extern const double h2s_Tref;
extern const cat_entry_t h2s_cat[];
extern const int h2s_num_lines;

extern const double h2s_Qtab[];
extern const int h2s_Qtab_cols;
extern const int h2s_Qtab_rows;

#endif /* AM_H2S_H */
