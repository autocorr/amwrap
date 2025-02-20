/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* hbr.h                         S. Paine rev. 2024 August 06
*
* External declarations for catalog data defined in hbr.c.
************************************************************/

#ifndef AM_HBr_H
#define AM_HBr_H

#include "am_types.h"

extern const double hbr_abundance_tab[];
extern const double hbr_mass_tab[];

extern const double hbr_Tref;
extern const cat_entry_t hbr_cat[];
extern const int hbr_num_lines;

extern const double hbr_Qtab[];
extern const int hbr_Qtab_cols;
extern const int hbr_Qtab_rows;

#endif /* AM_HBr_H */
