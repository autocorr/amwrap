/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* hno3.h                        S. Paine rev. 2024 August 06
*
* External declarations for catalog data defined in hno3.c.
************************************************************/

#ifndef AM_HNO3_H
#define AM_HNO3_H

#include "am_types.h"

extern const double hno3_abundance_tab[];
extern const double hno3_mass_tab[];

extern const double hno3_Tref;
extern const cat_entry_t hno3_cat[];
extern const int hno3_num_lines;

extern const double hno3_Qtab[];
extern const int hno3_Qtab_cols;
extern const int hno3_Qtab_rows;

#endif /* AM_HNO3_H */
