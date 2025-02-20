/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* ch3oh.h                       S. Paine rev. 2024 August 06
*
* External declarations for catalog data defined in ch3oh.c.
************************************************************/

#ifndef AM_CH3OH_H
#define AM_CH3OH_H

#include "am_types.h"

extern const double ch3oh_abundance_tab[];
extern const double ch3oh_mass_tab[];

extern const double ch3oh_Tref;
extern const cat_entry_t ch3oh_cat[];
extern const int ch3oh_num_lines;

extern const double ch3oh_Qtab[];
extern const int ch3oh_Qtab_cols;
extern const int ch3oh_Qtab_rows;

#endif /* AM_CH3OH_H */
