/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* ch3cn.h                       S. Paine rev. 2024 August 06
*
* External declarations for catalog data defined in ch3cn.c.
************************************************************/

#ifndef AM_CH3CN_H
#define AM_CH3CN_H

#include "am_types.h"

extern const double ch3cn_abundance_tab[];
extern const double ch3cn_mass_tab[];

extern const double ch3cn_Tref;
extern const cat_entry_t ch3cn_cat[];
extern const int ch3cn_num_lines;

extern const double ch3cn_Qtab[];
extern const int ch3cn_Qtab_cols;
extern const int ch3cn_Qtab_rows;

#endif /* AM_CH3CN_H */
