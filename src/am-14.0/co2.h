/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* co2.h                         S. Paine rev. 2024 August 06
*
* External declarations for catalog data defined in co2.c.
************************************************************/

#ifndef AM_CO2_H
#define AM_CO2_H

#include "am_types.h"

extern const double co2_abundance_tab[];
extern const double co2_mass_tab[];

extern const double co2_Tref;
extern const cat_entry_t co2_cat[];
extern const int co2_num_lines;

extern const double co2_Qtab[];
extern const int co2_Qtab_cols;
extern const int co2_Qtab_rows;

#endif /* AM_CO2_H */
