/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* hf.h                          S. Paine rev. 2024 August 06
*
* External declarations for catalog data defined in hf.c.
************************************************************/

#ifndef AM_HF_H
#define AM_HF_H

#include "am_types.h"

extern const double hf_abundance_tab[];
extern const double hf_mass_tab[];

extern const double hf_Tref;
extern const cat_entry_t hf_cat[];
extern const int hf_num_lines;

extern const double hf_Qtab[];
extern const int hf_Qtab_cols;
extern const int hf_Qtab_rows;

#endif /* AM_HF_H */
