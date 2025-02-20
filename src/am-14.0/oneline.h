/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* oneline.h                   S. Paine rev. 2014 September 2
*
* External declarations for catalog data defined in oneline.c
************************************************************/

#include "am_types.h"

extern const double oneline_abundance_tab[];
extern const double oneline_mass_tab[];
extern const int oneline_num_isotopologues;

extern const double oneline_Tref;
extern const cat_entry_t oneline_cat[];
extern const int oneline_num_lines;

extern const double oneline_Qtab[];
extern const int oneline_Qtab_cols;
extern const int oneline_Qtab_rows;
