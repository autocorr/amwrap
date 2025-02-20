/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* ocs.h                         S. Paine rev. 2024 August 06
*
* External declarations for catalog data defined in ocs.c.
************************************************************/

#ifndef AM_OCS_H
#define AM_OCS_H

#include "am_types.h"

extern const double ocs_abundance_tab[];
extern const double ocs_mass_tab[];

extern const double ocs_Tref;
extern const cat_entry_t ocs_cat[];
extern const int ocs_num_lines;

extern const double ocs_Qtab[];
extern const int ocs_Qtab_cols;
extern const int ocs_Qtab_rows;

#endif /* AM_OCS_H */
