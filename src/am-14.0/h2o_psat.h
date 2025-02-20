/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* h2o_psat.h                    S. Paine rev. 2014 August 29
*
* Declarations for h2o_psat.c
************************************************************/

#ifndef AM_H2O_PSAT_H
#define AM_H2O_PSAT_H

extern double h2o_abundance_tab[];

double H2O_liquid_Psat(double, int);
double H2O_ice_Psat(double, int);

#endif /* AM_H2O_PSAT_H */
