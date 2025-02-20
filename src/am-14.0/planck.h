/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* planck.h                      S. Paine rev. 2019 January 8
*
* Declarations for rt.c
************************************************************/

#ifndef AM_PLANCK_H
#define AM_PLANCK_H

void   B_Planck(double*, double*, double, gridsize_t, double);
double T_Planck(double, double);
double T_Rayleigh_Jeans(double, double);
void   planck_benchmarks(void);

#endif /* AM_PLANCK_H */
