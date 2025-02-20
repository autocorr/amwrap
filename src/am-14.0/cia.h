/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* cia.h                         S. Paine rev. 2019 January 9
*
* Constants relating to collision-induced absorption, and
* declarations for cia.c
************************************************************/

#ifndef AM_CIA_H
#define AM_CIA_H

/*
 * Integrated intensity ratios used for N2-air and O2-air CIA, from
 * Table 2 of J. Boissoles, C. Boulet, R.H. Tipping, Alex Brown,
 * and Q. Ma, JQSRT 82:505 (2003).
 */
#define S_N2O2_ON_S_N2N2 1.143
#define S_O2N2_ON_S_O2O2 0.822

int N2N2_cia(double*, const double*, const gridsize_t, const double);
int O2O2_cia(double*, const double*, const gridsize_t, const double);

#endif /* AM_CIA_H */
