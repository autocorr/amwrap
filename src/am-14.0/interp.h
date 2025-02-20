/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* interp.h                         S. Paine rev. 2021 Mar 30
*
* Program version number and declarations for interp.c
************************************************************/

#ifndef AM_INTERP_H
#define AM_INTERP_H

double lin_interp(
        const double, const double, const double, const double, const double);
double log_x_interp(
        const double, const double, const double, const double, const double);
double log_y_interp(
        const double, const double, const double, const double, const double);

#endif /* AM_INTERP_H */
