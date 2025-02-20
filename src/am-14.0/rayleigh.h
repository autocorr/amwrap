/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* rayleigh.h                   S. Paine rev. 2017 January 13
*
* Declarations for rayleigh.c
************************************************************/

#ifndef AM_RAYLEIGH_H
#define AM_RAYLEIGH_H

int Rayleigh_mol_abscoeff(
    double*,
    const double*,
    const gridsize_t,
    const double*,
    const double*,
    const double);

#endif /* AM_RAYLEIGH_H */
