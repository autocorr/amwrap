/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* mapping.h                     S. Paine rev. 2014 August 20
*
* Constants and declarations for mapping.c
************************************************************/

#ifndef AM_MAPPING_H
#define AM_MAPPING_H

enum {
    MAPPING_NONE,
    MAPPING_VMR,
    MAPPING_EXP
};

double map_differential(double, double, int);
double map_variable(double, int);
double unmap_differential(double, double, int);
double unmap_variable(double, int);

#endif /* AM_MAPPING_H */
