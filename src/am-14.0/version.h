/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* version.h                     S. Paine rev. 2024 August 12
*
* Program version number and declarations for version.c
************************************************************/

#ifndef AM_VERSION_H
#define AM_VERSION_H


/*
 * Cached absorption coefficient files are labeled by program
 * version number.  When line catalog or other spectral data are
 * updated, the program version number must be changed to prevent
 * reuse of computations based on old data.
 *
 * Normally, the user manual version and program version number
 * will be the same, unless the manual falls behind.  (The
 * location of the manual is kept in the string constant
 * AM_DOC_URL[] in doc.c.)
 *
 * Because major and minor version numbers are combined into a
 * single unsigned int in dcache file headers, they are limited
 * to the range 0 - 65535.
 */
enum {
    AM_VERSION_MAJOR     = 14,
    AM_VERSION_MINOR     = 0,
    AM_DOC_VERSION_MAJOR = 14,
    AM_DOC_VERSION_MINOR = 0,
    AM_DOC_YEAR          = 2024
};


void version(FILE*);
void write_version_line(FILE*);

#endif /* AM_VERSION_H */
