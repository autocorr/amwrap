/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* fileops.h                  S. Paine rev. 2019 September 25
*
* Constants and declarations for fileops.c
************************************************************/

#ifndef AM_FILEOPS_H
#define AM_FILEOPS_H

#include <stdio.h>

enum {
    AM_MAX_DIRPATH   = 256,
    AM_MAX_FNAMESIZE = 256,
    AM_FIT_EXTSIZE   = 5 /* max extension length on fit output filenames */
};

enum {
    AM_MAX_PATH = (AM_MAX_DIRPATH + AM_MAX_FNAMESIZE)
};

FILE *am_tmpfile(const char*, char*);
int check_for_dir_separator(char*);
int remove_with_retry(const char*);
int rename_with_retry(const char*, const char*);

#endif /* AM_FILEOPS_H */
