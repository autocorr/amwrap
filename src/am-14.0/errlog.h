/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* errlog.h                      S. Paine rev. 2014 August 28
*
* Declarations for errlog.c
************************************************************/

#ifndef AM_ERRLOG_H
#define AM_ERRLOG_H

int errlog(const int, const int);
int errstat(void);
int errtest(const int);
int print_errlog(void);

#endif /* AM_ERRLOG_H */
