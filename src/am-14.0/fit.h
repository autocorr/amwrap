/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* fit.h                           S. Paine rev. 2024 July 23
*
* Declarations and constants for fit.c
************************************************************/

#ifndef AM_FIT_H
#define AM_FIT_H

#include "am_types.h"

/*
 * Fit estimator types
 */
extern struct fit_estimator_type_tabentry {
    const char *name; /* text name of type */
} fit_estimator_type[]; /* defined in fit.c */

enum {
    FIT_ESTIMATOR_NONE,
    FIT_ESTIMATOR_ABSOLUTE_RESIDUALS,
    FIT_ESTIMATOR_SQUARED_RESIDUALS,
    FIT_ESTIMATOR_LOGARITHMIC,
    FIT_ESTIMATOR_END_OF_TABLE
};

/*
 * Fit output mode control bits.
 */
enum {
    FIT_OUTPUT_CONFIG   = 0x1,
    FIT_OUTPUT_SPECTRUM = 0x2,
    FIT_OUTPUT_RESIDUAL = 0x4,
    FIT_OUTPUT_VERBOSE  = 0x8
};


int  fit(model_t*, model_t*, fit_data_t*, simplex_t*);
char *fit_data_delimiters(const char*);
char *fit_data_format(const char*);
void report_fit_env_info(FILE*);
void reset_estimated_residuals(fit_data_t*);

#endif /* AM_FIT_H */
