/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* output.h                       S. Paine rev. 2024 April 26
*
* Functions, constants, and data structures related to
* program output.
************************************************************/

#ifndef AM_OUTPUT_H
#define AM_OUTPUT_H

#include <stdio.h>

#include "am_types.h"

/*
 * The global table output[] is an array of the following structures
 * defined in output.c.
 *
 * For each output type enumerated below,
 *   spectrum[i] is the value of the spectrum at frequency f[i].
 *   jacobian[j][i] is the derivative at f[i] with respect to the jth
 *     user-defined differentiation variable.
 *   jacobian_err[j][i] is a corresponding error estimate for the
 *     finite difference approximation for the derivative.
 */
extern struct output_tabentry {
    const char *name;               /* config file name for output array   */
    double  *spectrum;              /* pointer to spectral data array      */
    double **jacobian;              /* Jacobian                            */
    double **jacobian_round_err;    /* Jacobian estimated rounding error   */
    double **jacobian_trunc_err;    /* Jacobian estimated truncation error */
    double **jacobian_sorted_err;   /* Jacobian sorted total err for cdf   */
    int default_unitnum;            /* index number of default unit        */
    int unitnum;                    /* index number of output unit         */
    int k_typenum;                  /* typenum for k output header         */
    int flags;                      /* flag bits defined below             */ 
    int format;                     /* format                              */ 
} output[];


enum {
    OUTPUT_NONE,
    OUTPUT_FREQUENCY,
    OUTPUT_OPACITY,
    OUTPUT_TRANSMITTANCE,
    OUTPUT_RADIANCE,
    OUTPUT_RADIANCE_DIFF,
    OUTPUT_TB_PLANCK,
    OUTPUT_TB_RAYLEIGH_JEANS,
    OUTPUT_TSYS,
    OUTPUT_Y,
    OUTPUT_DELAY,
    OUTPUT_FREE_SPACE_LOSS,
    OUTPUT_K,
    OUTPUT_END_OF_TABLE
};

/*
 * The last table entry is used to indicate states and actions
 * that affect all outputs.
 */
enum {ALL_OUTPUTS = OUTPUT_END_OF_TABLE};

/*
 * Output flag bits.
 */
enum {
    OUTPUT_USER          =   0x1, /* spectrum is user-requested output      */
    OUTPUT_FITTED        =   0x2, /* spectrum needed for fit                */
    OUTPUT_ACTIVE        =   0x4, /* computation of this spectrum enabled   */
    ILS_ALLOWED          =   0x8, /* spectrum can be convolved with ILS     */
    ILS_APPLIED          =  0x10, /* spectrum is convolved with ILS, if any */
    JACOBIAN_ALLOWED     =  0x20, /* Jacobian can be computed               */
    OUTPUT_JACOBIAN      =  0x40, /* Jacobian is computed                   */
    OUTPUT_JACOBIAN_ERRS =  0x80, /* output error estimates with Jacobians  */ 
    OUTPUT_AM_ONLY       = 0x100, /* atmospheric model only, no spectra     */
    OUTPUT_HEADERS       = 0x200, /* include headers in output spectra      */
};

/*
 * Output formats.
 */
extern const char* output_format_name[];

enum {
    OUTPUT_FORMAT_TEXT,
    OUTPUT_FORMAT_CSV,
    OUTPUT_FORMAT_NPY,
    OUTPUT_FORMAT_END_OF_TABLE
};

/*
 * outcol[] is an index table, defined in output.c, which maps
 * from output column numbers to entry numbers in the output[]
 * table.  For example, if the first column of program output is
 * frequency, then outcol[0] will be equal to OUTPUT_FREQUENCY.
 */
extern int outcol[];


int  set_active_outputs(const int);
void write_model_config_data(FILE*, model_t*, fit_data_t*, simplex_t*);
void write_PTmode(FILE*, int);
void write_model_spectra(FILE*, model_t*, simplex_t*);
void write_fit_residuals(FILE*, fit_data_t*);


#endif /* AM_OUTPUT_H */
