/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* config.c                        S. Paine rev. 2024 July 22
*
* model configuration
************************************************************/

#include <ctype.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

#include "abscoeff.h"
#include "am_alloc.h"
#include "am_types.h"
#include "column.h"
#include "config.h"
#include "errlog.h"
#include "fileops.h"
#include "fit.h"
#include "ils.h"
#include "layer.h"
#include "linesum.h"
#include "mapping.h"
#include "math_const.h"
#include "model.h"
#include "nscale.h"
#include "output.h"
#include "phys_const.h"
#include "simplex.h"
#include "tags.h"
#include "units.h"

/*
 * This is the maximum grid frequency for am, corresponding to the
 * upper limit of internal line catalogs and continuum models.
 */
#define AM_MAX_FREQ     15000.0 /* GHz */

enum {
    LINEBUFSIZE     = 1024, /* size of config file line buffer               */
    VARNAME_BUFSIZE = 1024, /* sizeof buffer for constructing variable names */
    MAX_NTOK        = 32    /* max number of tokens on a config file line    */
};

/*
 * Configuration file keywords.
 */
enum kw_enum {
    KW_COL,
    KW_COLUMN,
    KW_DG_DZ,
    KW_DP,
    KW_DSB_UTOL_RATIO,
    KW_F,
    KW_FIF,
    KW_FIT,
    KW_FITS,
    KW_FIT_DATA_COLS,
    KW_FIT_DATA_COLUMNS,
    KW_FIT_DATA_DELIMITERS,
    KW_FIT_DATA_FORMAT,
    KW_FIT_DATA_UNITS,
    KW_FIT_ESTIMATOR,
    KW_FIT_ITER,
    KW_FIT_OUTPUT,
    KW_FIT_REINIT,
    KW_FIT_TOL,
    KW_FIT_TRACK_RESIDUALS,
    KW_FIT_VERBOSE,
    KW_FOUT,
    KW_FREQ,
    KW_FREQUENCY,
    KW_G,
    KW_GEOMETRY,
    KW_H,
    KW_HEADERS,
    KW_IFSPEC,
    KW_ILS,
    KW_ILSMODE,
    KW_JACOBIAN,
    KW_KCACHE,
    KW_LAYER,
    KW_LINESHAPE,
    KW_MAIR,
    KW_NSCALE,
    KW_OUTPUT,
    KW_P,
    KW_PBASE,
    KW_POBS,
    KW_PSOURCE,
    KW_PTAN,
    KW_PTMODE,
    KW_R0,
    KW_REFRACT,
    KW_REVERSE,
    KW_RH_OFFSET,
    KW_RH_SCALE,
    KW_RUNTIME,
    KW_RX_GAIN_FACTOR,
    KW_SEC_ZA,
    KW_SELFBROAD_VMR_TOL,
    KW_SIMPLEX_LOG,
    KW_T,
    KW_T0,
    KW_TBASE,
    KW_TREF,
    KW_TRX,
    KW_TOL,
    KW_TOLERANCE,
    KW_Z0,
    KW_ZA,
    KW_ZOBS,
    KW_ZSOURCE,
    KW_ZTAN,
    KW_END_OF_TABLE
};

static struct kw_table_entry {
    const char *name;
} const kw_table[] = {
    {"col"},
    {"column"},
    {"dg_dz"},
    {"dP"},
    {"dsb_utol_ratio"},
    {"f"},
    {"fif"},
    {"fit"},
    {"fits"},
    {"fit_data_cols"},
    {"fit_data_columns"},
    {"fit_data_delimiters"},
    {"fit_data_format"},
    {"fit_data_units"},
    {"fit_estimator"},
    {"fit_iter"},
    {"fit_output"},
    {"fit_reinit"},
    {"fit_tol"},
    {"fit_track_residuals"},
    {"fit_verbose"},
    {"fout"},
    {"freq"},
    {"frequency"},
    {"g"},
    {"geometry"},
    {"h"},
    {"headers"},
    {"ifspec"},
    {"ils"},
    {"ilsmode"},
    {"jacobian"},
    {"kcache"},
    {"layer"},
    {"lineshape"},
    {"Mair"},
    {"Nscale"},
    {"output"},
    {"P"},
    {"Pbase"},
    {"Pobs"},
    {"Psource"},
    {"Ptan"},
    {"PTmode"},
    {"R0"},
    {"refract"},
    {"reverse"},
    {"RH_offset"},
    {"RH_scale"},
    {"runtime"},
    {"rx_gain_factor"},
    {"sec_za"},
    {"selfbroad_vmr_tol"},
    {"simplex_log"},
    {"T"},
    {"T0"},
    {"Tbase"},
    {"Tref"},
    {"Trx"},
    {"tol"},
    {"tolerance"},
    {"z0"},
    {"za"},
    {"zobs"},
    {"zsource"},
    {"ztan"},
    {""}
};

/*
 * On a configuration file line, the occurrence of any of the
 * characters in COMMENTSTR as the first character of a non-quoted
 * token will cause that character and the rest of the line to be
 * ignored.
 *
 * HELP_CHAR is a special comment character.  If it is the first
 * character on a line, it flags that line as help text.  Help text
 * is echoed to stderr if the user fails to supply all the command
 * line parameters required by a config file.
 */
static const char COMMENTSTR[] = "~!@#$^&*()|{}[]:;?<>/";
static const char HELP_CHAR    = '?';

/*
 * The following buffer is used for constructing variable names in the
 * various keyword processing functions.
 */
static char VARNAME[VARNAME_BUFSIZE];


static int   CONFIG_FILE_LINENUM;
static char* CONFIG_FILE_NAME;

static int  ci_strcmp(const char*, const char*);
static int  get_bool_val(int*, const char*);
static int  get_col_typenum(const char*);
static int  get_dbl_val(
        double*, const char*, const char*, const int, const int);
static int  get_keyword_num(const char*);
static int  get_nonneg_dbl_val(
        double*, const char*, const char*, const int, const int);
static int  get_nonneg_int_val(int*, const char*);
static int  get_pos_dbl_val(
        double*, const char*, const char*, const int, const int);
static int  get_RH_value(double*, const char*);
static void init_freq_grid(model_t*);
static int  is_numeric_token(const char *tok);
enum load_fit_file_list_retval {
    FILELIST_OK,
    FILELIST_CANNOT_OPEN,
    FILELIST_ERROR
};
static void list_allowed_keywords(const char**, const int);
static void list_allowed_output_array_names(void);
static void list_allowed_units(const int);
static void list_line_by_line_abscoeffs(void);
static int  load_fit_file_list(const char*, fit_data_t*);

/*
 * parse_config_line() returns the number of tokens found, or
 * a negative value on error.  These are the error values.
 */
enum parse_config_line_retval {
    SUBST_PARAMETER_MISSING      = -1,
    SUBST_PARAMETER_SYNTAX_ERROR = -2
};
static int  parse_config_line(char*, int, char**, char**);
static void parse_error(const char*, ...);
static void post_config_hints(model_t*);

/*
 * The functions below process a configuration line after it has
 * been split into tokens.  The first token of each line is a
 * keyword which selects the function used to process the
 * remainder of the line.  These functions return one of the
 * values enumerated below.
 *
 * KWFUNC_SYNTAX_ERROR_CONTINUE signals that the line contains
 * either a true syntax error or an out-of-bounds parameter.  In
 * such cases, parsing of the config file continues, checking for
 * possible additional errors, but the model will not run.
 *
 * KWFUNC_SYNTAX_ERROR_STOP is for syntax errors that may cause a
 * large number of subsequent errors.  To avoid generating a
 * confusing slew of error messages, config file processing is
 * terminated.
 *
 * KWFUNC_FATAL_ERROR signals more serious problems, like failed
 * memory allocation.  In this case, config file processing is
 * terminated and the program will shut down.
 */
enum kwfunc_retval {
    KWFUNC_SUCCESS,
    KWFUNC_SYNTAX_ERROR_CONTINUE,
    KWFUNC_SYNTAX_ERROR_STOP,
    KWFUNC_FATAL_ERROR
};
static int kwfunc_column(char*[], const int, model_t*, simplex_t*);
static int kwfunc_dg_dz(char*[], const int, model_t*);
static int kwfunc_dP(char*[], const int, model_t*, simplex_t*, int);
static int kwfunc_dsb_utol_ratio(char*[], const int, model_t*, simplex_t*);
static int kwfunc_fif(char*[], const int, model_t*);
static int kwfunc_fit(char*[], const int, fit_data_t*);
static int kwfunc_fits(char*[], const int, fit_data_t*);
static int kwfunc_fit_data_columns(char*[], const int, fit_data_t*);
static int kwfunc_fit_data_delimiters(char*[], const int);
static int kwfunc_fit_data_format(char*[], const int);
static int kwfunc_fit_data_units(char*[], const int, fit_data_t*);
static int kwfunc_fit_estimator(char*[], const int, fit_data_t*);
static int kwfunc_fit_iter(char*[], const int, fit_data_t*);
static int kwfunc_fit_output(char*[], const int, fit_data_t*);
static int kwfunc_fit_reinit(char*[], const int, simplex_t*);
static int kwfunc_fit_tol(char*[], const int, fit_data_t*, simplex_t*);
static int kwfunc_fit_track_residuals(char*[], const int, fit_data_t*);
static int kwfunc_fit_verbose(char*[], const int, fit_data_t*);
static int kwfunc_fout(char*[], const int, model_t*);
static int kwfunc_frequency(char*[], const int, model_t*);
static int kwfunc_g(char*[], const int, model_t*);
static int kwfunc_geometry(char*[], const int, model_t*);
static int kwfunc_h(char*[], const int, model_t*, simplex_t*, int);
static int kwfunc_headers(char*[], const int);
static int kwfunc_ils(char*[], const int, model_t*, simplex_t*);
static int kwfunc_ifspec(char*[], const int, model_t*, simplex_t*);
static int kwfunc_ilsmode(char*[], const int, model_t*, simplex_t*);
static int kwfunc_jacobian(char*[], const int);
static int kwfunc_kcache(char*[], const int, model_t*);
static int kwfunc_layer(char*[], const int, model_t*);
static int kwfunc_lineshape(char*[], const int, model_t*);
static int kwfunc_Mair(char*[], const int, model_t*, int);
static int kwfunc_Nscale(char*[], const int, simplex_t*);
static int kwfunc_output(char*[], const int);
static int kwfunc_P(char*[], const int, model_t*, simplex_t*, int);
static int kwfunc_Pbase(char*[], const int, model_t*, int);
static int kwfunc_Pobs(char*[], const int, model_t*, simplex_t*);
static int kwfunc_Psource(char*[], const int, model_t*, simplex_t*);
static int kwfunc_Ptan(char*[], const int, model_t*, simplex_t*);
static int kwfunc_PTmode(char*[], const int, model_t*);
static int kwfunc_R0(char*[], const int, model_t*);
static int kwfunc_refract(char*[], const int, model_t*);
static int kwfunc_reverse(char*[], const int, model_t*);
static int kwfunc_RH_offset(char*[], const int, model_t*, simplex_t*);
static int kwfunc_RH_scale(char*[], const int, model_t*, simplex_t*);
static int kwfunc_runtime(char*[], const int, model_t*);
static int kwfunc_rx_gain_factor(char*[], const int, model_t*, simplex_t*);
static int kwfunc_sec_za(char*[], const int, model_t*, simplex_t*);
static int kwfunc_selfbroad_vmr_tol(char*[], const int, model_t*);
static int kwfunc_simplex_log(char*[], const int, simplex_t*);
static int kwfunc_T(char*[], const int, model_t*, simplex_t*, int);
static int kwfunc_T0(char*[], const int, model_t*, simplex_t*);
static int kwfunc_Tbase(char*[], const int, model_t*, simplex_t*, int);
static int kwfunc_Tref(char*[], const int, model_t*, simplex_t*);
static int kwfunc_Trx(char*[], const int, model_t*, simplex_t*);
static int kwfunc_tol(char*[], const int, model_t*);
static int kwfunc_z0(char*[], const int, model_t*);
static int kwfunc_za(char*[], const int, model_t*, simplex_t*);
static int kwfunc_zobs(char*[], const int, model_t*, simplex_t*);
static int kwfunc_zsource(char*[], const int, model_t*, simplex_t*);
static int kwfunc_ztan(char*[], const int, model_t*, simplex_t*);


/***********************************************************
* int parse_config_file(
*         int argc,
*         char **argv,
*         model_t *model,
*         fit_data_t *fit_data,
*         simplex_t *simplex)
*
* Purpose:
*   Parses the model configuration file.  Sets up the
*   model, fit_data, and simplex data structures.
*
* Arguments:
*   int        argc      - command line argument count
*   char       **argv    - command line argument list
*   model_t    *model    - pointer to model structure
*   fit_data_t *fit_data - pointer to fit data structure
*   simplex_t  *simplex  - pointer to simplex data structure
*
* Return:
*   0 on success, 1 otherwise
************************************************************/

int parse_config_file(
        int argc,
        char **argv,
        model_t *model,
        fit_data_t *fit_data,
        simplex_t *simplex)
{
    FILE *cfgfile;
    char buf[LINEBUFSIZE];
    char *tok[MAX_NTOK];
    int  ntok;
    int  status;
    int  syntax_errflag = 0;
    int  param_errflag  = 0;

    /*
     * Open the configuration file.
     */
    if ((strcmp(argv[1], "stdin") == 0) || (strcmp(argv[1], "-") == 0)) {
        cfgfile = stdin;
    } else if ((cfgfile = fopen(argv[1], "r")) == NULL) {
        fprintf(stderr, "am : cannot open file %s\n", argv[1]);
        errlog(19, 0);
        return 1;
    }
    CONFIG_FILE_NAME = argv[1];
    /*
     * Get lines from the config file until EOF, or until a line with
     * a missing substitution parameter is encountered.
     */
    CONFIG_FILE_LINENUM = 0;
    while (fgets(buf, (int)sizeof(buf), cfgfile) != NULL) {
        ++CONFIG_FILE_LINENUM;
        if (strlen(buf) == 0)
            continue;
        ntok = parse_config_line(buf, argc, argv, tok);
        /*
         * If any of the following is true, the line is not processed
         * further.
         */
        if (ntok == 0) {
            /* No tokens found (blank line) */
            continue;
        } else if (strchr(COMMENTSTR, tok[0][0]) != NULL) {
            /* First token starts with a comment character */
            continue;
        } else if (ntok == MAX_NTOK) {
            /* Too many tokens */
            parse_error("The line contains too many tokens.\n");
            syntax_errflag = 1;
            continue;
        } else if (ntok < 0) {
            /*
             * Substitution parameter errors typically cause
             * lots of errors on subsequent lines, so break out
             * of the line parsing loop.  If the error was
             * a missing command-line argument, set a flag to
             * trigger printing any embedded usage message.
             */
            parse_error("The file was not checked beyond this line.\n");
            syntax_errflag = 1;
            if (ntok == SUBST_PARAMETER_MISSING)
                param_errflag  = 1;
            break;
        } else if (ntok == SUBST_PARAMETER_SYNTAX_ERROR) {
            parse_error("The file was not checked beyond this line.\n");
        }
        /*
         * Process the line according to the keyword matching the
         * first token, or report an error if no keyword match is
         * found.
         */
        switch (get_keyword_num(tok[0])) {
            case KW_COL:
            case KW_COLUMN:
                status = kwfunc_column(tok, ntok, model, simplex);
                break;
            case KW_DG_DZ:
                status = kwfunc_dg_dz(tok, ntok, model);
                break;
            case KW_DP:
                status = kwfunc_dP(
                        tok, ntok, model, simplex, model->nlayers - 1);
                break;
            case KW_DSB_UTOL_RATIO:
                status = kwfunc_dsb_utol_ratio(tok, ntok, model, simplex);
                break;
            case KW_F:
            case KW_FREQ:
            case KW_FREQUENCY:
                status = kwfunc_frequency(tok, ntok, model);
                break;
            case KW_FIF:
                status = kwfunc_fif(tok, ntok, model);
                break;
            case KW_FIT:
                status = kwfunc_fit(tok, ntok, fit_data);
                break;
            case KW_FITS:
                status = kwfunc_fits(tok, ntok, fit_data);
                break;
            case KW_FIT_DATA_COLS:
            case KW_FIT_DATA_COLUMNS:
                status = kwfunc_fit_data_columns(tok, ntok, fit_data);
                break;
            case KW_FIT_DATA_DELIMITERS:
                status = kwfunc_fit_data_delimiters(tok, ntok);
                break;
            case KW_FIT_DATA_FORMAT:
                status = kwfunc_fit_data_format(tok, ntok);
                break;
            case KW_FIT_DATA_UNITS:
                status = kwfunc_fit_data_units(tok, ntok, fit_data);
                break;
            case KW_FIT_ESTIMATOR:
                status = kwfunc_fit_estimator(tok, ntok, fit_data);
                break;
            case KW_FIT_ITER:
                status = kwfunc_fit_iter(tok, ntok, fit_data);
                break;
            case KW_FIT_OUTPUT:
                status = kwfunc_fit_output(tok, ntok, fit_data);
                break;
            case KW_FIT_REINIT:
                status = kwfunc_fit_reinit(tok, ntok, simplex);
                break;
            case KW_FIT_TOL:
                status = kwfunc_fit_tol(tok, ntok, fit_data, simplex);
                break;
            case KW_FIT_TRACK_RESIDUALS:
                status = kwfunc_fit_track_residuals(tok, ntok, fit_data);
                break;
            case KW_FIT_VERBOSE:
                status = kwfunc_fit_verbose(tok, ntok, fit_data);
                break;
            case KW_FOUT:
                status = kwfunc_fout(tok, ntok, model);
                break;
            case KW_G:
                status = kwfunc_g(tok, ntok, model);
                break;
            case KW_GEOMETRY:
                status = kwfunc_geometry(tok, ntok, model);
                break;
            case KW_H:
                status = kwfunc_h(
                        tok, ntok, model, simplex, model->nlayers - 1);
                break;
            case KW_HEADERS:
                status = kwfunc_headers(tok, ntok);
                break;
            case KW_IFSPEC:
                status = kwfunc_ifspec(tok, ntok, model, simplex);
                break;
            case KW_ILS:
                status = kwfunc_ils(tok, ntok, model, simplex);
                break;
            case KW_ILSMODE:
                status = kwfunc_ilsmode(tok, ntok, model, simplex);
                break;
            case KW_JACOBIAN:
                status = kwfunc_jacobian(tok, ntok);
                break;
            case KW_KCACHE:
                status = kwfunc_kcache(tok, ntok, model);
                break;
            case KW_LAYER:
                status = kwfunc_layer(tok, ntok, model);
                break;
            case KW_LINESHAPE:
                status = kwfunc_lineshape(tok, ntok, model);
                break;
            case KW_MAIR:
                status = kwfunc_Mair(tok, ntok, model, model->nlayers - 1);
                break;
            case KW_NSCALE:
                status = kwfunc_Nscale(tok, ntok, simplex);
                break;
            case KW_OUTPUT:
                status = kwfunc_output(tok, ntok);
                break;
            case KW_P:
                status = kwfunc_P(
                        tok, ntok, model, simplex, model->nlayers - 1);
                break;
            case KW_PBASE:
                status = kwfunc_Pbase(tok, ntok, model, model->nlayers - 1);
                break;
            case KW_POBS:
                status = kwfunc_Pobs(tok, ntok, model, simplex);
                break;
            case KW_PSOURCE:
                status = kwfunc_Psource(tok, ntok, model, simplex);
                break;
            case KW_PTAN:
                status = kwfunc_Ptan(tok, ntok, model, simplex);
                break;
            case KW_PTMODE:
                status = kwfunc_PTmode(tok, ntok, model);
                break;
            case KW_R0:
                status = kwfunc_R0(tok, ntok, model);
                break;
            case KW_REFRACT:
                status = kwfunc_refract(tok, ntok, model);
                break;
            case KW_REVERSE:
                status = kwfunc_reverse(tok, ntok, model);
                break;
            case KW_RH_OFFSET:
                status = kwfunc_RH_offset(tok, ntok, model, simplex);
                break;
            case KW_RH_SCALE:
                status = kwfunc_RH_scale(tok, ntok, model, simplex);
                break;
            case KW_RUNTIME:
                status = kwfunc_runtime(tok, ntok, model);
                break;
            case KW_RX_GAIN_FACTOR:
                status = kwfunc_rx_gain_factor(tok, ntok, model, simplex);
                break;
            case KW_SEC_ZA:
                status = kwfunc_sec_za(tok, ntok, model, simplex);
                break;
            case KW_SELFBROAD_VMR_TOL:
                status = kwfunc_selfbroad_vmr_tol(tok, ntok, model);
                break;
            case KW_SIMPLEX_LOG:
                status = kwfunc_simplex_log(tok, ntok, simplex);
                break;
            case KW_T:
                status = kwfunc_T(
                        tok, ntok, model, simplex, model->nlayers - 1);
                break;
            case KW_T0:
                status = kwfunc_T0(tok, ntok, model, simplex);
                break;
            case KW_TBASE:
                status =
                    kwfunc_Tbase(tok, ntok, model, simplex, model->nlayers - 1);
                break;
            case KW_TREF:
                status = kwfunc_Tref(tok, ntok, model, simplex);
                break;
            case KW_TRX:
                status = kwfunc_Trx(tok, ntok, model, simplex);
                break;
            case KW_TOL:
            case KW_TOLERANCE:
                status = kwfunc_tol(tok, ntok, model);
                break;
            case KW_Z0:
                status = kwfunc_z0(tok, ntok, model);
                break;
            case KW_ZA:
                status = kwfunc_za(tok, ntok, model, simplex);
                break;
            case KW_ZOBS:
                status = kwfunc_zobs(tok, ntok, model, simplex);
                break;
            case KW_ZSOURCE:
                status = kwfunc_zsource(tok, ntok, model, simplex);
                break;
            case KW_ZTAN:
                status = kwfunc_ztan(tok, ntok, model, simplex);
                break;
            case KW_END_OF_TABLE:
            default:
                status = KWFUNC_SYNTAX_ERROR_CONTINUE;
                parse_error("Unrecognized keyword \"%s\"\n", tok[0]);
                break;
        }
        /*
         * If there was a serious error (such as a memory allocation
         * problem), return error status immediately.
         */
        if (status == KWFUNC_FATAL_ERROR) {
            fclose(cfgfile);
            return 1;
        }
        /*
         * For syntax errors likely to generate lots of cascading
         * erors, break out of the line parsing loop immediately.
         */
        if (status == KWFUNC_SYNTAX_ERROR_STOP) {
            syntax_errflag = 1;
            break;
        }
        /*
         * For other syntax errors, set a flag and keep going.
         * There might be more errors worth reporting.
         */
        if (status == KWFUNC_SYNTAX_ERROR_CONTINUE)
            syntax_errflag = 1;
    }
    /*
     * If there were any missing substitution parameters, rewind to
     * the top of the file and print out any comment lines that start
     * with the help character.
     */
    if (param_errflag && (cfgfile != stdin)) {
        rewind(cfgfile);
        while (fgets(buf, (int)sizeof(buf), cfgfile) != NULL) {
            if (buf[0] == HELP_CHAR)
                fprintf(stderr, "%s", buf + 1);
        }
    }
    fclose(cfgfile);
    /*
     * If a syntax error occurred in the line processing loop,
     * return with an error.
     */
    if (syntax_errflag)
        return 1;
    /*
     * Otherwise, continue on to further checks.
     */
    if (model->nlayers > 0) {
        if (model->layer[model->nlayers - 1]->P_unitnum == UNIT_NONE) {
            parse_error(
                    "Missing or invalid pressure definition on last layer.\n");
            syntax_errflag = 1;
        }
        if (model->layer[model->nlayers - 1]->T_unitnum == UNIT_NONE) {
            parse_error(
                    "Missing or invalid temperature definition"
                    " on last layer.\n");
            syntax_errflag = 1;
        }
    }
    /*
     * Check for parameter definitions or geometry modes that are
     * not valid in non-hydrostatic models.
     */
    if (!(model->PTmode & PTMODE_HYDROSTATIC)) {
        if (model->geometry & GEOMETRY_OBS_LEVEL_USER_DEFINED) {
            parse_error(
                    "User-specified observing level (Pobs or zobs) is not\n");
            parse_error("\tapplicable to non-hydrostatic models.\n");
            syntax_errflag = 1;
        }
        if (model->geometry & GEOMETRY_SOURCE_LEVEL_USER_DEFINED) {
            parse_error(
                    "User-specified source level (Psource or zsource) is\n");
            parse_error("\tnot applicable to non-hydrostatic models.\n");
            syntax_errflag = 1;
        }
        if (model->geometry & GEOMETRY_TAN_LEVEL_USER_DEFINED) {
            parse_error(
                    "User-specified tangent level (Ptan or ztan) is not\n");
            parse_error("\tapplicable to non-hydrostatic models.\n");
            syntax_errflag = 1;
        }
        if (model->geometry & GEOMETRY_R0_USER_DEFINED) {
            parse_error("R0 is not applicable to non-hydrostatic models.\n");
            syntax_errflag = 1;
        }
        if (model->geometry & GEOMETRY_Z0_USER_DEFINED) {
            parse_error("z0 is not applicable to non-hydrostatic models.\n");
            syntax_errflag = 1;
        }
        if ((model->geometry & GEOMETRY_SPHERICAL) ||
                (model->geometry & GEOMETRY_LIMB)) {
            parse_error(
                    "Spherical and limb geometry modes can only be used in\n");
            parse_error("\thydrostatic models.\n");
            syntax_errflag = 1;
        }
    }
    /*
     * Geometry mode consistency checks.
     */
    if (model->geometry & GEOMETRY_TAN_LEVEL_USER_DEFINED) {
        if (model->geometry & GEOMETRY_PLANE_PARALLEL) {
            parse_error(
                    "User-specified tangent level (Ptan or ztan) is\n");
            parse_error("\tnot applicable to plane-parallel geometry.\n");
            syntax_errflag = 1;
        } 
        if (model->geometry & GEOMETRY_SPHERICAL) {
            parse_error(
                    "In spherical mode, the tangent level is derived from\n");
            parse_error(
                    "\tthe observing level.  To specify a path by tangent\n");
            parse_error(
                    "\tlevel, use limb mode.\n");
            syntax_errflag = 1;
        }
        if (model->geometry &
                (GEOMETRY_ZA_USER_DEFINED | GEOMETRY_SEC_ZA_USER_DEFINED)) {
            parse_error(
                    "User-specified zenith angle is not applicable in limb\n");
            parse_error(
                    "\tmode.  The path geometry is set by specifying the\n");
            parse_error("\ttangent level.\n");
            syntax_errflag = 1;
        }
    }
    if (model->geometry &
            (GEOMETRY_OBS_LEVEL_USER_DEFINED |
             GEOMETRY_SOURCE_LEVEL_USER_DEFINED)) {
        if (model->geometry & GEOMETRY_LIMB) {
            parse_error(
                    "Observing level or source level cannot be defined in\n");
            parse_error("\tlimb mode");
            syntax_errflag = 1;
        }
    }
    /*
     * Return error status if there were any syntax errors found
     * after the line parsing loop.
     */
    if (syntax_errflag)
        return 1;
    /*
     * Finish up with memory allocation and ILS initialization,
     * unless in atmospheric model only mode.
     */
    if (!(output[ALL_OUTPUTS].flags & OUTPUT_AM_ONLY)) {
        if (alloc_model_arrays(model))
            return 1;
        if (alloc_jacobians(model, simplex))
            return 1;
        init_freq_grid(model);
        if (model->ils_typenum != ILS_NONE) {
            if (model->npad == 0) {
                errlog(62, 0);
                return 1;
            }
            if (initialize_ils(model))
                return 1;
        }
        if (set_active_outputs(ALL_OUTPUTS))
            return 1;
    }
    /*
     * Provide hints to the user regarding questionable model
     * configurations.
     */
    post_config_hints(model);
    return 0;
}   /* parse_config_file() */


/***********************************************************
* int set_config_parameter(
*         char *cmdline,
*         char *filename,
*         int linenum,
*         model_t *model,
*         fit_data_t *fit_data,
*         simplex_t *simplex)
*
* Purpose:
*   Processes a single configuration-file-style line to set
*   a configuration parameter.  This facility is used in
*   multi-dataset fits to change certain configuration
*   parameters on the fly from lines embedded in the
*   data stream.
*
* Arguments:
*   char       *cmdline  - input line buffer
*   char       *filename - file name, for error reporting
*   int        linenum   - line number, for error reporting
*   model_t    *model    - pointer to model structure
*   fit_data_t *fit_data - pointer to fit data structure
*   simplex_t  *simplex  - pointer to fit simplex structure
*
* Return:
*   0 on success, 1 otherwise
************************************************************/

int set_config_parameter(
        char *cmdline,
        char *filename,
        int linenum,
        model_t *model,
        fit_data_t *fit_data,
        simplex_t *simplex)
{
    char *tok[MAX_NTOK];
    int i, lnum, n_def_layers;
    int ntok;
    int status;

    CONFIG_FILE_NAME    = filename;
    CONFIG_FILE_LINENUM = linenum;
    ntok = parse_config_line(cmdline, 0, NULL, tok);
    if (ntok <= 0) {
        parse_error("No parameter specified.\n");
        return 1;
    } else if (ntok == MAX_NTOK) {
        parse_error("Too many tokens\n");
        return 1;
    }
    switch (get_keyword_num(tok[0])) {
        /*
         * Allowed keywords correspond to parameters which can
         * be changed without changing the model or simplex
         * dimensions.  For model variables, NULL is passed to
         * the corresponding keyword function for the simplex
         * pointer, to protect against simplex dimensionality
         * changes.
         */
        case KW_DG_DZ:
            status = kwfunc_dg_dz(tok, ntok, model);
            break;
        case KW_DSB_UTOL_RATIO:
            status = kwfunc_dsb_utol_ratio(tok, ntok, model, NULL);
            break;
        case KW_FIF:
            status = kwfunc_fif(tok, ntok, model);
            break;
        case KW_FIT_DATA_COLS:
        case KW_FIT_DATA_COLUMNS:
            status = kwfunc_fit_data_columns(tok, ntok, fit_data);
            break;
        case KW_FIT_DATA_DELIMITERS:
            status = kwfunc_fit_data_delimiters(tok, ntok);
            break;
        case KW_FIT_DATA_FORMAT:
            status = kwfunc_fit_data_format(tok, ntok);
            break;
        case KW_FIT_DATA_UNITS:
            status = kwfunc_fit_data_units(tok, ntok, fit_data);
            break;
        case KW_FIT_ESTIMATOR:
            status = kwfunc_fit_estimator(tok, ntok, fit_data);
            break;
        case KW_FIT_ITER:
            status = kwfunc_fit_iter(tok, ntok, fit_data);
            break;
        case KW_FIT_OUTPUT:
            status = kwfunc_fit_output(tok, ntok, fit_data);
            break;
        case KW_FIT_REINIT:
            status = kwfunc_fit_reinit(tok, ntok, simplex);
            break;
        case KW_FIT_TOL:
            status = kwfunc_fit_tol(tok, ntok, fit_data, simplex);
            break;
        case KW_FIT_TRACK_RESIDUALS:
            status = kwfunc_fit_track_residuals(tok, ntok, fit_data);
            break;
        case KW_FIT_VERBOSE:
            status = kwfunc_fit_verbose(tok, ntok, fit_data);
            break;
        case KW_FOUT:
            status = kwfunc_fout(tok, ntok, model);
            break;
        case KW_G:
            status = kwfunc_g(tok, ntok, model);
            break;
        case KW_HEADERS:
            status = kwfunc_headers(tok, ntok);
            break;
        case KW_IFSPEC:
            status = kwfunc_ifspec(tok, ntok, model, NULL);
            break;
        case KW_ILS:
            status = kwfunc_ils(tok, ntok, model, NULL);
            break;
        case KW_ILSMODE:
            status = kwfunc_ilsmode(tok, ntok, model, NULL);
            break;
        case KW_LAYER:
            if (ntok < 2) {
                parse_error("Nothing found after %s keyword.\n",
                        kw_table[KW_LAYER].name);
                status = KWFUNC_SYNTAX_ERROR_STOP;
                break;
            }
            if (get_nonneg_int_val(&lnum, tok[1])) {
                parse_error("Invalid layer number %s\n", tok[1]);
                status = KWFUNC_SYNTAX_ERROR_STOP;
                break;
            }
            /*
             * "set layer" only applies to the originally-defined
             * model layers, not interpolated layers.  lnum is
             * checked against the number of model definition
             * layers, then adjusted to point to the correct
             * layer.
             */
            n_def_layers = 0;
            for (i = 0; i < model->nlayers; ++i) {
                if (model->layer[i]->type == LAYER_TYPE_DEF)
                    ++n_def_layers;
            }
            if (lnum < 0 || lnum >= n_def_layers) {
                parse_error(
                        "The specified layer number is outside the range"
                        " [0 .. nlayers-1]\n");
                status = KWFUNC_SYNTAX_ERROR_STOP;
                break;
            }
            i = 0;
            while (i < lnum) {
                if (model->layer[i++]->type != LAYER_TYPE_DEF)
                    ++lnum;
            }
            if (ntok < 3) {
                parse_error("Empty layer configuration statement.\n");
                status = KWFUNC_SYNTAX_ERROR_STOP;
                break;
            }
            switch (get_keyword_num(tok[2])) {
                case KW_H:
                    status = kwfunc_h(tok + 2, ntok - 2, model, NULL, lnum);
                    break;
                case KW_DP:
                    status = kwfunc_dP(tok + 2, ntok - 2, model, NULL, lnum);
                    break;
                case KW_MAIR:
                    status = kwfunc_Mair(tok + 2, ntok - 2, model, lnum);
                    break;
                case KW_P:
                    status = kwfunc_P(tok + 2, ntok - 2, model, NULL, lnum);
                    break;
                case KW_PBASE:
                    status = kwfunc_Pbase(tok + 2, ntok - 2, model, lnum);
                    break;
                case KW_T:
                    status = kwfunc_T(tok + 2, ntok - 2, model, NULL, lnum);
                    break;
                case KW_TBASE:
                    status = kwfunc_Tbase(tok + 2, ntok - 2, model, NULL, lnum);
                    break;
                case KW_END_OF_TABLE:
                    status = KWFUNC_SYNTAX_ERROR_STOP;
                    parse_error("Unrecognized keyword \"%s\"\n", tok[2]);
                    break;
                default:
                    status = KWFUNC_SYNTAX_ERROR_STOP;
                    parse_error(
                            "The keyword \"%s\" may not be used here\n",
                            tok[2]);
                    break;
            }
            break;
        case KW_NSCALE:
            status = kwfunc_Nscale(tok, ntok, NULL);
            break;
        case KW_POBS:
            status = kwfunc_Pobs(tok, ntok, model, NULL);
            break;
        case KW_PSOURCE:
            status = kwfunc_Psource(tok, ntok, model, NULL);
            break;
        case KW_PTAN:
            status = kwfunc_Ptan(tok, ntok, model, NULL);
            break;
        case KW_R0:
            status = kwfunc_R0(tok, ntok, model);
            break;
        case KW_REVERSE:
            status = kwfunc_reverse(tok, ntok, model);
            break;
        case KW_RH_OFFSET:
            status = kwfunc_RH_offset(tok, ntok, model, NULL);
            break;
        case KW_RH_SCALE:
            status = kwfunc_RH_scale(tok, ntok, model, NULL);
            break;
        case KW_RUNTIME:
            status = kwfunc_runtime(tok, ntok, model);
            break;
        case KW_RX_GAIN_FACTOR:
            status = kwfunc_rx_gain_factor(tok, ntok, model, NULL);
            break;
        case KW_SEC_ZA:
            status = kwfunc_sec_za(tok, ntok, model, NULL);
            break;
        case KW_T0:
            status = kwfunc_T0(tok, ntok, model, NULL);
            break;
        case KW_TREF:
            status = kwfunc_Tref(tok, ntok, model, NULL);
            break;
        case KW_TRX:
            status = kwfunc_Trx(tok, ntok, model, NULL);
            break;
        case KW_Z0:
            status = kwfunc_z0(tok, ntok, model);
            break;
        case KW_ZA:
            status = kwfunc_za(tok, ntok, model, NULL);
            break;
        case KW_ZOBS:
            status = kwfunc_zobs(tok, ntok, model, NULL);
            break;
        case KW_ZSOURCE:
            status = kwfunc_zsource(tok, ntok, model, NULL);
            break;
        case KW_ZTAN:
            status = kwfunc_ztan(tok, ntok, model, NULL);
            break;
            /*
             * Keywords which should have been preceded by "layer
             * (layer number)"
             */
        case KW_H:
        case KW_DP:
        case KW_MAIR:
        case KW_P:
        case KW_PBASE:
        case KW_T:
        case KW_TBASE:
            status = KWFUNC_SYNTAX_ERROR_STOP;
            parse_error(
                    "\"%s\" must be preceded by \"layer (layer number)\".\n",
                    tok[0]);
            break;
            /*
             * Disallowed keywords.
             */
        case KW_COL:
        case KW_COLUMN:
        case KW_F:
        case KW_FREQ:
        case KW_FREQUENCY:
        case KW_FIT:
        case KW_FITS:
        case KW_GEOMETRY:
        case KW_JACOBIAN:
        case KW_KCACHE:
        case KW_LINESHAPE:
        case KW_OUTPUT:
        case KW_PTMODE:
        case KW_REFRACT:
        case KW_SELFBROAD_VMR_TOL:
        case KW_SIMPLEX_LOG:
        case KW_TOL:
        case KW_TOLERANCE:
            status = KWFUNC_SYNTAX_ERROR_STOP;
            parse_error(
                    "\"%s\" may only be used for initial model"
                    " configuration.\n", tok[0]);
            break;
        case KW_END_OF_TABLE:
        default:
            status = KWFUNC_SYNTAX_ERROR_STOP;
            parse_error("Unrecognized keyword \"%s\"\n", tok[0]);
            break;
    }
    if (status == KWFUNC_FATAL_ERROR ||
            status == KWFUNC_SYNTAX_ERROR_CONTINUE ||
            status == KWFUNC_SYNTAX_ERROR_STOP)
        return 1;
    return 0;
} /* set_config_parameter() */


/***********************************************************
* static int ci_strcmp(const char *cs, const char *ct)
*
* Purpose:
*   Does a case-insensitive string comparison.
*
* Arguments:
*   char *cs, *ct - strings to compare
*
* Return:
*   -1 if cs <  ct
*    0 if cs == ct
*    1 if cs >  ct
************************************************************/

static int ci_strcmp(const char *cs, const char *ct)
{
    int diff;

    /*
     * Scan until different characters are encountered.  Note that
     * tolower() takes an int.  Also, for consistency between
     * machines, force the char to unsigned.
     */
    while((diff = ( tolower((int)*(unsigned char*)cs) -
                    tolower((int)*(unsigned char*)ct) )) == 0) {
        /*
         * diff == 0, so hitting the null terminator at the end
         * of cs means the same for ct.  The strings are equal.
         */
        if (*cs == '\0')
            return 0;
        ++cs;
        ++ct;
    }
    /* diff != 0, sign determines lexical order */
    return ((diff < 0) ? -1 : 1);
}   /* ci_strcmp() */


/***********************************************************
* static int get_bool_val(int *val, const char *numstr)
*
* Purpose:
*   Converts a string to an int, constrained to be 0 or 1.
*
* Arguments:
*   int        *val    - address of converted value
*   const char *numstr - numeric string to be converted
*
* Return:
*   0 if OK
*   1 on error
************************************************************/

static int get_bool_val(int *val, const char *numstr)
{
    char *endp;

    *val = (int)strtol(numstr, &endp, 10);
    if (((*val != 0) && (*val != 1)) || (endp[0] != '\0')) {
        parse_error("Expected either 0 or 1 (%s?)\n", numstr);
        return 1;
    }
    return 0;
}   /* get_bool_val() */


/***********************************************************
* static int get_col_typenum(const int *name)
*
* Purpose:
*   Gets the column number associated with a column name.
*
* Arguments:
*   const char *name - column type name
*
* Return:
*   column number if found, otherwise COL_TYPE_NONE
************************************************************/

static int get_col_typenum(const char *name)
{
    int i;

    for (i = 0; i < COL_TYPE_END_OF_TABLE; ++i) {
        if (!ci_strcmp(name, col_type[i].name))
            break;
    }
    if (i == COL_TYPE_END_OF_TABLE)
        return COL_TYPE_NONE;
    else
        return i;
}   /* get_col_typenum() */


/***********************************************************
* static int get_dbl_val(
*         double     *val,
*         const char *numstr,
*         const char *unitstr,
*         const int  ugroup,
*         const int  diff_flag)
*
* Purpose:
*   Converts a string to a double, and optionally converts
*   the value thus obtained to am native units.  The source
*   unit must be in the unit group ugroup.  If diff_flag is
*   non-zero, *val is a differential quantity, and no offset
*   (as for degrees C) is applied in the unit conversion.
*
*   If the string is not successfully converted, or if the
*   source unit is not in the unit group ugroup, *val is
*   left unchanged, and the return value is set to 1.
*
* Arguments:
*   double     *val     - address of converted value
*   const char *numstr  - numeric string to be converted
*   const char *unitstr -
*       name of user unit (NULL if val is dimensionless)
*   const int  ugroup   - target unit group
*   const int diff_flag - flag for differential quantity
*
* Return:
*   0 if OK
*   1 on error
************************************************************/

static int get_dbl_val(
        double     *val,
        const char *numstr,
        const char *unitstr,
        const int  ugroup,
        const int  diff_flag)
{
    double x;
    char   *endp;
    int    unitnum;

    x = strtod(numstr, &endp);
    if (!(x == x)) {
        parse_error("NaN encountered in input.\n");
        return 1;
    }
    if (endp[0] != '\0') {
        parse_error("Number conversion failed. (%s?)\n", numstr);
        return 1;
    }
    if (unitstr != NULL) {
        if ((unitnum = get_unitnum(unitstr)) >= UNIT_END_OF_TABLE) {
            parse_error("Unrecognized unit (%s?)\n", unitstr);
            return 1;
        }
        if (unit_tab[unitnum].group != ugroup) {
            parse_error(
                    "Inappropriate unit (%s). Valid choices here are:\n",
                    unitstr);
            list_allowed_units(ugroup);
            return 1;
        }
        if (convert_to_am_unit(&x, unitnum, diff_flag)) {
            parse_error("Unit conversion failed. (%s?)\n", unitstr);
            return 1;
        }
    }
    *val = x;
    return 0;
}   /* get_dbl_val() */


/***********************************************************
* static int get_keyword_num(const char *kwname)
*
* Purpose:
*   Gets the index number associated with a configuration
*   file keyword name.  The comparison is case-insensitive.
*
* Arguments:
*   const char *kwname - pointer to string holding keyword
*
* Return:
*   keyword index number, as an int, or KW_END_OF_TABLE if
*   the keyword is not found.
************************************************************/

static int get_keyword_num(const char *kwname)
{
    int kwnum;

    for (kwnum = 0; kwnum < KW_END_OF_TABLE; ++kwnum) {
        if (ci_strcmp(kwname, kw_table[kwnum].name) == 0)
            break;
    }
    return kwnum;
}   /* get_keyword_num() */


/***********************************************************
* static int get_nonneg_dbl_val(
*         double *val,
*         const char *numstr,
*         const char *unitstr,
*         const int  ugroup,
*         const int  diff_flag)
*
* Purpose:
*   Calls get_dbl_val() with an additional check that the
*   converted value is greater than or equal to 0.0.
*
*   If the value read is not nonnegative, then *val is left
*   unchanged, and the return value is set to 1.
*
* Arguments:
*   double     *val     - address of converted value
*   const char *numstr  - numeric string to be converted
*   const char *unitstr -
*       name of user unit (NULL if val is dimensionless)
*   const int  ugroup   - target unit group
*   const int diff_flag - flag for differential quantity
*
* Return:
*   0 if OK
*   1 on error
************************************************************/

static int get_nonneg_dbl_val(
        double *val,
        const char *numstr,
        const char *unitstr,
        const int  ugroup,
        const int  diff_flag)
{
    double x;

    if (get_dbl_val(
                &x,
                numstr,
                unitstr,
                ugroup,
                diff_flag)) {
        return 1;
    }
    if (x < 0.0) {
        parse_error("Non-negative value expected, found %s\n", numstr);
        return 1;
    }
    *val = x;
    return 0;
}   /* get_nonneg_dbl_val() */


/***********************************************************
* static int get_nonneg_int_val(int *val, const char *numstr)
*
* Purpose:
*   Converts a string to an int.  The resulting value must
*   be greater than or equal to 0.
*
*   If the value read is not positive, then *val is left
*   unchanged, and the return value is set to 1.
*
* Arguments:
*   int        *val    - address of converted value
*   const char *numstr - numeric string to be converted
*
* Return:
*   0 if OK
*   1 on error
************************************************************/

static int get_nonneg_int_val(int *val, const char *numstr)
{
    int n;
    char *endp;

    n = (int)strtol(numstr, &endp, 10);
    if (endp[0] != '\0') {
        parse_error("Integer conversion failed. (%s?)\n", numstr);
        return 1;
    }
    if (n < 0) {
        parse_error("Positive value expected, found %s\n", numstr);
        return 1;
    }
    *val = n;
    return 0;
}   /* get_nonneg_int_val() */


/***********************************************************
* static int get_pos_dbl_val(
*         double *val,
*         const char *numstr,
*         const char *unitstr,
*         const int  ugroup,
*         const int  diff_flag)
*
* Purpose:
*   Calls get_dbl_val() with an additional check that the
*   converted value is greater than 0.0.
*
*   If the value read is not positive, then *val is left
*   unchanged, and the return value is set to 1.
*
* Arguments:
*   double     *val     - address of converted value
*   const char *numstr  - numeric string to be converted
*   const char *unitstr -
*       name of user unit (NULL if val is dimensionless)
*   const int  ugroup   - target unit group
*   const int diff_flag - flag for differential quantity
*
* Return:
*   0 if OK
*   1 on error
************************************************************/

static int get_pos_dbl_val(
        double *val,
        const char *numstr,
        const char *unitstr,
        const int  ugroup,
        const int  diff_flag)
{
    double x;

    if (get_dbl_val(
                &x,
                numstr,
                unitstr,
                ugroup,
                diff_flag)) {
        return 1;
    }
    if (x <= 0.0) {
        parse_error("Positive value expected, found %s\n", numstr);
        return 1;
    }
    *val = x;
    return 0;
}   /* get_pos_dbl_val() */


/***********************************************************
* static int get_RH_value(double *val, const char *numstr)
*
* Purpose:
*   Converts a string to a relative humidity value in %.
*
*   The string is required to be in the form of a number
*   ending in %.  If it is not in this form, then *val
*   is left unchanged, and the return value is set to 1.
*
* Arguments:
*   double     *val    - address of converted value
*   const char *numstr - numeric string to be converted
*
* Return:
*   0 if OK
*   1 on error
************************************************************/

static int get_RH_value(double *val, const char *numstr)
{
    double x;
    char *endp;

    x = strtod(numstr, &endp);
    if (!(x == x)) {
        parse_error("NaN encountered in input.\n");
        return 1;
    }
    if (endp[0] == '\0') {
        parse_error("RH values must end in %%, with no preceding space.\n");
        parse_error("\tExample: column h2o RH 45%%\n");
        return 1;
    }
    if (strlen(endp) != 1) {
        parse_error("Extra characters in RH value. (%s?)\n", numstr);
        return 1;
    }
    if (endp[0] != '%') {
        parse_error("RH values must end in %%\n");
        parse_error("\tExample: column h2o RH 45%%\n");
        return 1;
    }
    *val = x;
    return 0;
}   /* get_RH_value() */


/***********************************************************
* static void init_freq_grid(model_t *model)
*
* Purpose:
*   initialize the frequency and squared frequency arrays.
*
* Arguments:
*   model_t *model - pointer to model structure
*
* Return:
*   0 if OK
*   1 on error
************************************************************/

static void init_freq_grid(model_t *model)
{
    double f0;
    gridsize_t i;

    /*
     * The frequency grid is aligned to the origin.  That is, it
     * consists of all points i * df, where i is an integer, on the
     * closed interval f_min <= f <= f_max.  This alignment is
     * required for computation of delay spectra by discrete Hilbert
     * transformation.  It also facilitates caching of absorption
     * coefficients to disk, since all grids with the same df are
     * aligned to one another.
     */
    f0 = model->df * ceil((model->fmin * (1.0 - DBL_EPSILON)) / model->df);
    for (i = 0; i < model->ngrid; ++i) {
        model->f[i]  = f0 + (double)i * model->df;
        model->f2[i] = model->f[i] * model->f[i];
    }
    return;
}   /* init_freq_grid() */


/***********************************************************
* static int is_numeric_token(const char *tok)
*
* Purpose:
*   Checks a token *tok to see if it looks like a number, by
*   attempting to convert it to a double.
*
* Arguments:
*   double *tok - token to be checked
*
* Return:
*   1 if tok is a number
*   0 otherwise
************************************************************/

static int is_numeric_token(const char *tok)
{
    double x;
    char *endp;

    if (tok == NULL || tok[0] == '\0')
        return 0;
    x = strtod(tok, &endp);
    if (!(x == x))
        return 0;   /* NaN is not a number */
    if (endp[0] != '\0')
        return 0;
    else
        return 1;
}   /* is_numeric_token() */


/***********************************************************
* static void list_allowed_keywords(
*         const char** kw, const int nkw)
*
* Purpose:
*   For use in error messages, prints a list of allowed
*   keyword strings from an array kw of nkw strings.
************************************************************/

static void list_allowed_keywords(
        const char** kw, const int nkw)
{
    int i, nchar = 0;

    for (i = 0; i < nkw; ++i) {
        if (nchar > 40) {
            fprintf(stderr, "\n");
            nchar = 0;
        }
        if (nchar == 0)
            parse_error("\t");
        else 
            nchar += fprintf(stderr, "  ");
        nchar += fprintf(stderr, "%s", kw[i]);
    }
    fprintf(stderr, "\n");
    return;
}   /* list_allowed_keywords() */


/***********************************************************
* static void list_allowed_output_array_names(void)
*
* Purpose:
*   For use in error messages, prints a list of allowed
*   output array names.
************************************************************/

static void list_allowed_output_array_names(void)
{
    int i, nchar = 0;

    for (i = 1; i < OUTPUT_END_OF_TABLE; ++i) {
        if (nchar > 40) {
            fprintf(stderr, "\n");
            nchar = 0;
        }
        if (nchar == 0)
            parse_error("\t");
        else 
            nchar += fprintf(stderr, "  ");
        nchar += fprintf(stderr, "%s", output[i].name);
    }
    fprintf(stderr, "\n");
    return;
}   /* list_allowed_output_array_names() */


/***********************************************************
* static void list_allowed_units(const int ugroup)
*
* Purpose:
*   For use in error messages, prints a list of allowed
*   units belonging to the unit group ugroup.
*
* Arguments:
*   const int ugroup - unit group number
************************************************************/

static void list_allowed_units(const int ugroup)
{
    int i, nchar = 0;

    for (i = 0; i < UNIT_END_OF_TABLE; ++i) {
        if (unit_tab[i].group == ugroup) {
            if (nchar > 40) {
                fprintf(stderr, "\n");
                nchar = 0;
            }
            if (nchar == 0)
                parse_error("\t");
            else 
                nchar += fprintf(stderr, "  ");
            nchar += fprintf(stderr, "%s", unit_tab[i].name);
        }
    }
    fprintf(stderr, "\n");
    return;
}   /* list_allowed_units() */


/***********************************************************
* static void list_line_by_line_abscoeffs(void)
*
* Purpose:
*   For use in error messages, prints a list of line-by-line
*   absorption coefficient types.
************************************************************/

static void list_line_by_line_abscoeffs(void)
{
    int i, nchar = 0;

    for (i = 0; i < K_TYPE_END_OF_TABLE; ++i) {
        if (k_type[i].dep_flags & DEP_ON_LSHAPE) {
            if (nchar > 40) {
                fprintf(stderr, "\n");
                nchar = 0;
            }
            if (nchar == 0)
                parse_error("\t");
            else 
                nchar += fprintf(stderr, "  ");
            nchar += fprintf(stderr, "%s", k_type[i].name);
        }
    }
    fprintf(stderr, "\n");
    return;
}   /* list_line_by_line_abscoeffs() */


/***********************************************************
* static int load_fit_file_list(
*         const char *file_list,
*         fit_data_t *fit_data)
*
* Purpose:
*   loads the list of fit data file names into the fit data
*   structure.
*
* Arguments:
*   char       *file_list -
*       string containing path to list of fit data files
*   fit_data_t *fit_data  - fit data structure
*
* Return:
*   0 on success, 1 otherwise
************************************************************/

static int load_fit_file_list(
        const char *file_list,
        fit_data_t *fit_data)
{
    char nbuf[AM_MAX_FNAMESIZE - AM_FIT_EXTSIZE];
    char **tptr;
    int i;
    FILE *fp;

    if ((fp = fopen(file_list, "r")) == NULL)
        return FILELIST_CANNOT_OPEN;
    while (fgets(nbuf, (int)sizeof(nbuf), fp) != NULL) {
        ++fit_data->nfiles;
        if ((tptr = (char**)realloc(fit_data->filename,
                        fit_data->nfiles * sizeof(char*))) == NULL) {
            errlog(33, 0);
            return FILELIST_ERROR;
        }
        fit_data->filename = tptr;
        if ((fit_data->filename[fit_data->nfiles - 1] =
                    (char*)malloc((1 + strlen(nbuf)) * sizeof(char)))
                == NULL) {
            errlog(33, 0);
            return FILELIST_ERROR;
        }
        /*
         * Chop any whitespace off the end of the input line (including
         * CR and LF).
         */
        i = (int)strlen(nbuf) - 1;
        while (isspace(nbuf[i]))
            nbuf[i--] = '\0';
        strcpy(fit_data->filename[fit_data->nfiles - 1], nbuf);
    }
    return FILELIST_OK;
}   /* load_fit_file_list() */


/***********************************************************
* static int parse_config_line(char *buf, int argc, char **argv, char **tok)
*
* Purpose:
*   Breaks a configuration line into tokens, according to
*   the following rules:
*   - Quoted strings become a single token, sans quotes.
*   - Tokens beginning with a comment character are ignored,
*     along with the rest of the line.
*   - Tokens such as %1, %2,... are replaceable parameters,
*     and are replaced with argv[i], argv[i+1],... where
*     argv[i] is the first command line argument remaining
*     after processing of command line switches and the
*     config file name.
*   - No replacements are made in quoted strings, such as
*     format strings and quoted filenames.
*
* Arguments:
*   char *buf   - line buffer
*   int  argc   - command line argument count
*   char **argv - command line argument list
*   char **tok  - array to receive tokens
*
* Return:
*   Number of tokens found, or MAX_NTOK if there are too
*   many tokens.
*
*   Negative values are returned for errors involving
*   replaceable paramters.  These are:
*   SUBST_PARAMETER_MISSING      - A replaceable parameter
*           has no corresponding command-line argument.
*   SUBST_PARAMETER_SYNTAX_ERROR - Encountered a token
*           starting with '%' that was not followed by an
*           integer index number.
************************************************************/

static int parse_config_line(char *buf, int argc, char **argv, char **tok)
{
    char *bufptr;
    enum {QUOTE, SPACE} delim;
    int  ntok;

    bufptr = buf;
    delim  = SPACE;
    for (ntok = 0; ntok < MAX_NTOK; ++ntok) {
        while (isspace(*bufptr)) {
            delim = SPACE;
            ++bufptr;
        }
        while (*bufptr == '\"') {
            delim = QUOTE;
            ++bufptr;
        }
        if ((*bufptr == '\0') || (*bufptr == '\n') || (*bufptr == '\r')) {
            break;
        } else {
            tok[ntok] = bufptr++;
        }
        if (delim == SPACE) {
            while (!isspace(*bufptr) && (*bufptr != '\0'))
                ++bufptr;
        } else if (delim == QUOTE) {
            while ((*bufptr != '\"') && (*bufptr != '\0')
                    && (*bufptr != '\n') && (*bufptr != '\r'))
                ++bufptr;
        }
        if (*bufptr != '\0')
            *bufptr++ = '\0';
        if (strchr(COMMENTSTR, tok[ntok][0]) != NULL) {
            break;
        } else if ((delim == SPACE) && tok[ntok][0] == '%') {
            int i;
            if (tok[ntok][1] == '\0') {
                parse_error(
                        "Missing replaceable parameter index after '%%'\n");
                ntok = SUBST_PARAMETER_SYNTAX_ERROR;
                break;
            }
            if (get_nonneg_int_val(&i, tok[ntok] + 1)) {
                parse_error("Error in replaceable parameter format.\n");
                ntok = SUBST_PARAMETER_SYNTAX_ERROR;
                break;
            }
            ++i;
            if (i <= 1 || i >= argc) {
                parse_error(
                        "The replaceable parameter %s has no matching"
                        " command-line\n", tok[ntok]);
                parse_error("\tparameter.\n");
                ntok = SUBST_PARAMETER_MISSING;
                break;
            } else {
                tok[ntok] = argv[i];
            }
        }
    }
    return ntok;
}   /* parse_config_line() */


/***********************************************************
* static void parse_error(const char* format, ...)
*
* Purpose:
*   Prints a configuration file parse error message,
*   identified by file name and line number.  If format
*   begins with a tab character, a consistently-indented
*   continuation line is printed without the file name and
*   line number.  The format string should include any
*   needed newline character.
*
*   Output (to stderr) is of the form
*
*   filename(linenum) : error message
*                       continuation line
*
*   or 
*
*   longer_filename(linenum) :
*       error message
*       continuation line
*
*   The parse error count in the global errlog is
*   incremented when a non-continuation line is printed.
*
* Arguments:
*   char *format, ... - printf-style format string and
*       argument list
************************************************************/

static void parse_error(const char* format, ...)
{
    enum {
        LINE_ID_BUFSIZE = 256,
        MAX_INDENT      = 20,
        FALLBACK_INDENT = 4
    };
    /*
     * Here, buf[] must be large enough to hold any
     * string printed with LINE_ID_FMT[].
     */
    char buf[LINE_ID_BUFSIZE];
    const char LINE_ID_FMT[] = "%.200s(%d) : ";
    int indent;
    va_list args;

    va_start(args, format);
    indent = snprintf(buf, sizeof(buf),
            LINE_ID_FMT, CONFIG_FILE_NAME, CONFIG_FILE_LINENUM);
    if (format[0] == '\t') {
        ++format; /* strip tab from format */
        if (indent > MAX_INDENT)
            fprintf(stderr, "%*s", FALLBACK_INDENT, "");
        else
            fprintf(stderr, "%*s", indent, "");
    } else {
        errlog(20, 0); /* increment parse error count in errlog */
        if (indent > MAX_INDENT)
            fprintf(stderr, "%s\n%*s", buf, FALLBACK_INDENT, "");
        else
            fprintf(stderr, "%s", buf);
    }
    vfprintf(stderr, format, args);
    va_end(args);
    return;
}   /* parse_error() */


/***********************************************************
* static void post_config_hints(model_t *model)
*
* Purpose:
*   Checks for unusual model configurations that might not
*   have been what the user intended, and posts advisory
*   warnings to the error log.
*
* Arguments:
*   model_t *model - pointer to model structure
************************************************************/

static void post_config_hints(model_t *model)
{
    /*
     * If opacity is a requested output, and a DSB spectrum
     * is being computed, warn the user that this is only
     * reasonable in the optically-thin limit.
     */
    if ((model->ifmode & IFMODE_DSB) &&
            (output[OUTPUT_OPACITY].flags &
            (OUTPUT_USER | OUTPUT_FITTED | OUTPUT_JACOBIAN)))
        errlog(183, 0);
    /*
     * If no spectra were explicitly specified in an ILS statement,
     * warn if none or some of the spectra were ineligible for
     * ILS convolution.
     */
    if (model->ils_typenum != ILS_NONE &&
            !(output[ALL_OUTPUTS].flags & ILS_APPLIED)) {
        int i, j = 0;
        for (i = 0; outcol[i] != 0; ++i) {
            if (output[outcol[i]].flags & ILS_ALLOWED)
                ++j;
        }
        if (j == 0)
            errlog(184, 0);
        else if (j < i - 1)
            errlog(185, 0);
    }
    /*
     * Spherical or limb geometry with refraction turned
     * off by default.
     */
    if ((model->geometry & (GEOMETRY_SPHERICAL | GEOMETRY_LIMB)) &&
            !(model->geometry & GEOMETRY_REFRACT_USER_DEFINED) &&
            (model->geometry & GEOMETRY_REFRACT_NONE))
        errlog(174, 0);
    /*
     * R0 defined in a plane-parallel model
     */
    if ((model->geometry & GEOMETRY_PLANE_PARALLEL) &&
            (model->geometry & GEOMETRY_R0_USER_DEFINED))
        errlog(175, 0);
    /*
     * dg/dz was given a positive value.  Under normal
     * circumstances, it should be negative or zero.
     */
    if (model->dg_dz > 0.0)
        errlog(176, 0);
    /*
     * A plume or instrument layer was defined in Tbase
     * mode, and is interacting with the temperature of
     * adjacent nonempty layer(s).  The user might have
     * assume the plume or instrument layer temperature
     * would be decoupled from the surrounding atmosphere.
     * The fix is to use dummy layers.
     */
    if (model->PTmode & PTMODE_TBASE) {
        int i;
        for (i = 0; i < model->nlayers; ++i) {
            layer_t *layer = model->layer[i];
            if (layer->h > 0) {
                if (i > 0) {
                    layer_t *layer_above = model->layer[i-1];
                    if (layer_above->ncols &&
                            layer_above->Tbase != layer->Tbase)
                        errlog(207, i);
                }
                if (i < model->nlayers - 2) {
                    layer_t *layer_below = model->layer[i+1];
                    if (layer_below->ncols &&
                            layer_below->Tbase != layer->Tbase)
                        errlog(208, i);
                }
            }
        }
    }
    /*
     * A PTmode statement specified a non-isothermal base
     * extrapolation mode for a model with layers defined by
     * midpoint temperature.  For such models, extrapolation is
     * always isothermal.
     */
    if (model->PTmode & PTMODE_T &&
            !(model->PTmode & PTMODE_EXTEND_ISOTHERMAL))
        errlog(209, 0);
    return;
}   /* post_config_hints() */


/***********************************************************
* static int kwfunc_column(
*   char *tok[],
*   const int ntok,
*   model_t *model,
*   simplex_t *simplex)
*
* keywords:
*   col or column - set column type and column density
*
* syntax:
*   column column_type col_density unit [dcol_density unit]
*   column column_type vmr [vmr [dvmr]]
*   column column_type hydrostatic [vmr [dvmr]]
*   column h2o_column_type RH RH% [dRH%]
************************************************************/

static int kwfunc_column(
        char *tok[],
        const int ntok,
        model_t *model,
        simplex_t *simplex)
{
    int m, n;
    int col_typenum;
    double dN;
    if (model->nlayers == 0) {
        parse_error("Attempted to create column in non-existent layer\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (ntok < 3 || ntok > 6) {
        parse_error("Expected\n");
        parse_error(
                "\t  \"column column_type col_density unit"
                " [dcol_density unit]\"\n");
        parse_error("\tor\n");
        parse_error(
                "\t  \"column column_type (hydrostatic | vmr)"
                " [vmr [dvmr]]\"\n");
        parse_error("\t    where vmr = volume mixing ratio"
                " (none = use default)\n");
        parse_error("\tor\n");
        parse_error("\t  \"column column_type RH RH%% [dRH%%]\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if ((col_typenum = get_col_typenum(tok[1])) == COL_TYPE_NONE) {
        parse_error("Unrecognized column type \"%s\"\n", tok[1]);
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    /*
     * Add the column to the current layer
     */
    n = model->nlayers - 1;
    if (add_column(model->layer[n], col_typenum))
        return KWFUNC_FATAL_ERROR;
    m = model->layer[n]->ncols - 1;
    if (!ci_strcmp(tok[2], "hydrostatic") || 
            !ci_strcmp(tok[2], "vmr")) {
        /*
         * Column density specified by volume mixing ratio.
         */
        double vmr;
        if (col_type[col_typenum].flags & COL_PARAMETRIC) {
            parse_error(
                    "Computation of column density by mixing ratio is not\n");
            parse_error("\tapplicable to parametric column types.\n");
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
        if (model->PTmode & PTMODE_HYDROSTATIC) {
            /*
             * For hydrostatic models, "hydrostatic" and "vmr"
             * are both accepted on input for setting hydrostatic
             * column densities N, though "hydrostatic" is used
             * on output.  However, in embedded plume or
             * instrument layers, N is set by distance and the
             * keyword "vmr" will be used in output.
             */
            if (model->layer[n]->h >= 0.0)
                model->layer[n]->column[m]->N_mode = N_BY_DISTANCE;
            else
                model->layer[n]->column[m]->N_mode = N_HYDROSTATIC;
        } else if (!ci_strcmp(tok[2], "hydrostatic")) {
            /*
             * non-hydrostatic model, and tok[2] was "hydrostatic",
             * which is not allowed.
             */
            parse_error("Hydrostatic column density is not defined in ");
            write_PTmode(stderr, model->PTmode);
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        } else {
            /*
             * non-hydrostatic model, and tok[2] was "vmr"
             */
            model->layer[n]->column[m]->N_mode = N_BY_DISTANCE;
            if (model->layer[n]->h < 0.0) {
                parse_error(
                        "The layer thickness h has not been defined,"
                        " but is needed\n");
                parse_error(
                        "\tfor computation of the column density from"
                        " a mixing ratio\n");
                parse_error("\tin non-hydrostatic PTmodes.\n");
                return KWFUNC_SYNTAX_ERROR_CONTINUE;
            }
        }
        if (ntok == 3) {
            /*
             * No mixing ratio specified.  Check for a default.
             */
            if ((vmr = col_type[col_typenum].default_vmr) == 0.0) {
                parse_error(
                        "For %s, there is no default mixing ratio."
                        " Specify a mixing ratio\n",
                        col_type[col_typenum].name);
                parse_error("after \"hydrostatic\" or \"vmr\" keyword.\n");
                return KWFUNC_SYNTAX_ERROR_CONTINUE;
            }
        } else {
            if (get_dbl_val(
                        &vmr,
                        tok[3],
                        NULL,
                        UGROUP_NONE,
                        0)) {
                return KWFUNC_SYNTAX_ERROR_CONTINUE;
            }
            if (vmr < 0.0 || vmr > 1.0) {
                parse_error(
                        "The volume mixing ratio must be"
                        " in the range 0 to 1.\n");
                return KWFUNC_SYNTAX_ERROR_CONTINUE;
            }
            model->layer[n]->column[m]->vmr_stat |= VMR_USER_DEFINED;
        }
        /*
         * To facilitate fits, mixing ratios are mapped from the
         * interval [0, 1] to [0, 1 / FLT_EPSILON].  The mapping is
         * xvmr = vmr / (1 + FLT_EPSILON - vmr).
         */
        model->layer[n]->column[m]->xvmr = map_variable(vmr, MAPPING_VMR);
        if (ntok == 5) {
            double dvmr, dxvmr;
            if (get_pos_dbl_val(
                        &dvmr,
                        tok[4],
                        NULL,
                        UGROUP_NONE,
                        1)) {
                return KWFUNC_SYNTAX_ERROR_CONTINUE;
            }
            dxvmr = map_differential(vmr, dvmr, MAPPING_VMR);
            snprintf(VARNAME, sizeof(VARNAME), "%s vmr (layer %d)",
                    col_type[col_typenum].name, n);
            if (add_var_to_simplex(
                        simplex,
                        VARNAME,
                        &(model->layer[n]->column[m]->xvmr),
                        model->layer[n]->column[m]->xvmr,
                        dxvmr,
                        UNIT_NONE,
                        MAPPING_VMR)) {
                return KWFUNC_FATAL_ERROR;
            }
        }
    } else if (!ci_strcmp(tok[2], "RH") || !ci_strcmp(tok[2], "RHi")) {
        /*
         * Column density computed by relative humidity.
         */
        if (!(col_type[col_typenum].flags & COL_RH_ALLOWED)) {
            parse_error(
                    "RH (relative humidity) is not applicable to column\n");
            parse_error("\ttype %s\n", col_type[col_typenum].name);
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
        if (model->PTmode & PTMODE_HYDROSTATIC) {
            model->layer[n]->column[m]->N_mode = N_HYDROSTATIC;
        } else {
            model->layer[n]->column[m]->N_mode = N_BY_DISTANCE;
            if (model->layer[n]->h < 0.0) {
                parse_error(
                        "The layer thickness h has not been defined,"
                        " but is needed\n");
                parse_error(
                        "\tto compute column density from relative"
                        " humidity in\n");
                parse_error("\tnon-hydrostatic PTmodes.\n");
                return KWFUNC_SYNTAX_ERROR_CONTINUE;
            }
        }
        model->layer[n]->column[m]->vmr_stat |= VMR_USER_DEFINED;
        if (!ci_strcmp(tok[2], "RH"))
            model->layer[n]->column[m]->vmr_stat |= VMR_BY_RH_LIQUID;
        else
            model->layer[n]->column[m]->vmr_stat |= VMR_BY_RH_ICE;
        if (get_RH_value(&(model->layer[n]->column[m]->RH), tok[3]))
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        if (model->layer[n]->column[m]->RH < 0.0) {
            parse_error(
                    "Relative humidity must be greater than"
                    " or equal to 0%%.\n");
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
        if (ntok == 5) {
            double dRH;
            if (get_RH_value(&dRH, tok[4]))
                return KWFUNC_SYNTAX_ERROR_CONTINUE;
            if (dRH <= 0.0) {
                parse_error("dRH must be greater than 0%.\n");
                return KWFUNC_SYNTAX_ERROR_CONTINUE;
            }
            snprintf(VARNAME, sizeof(VARNAME),
                    "%s RH (layer %d)", col_type[col_typenum].name, n);
            if (add_var_to_simplex(
                        simplex,
                        VARNAME,
                        &(model->layer[n]->column[m]->RH),
                        model->layer[n]->column[m]->RH,
                        dRH,
                        UNIT_NONE,
                        MAPPING_NONE)) {
                return KWFUNC_FATAL_ERROR;
            }
        }
    } else {
        /*
         * Numerically-specified column density, or parameter for
         * parametric column types.
         */
        int ugroup = unit_tab[col_type[col_typenum].common_unit].group;
        model->layer[n]->column[m]->N_mode = N_EXPLICIT;
        if (ntok != 4 && ntok != 6) {
            parse_error(
                    "Expected \"column column_type col_density unit"
                    " [dcol_density unit]\"\n");
            parse_error(
                    "\t      or \"column column_type (hydrostatic | vmr)"
                    " [vmr [dvmr]]\"\n");
            parse_error(
                    "\t      where vmr = volume mixing ratio"
                    " (none = use default)\n");
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
        /*
         * If a model definition layer is subsequently split by an
         * interpolated level or levels, an explicit column density
         * will be split proportionately as well.  In a column struct,
         * N_def is the defining value, and N is the (possibly) split
         * value set upon each model computation.
         */
        if (get_nonneg_dbl_val(
                    &(model->layer[n]->column[m]->N_def),
                    tok[2],
                    tok[3],
                    ugroup,
                    0)) {
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
        model->layer[n]->column[m]->N_unitnum = get_unitnum(tok[3]);
        if (ntok == 6) {
            if (get_pos_dbl_val(
                        &dN,
                        tok[4],
                        tok[5],
                        ugroup,
                        1)) {
                return KWFUNC_SYNTAX_ERROR_CONTINUE;
            }
            snprintf(VARNAME, sizeof(VARNAME),
                    "N(%s) (layer %d)", col_type[col_typenum].name, n);
            if (add_var_to_simplex(
                        simplex,
                        VARNAME,
                        &(model->layer[n]->column[m]->N_def),
                        model->layer[n]->column[m]->N_def,
                        dN,
                        model->layer[n]->column[m]->N_unitnum,
                        MAPPING_NONE)) {
                return KWFUNC_FATAL_ERROR;
            }
        }
    }
    return KWFUNC_SUCCESS;
} /* kwfunc_column() */


/***********************************************************
* static int kwfunc_dg_dz(
*   char *tok[],
*   const int ntok,
*   model_t *model)
*
* keyword:
*   dg_dz - set vertical gradient of gravitational
*           acceleration
*
* syntax:
*   dg_dz vertical_gradient unit
************************************************************/

static int kwfunc_dg_dz(
        char *tok[],
        const int ntok,
        model_t *model)
{
    if (ntok < 3) {
        parse_error("Expected \"dg_dz vertical_gradient unit\"\n");
        parse_error("\texample: \"dg_dz -3.086e-6 s^-2\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_dbl_val(
                &(model->dg_dz),
                tok[1],
                tok[2],
                UGROUP_ACCEL_GRADIENT,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    model->dg_dz_unitnum = get_unitnum(tok[2]);
    return KWFUNC_SUCCESS;
} /* kwfunc_dg_dz() */


/***********************************************************
* static int kwfunc_dP(
*   char *tok[],
*   const int ntok,
*   model_t *model,
*   simplex_t *simplex,
*   int lnum)
*
* keyword:
*   dP - set pressure drop across layer
*
* syntax:
*   dP pressure unit [dpressure unit]
************************************************************/

static int kwfunc_dP(
        char *tok[],
        const int ntok,
        model_t *model,
        simplex_t *simplex,
        int lnum)
{
    if (lnum < 0 || lnum >= model->nlayers) {
        parse_error("Attempted to set dP for non-existent layer.\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if ((model->PTmode & PTMODE_PMODES) && !(model->PTmode & PTMODE_DP)) {
        parse_error(
                "dP is inconsistent with a prior layer or PTmode setting.\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (ntok != 3 && ntok != 5) {
        parse_error("Expected \"dP pressure unit [dpressure unit]\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    model->PTmode |= PTMODE_DP;
    /*
     * If a model definition layer is subsequently split by an
     * interpolated level or levels, the pressure drop across the
     * level will be split proportionately as well.  In a column
     * struct, dP_def is the defining value, and dP is the (possibly)
     * split value set upon each model computation.
     */
    if (get_nonneg_dbl_val(
                &(model->layer[lnum]->dP_def),
                tok[1],
                tok[2],
                UGROUP_PRESSURE,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    model->layer[lnum]->P_unitnum = get_unitnum(tok[2]);
    if (ntok == 5 && simplex != NULL) {
        double ddP;
        if (get_pos_dbl_val(
                    &ddP,
                    tok[3],
                    tok[4],
                    UGROUP_PRESSURE,
                    1)) {
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
        snprintf(VARNAME, sizeof(VARNAME), "dP (layer %d)", lnum);
        if (add_var_to_simplex(
                    simplex,
                    VARNAME,
                    &(model->layer[lnum]->dP_def),
                    model->layer[lnum]->dP_def,
                    ddP,
                    model->layer[lnum]->P_unitnum,
                    MAPPING_NONE)) {
            return KWFUNC_FATAL_ERROR;
        }
    }
    /*
     * In hydrostatic models, h can only be defined on layers with
     * zero pressure drop.  Since dP_def is allowed to be a fit
     * variable, check for that too.
     */
    if (model->layer[lnum]->h >= 0.0) {
        if (model->layer[lnum]->dP_def > 0.0 ||
                isvar(simplex, &(model->layer[lnum]->dP_def))) {
            parse_error(
                    "In hydrostatic models, layer thickness h can only be\n");
            parse_error(
                    "\tassigned if there is no pressure drop across the\n");
            parse_error(
                    "\tlayer.\n");
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
    }
    return KWFUNC_SUCCESS;
} /* kwfunc_dP() */


/***********************************************************
* static int kwfunc_dsb_utol_ratio(
*   char *tok[],
*   const int ntok,
*   model_t *model,
*   simplex_t *simplex)
*
* keyword:
*   dsb_utol_ratio - set USB/LSB gain ratio for DSB ILS or
*       DSB IF spectrum.
*
* syntax:
*   dsb_utol_ratio ratio [dratio]
************************************************************/

static int kwfunc_dsb_utol_ratio(
        char *tok[],
        const int ntok,
        model_t *model,
        simplex_t *simplex)
{
    double dratio;
    if (ntok != 2 && ntok != 3) {
        parse_error("Expected \"dsb_utol_ratio ratio [dratio]\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_nonneg_dbl_val(
                &(model->dsb_utol_ratio),
                tok[1],
                NULL,
                UGROUP_NONE,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (ntok == 3 && simplex != NULL) {
        if (get_pos_dbl_val(
                    &dratio,
                    tok[2],
                    NULL,
                    UGROUP_NONE,
                    1)) {
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
        if (add_var_to_simplex(
                    simplex,
                    "dsb_utol_ratio",
                    &(model->dsb_utol_ratio),
                    model->dsb_utol_ratio,
                    dratio,
                    UNIT_NONE,
                    MAPPING_NONE)) {
            return KWFUNC_FATAL_ERROR;
        }
    }
    return KWFUNC_SUCCESS;
} /* kwfunc_dsb_utol_ratio() */


/***********************************************************
* static int kwfunc_fif(
*   char *tok[],
*   const int ntok,
*   model_t *model)
*
* keyword:
*   fif - restrict IF spectrum frequency range
*
* syntax:
*   fif fmin unit fmax unit
************************************************************/

static int kwfunc_fif(
        char *tok[],
        const int ntok,
        model_t *model)
{
    double fif_min, fif_max;
    if (ntok != 5) {
        parse_error("Expected \"fif fmin unit fmax unit\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_nonneg_dbl_val(
                &fif_min,
                tok[1],
                tok[2],
                UGROUP_FREQUENCY,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_nonneg_dbl_val(
                &fif_max,
                tok[3],
                tok[4],
                UGROUP_FREQUENCY,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (fif_min > fif_max) {
        parse_error("The IF frequency range appears to be out of order.\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    model->fif_min = fif_min;
    model->fif_min_unitnum = get_unitnum(tok[4]);
    model->fif_max = fif_max;
    model->fif_max_unitnum = get_unitnum(tok[2]);
    return KWFUNC_SUCCESS;
} /* kwfunc_fif() */


/***********************************************************
* static int kwfunc_fit(
*   char *tok[],
*   const int ntok,
*   fit_data_t *fit_data)
*
* keyword:
*   fit - fit a single data file
*
* syntax:
*   fit data_type filename
************************************************************/

static int kwfunc_fit(
        char *tok[],
        const int ntok,
        fit_data_t *fit_data)
{
    int i;
    if (fit_data->filename != NULL) {
        parse_error("Redundant fit definition.\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (ntok != 3) {
        parse_error("Expected \"fit data_type filename\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    for (i = 0; i < OUTPUT_END_OF_TABLE; ++i) {
        if (!ci_strcmp(tok[1], output[i].name))
            break;
    }
    if (i == OUTPUT_END_OF_TABLE) {
        parse_error("Unrecognized fit data type \"%s\"\n", tok[1]);
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    fit_data->data_type = i;
    /*
     * When an output is flagged for fitting, set the flag on
     * output[ALL_OUTPUTS] also.
     */
    output[ALL_OUTPUTS].flags |= OUTPUT_FITTED;
    output[i].flags |= OUTPUT_FITTED;
    fit_data->nfiles = 1;
    if ((fit_data->filename = (char**)malloc(sizeof(char*))) == NULL) {
        errlog(33, 0);
        return KWFUNC_FATAL_ERROR;
    }
    if ((fit_data->filename[0]
                = (char*)malloc((1 + strlen(tok[2])) * sizeof(char)))
            == NULL) {
        errlog(33, 0);
        return KWFUNC_FATAL_ERROR;
    }
    strcpy(fit_data->filename[0], tok[2]);
    /*
     * Set the default fit spectrum unit based on the fit data type,
     * unless it has already been set.
     */
    if (fit_data->s_unitnum == UNIT_NONE)
        fit_data->s_unitnum = output[i].default_unitnum;
    return KWFUNC_SUCCESS;
} /* kwfunc_fit() */


/***********************************************************
* static int kwfunc_fits(
*   char *tok[],
*   const int ntok,
*   fit_data_t *fit_data)
*
* keyword:
*   fits - fit multiple files sequentially, using a list of
*   data file names taken from a file.
*
* syntax:
*   fits data_type filename
************************************************************/

static int kwfunc_fits(
        char *tok[],
        const int ntok,
        fit_data_t *fit_data)
{
    int i, status;
    if (fit_data->filename != NULL) {
        parse_error("Redundant fit definition.\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (ntok != 3) {
        parse_error("Expected \"fits data_type filename\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    for (i = 0; i < OUTPUT_END_OF_TABLE; ++i) {
        if (!ci_strcmp(tok[1], output[i].name))
            break;
    }
    if (i == OUTPUT_END_OF_TABLE) {
        parse_error("Unrecognized fit data type \"%s\"\n", tok[1]);
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    fit_data->data_type = i;
    /*
     * When an output is flagged for fitting, set the flag on
     * output[ALL_OUTPUTS] also.
     */
    output[ALL_OUTPUTS].flags |= OUTPUT_FITTED;
    output[i].flags |= OUTPUT_FITTED;
    status = load_fit_file_list(tok[2], fit_data);
    if (status == FILELIST_CANNOT_OPEN) {
        parse_error("Error loading fit data file list \"%s\"\n", tok[2]);
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    } else if (status == FILELIST_ERROR) {
        return KWFUNC_FATAL_ERROR;
    }
    /*
     * Set the default fit spectrum unit based on the fit data type,
     * unless it has already been set.
     */
    if (fit_data->s_unitnum == UNIT_NONE)
        fit_data->s_unitnum = output[i].default_unitnum;
    return KWFUNC_SUCCESS;
} /* kwfunc_fits() */


/***********************************************************
* static int kwfunc_fit_data_columns(
*   char *tok[],
*   const int ntok,
*   fit_data_t *fit_data)
*
* keyword:
*   fit_data_cols or
*   fit_data_columns - set columns for fit input data.
*
* syntax:
*   fit_data_columns c1 c2 [c3 [c4]]
************************************************************/

static int kwfunc_fit_data_columns(
        char *tok[],
        const int ntok,
        fit_data_t *fit_data)
{
    if (ntok < 3 || ntok > 5) {
        parse_error("Expected \"fit_data_columns c1 c2 [c3 [c4]]\"\n");
        parse_error("\twhere the column numbers c1 - c4 are:\n");
        parse_error("\t  c1 = column for frequency data    (default 1)\n");
        parse_error("\t  c2 = column for spectrum data     (default 2)\n");
        parse_error("\t  c3 = column for channel bandwidth (default 0)\n");
        parse_error("\t  c4 = column for weight factor     (default 0)\n");
        parse_error(
                "\tA column number of 0 indicates that default values\n");
        parse_error(
                "\twill be assigned to the channel bandwidth (0.0) or\n");
        parse_error(
                "\tweight factor (1.0)\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_nonneg_int_val(&fit_data->f_col, tok[1]))
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    if (fit_data->f_col == 0) {
        parse_error(
                "A non-zero column must be specified for frequency data.\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_nonneg_int_val(&fit_data->s_col, tok[2]))
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    if (fit_data->s_col == 0) {
        parse_error(
                "A non-zero column must be specified for spectral data.\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (ntok > 3) {
        if (get_nonneg_int_val(&fit_data->b_col, tok[3]))
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (ntok > 4) {
        if (get_nonneg_int_val(&fit_data->w_col, tok[4]))
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    return KWFUNC_SUCCESS;
} /* kwfunc_fit_data_columns() */


/***********************************************************
* static int kwfunc_fit_data_delimiters(
*   char *tok[],
*   const int ntok)
*
* keyword:
*   fit_data_delimiters - set fit data column delimiters
*
* syntax:
*   fit_data_delimiters string
************************************************************/

static int kwfunc_fit_data_delimiters(
        char *tok[],
        const int ntok)
{
    if (ntok < 2) {
        parse_error("Expected \"fit_data_delimiters string\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (fit_data_delimiters(tok[1]) == NULL) {
        parse_error("Bad delimiter string\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    return KWFUNC_SUCCESS;
} /* kwfunc_fit_data_delimiters() */


/***********************************************************
* static int kwfunc_fit_data_format(
*   char *tok[],
*   const int ntok)
*
* keyword:
*   fit_data_format - set fit data file format
*
* syntax:
*   fit_data_format string
************************************************************/

static int kwfunc_fit_data_format(
        char *tok[],
        const int ntok)
{
    if (ntok < 2) {
        parse_error("Expected \"fit_data_format string\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (fit_data_format(tok[1]) == NULL) {
        parse_error("Bad fit data format string\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    return KWFUNC_SUCCESS;
} /* kwfunc_fit_data_format() */


/***********************************************************
* static int kwfunc_fit_data_units(
*   char *tok[],
*   const int ntok,
*   fit_data_t *fit_data)
*
* keyword:
*   fit_data_units - set fit data units
*
* syntax:
*   fit_data_units freq_unit spectrum_unit
************************************************************/

static int kwfunc_fit_data_units(
        char *tok[],
        const int ntok,
        fit_data_t *fit_data)
{
    if (ntok < 3) {
        parse_error("Expected \"fit_data_units freq_unit spectrum_unit\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    fit_data->f_unitnum = get_unitnum(tok[1]);
    if (fit_data->f_unitnum == UNIT_END_OF_TABLE) {
        parse_error("Unrecognized unit: \"%s\"\n", tok[1]);
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    fit_data->s_unitnum = get_unitnum(tok[2]);
    if (fit_data->s_unitnum == UNIT_END_OF_TABLE) {
        parse_error("Unrecognized unit: \"%s\"\n", tok[2]);
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    return KWFUNC_SUCCESS;
} /* kwfunc_fit_data_units() */


/***********************************************************
* static int kwfunc_fit_estimator(
*   char *tok[],
*   const int ntok,
*   fit_data_t *fit_data)
*
* keyword:
*   fit_estimator - set fit estimator type
*
* syntax:
*   fit_estimator estimator_type
************************************************************/

static int kwfunc_fit_estimator(
        char *tok[],
        const int ntok,
        fit_data_t *fit_data)
{
    int i;
    if (ntok < 2) {
        parse_error("Expected \"fit_estimator estimator_type\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    for (i = 0; i < FIT_ESTIMATOR_END_OF_TABLE; ++i) {
        if (!ci_strcmp(tok[1], fit_estimator_type[i].name))
            break;
    }
    if (i == FIT_ESTIMATOR_END_OF_TABLE) {
        parse_error("Unrecognized fit estimator type \"%s\"\n", tok[1]);
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    fit_data->estimator_type = i;
    return KWFUNC_SUCCESS;
} /* kwfunc_fit_estimator() */


/***********************************************************
* static int kwfunc_fit_iter(
*   char *tok[],
*   const int ntok,
*   fit_data_t *fit_data)
*
* keyword:
*   fit_iter - set fit maximum iteration count
*
* syntax:
*   fit_iter n_iter
************************************************************/

static int kwfunc_fit_iter(
        char *tok[],
        const int ntok,
        fit_data_t *fit_data)
{
    if (ntok < 2) {
        parse_error("Expected \"fit_iter n_iter\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_nonneg_int_val(&(fit_data->max_iter), tok[1])) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    return KWFUNC_SUCCESS;
} /* kwfunc_fit_iter() */


/***********************************************************
* static int kwfunc_fit_output(
*   char *tok[],
*   const int ntok,
*   fit_data_t *fit_data)
*
* keyword:
*   fit_output - enable or disable various fit outputs
*
* syntax:
*   fit_output (config | spectrum | residual) 0|1
************************************************************/

static int kwfunc_fit_output(
        char *tok[],
        const int ntok,
        fit_data_t *fit_data)
{
    int mode, mask;
    if (ntok < 2) {
        parse_error(
                "Expected \"fit_output"
                " (config | spectrum | residual) 0|1\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (!ci_strcmp(tok[1], "config")) {
        mask = FIT_OUTPUT_CONFIG;
    } else if (!ci_strcmp(tok[1], "spectrum")) {
        mask = FIT_OUTPUT_SPECTRUM;
    } else if (!ci_strcmp(tok[1], "residual")) {
        mask = FIT_OUTPUT_RESIDUAL;
    } else {
        parse_error(
                "Expected \"fit_output"
                " (config | spectrum | residual) 0|1\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_bool_val(&mode, tok[2]))
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    if (mode)
        fit_data->output_mode |= mask;
    else
        fit_data->output_mode &= ~mask;
    return KWFUNC_SUCCESS;
} /* kwfunc_fit_output() */


/***********************************************************
* static int kwfunc_fit_reinit(
*   char *tok[],
*   const int ntok,
*   simplex_t *simplex)
*
* keyword:
*   fit_reinit - set simplex reinit mode
*
* syntax:
*   fit_reinit 0|1
************************************************************/

static int kwfunc_fit_reinit(
        char *tok[],
        const int ntok,
        simplex_t *simplex)
{
    if (ntok < 2) {
        parse_error("Expected \"fit_reinit 0|1\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_bool_val(&(simplex->reinit), tok[1]))
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    return KWFUNC_SUCCESS;
} /* kwfunc_fit_reinit() */


/***********************************************************
* static int kwfunc_fit_tol(
*   char *tok[],
*   const int ntok,
*   fit_data_t *fit_data,
*   simplex_t *simplex)
*
* keyword:
*   fit_tol - set fit convergence tolerance and maximum number
*   of restarts of a converged fit.
*
* syntax:
*   fit_tol tol [max_restarts]
************************************************************/

static int kwfunc_fit_tol(
        char *tok[],
        const int ntok,
        fit_data_t *fit_data,
        simplex_t *simplex)
{
    if (ntok < 2) {
        parse_error("Expected \"fit_tol tol [max_restarts]\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_nonneg_dbl_val(
                &(simplex->tol),
                tok[1],
                NULL,
                UGROUP_NONE,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (ntok == 3) {
        if (get_nonneg_int_val(&(fit_data->max_restarts), tok[2])) {
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
    }
    return KWFUNC_SUCCESS;
} /* kwfunc_fit_tol() */


/***********************************************************
* static int kwfunc_fit_track_residuals(
*   char *tok[],
*   const int ntok,
*   fit_data_t *fit_data)
*
* keyword:
*   fit_track_residuals - set filter gain for residual
*   tracking, or reset residual tracking.
*
* syntax:
*   fit_track_residuals gain | reset | off
************************************************************/

static int kwfunc_fit_track_residuals(
        char *tok[],
        const int ntok,
        fit_data_t *fit_data)
{
    if (ntok < 2) {
        parse_error(
                "Expected \"fit_track_residuals gain | reset | off\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (!ci_strcmp(tok[1], "reset")) {
        reset_estimated_residuals(fit_data);
        return KWFUNC_SUCCESS;
    }
    if (!ci_strcmp(tok[1], "off")) {
        reset_estimated_residuals(fit_data);
        fit_data->res_track_gain = -1.0;
        return KWFUNC_SUCCESS;
    }
    if (get_nonneg_dbl_val(
                &(fit_data->res_track_gain),
                tok[1],
                NULL,
                UGROUP_NONE,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (fit_data->res_track_gain > 1.0) {
        parse_error("Filter gain cannot exceed unity.\n");
        fit_data->res_track_gain = 0.0;
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    return KWFUNC_SUCCESS;
} /* kwfunc_fit_track_residuals() */


/***********************************************************
* static int kwfunc_fit_verbose(
*   char *tok[],
*   const int ntok,
*   fit_data_t *fit_data)
*
* keyword:
*   fit_verbose - enable fit verbose mode
*
* syntax:
*   fit_verbose 0|1
************************************************************/

static int kwfunc_fit_verbose(
        char *tok[],
        const int ntok,
        fit_data_t *fit_data)
{
    int mode;
    if (ntok < 2) {
        parse_error("Expected \"fit_verbose 0|1\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_bool_val(&mode, tok[1]))
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    if (mode)
        fit_data->output_mode |= FIT_OUTPUT_VERBOSE;
    else
        fit_data->output_mode &= ~FIT_OUTPUT_VERBOSE;
    return KWFUNC_SUCCESS;
} /* kwfunc_fit_verbose() */


/***********************************************************
* static int kwfunc_fout(
*   char *tok[],
*   const int ntok,
*   model_t *model)
*
* keyword:
*   fout - restrict output frequency range
*
* syntax:
*   fout fmin unit fmax unit
************************************************************/

static int kwfunc_fout(
        char *tok[],
        const int ntok,
        model_t *model)
{
    double fout_min, fout_max;
    if (ntok != 5) {
        parse_error("Expected \"fout fmin unit fmax unit\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_nonneg_dbl_val(
                &fout_min,
                tok[1],
                tok[2],
                UGROUP_FREQUENCY,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_nonneg_dbl_val(
                &fout_max,
                tok[3],
                tok[4],
                UGROUP_FREQUENCY,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (fout_min > fout_max) {
        parse_error(
                "The output frequency range appears to be out of order.\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    model->fout_min = fout_min;
    model->fout_min_unitnum = get_unitnum(tok[4]);
    model->fout_max = fout_max;
    model->fout_max_unitnum = get_unitnum(tok[2]);
    return KWFUNC_SUCCESS;
} /* kwfunc_fout() */


/***********************************************************
* static int kwfunc_frequency(
*   char *tok[],
*   const int ntok,
*   model_t *model)
*
* keywords:
*   f, freq, frequency - set frequency grid
*
* syntax:
*   f fmin unit fmax unit df unit
************************************************************/

static int kwfunc_frequency(
        char *tok[],
        const int ntok,
        model_t *model)
{
    double f_min, f_max, df;
    double dbl_imin, dbl_imax;
    double gridsize, gridsize_t_limit, fseek_limit;
    int n;
    if (model->ngrid != 0) {
        parse_error("Redundant definition of frequency grid.\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (ntok != 7) {
        parse_error("Expected \"f fmin unit fmax unit df unit\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_nonneg_dbl_val(
                &f_min,
                tok[1],
                tok[2],
                UGROUP_FREQUENCY,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_nonneg_dbl_val(
                &f_max,
                tok[3],
                tok[4],
                UGROUP_FREQUENCY,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_nonneg_dbl_val(
                &df,
                tok[5],
                tok[6],
                UGROUP_FREQUENCY,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    /*
     * Save fmin, fmax, df, and their unit numbers in the model
     * structure.
     */
    model->fmin = f_min;
    model->fmin_unitnum = get_unitnum(tok[2]);
    model->fmax = f_max;
    model->fmax_unitnum = get_unitnum(tok[4]);
    model->df = df;
    model->df_unitnum = get_unitnum(tok[6]);
    /*
     * Bypass remaining frequency grid range checks and length
     * computations if we're in atmospheric model only (-a) mode.
     * This prevents checking for frequency grid errors when it
     * doesn't matter, so the user can just do
     *
     *   f 0 GHz 0 GHz 0 GHz
     *
     * for example.  In -a mode, it's also OK to have no
     * frequency grid statement at all, but if a grid was defined
     * in the original config file, we do want to preserve it for
     * the stderr output.
     */
    if (output[ALL_OUTPUTS].flags & OUTPUT_AM_ONLY)
        return KWFUNC_SUCCESS;
    /*
     * Frequency grid range checks
     */
    if (f_max < f_min) {
        parse_error("The frequency grid appears to be out of order.\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (df <= DBL_EPSILON * (f_max - f_min)) {
        parse_error("The frequency grid interval is too small.\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if ((f_max - AM_MAX_FREQ) > (DBL_EPSILON * AM_MAX_FREQ)) {
        parse_error("The maximum grid frequency is %g GHz, which is the\n",
                AM_MAX_FREQ);
        parse_error("\tupper limit of am's internal spectroscopic data.\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    /*
     * Calculate the size of the frequency grid, and the padded
     * array sizes for Fourier and Hilbert transforms.  Note that
     * the frequency grid is aligned to the origin; it consists
     * of all points i * df, where i is an integer, on the closed
     * interval fmin <= f <= fmax.
     *
     * There are two limits on the size of the frequency grid;
     * these limits also apply to the padded grids, but will only
     * matter if padded arrays are actually needed.
     *
     * (1) Arrays are indexed with a signed integer of type
     * gridsize_t, defined in am_types.h.  For compatibility with
     * OpenMP's parallel loop construct, which requires a signed
     * integral loop index, gridsize_t is of signed integral
     * type.  Normally, gridsize_t will be a 32-bit int in both
     * 32- and 64-bit environments.
     *
     * (2)  The number of bytes in an array cannot exceed the
     * maximum signed long int.  Note that this is a factor of 2
     * smaller than the maximum pointer size on practically any
     * machine.  The reason for the more restrictive limit is
     * that fseek(), which is used to read arrays cached to disk,
     * uses a signed long int for the file offset.
     *
     * In 32-bit environments, (2) will set the upper limit on
     * array size.  In 64-bit environments (except Win64, which
     * has 32-bit longs), (1) will be the limit if gridsize_t is
     * defined to be a 32-bit int.
     */
    gridsize_t_limit = (double)GRIDSIZE_T_MAX;
    fseek_limit = (double)LONG_MAX / (double)sizeof(double);
    dbl_imin = ceil((f_min * (1.0 - DBL_EPSILON)) / df);
    dbl_imax = floor((f_max * (1.0 + DBL_EPSILON)) / df);
    gridsize = 1.0 + dbl_imax - dbl_imin;
    if (gridsize > gridsize_t_limit || gridsize > fseek_limit) {
        parse_error("The defined frequency grid has too many points.\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    model->ngrid = (gridsize_t)(gridsize + 0.5);
    if (model->ngrid < 1) {
        parse_error("The defined frequency grid contains no points.\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    /*
     * model->npad is the padded array size for convolution of
     * spectra with an instrumental line shape.  It is the
     * smallest power of 2 which is at least three times the
     * model frequency grid size.  See ils.c for more details on
     * the array layouts.  If the array size would be too large,
     * signal this by setting model->npad to zero.  This won't
     * generate an error unless an ils is defined.
     */
    (void)frexp(3.0 * (double)model->ngrid, &n);
    gridsize = ldexp(1.0, n);
    if (gridsize > gridsize_t_limit || gridsize > fseek_limit) {
        model->npad = 0;
    } else {
        model->npad = 1 << n;
    }
    /*
     * model->Lpad is the padded array size for the Hilbert
     * transform to compute the delay spectrum.  It is the
     * smallest power of 2 which can hold the two-sided spectrum
     * (including any internal zero padding from -f_min to f_min,
     * if f_min != 0) further zero- padded by a factor of 4.  If
     * the array size would be too large, set model->Lpad to
     * zero.  This won't generate an error unless L is going to
     * be computed.
     */
    (void)frexp(dbl_imax, &n);
    n += 3;
    gridsize = ldexp(1.0, n);
    if (gridsize > gridsize_t_limit || gridsize > fseek_limit) {
        model->nLpad = 0;
    } else {
        model->nLpad = 1 << n;
    }
    return KWFUNC_SUCCESS;
} /* kwfunc_frequency() */


/***********************************************************
* static int kwfunc_g(
*   char *tok[],
*   const int ntok,
*   model_t *model)
*
* keyword:
*   g - set gravitational acceleration
*
* syntax:
*   g acceleration unit
************************************************************/

static int kwfunc_g(
        char *tok[],
        const int ntok,
        model_t *model)
{
    if (ntok < 3) {
        parse_error("Expected \"g acceleration unit\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_pos_dbl_val(
                &(model->g),
                tok[1],
                tok[2],
                UGROUP_ACCEL,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    model->g_unitnum = get_unitnum(tok[2]);
    return KWFUNC_SUCCESS;
} /* kwfunc_g() */


/***********************************************************
* static int kwfunc_geometry(
*   char *tok[],
*   const int ntok,
*   model_t *model)
*
* keyword:
*   geometry - set geometry mode
*
* syntax:
*   geometry plane-parallel | spherical | limb
************************************************************/

static int kwfunc_geometry(
        char *tok[],
        const int ntok,
        model_t *model)
{
    if (ntok < 2) {
        parse_error(
                "Expected \"geometry plane-parallel | spherical | limb\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (!ci_strcmp(tok[1], "plane-parallel")) {
        model->geometry &= ~GEOMETRY_MODE_BITS;
        model->geometry |= GEOMETRY_PLANE_PARALLEL;
    } else if (!ci_strcmp(tok[1], "spherical")) {
        model->geometry &= ~GEOMETRY_MODE_BITS;
        model->geometry |= GEOMETRY_SPHERICAL;
    } else if (!ci_strcmp(tok[1], "limb")) {
        model->geometry &= ~GEOMETRY_MODE_BITS;
        model->geometry |= GEOMETRY_LIMB;
    } else {
        parse_error(
                "Expected \"geometry plane-parallel | spherical | limb\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    return KWFUNC_SUCCESS;
}   /* kwfunc_geometry() */


/***********************************************************
* static int kwfunc_h(
*   char *tok[],
*   const int ntok,
*   model_t *model,
*   simplex_t *simplex,
*   int lnum)
*
* keyword:
*   h - set layer thickness
*
* syntax:
*   h length unit [dlength unit]
************************************************************/

static int kwfunc_h(
        char *tok[],
        const int ntok,
        model_t *model,
        simplex_t *simplex,
        int lnum)
{
    if (lnum < 0 || lnum >= model->nlayers) {
        parse_error("Attempted to set h for non-existent layer.\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    /*
     * If this is a hydrostatic model, assignment of a layer
     * thickness h is disallowed if the layer pressure drop is
     * already known to be non-zero.
     */
    if (model->PTmode & PTMODE_HYDROSTATIC) {
        int errflag = 0;
        if (model->PTmode & PTMODE_PBASE) {
            if (lnum == 0) {
                errflag = (model->layer[lnum]->Pbase > 0.0);
            } else {
                errflag = (fabs(model->layer[lnum]->Pbase -
                            model->layer[lnum]->Pbase) > 
                        (DBL_EPSILON * model->layer[lnum]->Pbase));
            }
        } else if (model->PTmode & PTMODE_DP) {
            errflag = model->layer[lnum]->dP_def > 0.0 ||
                isvar(simplex, &(model->layer[lnum]->dP_def));
        }
        if (errflag) {
            parse_error(
                    "In hydrostatic models, layer thickness h can only be\n");
            parse_error(
                    "\tassigned if there is no pressure drop across the\n");
            parse_error(
                    "\tlayer.\n");
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
    }
    if (ntok != 3 && ntok != 5) {
        parse_error("Expected \"h length unit [dlength unit]\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_nonneg_dbl_val(
                &(model->layer[lnum]->h),
                tok[1],
                tok[2],
                UGROUP_DELAY_DIST,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    model->layer[lnum]->h_unitnum = get_unitnum(tok[2]);
    if (ntok == 5 && simplex != NULL) {
        double dh;
        if (get_pos_dbl_val(
                    &dh,
                    tok[3],
                    tok[4],
                    UGROUP_DELAY_DIST,
                    1)) {
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
        snprintf(VARNAME, sizeof(VARNAME), "h (layer %d)", lnum);
        if (add_var_to_simplex(
                    simplex,
                    VARNAME,
                    &(model->layer[lnum]->h),
                    model->layer[lnum]->h,
                    dh,
                    model->layer[lnum]->h_unitnum,
                    MAPPING_NONE)) {
            return KWFUNC_FATAL_ERROR;
        }
    }
    return KWFUNC_SUCCESS;
} /* kwfunc_h() */


/***********************************************************
* static int kwfunc_headers(
*   char *tok[],
*   const int ntok)
*
* keyword:
*   headers - print output column headers
*
* syntax:
*   headers 0|1
************************************************************/

static int kwfunc_headers(
        char *tok[],
        const int ntok)
{
    int headers;
    if (ntok < 2) {
        parse_error("Expected \"headers 0|1\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_bool_val(&headers, tok[1]))
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    if (headers)
        output[ALL_OUTPUTS].flags |= OUTPUT_HEADERS;
    else
        output[ALL_OUTPUTS].flags &= ~OUTPUT_HEADERS;
    return KWFUNC_SUCCESS;
} /* kwfunc_headers() */


/***********************************************************
* static int kwfunc_ifspec(
*   char *tok[],
*   const int ntok,
*   model_t *model,
*   simplex_t *simplex)
*
* keyword:
*   ifspec - set IF spectrum mode
*
* syntax:
*   ifspec (dsb | usb | lsb) flo unit [dflo unit])
************************************************************/

static int kwfunc_ifspec(
        char *tok[],
        const int ntok,
        model_t *model,
        simplex_t *simplex)
{
    if (ntok != 4 && ntok != 6) {
        parse_error(
                "Expected \"ifspec (dsb | usb | lsb)"
                " flo unit [dflo unit]\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (!ci_strcmp(tok[1], "dsb")) {
        model->ifmode = IFMODE_DSB;
    } else if (!ci_strcmp(tok[1], "usb")) {
        model->ifmode = IFMODE_USB;
    } else if (!ci_strcmp(tok[1], "lsb")) {
        model->ifmode = IFMODE_LSB;
    } else {
        parse_error(
                "Expected \"ifspec (dsb | usb | lsb)"
                " flo unit [dflo unit]\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_dbl_val(
                &(model->flo),
                tok[2],
                tok[3],
                UGROUP_FREQUENCY,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    model->flo_unitnum = get_unitnum(tok[3]);
    if (ntok == 6 && simplex != NULL) {
        double dflo;
        if (get_pos_dbl_val(
                    &dflo,
                    tok[4],
                    tok[5],
                    UGROUP_FREQUENCY,
                    1)) {
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
        if (add_var_to_simplex(
                    simplex,
                    "flo",
                    &(model->flo),
                    model->flo,
                    dflo,
                    model->flo_unitnum,
                    MAPPING_NONE)) {
            return KWFUNC_FATAL_ERROR;
        }
    }
    return KWFUNC_SUCCESS;
} /* kwfunc_ifspec() */


/***********************************************************
* static int kwfunc_ils(
*   char *tok[],
*   const int ntok,
*   model_t *model,
*   simplex_t *simplex)
*
* keyword:
*   ils - instrumental line shape
*
* syntax:
*   ils type_name fwhm unit [dfwhm unit] {array_name}*
************************************************************/

static int kwfunc_ils(
        char *tok[],
        const int ntok,
        model_t *model,
        simplex_t *simplex)
{
    int i, itok;
    int status;
    if (ntok < 4) {
        parse_error(
                "Expected \"ils type_name fwhm unit [dfwhm unit]"
                " {array_name}*\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    for (i = 0; i < ILS_END_OF_TABLE; ++i) {
        if (!ci_strcmp(tok[1], ils_type[i].name))
            break;
    }
    if (i == ILS_END_OF_TABLE) {
        parse_error("Unrecognized ILS type \"%s\"\n", tok[1]);
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (i == ILS_NONE)
        return KWFUNC_SUCCESS;
    model->ils_typenum = i;
    if (get_nonneg_dbl_val(
                &(model->ils_fwhm),
                tok[2],
                tok[3],
                UGROUP_FREQUENCY,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    model->ils_fwhm_unitnum = get_unitnum(tok[3]);
    /*
     * If there are at least 6 tokens, and tok[5] is a unit, then
     * the ils fwhm is a fit variable.
     */
    if (ntok >= 6 && is_unit(tok[5]) && simplex != NULL) {
        double dfwhm;
        if (get_pos_dbl_val(
                    &dfwhm,
                    tok[4],
                    tok[5],
                    UGROUP_FREQUENCY,
                    1)) {
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
        if (add_var_to_simplex(
                    simplex,
                    "ils fwhm",
                    &(model->ils_fwhm),
                    model->ils_fwhm,
                    dfwhm,
                    model->ils_fwhm_unitnum,
                    MAPPING_NONE)) {
            return KWFUNC_FATAL_ERROR;
        }
        itok = 6;
    } else {
        itok = 4;
    }
    /*
     * If there are any named arrays following the ils
     * specification, these are the user-specified arrays to be
     * convolved with the ils.  Clear the default flags from the
     * global output table, and put in the user-specified ones.
     * Also, set the ILS_APPLIED flag in output[ALL_OUTPUTS] to
     * indicate user specification of arrays convolved with the
     * ils.
     */
    if (itok < ntok) {
        for (i = 0; i < OUTPUT_END_OF_TABLE; ++i)
            output[i].flags &= ~(ILS_APPLIED);
        output[ALL_OUTPUTS].flags |= ILS_APPLIED;
    }
    status = KWFUNC_SUCCESS;
    while (itok < ntok) {
        for (i = 0; i < OUTPUT_END_OF_TABLE; ++i) {
            if (!ci_strcmp(tok[itok], output[i].name))
                break;
        }
        if (i == OUTPUT_END_OF_TABLE) {
            parse_error(
                    "Expected \"ils type_name num unit [num unit]"
                    " {array_name}*\"\n");
            parse_error("\t(\"%s\" is not a valid array name.)\n", tok[itok]);
            status = KWFUNC_SYNTAX_ERROR_CONTINUE;
        } else if (!(output[i].flags & ILS_ALLOWED)) {
            parse_error(
                    "The array \"%s\" may not be convolved with an ILS.\n",
                    output[i].name);
            status = KWFUNC_SYNTAX_ERROR_CONTINUE;
        } else {
            output[i].flags |= ILS_APPLIED;
        }
        ++itok;
    }
    return status;
} /* kwfunc_ils() */


/***********************************************************
* static int kwfunc_ilsmode(
*   char *tok[],
*   const int ntok,
*   model_t *model,
*   simplex_t *simplex)
*
* keyword:
*   ilsmode - instrumental lineshape mode
*
* syntax:
*   ilsmode normal | ((dsb | usb | lsb) fif unit [dfif unit])
************************************************************/

static int kwfunc_ilsmode(
        char *tok[],
        const int ntok,
        model_t *model,
        simplex_t *simplex)
{
    if (ntok < 2) {
        parse_error(
                "Expected \"ilsmode normal | ((dsb | usb | lsb) fif unit"
                " [dfif unit])\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (!ci_strcmp(tok[1], "normal")) {
        model->ilsmode = ILSMODE_NORMAL;
    } else if (!ci_strcmp(tok[1], "dsb")) {
        model->ilsmode = ILSMODE_DSB;
    } else if (!ci_strcmp(tok[1], "usb")) {
        model->ilsmode = ILSMODE_USB;
    } else if (!ci_strcmp(tok[1], "lsb")) {
        model->ilsmode = ILSMODE_LSB;
    } else {
        parse_error(
                "Expected \"ilsmode normal | ((dsb | usb | lsb) fif unit"
                " [dfif unit])\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if ((model->ilsmode == ILSMODE_NORMAL && ntok != 2) ||
            (model->ilsmode != ILSMODE_NORMAL && (ntok != 4 && ntok != 6))) {
        parse_error(
                "Expected \"ilsmode normal | ((dsb | usb | lsb) fif unit"
                " [dfif unit])\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_dbl_val(
                &(model->ils_fif),
                tok[2],
                tok[3],
                UGROUP_FREQUENCY,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    model->ils_fif_unitnum = get_unitnum(tok[3]);
    if (ntok == 6 && simplex != NULL) {
        double dfif;
        if (get_pos_dbl_val(
                    &dfif,
                    tok[4],
                    tok[5],
                    UGROUP_FREQUENCY,
                    1)) {
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
        if (add_var_to_simplex(
                    simplex,
                    "ils fif",
                    &(model->ils_fif),
                    model->ils_fif,
                    dfif,
                    model->ils_fif_unitnum,
                    MAPPING_NONE)) {
            return KWFUNC_FATAL_ERROR;
        }
    }
    return KWFUNC_SUCCESS;
} /* kwfunc_ilsmode() */


/***********************************************************
* static int kwfunc_jacobian(
*   char *tok[],
*   const int ntok)
*
* keyword:
*   jacobian - flag selected outputs for Jacobian computation
*
* syntax:
*   jacobian {array_name}+ [include_estimated_errors]
************************************************************/

static int kwfunc_jacobian(
        char *tok[],
        const int ntok)
{
    int i, itok;
    int status = KWFUNC_SUCCESS;
    if (ntok < 2) {
        parse_error(
                "Expected \"jacobian {array_name}+"
                " [include_estimated_errors]\"\n");
        status = KWFUNC_SYNTAX_ERROR_CONTINUE;
        return status;
    }
    for (itok = 1; itok < ntok; ++itok) {
        /*
         * Look for "include_estimated_errors" keyword.
         */
        if (!ci_strcmp(tok[itok], "include_estimated_errors")) {
            output[ALL_OUTPUTS].flags |= OUTPUT_JACOBIAN_ERRS;
            continue;
        }
        /*
         * Look for a match to an array name.
         */
        for (i = 0; i < OUTPUT_END_OF_TABLE; ++i) {
            if (!ci_strcmp(tok[itok], output[i].name))
                break;
        }
        if (i == OUTPUT_END_OF_TABLE) {
            parse_error("Unrecognized array name (%s?)\n", tok[itok]);
            status = KWFUNC_SYNTAX_ERROR_CONTINUE;
        } else if (output[i].flags & OUTPUT_JACOBIAN) {
            parse_error(
                    "redundant Jacobian specification (%s)\n", tok[itok]);
            status = KWFUNC_SYNTAX_ERROR_CONTINUE;
        } else if (!(output[i].flags & JACOBIAN_ALLOWED)) {
            parse_error(
                    "Jacobian cannot be computed for %s.\n", tok[itok]);
            status = KWFUNC_SYNTAX_ERROR_CONTINUE;
        } else {
            /*
             * When any output is flagged for Jacobians, set the flag
             * on output[ALL_OUTPUTS] also.
             */
            output[ALL_OUTPUTS].flags |= OUTPUT_JACOBIAN;
            output[i].flags |= OUTPUT_JACOBIAN;
            /*
             * When a Jacobian is being computed, the corresponding
             * reference state is also included as output.
             */
            output[i].flags |= OUTPUT_USER;
        }
    }
    return status;
} /* kwfunc_jacobian() */


/***********************************************************
* static int kwfunc_kcache(
*   char *tok[],
*   const int ntok,
*   model_t *model)
*
* keyword:
*   kcache - control absorption coefficient caching
*
* syntax:
*   kcache off | Tmin unit Tmax unit dT unit
************************************************************/

static int kwfunc_kcache(
        char *tok[],
        const int ntok,
        model_t *model)
{
    if (ntok != 2 && ntok != 7) {
        parse_error(
                "Expected \"kcache off | Tmin unit Tmax unit dT unit\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (ntok == 2) {
        if (!ci_strcmp(tok[1], "off")) {
            model->nkcache = 0;
            return KWFUNC_SUCCESS;
        } else {
            parse_error(
                    "Expected \"kcache off | Tmin unit Tmax unit dT unit\"\n");
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
    }
    if (get_nonneg_dbl_val(
                &(model->kcache_Tmin),
                tok[1],
                tok[2],
                UGROUP_TEMPERATURE,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_nonneg_dbl_val(
                &(model->kcache_Tmax),
                tok[3],
                tok[4],
                UGROUP_TEMPERATURE,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_nonneg_dbl_val(
                &(model->kcache_dT),
                tok[5],
                tok[6],
                UGROUP_TEMPERATURE,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if ((model->kcache_Tmin < 0.) ||
            (model->kcache_Tmin >= model->kcache_Tmax) ||
            (model->kcache_dT <= 0.) ||
            (model->kcache_dT > (model->kcache_Tmax - model->kcache_Tmin))
       ) {
        parse_error("kcache range appears to be out of order.\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    model->nkcache = 1 +
        (int)(0.5 + ((model->kcache_Tmax - model->kcache_Tmin)
                    / model->kcache_dT));
    return KWFUNC_SUCCESS;
} /* kwfunc_kcache() */


/***********************************************************
* static int kwfunc_layer(
*   model_t *model)
*
* keyword:
*   layer - add layer to model
*
* syntax:
*   layer [tag]
************************************************************/

static int kwfunc_layer(
        char *tok[],
        const int ntok,
        model_t *model)
{
    if (model->nlayers > 0) {
        int status = KWFUNC_SUCCESS;
        if (model->layer[model->nlayers - 1]->P_unitnum == UNIT_NONE) {
            parse_error(
                    "Missing or invalid pressure definition"
                    " on prior layer.\n");
            status = KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
        if (model->layer[model->nlayers - 1]->T_unitnum == UNIT_NONE) {
            parse_error(
                    "Missing or invalid temperature definition"
                    " on prior layer.\n");
            status = KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
        if (status != KWFUNC_SUCCESS)
            return status;
    }
    if (add_layer(model, LAYER_TYPE_DEF))
        return KWFUNC_FATAL_ERROR;
    if (ntok > 1) {
        int n;
        if (get_col_typenum(tok[1]) != COL_TYPE_NONE) {
            parse_error(
                    "Column type names may not be used as layer tags.\n");
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
        n = tag_string_index(tok[1]);
        if (n < 0)
            return KWFUNC_FATAL_ERROR;
        model->layer[model->nlayers - 1]->tagnum = n;
    }
    return KWFUNC_SUCCESS;
} /* kwfunc_layer() */


/***********************************************************
* static int kwfunc_lineshape(
*   char *tok[],
*   const int ntok,
*   model_t *model)
*
* keyword:
*   lineshape - set lineshape types and self-broadening mode
*   for the current layer
*
* syntax:
*   lineshape type_name [strict_self_broadening] {k_type}*
************************************************************/

static int kwfunc_lineshape(
        char *tok[],
        const int ntok,
        model_t *model)
{
    layer_t *layer;
    int lineshape, strict_selfbroad;
    int itok;
    if (model->nlayers == 0) {
        parse_error("Attempted to set lineshape on a non-existent layer.\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    layer = model->layer[model->nlayers - 1];
    if (ntok < 2) {
        parse_error(
                "Expected \"lineshape type_name [strict_self_broadening]"
                " {k_type}*\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    for (lineshape = 0; lineshape < LINESHAPE_END_OF_TABLE; ++lineshape) {
        if (!ci_strcmp(tok[1], lineshape_type[lineshape].name))
            break;
    }
    if (lineshape == LINESHAPE_END_OF_TABLE) {
        parse_error("Unrecognized lineshape type \"%s\"\n", tok[1]);
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    strict_selfbroad = 0;
    itok = 2;
    if (itok < ntok) {
        if (!ci_strcmp(tok[2], "strict_self_broadening")) {
            strict_selfbroad = 1;
            ++itok;
        }
    }
    /*
     * If no absorption coefficient types are explicitly named,
     * apply the lineshape type and strict self broadening flag
     * to all line-by-line absorption coefficient types on this
     * layer.  This behavior is backward-compatible with older
     * versions of am.
     *
     * Otherwise, apply to the named absorption coefficient
     * types.
     */
    if (itok == ntok) {
        int knum;
        for (knum = 0; knum < K_TYPE_END_OF_TABLE; ++knum) {
            if (k_type[knum].dep_flags & DEP_ON_LSHAPE)
                layer->lineshape[knum] = lineshape;
            if (k_type[knum].dep_flags & DEP_ON_VMR_SELFBROAD)
                layer->strict_selfbroad[knum] = strict_selfbroad;
        }
    } else {
        while (itok < ntok) {
            int knum;
            for (knum = 0; knum < K_TYPE_END_OF_TABLE; ++knum) {
                if (!ci_strcmp(tok[itok], k_type[knum].name))
                    break;
            }
            if (k_type[knum].dep_flags & DEP_ON_LSHAPE) {
                layer->lineshape[knum] = lineshape;
            } else { /* not line-by-line, or not found in k_type table */
                parse_error(
                        "%s is not an absorption coefficient type with a\n",
                        tok[itok]);
                parse_error(
                        "\tuser-modifiable lineshape.  Valid choices are:\n");
                list_line_by_line_abscoeffs();
                return KWFUNC_SYNTAX_ERROR_CONTINUE;
            }
            if (k_type[knum].dep_flags & DEP_ON_VMR_SELFBROAD)
                layer->strict_selfbroad[knum] = strict_selfbroad;
            ++itok;
        }
    }
    return KWFUNC_SUCCESS;
} /* kwfunc_lineshape() */


/***********************************************************
* static int kwfunc_Mair(
*   char *tok[],
*   const int ntok,
*   model_t *model,
*   int lnum)
*
* keyword:
*   Mair - set relative molar mass for air on a layer
*
* syntax:
*   Mair [M0] {column_type}*
************************************************************/

static int kwfunc_Mair(
        char *tok[],
        const int ntok,
        model_t *model,
        int lnum)
{
    int i;
    layer_t *layer;
    if (lnum < 0 || lnum >= model->nlayers) {
        parse_error("Attempted to set Mair for non-existent layer.\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (ntok < 2) {
        parse_error("Expected \"Mair [M0] {column_type}*\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    layer = model->layer[lnum];
    for (i = 0; i < COL_TYPE_END_OF_TABLE; ++i)
        layer->Mair_flag[i] = 0;
    /*
     * Mair_flag[0] will serve as a flag indicating whether some other
     * element of Mair_flag[] has been set.
     */
    layer->Mair_flag[0] = 0;
    if (is_numeric_token(tok[1])) {
        if (get_nonneg_dbl_val(
                    &(layer->M0),
                    tok[1],
                    NULL,
                    UGROUP_NONE,
                    0)) {
            parse_error("Expected \"Mair [M0] {column_type}*\"\n");
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
        i = 2;
    } else {
        i = 1;
    }
    while (i < ntok) {
        int col_typenum = get_col_typenum(tok[i++]);
        if (col_typenum == COL_TYPE_NONE) {
            parse_error("Expected \"Mair [M0] {column_type}*\"\n");
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
        layer->Mair_flag[col_typenum] = 1;
        layer->Mair_flag[0] = 1;
    }
    return KWFUNC_SUCCESS;
} /* kwfunc_Mair() */


/***********************************************************
* static int kwfunc_Nscale(
*   char *tok[],
*   const int ntok,
*   simplex_t *simplex)
*
* keyword:
*   Nscale - define a column density scale factor for a column type
*
* syntax:
*   Nscale column_type scale_factor [dscale_factor]
************************************************************/

static int kwfunc_Nscale(
        char *tok[],
        const int ntok,
        simplex_t *simplex)
{
    double Nscale;
    Nscale_list_t *list_entry;
    int i, col_typenum, tagnum;
    if (ntok < 3 || ntok > 5) {
        parse_error(
                "Expected \"Nscale [tag] column_type"
                " scale_factor [dscale_factor]\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    /*
     * If the second token is a recognized column type, then
     * there is no associated tag, and tok[2] should be the scale
     * factor.
     */
    if ((col_typenum = get_col_typenum(tok[1])) != COL_TYPE_NONE) {
        tagnum = 0;
        i = 2;
    /*
     * If there were only three tokens, but the first one wasn't
     * recognized as a column type, write an error which depends
     * on whether the third token looks like a number.  If it
     * does, then the column type name probably contains a typo,
     * so tell the user it was not recognized.  If the third
     * token is not a number, then give the user the full Nscale
     * statement syntax.
     */
    } else if (ntok == 3) {
        if (is_numeric_token(tok[2])) {
            parse_error("Unrecognized column type \"%s\"\n", tok[1]);
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        } else {
            parse_error(
                    "Expected \"Nscale [tag] column_type"
                    " scale_factor [dscale_factor]\"\n");
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
    /*
     * If the third token is a recognized column type, and there
     * are at least four tokens, then the first token is a layer
     * tag, and tok[3] should be the scale factor.
     */
    } else if ((col_typenum = get_col_typenum(tok[2])) != COL_TYPE_NONE) {
        tagnum = tag_string_index(tok[1]);
        i = 3;
    /*
     * If the third token should have been a column type, but
     * wasn't, write an error message, again depending on whether
     * it looked like a number.
     */
    } else if (is_numeric_token(tok[2])) {
        parse_error(
                "Expected \"Nscale [tag] column_type"
                " scale_factor [dscale_factor]\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    } else {
        parse_error("Unrecognized column type \"%s\"\n", tok[2]);
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_nonneg_dbl_val(
                &Nscale,
                tok[i],
                NULL,
                UGROUP_NONE,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    list_entry = create_Nscale_list_entry(col_typenum, tagnum, Nscale);
    if (list_entry == NULL)
        return KWFUNC_FATAL_ERROR;
    if (++i < ntok && simplex != NULL) {
        double dNscale;
        if (get_pos_dbl_val(
                    &dNscale,
                    tok[i],
                    NULL,
                    UGROUP_NONE,
                    1)) {
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
        if (list_entry->tagnum == 0) { /* no tag */
            snprintf(VARNAME, sizeof(VARNAME),
                    "Nscale %s",
                    col_type[list_entry->col_typenum].name);
        } else {
            snprintf(VARNAME, sizeof(VARNAME),
                    "Nscale %s %s",
                    tag_string(list_entry->tagnum),
                    col_type[list_entry->col_typenum].name);
        }
        if (add_var_to_simplex(
                    simplex,
                    VARNAME,
                    &(list_entry->Nscale),
                    list_entry->Nscale,
                    dNscale,
                    UNIT_NONE,
                    MAPPING_NONE)) {
            return KWFUNC_FATAL_ERROR;
        }
    }
    return KWFUNC_SUCCESS;
} /* kwfunc_Nscale() */


/***********************************************************
* static int kwfunc_output(
*   char *tok[],
*   const int ntok)
*
* keyword:
*   output - define output columns and units
*
* syntax:
*   output {array_name [unit]}+ |
*   output format {array_name [unit]}*
************************************************************/

static int kwfunc_output(
        char *tok[],
        const int ntok)
{
    int i, itok, colnum, unitnum;
    int status;
    if (ntok < 2) {
        parse_error("Expected\n");
        parse_error("\t  \"output {array_name [unit]}+\"\n");
        parse_error("\tor\n");
        parse_error("\t  \"output format {array_name [unit]}*\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    /*
     * If the first token is an output format name, set the
     * output format.
     */
    itok = 1;
    for (i = 0; i < OUTPUT_FORMAT_END_OF_TABLE; ++i) {
        if (!ci_strcmp(tok[1], output_format_name[i])) {
            output[ALL_OUTPUTS].format = i;
            ++itok;
            break;
        }
    }
    if (itok >= ntok)
        return KWFUNC_SUCCESS;
    /*
     * Clear the defaults from the global output table and column
     * index table.
     */
    for (i = 0; i < OUTPUT_END_OF_TABLE; ++i) {
        output[i].flags &= ~(OUTPUT_USER | OUTPUT_ACTIVE);
        outcol[i] = 0;
    }
    colnum = 0;
    status = KWFUNC_SUCCESS;
    while (itok < ntok) {
        /*
         * Look for a match to the output data array name.
         */
        for (i = 0; i < OUTPUT_END_OF_TABLE; ++i) {
            if (!ci_strcmp(tok[itok], output[i].name))
                break;
        }
        if (i == OUTPUT_END_OF_TABLE) {
            parse_error("Expected\n");
            parse_error("\t  \"output {array_name [unit]}+\"\n");
            parse_error("\tor\n");
            parse_error("\t  \"output format {array_name [unit]}*\"\n");
            if (itok == 1) {
                parse_error(
                        "Unrecognized format, array_name, or unit (%s?)\n",
                        tok[itok]);
                parse_error(
                        "Valid formats are:\n");
                list_allowed_keywords(
                        output_format_name, OUTPUT_FORMAT_END_OF_TABLE);
            } else {
                parse_error(
                        "Unrecognized array_name or unit (%s?)\n",
                        tok[itok]);
            }
            parse_error("Valid array_names are:\n");
            list_allowed_output_array_names();
            status = KWFUNC_SYNTAX_ERROR_CONTINUE;
        } else if (output[i].flags & OUTPUT_USER) {
            parse_error(
                    "Redundant output array specification (%s)\n", tok[itok]);
            status = KWFUNC_SYNTAX_ERROR_CONTINUE;
        } else {
            output[i].flags |= OUTPUT_USER;
            outcol[colnum] = i;
            ++colnum;
        }
        /*
         * If there are no more tokens, then this is the last
         * output column, and no user unit has been defined.
         */
        if (++itok >= ntok)
            break;
        /*
         * See if the next token is a unit.  If so, it replaces
         * the default output unit.  Otherwise, it should be the
         * next array name.
         */
        if ((unitnum = get_unitnum(tok[itok])) != UNIT_END_OF_TABLE) {
            int old_ugroup = unit_tab[output[i].unitnum].group;
            int new_ugroup = unit_tab[unitnum].group;
            if (old_ugroup != UNIT_AUTO && old_ugroup != new_ugroup) {
                parse_error(
                        "The unit \"%s\" cannot be applied to the output\n",
                        tok[itok]);
                parse_error(
                        "array \"%s\".  Valid choices are:\n",
                        output[i].name);
                list_allowed_units(old_ugroup);
                status = KWFUNC_SYNTAX_ERROR_CONTINUE;
            } else {
                output[i].unitnum = unitnum;
            }
            ++itok;
        }
    }
    return status;
} /* kwfunc_output() */


/***********************************************************
* static int kwfunc_P(
*   char *tok[],
*   const int ntok,
*   model_t *model,
*   simplex_t *simplex,
*   int lnum)
*
* keyword:
*   P - set layer midpoint pressure
*
* syntax:
*   P pressure unit [dpressure unit]
************************************************************/

static int kwfunc_P(
        char *tok[],
        const int ntok,
        model_t *model,
        simplex_t *simplex,
        int lnum)
{
    if (lnum < 0 || lnum >= model->nlayers) {
        parse_error("Attempted to set P for non-existent layer.\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if ((model->PTmode & PTMODE_PMODES) && !(model->PTmode & PTMODE_P)) {
        parse_error(
                "P is inconsistent with a prior layer or PTmode setting.\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (ntok != 3 && ntok != 5) {
        parse_error("Expected \"P pressure unit [dpressure unit]\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    model->PTmode |= PTMODE_P;
    if (get_nonneg_dbl_val(
                &(model->layer[lnum]->P),
                tok[1],
                tok[2],
                UGROUP_PRESSURE,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    model->layer[lnum]->P_unitnum = get_unitnum(tok[2]);
    if (ntok == 5 && simplex != NULL) {
        double dP;
        if (get_pos_dbl_val(
                    &dP,
                    tok[3],
                    tok[4],
                    UGROUP_PRESSURE,
                    1)) {
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
        snprintf(VARNAME, sizeof(VARNAME), "P (layer %d)", lnum);
        if (add_var_to_simplex(
                    simplex,
                    VARNAME,
                    &(model->layer[lnum]->P),
                    model->layer[lnum]->P,
                    dP,
                    model->layer[lnum]->P_unitnum,
                    MAPPING_NONE)) {
            return KWFUNC_FATAL_ERROR;
        }
    }
    return KWFUNC_SUCCESS;
} /* kwfunc_P() */


/***********************************************************
* static int kwfunc_Pbase(
*   char *tok[],
*   const int ntok,
*   model_t *model,
*   int lnum)
*
* keyword:
*   Pbase - set layer base pressure
*
* syntax:
*   Pbase pressure unit
*
*   Note that Pbase is not available as a fit variable, to
*   avoid the possibility of crossing base pressures.
************************************************************/

static int kwfunc_Pbase(
        char *tok[],
        const int ntok,
        model_t *model,
        int lnum)
{
    double Pbase;

    if (lnum < 0 || lnum >= model->nlayers) {
        parse_error("Attempted to set Pbase for non-existent layer.\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if ((model->PTmode & PTMODE_PMODES) && !(model->PTmode & PTMODE_PBASE)) {
        parse_error(
                "Pbase is inconsistent with a prior layer"
                " or PTmode setting.\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (ntok != 3) {
        parse_error("Expected \"Pbase pressure unit\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    model->PTmode |= PTMODE_PBASE;
    if (get_nonneg_dbl_val(
                &Pbase,
                tok[1],
                tok[2],
                UGROUP_PRESSURE,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if ((lnum > 0) && (Pbase < model->layer[lnum - 1]->Pbase)) {
        parse_error("Layer base pressures would be out of order.\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if ((lnum < model->nlayers - 1) &&
            (Pbase > model->layer[lnum + 1]->Pbase)) {
        parse_error("Layer base pressures would be out of order.\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    model->layer[lnum]->Pbase = Pbase;
    model->layer[lnum]->P_unitnum = get_unitnum(tok[2]);
    /*
     * In hydrostatic models, h can be assigned only on layers
     * with zero pressure drop.
     */
    if (model->layer[lnum]->h >= 0.0) {
        int errflag = 0;
        if (lnum == 0) {
            errflag = (model->layer[lnum]->Pbase > 0.0);
        } else {
            errflag = (fabs(model->layer[lnum]->Pbase -
                        model->layer[lnum]->Pbase) > 
                    (DBL_EPSILON * model->layer[lnum]->Pbase));
        }
        if (errflag) {
            parse_error(
                    "In hydrostatic models, layer thickness h can only be\n");
            parse_error(
                    "\tassigned if there is no pressure drop across the\n");
            parse_error(
                    "\tlayer.\n");
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
    }
    return KWFUNC_SUCCESS;
} /* kwfunc_Pbase() */


/***********************************************************
* static int kwfunc_Pobs(
*   char *tok[],
*   const int ntok,
*   model_t *model,
*   simplex_t *simplex)
*
* keyword:
*   Pobs - set observing level by pressure.
*
* syntax:
*   Pobs pressure unit [dpressure unit]
************************************************************/

static int kwfunc_Pobs(
        char *tok[],
        const int ntok,
        model_t *model,
        simplex_t *simplex)
{
    if (ntok != 3 && ntok != 5) {
        parse_error("Expected \"Pobs pressure unit [dpressure unit]\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (model->geometry & GEOMETRY_ZOBS_USER_DEFINED) {
        parse_error("Pobs and zobs cannot both be user-defined.\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    /*
     * Zero is accepted as input here, but only so we can give
     * an informative error message for that specific case.
     */
    if (get_nonneg_dbl_val(
                &(model->Pobs),
                tok[1],
                tok[2],
                UGROUP_PRESSURE,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (model->Pobs == 0.0) {
        parse_error("Cannot set Pobs to zero.  Instead, use zobs to\n");
        parse_error("\tlocate the observing level outside the atmosphere\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    model->Pobs_unitnum = get_unitnum(tok[2]);
    model->geometry |= GEOMETRY_POBS_USER_DEFINED;
    if (ntok == 5 && simplex != NULL) {
        double dPobs;
        if (get_pos_dbl_val(
                    &dPobs,
                    tok[3],
                    tok[4],
                    UGROUP_PRESSURE,
                    1)) {
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
        if (add_var_to_simplex(
                    simplex,
                    "Pobs",
                    &(model->Pobs),
                    model->Pobs,
                    dPobs,
                    model->Pobs_unitnum,
                    MAPPING_NONE)) {
            return KWFUNC_FATAL_ERROR;
        }
    }
    return KWFUNC_SUCCESS;
} /* kwfunc_Pobs() */


/***********************************************************
* static int kwfunc_Psource(
*   char *tok[],
*   const int ntok,
*   model_t *model,
*   simplex_t *simplex)
*
* keyword:
*   Psource - set source level by pressure.
*
*   For paths below the astronomical horizon, the source
*   level is intersected on both the near and far sides
*   of the tangent point.  The optional "near" or "far"
*   arguments select which point should be the path
*   endpoint.  The default is "far".
*
* syntax:
*   Psource pressure unit [dpressure unit] [near | far]
************************************************************/

static int kwfunc_Psource(
        char *tok[],
        const int ntok,
        model_t *model,
        simplex_t *simplex)
{
    if (ntok < 3 || ntok > 6 || (ntok == 5 && !is_unit(tok[4]))) {
        parse_error("Expected\n");
        parse_error("\t  "
                "\"Psource pressure unit [dpressure unit] [near | far]\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (model->geometry & GEOMETRY_ZSOURCE_USER_DEFINED) {
        parse_error("Psource and zsource cannot both be user-defined.\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    /*
     * Zero is accepted as input here, but only so we can give
     * an informative error message for that specific case.
     */
    if (get_nonneg_dbl_val(
                &(model->Psource),
                tok[1],
                tok[2],
                UGROUP_PRESSURE,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (model->Psource == 0.0) {
        parse_error("Cannot set Psource to zero.  Instead, use zsource\n");
        parse_error("\tto locate the source level outside the atmosphere\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    model->Psource_unitnum = get_unitnum(tok[2]);
    model->geometry |= GEOMETRY_PSOURCE_USER_DEFINED;
    if (ntok >= 5 && is_unit(tok[4]) && simplex != NULL) {
        double dPsource;
        if (get_pos_dbl_val(
                    &dPsource,
                    tok[3],
                    tok[4],
                    UGROUP_PRESSURE,
                    1)) {
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
        if (add_var_to_simplex(
                    simplex,
                    "Psource",
                    &(model->Psource),
                    model->Psource,
                    dPsource,
                    model->Psource_unitnum,
                    MAPPING_NONE)) {
            return KWFUNC_FATAL_ERROR;
        }
    }
    if (ntok == 4 || ntok == 6) {
        if (!ci_strcmp(tok[ntok - 1], "near")) {
            model->geometry |= GEOMETRY_SOURCE_NEAR;
        } else if (!ci_strcmp(tok[ntok - 1], "far")) {
            model->geometry &= ~GEOMETRY_SOURCE_NEAR;
        } else {
            parse_error(
                    "Valid qualifiers for Psource are \"near\" or \"far\"\n");
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
    }
    return KWFUNC_SUCCESS;
} /* kwfunc_Psource() */


/***********************************************************
* static int kwfunc_Ptan(
*   char *tok[],
*   const int ntok,
*   model_t *model,
*   simplex_t *simplex)
*
* keyword:
*   Ptan - set tangent pressure level for limb geometry
*
* syntax:
*   Ptan pressure unit [dpressure unit]
************************************************************/

static int kwfunc_Ptan(
        char *tok[],
        const int ntok,
        model_t *model,
        simplex_t *simplex)
{
    if (ntok != 3 && ntok != 5) {
        parse_error("Expected \"Ptan pressure unit [dpressure unit]\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (model->geometry & GEOMETRY_ZTAN_USER_DEFINED) {
        parse_error("Ptan and ztan cannot both be user-defined.\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    /*
     * Zero is accepted as input here, but only so we can give
     * an informative error message for that specific case.
     */
    if (get_nonneg_dbl_val(
                &(model->Ptan),
                tok[1],
                tok[2],
                UGROUP_PRESSURE,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (model->Ptan == 0.0) {
        parse_error("Cannot set Ptan to zero.  Instead, use ztan to\n");
        parse_error("\tlocate the tangent level outside the atmosphere\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    model->Ptan_unitnum = get_unitnum(tok[2]);
    model->geometry |= GEOMETRY_PTAN_USER_DEFINED;
    if (ntok == 5 && simplex != NULL) {
        double dPtan;
        if (get_pos_dbl_val(
                    &dPtan,
                    tok[3],
                    tok[4],
                    UGROUP_PRESSURE,
                    1)) {
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
        if (add_var_to_simplex(
                    simplex,
                    "Ptan",
                    &(model->Ptan),
                    model->Ptan,
                    dPtan,
                    model->Ptan_unitnum,
                    MAPPING_NONE)) {
            return KWFUNC_FATAL_ERROR;
        }
    }
    return KWFUNC_SUCCESS;
} /* kwfunc_Ptan() */


/***********************************************************
* static int kwfunc_PTmode(
*   char *tok[],
*   const int ntok,
*   model_t *model)
*
* keyword:
*   PTmode - set P, T operating mode
*
* syntax:
*   PTmode midpoint | ((P | Pbase | dP) [T | Tbase]) or
*   PTmode extend (isothermal | dry_adiabatic | moist_adiabatic |
*   environmental)
************************************************************/

static int kwfunc_PTmode(
        char *tok[],
        const int ntok,
        model_t *model)
{
    int default_Tmode;
    if (ntok < 2) {
        parse_error("Expected\n");
        parse_error("\t  \"PTmode midpoint |"
                " ((P | Pbase | dP) [T | Tbase])\"\n");
        parse_error("\tor\n");
        parse_error("\t  \"PTmode extend (isothermal | dry_adiabatic\n");
        parse_error("\t      | moist_adiabatic | environmental)\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (!ci_strcmp(tok[1], "extend")) {
        if (ntok < 3) {
            parse_error("Expected \"PTmode extend"
                    " (isothermal | dry_adiabatic\n");
            parse_error("\t           | moist_adiabatic | environmental)\"\n");
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
        if (!ci_strcmp(tok[2], "isothermal")) {
            model->PTmode &= ~PTMODE_EXTEND_BITS;
            model->PTmode |=  PTMODE_EXTEND_ISOTHERMAL;
        } else if (!ci_strcmp(tok[2], "dry_adiabatic")) {
            model->PTmode &= ~PTMODE_EXTEND_BITS;
            model->PTmode |=  PTMODE_EXTEND_DRY_ADIABATIC;
        } else if (!ci_strcmp(tok[2], "moist_adiabatic")) {
            model->PTmode &= ~PTMODE_EXTEND_BITS;
            model->PTmode |=  PTMODE_EXTEND_MOIST_ADIABATIC;
        } else if (!ci_strcmp(tok[2], "environmental")) {
            model->PTmode &= ~PTMODE_EXTEND_BITS;
            model->PTmode |=  PTMODE_EXTEND_ENVIRONMENTAL;
        } else {
            parse_error("Expected \"PTmode extend"
                    " (isothermal | dry_adiabatic\n");
            parse_error("\t           | moist_adiabatic | environmental)\"\n");
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
        model->PTmode |= PTMODE_EXTEND_USER_DEFINED;
    } else {
        if (model->nlayers != 0) {
            parse_error(
                    "PTmode specification must precede"
                    " all layer specifications\n");
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
        if (model->PTmode & (PTMODE_PMODES | PTMODE_TMODES)) {
            parse_error("Redundant PTmode statement\n");
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
        if (!ci_strcmp(tok[1], "midpoint")) {
            model->PTmode |= PTMODE_MIDPOINT; /* P, T */
            return KWFUNC_SUCCESS;
        } else if (!ci_strcmp(tok[1], "P")) {
            model->PTmode |= PTMODE_P;
            default_Tmode  = PTMODE_T;
        } else if (!ci_strcmp(tok[1], "dp")) {
            model->PTmode |= PTMODE_DP;
            default_Tmode  = PTMODE_TBASE;
        } else if (!ci_strcmp(tok[1], "Pbase")) {
            model->PTmode |= PTMODE_PBASE;
            default_Tmode  = PTMODE_TBASE;
        } else {
            parse_error("Expected\n");
            parse_error("\t  \"PTmode midpoint |"
                    " ((P | Pbase | dP) [T | Tbase])\"\n");
            parse_error("\tor\n");
            parse_error("\t  \"PTmode extend (isothermal | dry_adiabatic)\"\n");
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
        if (ntok < 3) { /* no T mode specified, use default */
            model->PTmode |= default_Tmode;
        } else {
            if (!ci_strcmp(tok[2], "T")) {
                model->PTmode |= PTMODE_T;
            } else if (!ci_strcmp(tok[2], "Tbase")) {
                model->PTmode |= PTMODE_TBASE;
            } else {
                parse_error(
                        "Expected \"PTmode midpoint |"
                        " ((P | Pbase | dP) [T | Tbase])\"\n");
                parse_error( "\t  (%s?)\n", tok[1]);
                return KWFUNC_SYNTAX_ERROR_CONTINUE;
            }
        }
    }
    return KWFUNC_SUCCESS;
} /* kwfunc_PTmode() */


/***********************************************************
* static int kwfunc_R0(
*   char *tok[],
*   const int ntok,
*   model_t *model)
*
* keyword:
*   R0 - set planetary radius for spherical geometry
*
* syntax:
*   R0 radius unit
************************************************************/

static int kwfunc_R0(
        char *tok[],
        const int ntok,
        model_t *model)
{
    if (ntok < 3) {
        parse_error("Expected \"R0 radius unit\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_pos_dbl_val(
                &(model->R0),
                tok[1],
                tok[2],
                UGROUP_DELAY_DIST,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    model->R0_unitnum = get_unitnum(tok[2]);
    model->geometry |= GEOMETRY_R0_USER_DEFINED;
    return KWFUNC_SUCCESS;
} /* kwfunc_R0() */


/***********************************************************
* static int kwfunc_refract(
*   char *tok[],
*   const int ntok,
*   model_t *model)
*
* keyword:
*   refract - set refraction mode
*
* syntax:
*   refract none | radio | optical
************************************************************/

static int kwfunc_refract(
        char *tok[],
        const int ntok,
        model_t *model)
{
    if (ntok < 2) {
        parse_error("Expected \"refract none | radio | optical\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (!ci_strcmp(tok[1], "none")) {
        model->geometry &= ~GEOMETRY_REFRACT_BITS;
        model->geometry |=  GEOMETRY_REFRACT_NONE;
    } else if (!ci_strcmp(tok[1], "radio")) {
        model->geometry &= ~GEOMETRY_REFRACT_BITS;
        model->geometry |=  GEOMETRY_REFRACT_RADIO;
    } else if (!ci_strcmp(tok[1], "optical")) {
        model->geometry &= ~GEOMETRY_REFRACT_BITS;
        model->geometry |=  GEOMETRY_REFRACT_OPTICAL;
    } else {
        parse_error("Expected \"refract none | radio | optical\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    model->geometry |= GEOMETRY_REFRACT_USER_DEFINED;
    return KWFUNC_SUCCESS;
}   /* kwfunc_refract() */


/***********************************************************
* static int kwfunc_reverse(
*   char *tok[],
*   const int ntok,
*   model_t *model)
*
* keyword:
*   reverse - reverse order of propagation through layers
*
* syntax:
*   reverse 0|1
************************************************************/

static int kwfunc_reverse(
        char *tok[],
        const int ntok,
        model_t *model)
{
    int reverse;
    if (ntok < 2) {
        parse_error("Expected \"reverse 0|1\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_bool_val(&(reverse), tok[1]))
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    if (reverse)
        model->geometry |= GEOMETRY_REVERSE;
    else
        model->geometry &= ~GEOMETRY_REVERSE;
    return KWFUNC_SUCCESS;
} /* kwfunc_reverse() */


/***********************************************************
* static int kwfunc_RH_offset(
*   char *tok[],
*   const int ntok,
*   model_t *model,
*   simplex_t *simplex)
*
* keyword:
*   RH_offset - set RH_offset value.
*
* syntax:
*   RH_offset offset [doffset]
************************************************************/

static int kwfunc_RH_offset(
        char *tok[],
        const int ntok,
        model_t *model,
        simplex_t *simplex)
{
    double RH_offset;
    if (ntok != 2 && ntok != 3) {
        parse_error("Expected \"RH_offset offset [doffset]\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_RH_value(&RH_offset, tok[1]))
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    /*
     * For compatibility with the default log mapping of fit
     * variable space, RH_offset, which can be positive or
     * negative, is stored as exp(RH_offset).
     */
    if (RH_offset > log(DBL_MAX) || RH_offset < log(DBL_MIN)) {
        parse_error("RH_offset outside valid range.\n");
        return KWFUNC_FATAL_ERROR;
    }
    model->RH_offset_exp = map_variable(RH_offset, MAPPING_EXP);
    if (ntok == 3 && simplex != NULL) {
        double dRH_offset, dRH_offset_exp;
        if (get_RH_value(&dRH_offset, tok[2]))
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        if (dRH_offset <= 0.0) {
            parse_error("dRH_offset must be greater than 0%%\n");
            return KWFUNC_FATAL_ERROR;
        }
        dRH_offset_exp = map_differential(RH_offset, dRH_offset, MAPPING_EXP);
        if (add_var_to_simplex(
                    simplex,
                    "RH_offset",
                    &(model->RH_offset_exp),
                    model->RH_offset_exp,
                    dRH_offset_exp,
                    UNIT_NONE,
                    MAPPING_EXP)) {
            return KWFUNC_FATAL_ERROR;
        }
    }
    model->RH_adj_flag = 1;
    return KWFUNC_SUCCESS;
} /* kwfunc_RH_offset() */


/***********************************************************
* static int kwfunc_RH_scale(
*   char *tok[],
*   const int ntok,
*   model_t *model,
*   simplex_t *simplex)
*
* keyword:
*   RH_scale - set RH_scale value.
*
* syntax:
*   RH_scale scale [dscale]
************************************************************/

static int kwfunc_RH_scale(
        char *tok[],
        const int ntok,
        model_t *model,
        simplex_t *simplex)
{
    if (ntok != 2 && ntok != 3) {
        parse_error("Expected \"RH_scale scale [dscale]\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_nonneg_dbl_val(
                &(model->RH_scale),
                tok[1],
                NULL,
                UGROUP_NONE,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (ntok == 3 && simplex != NULL) {
        double dscale;
        if (get_nonneg_dbl_val(
                    &dscale,
                    tok[2],
                    NULL,
                    UGROUP_NONE,
                    1)) {
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
        if (add_var_to_simplex(
                    simplex,
                    "RH_scale",
                    &(model->RH_scale),
                    model->RH_scale,
                    dscale,
                    UNIT_NONE,
                    MAPPING_NONE)) {
            return KWFUNC_FATAL_ERROR;
        }
    }
    model->RH_adj_flag = 1;
    return KWFUNC_SUCCESS;
} /* kwfunc_RH_scale() */


/***********************************************************
* static int kwfunc_runtime(
*   char *tok[],
*   const int ntok,
*   model_t *model)
*
* keyword:
*   runtime - log run times and report in output
*
* syntax:
*   runtime 0|1
************************************************************/

static int kwfunc_runtime(
        char *tok[],
        const int ntok,
        model_t *model)
{
    if (ntok < 2) {
        parse_error("Expected \"runtime 0|1\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_bool_val(&(model->log_runtimes), tok[1]))
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    return KWFUNC_SUCCESS;
} /* kwfunc_runtime() */


/***********************************************************
* static int kwfunc_rx_gain_factor(
*   char *tok[],
*   const int ntok,
*   model_t *model,
*   simplex_t *simplex)
*
* keyword:
*   rx_gain_factor - set receiver gain correction factor
*
* syntax:
*   rx_gain_factor factor
************************************************************/

static int kwfunc_rx_gain_factor(
        char *tok[],
        const int ntok,
        model_t *model,
        simplex_t *simplex)
{
    if (ntok != 2 && ntok != 3) {
        parse_error("Expected \"rx_gain_factor factor [dfactor]\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_nonneg_dbl_val(
                &(model->rx_gain_factor),
                tok[1],
                NULL,
                UGROUP_NONE,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (ntok == 3 && simplex != NULL) {
        double dgain;
        if (get_pos_dbl_val(
                    &dgain,
                    tok[2],
                    NULL,
                    UGROUP_NONE,
                    1)) {
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
        if (add_var_to_simplex(
                    simplex,
                    "rx_gain_factor",
                    &(model->rx_gain_factor),
                    model->rx_gain_factor,
                    dgain,
                    UNIT_NONE,
                    MAPPING_NONE)) {
            return KWFUNC_FATAL_ERROR;
        }
    }
    return KWFUNC_SUCCESS;
} /* kwfunc_rx_gain_factor() */


/***********************************************************
* static int kwfunc_sec_za(
*   char *tok[],
*   const int ntok,
*   model_t *model,
*   simplex_t *simplex)
*
* keyword:
*   sec_za - set secant(za)
*   
* syntax:
*   sec_za value [dvalue]
************************************************************/

static int kwfunc_sec_za(
        char *tok[],
        const int ntok,
        model_t *model,
        simplex_t *simplex)
{
    if (ntok != 2 && ntok != 3) {
        parse_error("Expected \"sec_za value [dvalue]\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_nonneg_dbl_val(
                &(model->sec_za),
                tok[1],
                NULL,
                UGROUP_NONE,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    model->geometry &= ~GEOMETRY_ZA_USER_DEFINED;
    model->geometry |=  GEOMETRY_SEC_ZA_USER_DEFINED;
    if (ntok == 3 && simplex != NULL) {
        double dsec_za;
        if (get_pos_dbl_val(
                    &dsec_za,
                    tok[2],
                    NULL,
                    UGROUP_NONE,
                    1)) {
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
        if (add_var_to_simplex(
                    simplex,
                    "sec_za",
                    &(model->sec_za),
                    model->sec_za,
                    dsec_za,
                    UNIT_NONE,
                    MAPPING_NONE)) {
            return KWFUNC_FATAL_ERROR;
        }
    }
    return KWFUNC_SUCCESS;
} /* kwfunc_sec_za() */


/***********************************************************
* static int kwfunc_selfbroad_vmr_tol(
*   char *tok[],
*   const int ntok,
*   model_t *model)
*
* keywords:
*   selfbroad_vmr_tol - For a given molecule, if the actual
*     mixing ratio is within selfbroad_vmr_tol of the
*     default mixing ratio, self-broadening will be computed
*     at the default mixing ratio.  This behavior optimizes
*     dcache performance when self-broadening is negligible.
*
* syntax:
*   selfbroad_vmr_tol delta
************************************************************/

static int kwfunc_selfbroad_vmr_tol(
        char *tok[],
        const int ntok,
        model_t *model)
{
    if (ntok != 2) {
        parse_error("Expected \"selfbroad_vmr_tol delta\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_nonneg_dbl_val(
                &(model->selfbroad_vmr_tol),
                tok[1],
                NULL,
                UGROUP_NONE,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    return KWFUNC_SUCCESS;
} /* kwfunc_selfbroad_vmr_tol() */


/***********************************************************
* static int kwfunc_simplex_log(
*   char *tok[],
*   const int ntok,
*   simplex_t *simplex)
*
* keyword:
*   simplex_log - use logarithmic mapping of simplex
*     coordinates              
* syntax:
*   simplex_log 0|1
************************************************************/

static int kwfunc_simplex_log(
        char *tok[],
        const int ntok,
        simplex_t *simplex)
{
    if (ntok < 2) {
        parse_error("Expected \"simplex_log 0|1\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (simplex->n) {
        parse_error(
                "Any simplex_log statement must precede the definition of\n");
        parse_error(
                "\tany differentiation or fit variable.\n");
        return KWFUNC_SYNTAX_ERROR_STOP;
    }
    if (get_bool_val(&(simplex->logarithmic), tok[1]))
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    return KWFUNC_SUCCESS;
} /* kwfunc_simplex_log() */


/***********************************************************
* static int kwfunc_T(
*   char *tok[],
*   const int ntok,
*   model_t *model,
*   simplex_t *simplex,
*   int lnum)
*
* keyword:
*   T - set layer midpoint temperature
*
* syntax:
*   T temperature unit [dtemperature unit]
************************************************************/

static int kwfunc_T(
        char *tok[],
        const int ntok,
        model_t *model,
        simplex_t *simplex,
        int lnum)
{
    if (lnum < 0 || lnum >= model->nlayers) {
        parse_error("Attempted to set T for non-existent layer.\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if ((model->PTmode & PTMODE_TMODES) && !(model->PTmode & PTMODE_T)) {
        parse_error(
                "T is inconsistent with a prior layer or PTmode setting.\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (ntok != 3 && ntok != 5) {
        parse_error("Expected \"T temperature unit [dtemperature unit]\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    model->PTmode |= PTMODE_T;
    if (get_dbl_val(
                &(model->layer[lnum]->T),
                tok[1],
                tok[2],
                UGROUP_TEMPERATURE,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (model->layer[lnum]->T < 0.0) {
        parse_error(
                "negative absolute temperature: %.4g K\n",
                model->layer[lnum]->T);
    }
    model->layer[lnum]->T_unitnum = get_unitnum(tok[2]);
    if (ntok == 5 && simplex != NULL) {
        double dT;
        if (get_pos_dbl_val(
                    &dT,
                    tok[3],
                    tok[4],
                    UGROUP_TEMPERATURE,
                    1)) {
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
        snprintf(VARNAME, sizeof(VARNAME), "T (layer %d)", lnum);
        if (add_var_to_simplex(
                    simplex,
                    VARNAME,
                    &(model->layer[lnum]->T),
                    model->layer[lnum]->T,
                    dT,
                    model->layer[lnum]->T_unitnum,
                    MAPPING_NONE)) {
            return KWFUNC_FATAL_ERROR;
        }
    }
    return KWFUNC_SUCCESS;
} /* kwfunc_T() */


/***********************************************************
* static int kwfunc_T0(
*   char *tok[],
*   const int ntok,
*   model_t *model,
*   simplex_t *simplex)
*
* keyword:
*   T0 - set background temperature
*
* syntax:
*   T0 temperature unit [dtemperature unit]
************************************************************/

static int kwfunc_T0(
        char *tok[],
        const int ntok,
        model_t *model,
        simplex_t *simplex)
{
    double dT0;
    if (ntok != 3 && ntok != 5) {
        parse_error("Expected \"T0 temperature unit [dtemperature unit]\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_dbl_val(
                &(model->T0),
                tok[1],
                tok[2],
                UGROUP_TEMPERATURE,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (model->T0 < 0.0) {
        parse_error("negative absolute temperature: %.4g K\n", model->T0);
    }
    model->T0_unitnum = get_unitnum(tok[2]);
    if (ntok == 5 && simplex != NULL) {
        if (get_pos_dbl_val(
                    &dT0,
                    tok[3],
                    tok[4],
                    UGROUP_TEMPERATURE,
                    1)) {
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
        if (add_var_to_simplex(
                    simplex,
                    "T0",
                    &(model->T0),
                    model->T0,
                    dT0,
                    model->T0_unitnum,
                    MAPPING_NONE)) {
            return KWFUNC_FATAL_ERROR;
        }
    }
    return KWFUNC_SUCCESS;
} /* kwfunc_T0() */


/***********************************************************
* static int kwfunc_Tbase(
*   char *tok[],
*   const int ntok,
*   model_t *model,
*   simplex_t *simplex,
*   int lnum)
*
* keyword:
*   Tbase - set layer base temperature
*
* syntax:
*   Tbase temperature unit [dtemperature unit]
************************************************************/

static int kwfunc_Tbase(
        char *tok[],
        const int ntok,
        model_t *model,
        simplex_t *simplex,
        int lnum)
{
    if (lnum < 0 || lnum >= model->nlayers) {
        parse_error("Attempted to set Tbase for non-existent layer.\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if ((model->PTmode & PTMODE_TMODES) && !(model->PTmode & PTMODE_TBASE)) {
        parse_error(
                "Tbase is inconsistent with a prior"
                " layer or PTmode setting.\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    model->PTmode |= PTMODE_TBASE;
    if (ntok != 3 && ntok != 5) {
        parse_error(
                "Expected \"Tbase temperature unit [dtemperature unit]\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_dbl_val(
                &(model->layer[lnum]->Tbase),
                tok[1],
                tok[2],
                UGROUP_TEMPERATURE,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (model->layer[lnum]->Tbase < 0.0) {
        parse_error(
                "negative absolute temperature: %.4g K\n",
                model->layer[lnum]->Tbase);
    }
    model->layer[lnum]->T_unitnum = get_unitnum(tok[2]);
    if (ntok == 5 && simplex != NULL) {
        double dTbase;
        if (get_pos_dbl_val(
                    &dTbase,
                    tok[3],
                    tok[4],
                    UGROUP_TEMPERATURE,
                    1)) {
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
        snprintf(VARNAME, sizeof(VARNAME), "Tbase (layer %d)", lnum);
        if (add_var_to_simplex(
                    simplex,
                    VARNAME,
                    &(model->layer[lnum]->Tbase),
                    model->layer[lnum]->Tbase,
                    dTbase,
                    model->layer[lnum]->T_unitnum,
                    MAPPING_NONE)) {
            return KWFUNC_FATAL_ERROR;
        }
    }
    return KWFUNC_SUCCESS;
} /* kwfunc_Tbase() */


/***********************************************************
* static int kwfunc_Tref(
*   char *tok[],
*   const int ntok,
*   model_t *model,
*   simplex_t *simplex)
*
* keyword:
*   Tref - set reference temperature for Y factor
*
* syntax:
*   Tref temperature unit [dtemperature unit]
************************************************************/

static int kwfunc_Tref(
        char *tok[],
        const int ntok,
        model_t *model,
        simplex_t *simplex)
{
    double dTref;
    if (ntok != 3 && ntok != 5) {
        parse_error(
                "Expected \"Tref temperature unit [dtemperature unit]\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_dbl_val(
                &(model->Tref),
                tok[1],
                tok[2],
                UGROUP_TEMPERATURE,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (model->Tref < 0.0) {
        parse_error("negative absolute temperature: %.4g K\n", model->Tref);
    }
    model->Tref_unitnum = get_unitnum(tok[2]);
    if (ntok == 5 && simplex != NULL) {
        if (get_pos_dbl_val(
                    &dTref,
                    tok[3],
                    tok[4],
                    UGROUP_TEMPERATURE,
                    1)) {
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
        if (add_var_to_simplex(
                    simplex,
                    "Tref",
                    &(model->Tref),
                    model->Tref,
                    dTref,
                    model->Tref_unitnum,
                    MAPPING_NONE)) {
            return KWFUNC_FATAL_ERROR;
        }
    }
    return KWFUNC_SUCCESS;
} /* kwfunc_Tref() */


/***********************************************************
* static int kwfunc_Trx(
*   char *tok[],
*   const int ntok,
*   model_t *model,
*   simplex_t *simplex)
*
* keyword:
*   Trx - set receiver noise temperature
*
* syntax:
*   Trx temperature unit [dtemperature unit]
************************************************************/

static int kwfunc_Trx(
        char *tok[],
        const int ntok,
        model_t *model,
        simplex_t *simplex)
{
    double dTrx;
    if (ntok != 3 && ntok != 5) {
        parse_error(
                "Expected \"Trx temperature unit [dtemperature unit]\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_dbl_val(
                &(model->Trx),
                tok[1],
                tok[2],
                UGROUP_TEMPERATURE,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (model->Trx < 0.0) {
        parse_error("negative absolute temperature: %.4g K\n", model->Trx);
    }
    model->Trx_unitnum = get_unitnum(tok[2]);
    if (ntok == 5 && simplex != NULL) {
        if (get_pos_dbl_val(
                    &dTrx,
                    tok[3],
                    tok[4],
                    UGROUP_TEMPERATURE,
                    1)) {
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
        if (add_var_to_simplex(
                    simplex,
                    "Trx",
                    &(model->Trx),
                    model->Trx,
                    dTrx,
                    model->Trx_unitnum,
                    MAPPING_NONE)) {
            return KWFUNC_FATAL_ERROR;
        }
    }
    return KWFUNC_SUCCESS;
} /* kwfunc_Trx() */


/***********************************************************
* static int kwfunc_tol(
*   char *tok[],
*   const int ntok,
*   model_t *model)
*
* keywords:
*   tol, tolerance - set fast line sum tolerance
*
* syntax:
*   tol delta
************************************************************/

static int kwfunc_tol(
        char *tok[],
        const int ntok,
        model_t *model)
{
    if (ntok != 2) {
        parse_error("Expected \"tol delta\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_nonneg_dbl_val(
                &(model->tol),
                tok[1],
                NULL,
                UGROUP_NONE,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    return KWFUNC_SUCCESS;
} /* kwfunc_tol() */


/***********************************************************
* static int kwfunc_z0(
*   char *tok[],
*   const int ntok,
*   model_t *model)
*
* keyword:
*   z0 - set altitude of base layer
*
* syntax:
*   z0 altitude unit
************************************************************/

static int kwfunc_z0(
        char *tok[],
        const int ntok,
        model_t *model)
{
    if (ntok < 3) {
        parse_error("Expected \"z0 altitude unit\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_dbl_val(
                &(model->z0),
                tok[1],
                tok[2],
                UGROUP_DELAY_DIST,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    model->z0_unitnum = get_unitnum(tok[2]);
    model->geometry |= GEOMETRY_Z0_USER_DEFINED;
    return KWFUNC_SUCCESS;
} /* kwfunc_z0() */


/***********************************************************
* static int kwfunc_za(
*   char *tok[],
*   const int ntok,
*   model_t *model,
*   simplex_t *simplex)
*
* keyword:
*   za - set zenith angle
*
* syntax:
*   za angle deg|rad [dangle deg|rad]
************************************************************/

static int kwfunc_za(
        char *tok[],
        const int ntok,
        model_t *model,
        simplex_t *simplex)
{
    if (ntok != 3 && ntok != 5) {
        parse_error("Expected \"za angle deg|rad [dangle deg|rad]\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_nonneg_dbl_val(
                &(model->za),
                tok[1],
                tok[2],
                UGROUP_ANGLE,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    model->za_unitnum = get_unitnum(tok[2]);
    model->geometry &= ~GEOMETRY_SEC_ZA_USER_DEFINED;
    model->geometry |=  GEOMETRY_ZA_USER_DEFINED;
    if (ntok == 5 && simplex != NULL) {
        double dza;
        if (get_pos_dbl_val(
                    &dza,
                    tok[3],
                    tok[4],
                    UGROUP_ANGLE,
                    1)) {
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
        if (add_var_to_simplex(
                    simplex,
                    "za",
                    &(model->za),
                    model->za,
                    dza,
                    model->za_unitnum,
                    MAPPING_NONE)) {
            return KWFUNC_FATAL_ERROR;
        }
    }
    return KWFUNC_SUCCESS;
} /* kwfunc_za() */


/***********************************************************
* static int kwfunc_zobs(
*   char *tok[],
*   const int ntok,
*   model_t *model,
*   simplex_t *simplex)
*
* keyword:
*   zobs - set observing level by altitude.
*
* syntax:
*   zobs altitude unit [daltitude unit]
************************************************************/

static int kwfunc_zobs(
        char *tok[],
        const int ntok,
        model_t *model,
        simplex_t *simplex)
{
    if (ntok != 3 && ntok != 5) {
        parse_error("Expected \"zobs altitude unit [daltitude unit]\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (model->geometry & GEOMETRY_POBS_USER_DEFINED) {
        parse_error("Pobs and zobs cannot both be user-defined.\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_dbl_val(
                &(model->zobs),
                tok[1],
                tok[2],
                UGROUP_DELAY_DIST,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    model->zobs_unitnum = get_unitnum(tok[2]);
    model->geometry |= GEOMETRY_ZOBS_USER_DEFINED;
    if (ntok == 5 && simplex != NULL) {
        double dzobs;
        if (get_pos_dbl_val(
                    &dzobs,
                    tok[3],
                    tok[4],
                    UGROUP_DELAY_DIST,
                    1)) {
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
        if (add_var_to_simplex(
                    simplex,
                    "zobs",
                    &(model->zobs),
                    model->zobs,
                    dzobs,
                    model->zobs_unitnum,
                    MAPPING_NONE)) {
            return KWFUNC_FATAL_ERROR;
        }
    }
    return KWFUNC_SUCCESS;
} /* kwfunc_zobs() */


/***********************************************************
* static int kwfunc_zsource(
*   char *tok[],
*   const int ntok,
*   model_t *model,
*   simplex_t *simplex)
*
* keyword:
*   zsource - set source level by altitude.
*
*   For paths below the astronomical horizon, the source
*   level is intersected on both the near and far sides
*   of the tangent point.  The optional "near" or "far"
*   arguments select which point should be the path
*   endpoint.  The default is "far".
*
* syntax:
*   zsource altitude unit [daltitude unit] [near | far]
************************************************************/

static int kwfunc_zsource(
        char *tok[],
        const int ntok,
        model_t *model,
        simplex_t *simplex)
{
    if (ntok < 3 || ntok > 6 || (ntok == 5 && !is_unit(tok[4]))) {
        parse_error(
                "Expected \"zsource altitude unit [daltitude unit]"
                " [near | far]\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (model->geometry & GEOMETRY_PSOURCE_USER_DEFINED) {
        parse_error("Psource and zsource cannot both be user-defined.\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_dbl_val(
                &(model->zsource),
                tok[1],
                tok[2],
                UGROUP_DELAY_DIST,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    model->zsource_unitnum = get_unitnum(tok[2]);
    model->geometry |= GEOMETRY_ZSOURCE_USER_DEFINED;
    if (ntok >= 5 && is_unit(tok[4]) && simplex != NULL) {
        double dzsource;
        if (get_pos_dbl_val(
                    &dzsource,
                    tok[3],
                    tok[4],
                    UGROUP_DELAY_DIST,
                    1)) {
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
        if (add_var_to_simplex(
                    simplex,
                    "zsource",
                    &(model->zsource),
                    model->zsource,
                    dzsource,
                    model->zsource_unitnum,
                    MAPPING_NONE)) {
            return KWFUNC_FATAL_ERROR;
        }
    }
    if (ntok == 4 || ntok == 6) {
        if (!ci_strcmp(tok[ntok - 1], "near")) {
            model->geometry |= GEOMETRY_SOURCE_NEAR;
        } else if (!ci_strcmp(tok[ntok - 1], "far")) {
            model->geometry &= ~GEOMETRY_SOURCE_NEAR;
        } else {
            parse_error(
                    "Valid qualifiers for zsource are \"near\" or \"far\"\n");
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
    }
    return KWFUNC_SUCCESS;
} /* kwfunc_zsource() */


/***********************************************************
* static int kwfunc_ztan(
*   char *tok[],
*   const int ntok,
*   model_t *model,
*   simplex_t *simplex)
*
* keyword:
*   ztan - set tangent level altitude for limb geometry
*
* syntax:
*   ztan altitude unit [daltitude unit]
************************************************************/

static int kwfunc_ztan(
        char *tok[],
        const int ntok,
        model_t *model,
        simplex_t *simplex)
{
    if (ntok != 3 && ntok != 5) {
        parse_error("Expected \"ztan altitude unit [daltitude unit]\"\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (model->geometry & GEOMETRY_PTAN_USER_DEFINED) {
        parse_error("Ptan and ztan cannot both be user-defined.\n");
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    if (get_dbl_val(
                &(model->ztan),
                tok[1],
                tok[2],
                UGROUP_DELAY_DIST,
                0)) {
        return KWFUNC_SYNTAX_ERROR_CONTINUE;
    }
    model->ztan_unitnum = get_unitnum(tok[2]);
    model->geometry |= GEOMETRY_ZTAN_USER_DEFINED;
    if (ntok == 5 && simplex != NULL) {
        double dztan;
        if (get_pos_dbl_val(
                    &dztan,
                    tok[3],
                    tok[4],
                    UGROUP_DELAY_DIST,
                    1)) {
            return KWFUNC_SYNTAX_ERROR_CONTINUE;
        }
        if (add_var_to_simplex(
                    simplex,
                    "ztan",
                    &(model->ztan),
                    model->ztan,
                    dztan,
                    model->ztan_unitnum,
                    MAPPING_NONE)) {
            return KWFUNC_FATAL_ERROR;
        }
    }
    return KWFUNC_SUCCESS;
} /* kwfunc_ztan() */
