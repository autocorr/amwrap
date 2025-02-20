/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* fit.c                      S. Paine rev. 2024 September 29
*
* Model fits to spectra, by downhill simplex minimization of
* a user-specified goodness of fit estimator.
************************************************************/

#include <ctype.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "abscoeff.h"
#include "am_alloc.h"
#include "am_sysdep.h"
#include "am_types.h"
#include "config.h"
#include "dcache.h"
#include "errlog.h"
#include "fileops.h"
#include "fit.h"
#include "jacobian.h"
#include "kcache.h"
#include "layer.h"
#include "model.h"
#include "nscale.h"
#include "output.h"
#include "simplex.h"
#include "spectra.h"
#include "units.h"


enum {
    FIT_DATA_LBUFSIZE    = 1024, /* line buffer size for fit data files */
    FIT_MAX_DELIMSTRSIZE = 256,  /* max size of fit data col delim string */
    /*
     * If FIT_MAX_INPUT_FIELDS is changed, change the sscanf() args in
     * read_fit_data() also.
     */
    FIT_MAX_INPUT_FIELDS = 4,
    FIT_MAX_FMTSTRSIZE   = 256   /* max size of format string for reading data */
};

/*
 * Special strings that can be embedded in data files
 */
static const char EOB_STRING[]        = "end";     /* end of data block    */
static const char PARAM_CTRL_STRING[] = "set";     /* set config parameter */
static const char RESTART_STRING[]    = "restart"; /* new fit on same data */
/*
 * Environment variable names
 */
static const char FIT_INPUT_PATH_ENV_STRING[]  = "AM_FIT_INPUT_PATH";
static const char FIT_OUTPUT_PATH_ENV_STRING[] = "AM_FIT_OUTPUT_PATH";
/*
 * Fit output file extensions
 */
static const char MODEL_CONFIG_EXT[AM_FIT_EXTSIZE]      = ".amc";
static const char MODEL_SPECTRUM_EXT[AM_FIT_EXTSIZE]    = ".ams";
static const char RESIDUAL_SPECTRUM_EXT[AM_FIT_EXTSIZE] = ".amr";
/*
 * String to loop back to the top of a data file list and
 * restart.
 */
static const char RESTART_FNAME[] = "restart";

struct fit_estimator_type_tabentry fit_estimator_type[] = {
    {"none"},
    {"absolute_residuals"},
    {"squared_residuals"},
    {"logarithmic"},
    {"END"},
};

static int   close_fit_data_file(fit_data_t*);
static char *fit_indir(void);
static char *fit_outdir(void);
static int   init_fit_kcache(model_t*, simplex_t*);
static int   open_fit_data_file(fit_data_t*);
enum read_fit_data_retval {
    READ_FIT_DATA_EOF,
    READ_FIT_DATA_EOB,
    READ_FIT_DATA_RESTART_FILE_LIST,
    READ_FIT_DATA_ERROR
};
static int   read_fit_data(model_t*, fit_data_t*, simplex_t*);
static void  report_fit_status(FILE*, fit_data_t*, simplex_t*);
static void  compute_variance(fit_data_t*, simplex_t*);
static int   update_simplex_vertex(
        model_t*,
        model_t*,
        fit_data_t*,
        simplex_t*,
        double*,
        double*);
static int write_fit_output_files(model_t*, fit_data_t*, simplex_t*);
static void write_fit_summary_line(fit_data_t*, simplex_t*);
static void update_estimated_residuals(fit_data_t*);


/***********************************************************
* int fit(
*         model_t    *model,
*         model_t    *lmodel,
*         fit_data_t *fit_data,
*         simplex_t  *simplex)
*
* Purpose:
*   Performs model fits to spectra, by downhill simplex
*   minimization of a goodness of fit estimator.  Good
*   references on the algorithm are:
*
*   J. A. Nelder and R. Mead, "A simplex method for function
*   minimization" Computer Journal, 7:308 (1965).  The full
*   text is available at http://comjnl.oxfordjournals.org/.
*
*   Press, et al., Numerical Recipes in C, 2nd ed.,
*   section 10.4. (Cambridge University Press, 1992)
*
* Arguments:
*   model_t    *model    - pointer to model structure
*   model_t    *lmodel   - pointer to model structure
*                          containing scalar data from the
*                          prior computation of *model, if
*                          any.
*   fit_data_t *fit_data - pointer to fit_data structure
*   simplex_t  *simplex  - pointer to simplex data structure
*
* Return:
*   0 if successful, 1 otherwise
************************************************************/

int fit(
        model_t    *model,
        model_t    *lmodel,
        fit_data_t *fit_data,
        simplex_t  *simplex)
{
    double fit_tstart;

    /*
     * If no fit variables have been specified, this is a 
     * null fit.  The model and residuals will be computed,
     * but nothing will be varied.
     */
    if (simplex->n == 0) {
        create_null_simplex(simplex);
        fit_data->max_restarts = 0;
        errlog(34, 0);
    }
    if (model->nkcache != 0) {
        if (init_fit_kcache(model, simplex))
            return 1;
    }
    /*
     * OUTER LOOP OVER ALL DATA FILES/BLOCKS TO BE FIT.
     */
    fit_data->open_filenum = 0;
    fit_data->blocknum     = 0;
    while (fit_data->open_filenum < fit_data->nfiles) {
        gridsize_t npts_last;
        unsigned int i, j;
        int read_status;
        /*
         * Clear the cache logs, and start the timer for this
         * fit.
         */
        dcache_log(DCACHE_CLEAR_LOG);
        kcache_log(KCACHE_CLEAR_LOG);
        fit_tstart  = am_timer(0.0);
        npts_last   = fit_data->npts;
        read_status = read_fit_data(model, fit_data, simplex);
        if (read_status == READ_FIT_DATA_ERROR)
            return 1;
        if (read_status == READ_FIT_DATA_RESTART_FILE_LIST) {
            fit_data->open_filenum = 0;
            fit_data->blocknum = 0;
            continue;
        }
        if (fit_data->npts == 0) {
            if (read_status == READ_FIT_DATA_EOF) {
                ++fit_data->open_filenum;
                fit_data->blocknum = 0;
            } else {
                ++fit_data->blocknum;
            }
            continue;
        }
        /*
         * Reset the estimated residuals array if this is the
         * first fit data block, or if the size of the fit data
         * block changed from last time.  Otherwise, update the
         * estimated residuals array.
         */
        if (fit_data->npts != npts_last) {
            reset_estimated_residuals(fit_data);
        } else {
            if (fit_data->res_track_gain >= 0.0)
                update_estimated_residuals(fit_data);
        }
        /*
         * SIMPLEX INITIALIZATION.  Reset the simplex and initialize
         * the fit estimator E for each of the simplex vertices.
         */
        simplex->restart_count = 0;
        reset_simplex_vertices(simplex);
        set_active_outputs(fit_data->data_type);
        for (i = 0; i <= simplex->n; ++i) {
            if (update_simplex_vertex(
                        model,
                        lmodel,
                        fit_data,
                        simplex,
                        simplex->vertex[i],
                        &simplex->E[i])) {
                return 1;
            }
        }
        /*
         * DOWNHILL SIMPLEX MINIMIZATION LOOP.  Iterate until
         * convergence, or until the maximum iteration count is
         * reached.
         */
        for (simplex->iter = 0; simplex->iter < fit_data->max_iter;
            ++simplex->iter) {
            int nhigher;
            /*
             * Identify high and low simplex vertices (i.e. the
             * vertices with the higest and lowest fit
             * estimator.)
             */
            simplex->ilow = simplex->ihigh = 0;
            for (i = 1; i <= simplex->n; ++i) {
                if (simplex->E[i] > simplex->E[simplex->ihigh]) {
                    simplex->ihigh = i;
                } else if (simplex->E[i] < simplex->E[simplex->ilow]) {
                    simplex->ilow  = i;
                }
            }
            /*
             * If in verbose mode, update the run time and write
             * fit status to stderr.
             */
            if (fit_data->output_mode & FIT_OUTPUT_VERBOSE) {
                fit_data->runtime = am_timer(fit_tstart);
                report_fit_status(stderr, fit_data, simplex);
            }
            /*
             * CONVERGENCE TEST.  
             */
            if (simplex_scaled_diameter(simplex) < simplex->tol) {
                double d_c = 0.0;
                /*
                 * No restarts called for, so just break out of the fit.
                 */
                if (fit_data->max_restarts == 0)
                    break;
                /*
                 * If there has been at least one restart,
                 * computed the scaled distance from the last
                 * converged low vertex.
                 */
                if (simplex->restart_count >= 1) {
                    d_c = simplex_scaled_distance(
                            simplex,
                            simplex->vertex[simplex->ilow],
                            simplex->pc);
                    /*
                     * Convergence was consistent after restart,
                     * so break out of the fit.
                     */
                    if (d_c < simplex->tol) {
                        if (fit_data->output_mode & FIT_OUTPUT_VERBOSE) {
                            fprintf(stderr,
                                    "\nconverged after restart "
                                    "%d : dc %g (%g)\n",
                                    simplex->restart_count,
                                    d_c,
                                    simplex->tol);
                        }
                        break;
                    }
                }
                /*
                 * No consistent convergence yet, so increment
                 * the restart counter.  Break out of the fit if
                 * there have been too many restarts.
                 */
                if (++simplex->restart_count > fit_data->max_restarts) {
                    if (fit_data->output_mode & FIT_OUTPUT_VERBOSE) {
                        fprintf(stderr,
                                "\nfailed to converge after restart %d : "
                                "dc %g (%g)\n",
                                simplex->restart_count,
                                d_c,
                                simplex->tol);
                    }
                    break;
                }
                if (fit_data->output_mode & FIT_OUTPUT_VERBOSE) {
                    fprintf(stderr, "\nrestart %d", simplex->restart_count);
                    if (simplex->restart_count > 1) {
                        fprintf(stderr,
                                " : dc %g (%g)\n",
                                d_c,
                                simplex->tol);
                    } else {
                        fprintf(stderr, "\n");
                    }
                }
                /*
                 * Save the converged low vertex for comparison
                 * with the next.
                 */
                for (j = 0; j < simplex->n; ++j)
                    simplex->pc[j] = simplex->vertex[simplex->ilow][j];
                /*
                 * Swap the low vertex into position 0 and reset
                 * the simplex.
                 */
                if (simplex->ilow != 0) {
                    double tmp;
                    tmp = simplex->E[0];
                    simplex->E[0] = simplex->E[simplex->ilow];
                    simplex->E[simplex->ilow] = tmp;
                    for (j = 0; j < simplex->n; ++j) {
                        tmp = simplex->vertex[0][j];
                        simplex->vertex[0][j] =
                                simplex->vertex[simplex->ilow][j];
                        simplex->vertex[simplex->ilow][j] = tmp;
                    }
                    simplex->ilow = 0;
                }
                reset_simplex_vertices(simplex);
                /*
                 * Recompute the model and the fit estimator at
                 * each of the new vertices.
                 */
                for (i = 0; i <= simplex->n; ++i) {
                    if (update_simplex_vertex(
                            model,
                            lmodel,
                            fit_data,
                            simplex,
                            simplex->vertex[i],
                            &simplex->E[i])) {
                        return 1;
                    }
                }
                /*
                 * With the newly-reset simplex, resume the fit
                 * at the top of the fit loop.
                 */
                continue;
            }
            /*
             * Compute the centroid, pbar, of the face of the
             * simplex opposite the highest vertex.
             */
            for (j = 0; j < simplex->n; ++j)
                simplex->pbar[j] = 0.;
            for (i = 0; i <= simplex->n; ++i) {
                if (i == simplex->ihigh)
                    continue;
                for (j = 0; j < simplex->n; ++j)
                    simplex->pbar[j] += simplex->vertex[i][j];
            }
            for (j = 0; j < simplex->n; ++j)
                simplex->pbar[j] /= simplex->n;
            /*
             * Compute the trial vertex p1, found by reflecting
             * the highest vertex through pbar by a factor -1.
             * Also compute its fit estimator, E1.
             */
            for (j = 0; j < simplex->n; ++j)
                simplex->p1[j] = 2. * simplex->pbar[j] -
                        simplex->vertex[simplex->ihigh][j];
            if (update_simplex_vertex(
                    model,
                    lmodel,
                    fit_data,
                    simplex,
                    simplex->p1,
                    &simplex->E1)) {
                return 1;
            }
            /*
             * If p1 is lower the lowest vertex, try extending
             * the reflected point twice as far from pbar.  Call
             * this trial vertex p2, and its corresponding fit
             * estimator E2.
             */
            if (simplex->E1 < simplex->E[simplex->ilow]) {
                for (j = 0; j < simplex->n; ++j)
                    simplex->p2[j] = 3. * simplex->pbar[j] -
                            2. * simplex->vertex[simplex->ihigh][j];
                if (update_simplex_vertex(
                        model,
                        lmodel,
                        fit_data,
                        simplex,
                        simplex->p2,
                        &simplex->E2)) {
                    return 1;
                }
                /*
                 * Then, if p2 is also lower the lowest vertex,
                 * replace the highest vertex with p2.
                 * Otherwise, replace the highest vertex with p1.
                 * Continue to the next iteration.
                 */
                if (simplex->E2 < simplex->E[simplex->ilow]) {
                    for (j = 0; j < simplex->n; ++j)
                        simplex->vertex[simplex->ihigh][j] = simplex->p2[j];
                    simplex->E[simplex->ihigh] = simplex->E2;
                } else {
                    for (j = 0; j < simplex->n; ++j)
                        simplex->vertex[simplex->ihigh][j] = simplex->p1[j];
                    simplex->E[simplex->ihigh] = simplex->E1;
                }
                continue;
            }
            /*
             * p1 was not lower than the lowest vertex.  Count
             * the number of vertices which are higher than p1,
             * not including the highest vertex.
             */
            nhigher = 0;
            for (i = 0; i <= simplex->n; ++i) {
                if (i == simplex->ihigh)
                    continue;
                if (simplex->E[i] > simplex->E1)
                    ++nhigher;
            }
            /*
             * If, on replacing the the highest vertex with p1,
             * p1 would not be the highest vertex of the simplex,
             * then do this replacement and continue to the next
             * iteration.  (This condition avoids a situation
             * where the high vertex simply reflects back and
             * forth through the opposite face on successive
             * iterations.)
             */
            if (nhigher >= 1) {
                for (j = 0; j < simplex->n; ++j)
                    simplex->vertex[simplex->ihigh][j] = simplex->p1[j];
                simplex->E[simplex->ihigh] = simplex->E1;
                continue;
            }
            /*
             * Otherwise, try contracting the simplex, moving the
             * highest vertex halfway towards the centroid of the
             * opposite face.  Before doing this, replace the
             * high vertex with p1, if p1 is lower.
             */
            if (simplex->E[simplex->ihigh] > simplex->E1) {
                for (j = 0; j < simplex->n; ++j)
                    simplex->vertex[simplex->ihigh][j] = simplex->p1[j];
                simplex->E[simplex->ihigh] = simplex->E1;
            }
            /*
             * Now try the contraction, calling the contracted
             * point p2.
             */
            for (j = 0; j < simplex->n; ++j)
                simplex->p2[j] = 0.5 *
                        (simplex->pbar[j] + simplex->vertex[simplex->ihigh][j]);
            if (update_simplex_vertex(
                    model,
                    lmodel,
                    fit_data,
                    simplex,
                    simplex->p2,
                    &simplex->E2)) {
                return 1;
            }
            /*
             * If the contracted vertex p2 is lower than the highest vertex,
             * replace the highest vertex with p2 and continue to the next
             * iteration.
             */
            if (simplex->E2 < simplex->E[simplex->ihigh]) {
                for (j = 0; j < simplex->n; ++j)
                    simplex->vertex[simplex->ihigh][j] = simplex->p2[j];
                simplex->E[simplex->ihigh] = simplex->E2;
                continue;
            }
            /*
             * Otherwise, contract the simplex by moving all vertices, except
             * the lowest one, towards the lowest one.
             */
            for (i = 0; i <= simplex->n; ++i) {
                if (i == simplex->ilow)
                    continue;
                for (j = 0; j < simplex->n; ++j)
                    simplex->vertex[i][j] = 0.5 * (simplex->vertex[i][j] +
                        simplex->vertex[simplex->ilow][j]);
                if (update_simplex_vertex(
                        model,
                        lmodel,
                        fit_data,
                        simplex,
                        simplex->vertex[i],
                        &simplex->E[i])) {
                    return 1;
                }
            }
        }
        /*
         * END OF MINIMIZATION LOOP.  Either the fit converged,
         * or the maximum iteration count was reached.
         *
         * Swap the lowest vertex of the simplex with vertex 0.
         * (Vertex 0 will be the starting point for construction
         * of the initial simplex for the next fit.)
         */
        if (simplex->ilow != 0) {
            double tmp;
            tmp = simplex->E[0];
            simplex->E[0] = simplex->E[simplex->ilow];
            simplex->E[simplex->ilow] = tmp;
            for (j = 0; j < simplex->n; ++j) {
                tmp = simplex->vertex[0][j];
                simplex->vertex[0][j] = simplex->vertex[simplex->ilow][j];
                simplex->vertex[simplex->ilow][j] = tmp;
            }
            simplex->ilow = 0;
        }
        /*
         * Recompute the model and fit data at vertex 0, so that
         * all output data will reflect the parameter set at the
         * best vertex.  Before doing so, call
         * set_active_outputs() to restore computation of all
         * requested outputs.
         *
         * All, or nearly all, of the absorption coefficient
         * computations in this last step will be redundant, and
         * therefore take minimal time.  Because of this,
         * temporarily disable run time logging if it is enabled,
         * so that the reported runtimes will be those of the
         * model computation just prior to this one.  This will
         * make it easier to judge where the computational hot
         * spots are when trying to speed up a fit.
         *
         * If Jacobians are being computed, this is also done
         * while run time logging is turned off.
         */
        set_active_outputs(ALL_OUTPUTS);
        i = model->log_runtimes;
        model->log_runtimes = 0;
        if (compute_jacobians(model, lmodel, simplex))
                return 1;
        if (update_simplex_vertex(
                    model,
                    lmodel,
                    fit_data,
                    simplex,
                    simplex->vertex[0],
                    &simplex->E[0])) {
                return 1;
        }
        model->log_runtimes = i;
        compute_variance(fit_data, simplex);
        fit_data->runtime = am_timer(fit_tstart);
        if (write_fit_output_files(model, fit_data, simplex))
            return 1;
        /*
         * For multiple-file fits, and for fits to multiple data
         * blocks within a file, write a single line summary
         * output to stdout.
         */
        if ((fit_data->nfiles > 1) ||
            (read_status == READ_FIT_DATA_EOB) || 
            ((read_status == READ_FIT_DATA_EOF) && (fit_data->blocknum > 0))
            ) {
            write_fit_summary_line(fit_data, simplex);
        }
        /*
         * Increment filenum and blocknum, as needed, for next fit.
         */
        if (read_status == READ_FIT_DATA_EOF) {
            ++fit_data->open_filenum;
            fit_data->blocknum = 0;
        } else {
            ++fit_data->blocknum;
        }
    }
    /*
     * END OF LOOP OVER DATA FILES/BLOCKS
     */
    return 0;
}   /* fit() */


/***********************************************************
* static int close_fit_data_file(fit_data_t *fit_data)
*
* Purpose:
*   Closes the currently open fit data file.  If the fit
*   data is being read from stdin, stdin is left open, but
*   the fit data file pointer is set to NULL.
*
* Arguments:
*   fit_data_t *fit_data - fit_data structure
*
* Return:
*   0 if OK, EOF otherwise
************************************************************/

static int close_fit_data_file(fit_data_t *fit_data)
{
    int status;

    if ((strcmp(fit_data->filename[fit_data->open_filenum], "stdin") == 0) ||
            (strcmp(fit_data->filename[fit_data->open_filenum], "-") == 0)) {
        status = 0;
    } else {
        status = fclose(fit_data->fp);
    }
    fit_data->fp = NULL;
    return status;
} /* close_fit_data_file() */


/***********************************************************
* char *fit_data_delimiters(const char *delimstr)
*
* Purpose:
*   Sets the delimiter string for parsing fit data columns
*   to delimstr.
*
*   fit_data_delimiters(NULL) just returns a pointer to the
*   current delimiter string.
*
* Arguments:
*   char *delimstr - new fit data delimiter string, or NULL.
*
* Return:
*   pointer to fit_data_delimiters string, or NULL on error.
************************************************************/

char *fit_data_delimiters(const char *fmt)
{
    static char fit_data_delimstr[FIT_MAX_DELIMSTRSIZE] = " \t\n\r";
    /*
     * To handle the end of a string returned by fgets(), newline
     * and carriage return are always included in the delimiter
     * string.
     */
    const char  reqd_delims[] = "\n\r";

    if (fmt != NULL) {
        if (strlen(fmt) == 0) {
            return NULL;
        }
        if ((size_t)snprintf(
                    fit_data_delimstr,
                    sizeof(fit_data_delimstr),
                    "%s%s",
                    fmt,
                    reqd_delims
                    ) > sizeof(fit_data_delimstr)) {
            return NULL;
        }
    }
    return fit_data_delimstr;
}   /* fit_data_delimiters() */


/***********************************************************
* char *fit_data_format(const char *fmt)
*
* Purpose:
*   Sets the format string for reading fit data files equal
*   to fmt.
*
*   fit_data_format(NULL) just returns a pointer to the
*   current format string.
*
* Arguments:
*   char *format - new fit data format string, or NULL.
*
* Return:
*   pointer to fit_data_format string
************************************************************/

char *fit_data_format(const char *fmt)
{
    int len, n, i;
    static char fit_data_fmt[FIT_MAX_FMTSTRSIZE] = "";

    if (fmt != NULL) {
        len = (int)strlen(fmt);
        n = 0;
        i = 0;
        while (i < len - 1) {
            if (fmt[i] == '%') {
                if (fmt[i+1] == '*' || fmt[i+1] == '%') {
                    i += 2;
                    continue;
                }
                ++i;
                while ((i < len - 1) && isdigit((int)fmt[i]))
                    ++i;
                if (fmt[i] == 'l') {
                    if (strchr("efg", (int)fmt[i+1]) == NULL) {
                        errlog(91, 0);
                    } else {
                        ++n;
                        i += 2;
                    }
                } else {
                    errlog(91, 0);
                    return NULL;
                }
            } else {
                ++i;
            }
        }
        if (n < 2 || n > FIT_MAX_INPUT_FIELDS) {
            errlog(90, n);
            return NULL;
        }
        strncpy(fit_data_fmt, fmt, sizeof(fit_data_fmt) - 1);
        fit_data_fmt[sizeof(fit_data_fmt) - 1] = '\0';
    }
    if (strlen(fit_data_fmt) == 0) { /* not set */
        return NULL;
    }
    return (char *)fit_data_fmt;
}   /* fit_data_format() */


/***********************************************************
* static char *fit_indir(void)
*
* Purpose:
*   If a fit input data directory is named in the environment,
*   fit_indir() returns a pointer to a static string copied
*   from the environment.  If no directory separator character
*   is found at the end of the environment string, a default
*   one is appended to the copy.
*
*   If no string is named in the environment, a pointer to
*   an empty string is returned.
************************************************************/


static char *fit_indir(void)
{
    static char indir[AM_MAX_DIRPATH] = "";
    char *envstr;

    if ((envstr = getenv(FIT_INPUT_PATH_ENV_STRING)) != NULL) {
        if (strlen(envstr) < AM_MAX_DIRPATH) {
            strncpy(indir, envstr, sizeof(indir) - 1);
            indir[sizeof(indir) - 1] = '\0';
            check_for_dir_separator(indir);
        } else {
            errlog(31, 0);
        }
    }
    return (char *)indir;
}   /* fit_indir() */


/***********************************************************
* static char *fit_outdir(void)
*
* Purpose:
*   If a fit output directory is named in the environment,
*   fit_outdir() returns a pointer to a static string copied
*   from the environment.  If no directory separator character
*   is found at the end of the environment string, a default
*   one is appended to the copy.
*
*   If no string is named in the environment, a pointer to
*   an empty string is returned.
************************************************************/


static char *fit_outdir(void)
{
    static char outdir[AM_MAX_DIRPATH] = "";
    char *envstr;

    if ((envstr = getenv(FIT_OUTPUT_PATH_ENV_STRING)) != NULL) {
        if (strlen(envstr) < AM_MAX_DIRPATH) {
            strncpy(outdir, envstr, sizeof(outdir) - 1);
            outdir[sizeof(outdir) - 1] = '\0';
            check_for_dir_separator(outdir);
        } else {
            errlog(38, 0);
        }
    }
    return (char *)outdir;
}   /* fit_outdir() */


/***********************************************************
* static int init_fit_kcache(model_t *model, simplex_t *simplex)
*
* Purpose:
*   Initialize the index tables for the kcache, an in-memory
*   spectral absorption coefficient cache.  The kcache
*   functions on layers for for which temperature is a fit
*   variable, and is the only variable quantity affecting
*   column absorption coefficients on that layer.  In such
*   cases, a significant performance benefit is realized by
*   computing spectral absorption coefficients on a uniform
*   temperature grid, and interpolating between spectral
*   arrays to intermediate temperatures.
*
*   On layers which do not satisfy this condition, the cache
*   index tables are left uninitialized (NULL).  This is not
*   strictly needed to enforce correctness-- layer pressures
*   and column mixing ratios are checked for changes in
*   get_absorption_coefficient(), and kcache entries are
*   flushed there if necessary.  (Examples are when the
*   "set" facility is used to change a layer pressure or a
*   mixing ratio, or when an interpolated level moves to a
*   new layer.) However, anticipating these cases and not
*   initializing a kcache in the first place avoids kcache
*   thrashing that that hurts rather than helps performance.
*
* Arguments:
*   model_t *model     - pointer to model structure
*   simplex_t *simplex - pointer to simplex data structure
*
* Return:
*   1 on error, 0 otherwise
************************************************************/

static int init_fit_kcache(model_t *model, simplex_t *simplex)
{
    int lnum;
    
    /*
     * If any interpolated levels are fit variables, there
     * is a possibility that the variable level may move
     * from one parent layer to another.  A variable
     * interpolated level moving into a kcache-enabled level
     * would lead to thrashing of the kcache and dcache.
     * Because of this, if any variable interpolated levels
     * exist, the kcache is disabled on all levels.
     */
    if (model->geometry & GEOMETRY_POBS_USER_DEFINED &&
            isvar(simplex, &(model->Pobs)))
        return 0;
    if (model->geometry & GEOMETRY_ZOBS_USER_DEFINED &&
            isvar(simplex, &(model->zobs)))
        return 0;
    if (model->geometry & GEOMETRY_PSOURCE_USER_DEFINED &&
            isvar(simplex, &(model->Psource)))
        return 0;
    if (model->geometry & GEOMETRY_ZSOURCE_USER_DEFINED &&
            isvar(simplex, &(model->zsource)))
        return 0;
    if (model->geometry & GEOMETRY_PTAN_USER_DEFINED &&
            isvar(simplex, &(model->Ptan)))
        return 0;
    if (model->geometry & GEOMETRY_ZTAN_USER_DEFINED &&
            isvar(simplex, &(model->ztan)))
        return 0;
    for (lnum = 0; lnum < model->nlayers; ++lnum) {
        int cnum;
        layer_t *layer = model->layer[lnum];
        /*
         * If layer pressures are specified by dP, and dP_def is
         * variable on this layer, then k cannot be cached from
         * this layer onwards.  Break out of the loop.
         */
        if ((model->PTmode & PTMODE_DP) && isvar(simplex, &(layer->dP_def)))
            break;
        /*
         * If layer pressures are specified by P, and P is a fit
         * variable, then k cannot be cached on this layer.
         * (Note that Pbase is never a fit variable, and
         * does not need to be checked.)
         */
        if ((model->PTmode & PTMODE_P) && isvar(simplex, &(layer->P)))
            continue;
        /*
         * If layer temperatures are specified by T, and T is not a fit
         * variable, k does not need to be cached on this layer.
         */
        if ((model->PTmode & PTMODE_T) && !isvar(simplex, &(layer->T)))
            continue;
        /*
         * If layer temperatures are specified by Tbase, and Tbase does
         * not vary on this layer, or on the layer above, then k does not
         * need to be cached on this layer.
         */
        if ((model->PTmode & PTMODE_TBASE) &&
            !isvar(simplex, &(layer->Tbase))) {
            if (lnum == 0)
                continue;
            else if (!isvar(simplex, &(model->layer[lnum - 1]->Tbase)))
                continue;
        }
        /*
         * Skip this layer if it is empty.
         */
        if (layer->tau == NULL)
            continue;
        /*
         * k is cacheable in this layer.  Initialize the kcache table
         * for each column absorption coefficient in this layer for
         * which caching is possible.
         */
        for (cnum = 0; cnum < layer->ncols; ++cnum) {
            int knum;
            column_t *column = layer->column[cnum];
            /*
             * Skip this column if it is empty.
             */
            if (column->ztau == NULL)
                continue;
            for (knum = 0; knum < column->n_abscoeffs; ++knum) {
                abscoeff_t *abscoeff;
                Nscale_list_t *Nscale_list_entry_tagged;
                Nscale_list_t *Nscale_list_entry_global;
                Nscale_list_entry_tagged =
                    find_Nscale_list_entry(column->col_typenum, layer->tagnum);
                Nscale_list_entry_global =
                    find_Nscale_list_entry(column->col_typenum, 0);
                abscoeff = column->abscoeff[knum];
                /*
                 * If strict self broadening is in effect, then k cannot
                 * be cached unless the mixing ratio is user-defined and
                 * non-variable.  (Default mixing ratios are excluded
                 * because they are subject to automatic rescaling to
                 * accommodate changes in user-defined mixing ratios.)
                 */
                if (!layer->strict_selfbroad[abscoeff->k_typenum] ||
                    (column->vmr_stat & VMR_USER_DEFINED &&
                    !isvar(simplex, &(column->xvmr)) &&
                    !isvar(simplex, &(column->RH)) &&
                    !isvar(simplex, &(Nscale_list_entry_tagged->Nscale)) &&
                    !isvar(simplex, &(Nscale_list_entry_global->Nscale)))
                    ) {
                    if (init_kcache(abscoeff, model->nkcache)) {
                        errlog(52, lnum);
                        return 1;
                    }
                }
            }
        }
    }
    return 0;
}   /* init_fit_kcache() */


/***********************************************************
* static int open_fit_data_file(fit_data_t *fit_data)
*
* Purpose:
*   Opens the current fit data file and resets the data
*   block counter.  If the current fit data file name
*   is "stdin" or "-", the current fit data file pointer
*   is set to stdin.
*
* Arguments:
*   fit_data_t *fit_data - fit_data structure
*
* Return:
*   0 if OK, 1 otherwise
************************************************************/

static int open_fit_data_file(fit_data_t *fit_data)
{
    char path[AM_MAX_PATH];
    int status;

    if ((strcmp(fit_data->filename[fit_data->open_filenum], "stdin") == 0) ||
        (strcmp(fit_data->filename[fit_data->open_filenum], "-") == 0)) {
        fit_data->fp = stdin;
        status = 0;
    } else {
        snprintf(path,
                sizeof(path),
                "%s%s",
                fit_indir(),
                fit_data->filename[fit_data->open_filenum]);
        if ((fit_data->fp = fopen(path, "r")) == NULL) {
            fprintf(stderr, "Could not open file \"%s\"\n", path);
            status = 1;
        } else {
            status = 0;
        }
    }
    fit_data->blocknum = 0;
    return status;
} /* open_fit_data_file() */


/***********************************************************
* static int read_fit_data(
*   model_t *model,
*   fit_data_t *fit_data,
*   simplex_t *simplex)
*
* Purpose:
*   Reads a block of data from the current fit data file,
*   discarding data which fall outside the model frequency
*   grid range.  Also processes configuration statements
*   embedded in the data stream.
*
* Arguments:
*   model_t    *model    - pointer to model structure
*   fit_data_t *fit_data - pointer to fit_data structure
*   simplex_t  *simplex  - pointer to simplex data structure
*
* Return:
*   READ_FIT_DATA_EOF   - read to the end of the data file
*   READ_FIT_DATA_EOB   - read to the end of a data block
*                         within the data file
*   READ_FIT_DATA_ERROR - serious trouble, e.g. memory
*                         allocation failure or NaN in data
*   READ_FIT_DATA_RESTART_FILE_LIST
*                       - encountered the special file name
*                         "restart" in the file list.
************************************************************/

static int read_fit_data(
    model_t *model,
    fit_data_t *fit_data,
    simplex_t *simplex)
{
    double f_min = 0.0, f_max = 0.0;
    int status;
    int block_linenum;
    int in_data_block;
    static int file_linenum = 0;
    gridsize_t n_mod = 0;

    /*
     * Open the current fit data file, if it isn't open already.
     */
    if (fit_data->fp == NULL) {
        if (!strcmp(fit_data->filename[fit_data->open_filenum], RESTART_FNAME))
            return READ_FIT_DATA_RESTART_FILE_LIST;
        if (open_fit_data_file(fit_data)) {
            errlog(36,0);
            fit_data->npts = 0;
            return READ_FIT_DATA_EOF;
        }
        file_linenum = 0;
    }
    /*
     * Read through the data file a line at a time, until the end
     * of the file, or the end of a data block within the file,
     * is encountered.
     */
    block_linenum = 0;
    in_data_block = 0;
    for (;;) {
        char lbuf[FIT_DATA_LBUFSIZE];
        double f = DBL_MAX, s = DBL_MAX, b = 0.0 , w = 1.0;
        if (fgets(lbuf, (int)sizeof(lbuf), fit_data->fp) == NULL) {
            /*
             * If the end of a file (EOF) is encountered while in
             * a block of data, then this is the end of the data
             * block, and a fit will be run with the new data.
             * If EOF is encountered outside a data block, the
             * data block must already have been terminated with
             * an end of block string.  In this case, set the
             * number of points to zero to signal that there is
             * no new data.
             */
            if (!in_data_block) { 
                fit_data->npts = 0;
            } 
            status = READ_FIT_DATA_EOF; /* EOF, new data not yet fitted */
            break;
        }
        ++file_linenum;
        if (strncmp(lbuf, EOB_STRING, sizeof(EOB_STRING) - 1) == 0) {
            status = READ_FIT_DATA_EOB; /* end of data block in file */
            break;
        }
        ++block_linenum;
        if (strncmp(lbuf, PARAM_CTRL_STRING, sizeof(PARAM_CTRL_STRING) - 1)
            == 0) {
            /*
             * Embedded configuration statements are only allowed
             * before the first data point in a block, otherwise
             * they are ignored, and a warning is logged.
             */
            if (!in_data_block) {
                if (set_config_parameter(
                            lbuf + strlen(PARAM_CTRL_STRING),
                            fit_data->filename[fit_data->open_filenum],
                            file_linenum,
                            model,
                            fit_data,
                            simplex)) {
                    errlog(202, 0);
                    return READ_FIT_DATA_ERROR;
                }
                continue;
            } else {
                errlog(111, 0);
            }
        }
        if (strncmp(lbuf, RESTART_STRING, sizeof(RESTART_STRING) - 1) == 0) {
            /*
             * A restart is treated as a duplicate of the
             * previous data block.  If a restart string is
             * encountered within a data block, an error is
             * logged, and the restart is ignored.
             */
            if (!in_data_block) {
                status = READ_FIT_DATA_EOB;
                break;
            } else {
                errlog(115, 0);
            }
        }
        /*
         * Before processing the first data point, reset the
         * point count, and determine the frequency range covered
         * by the model.
         */
        if (!in_data_block) {
            fit_data->npts = 0;
            if (model->ifmode) {
                set_IF_spectrum_subgrid_ranges(model);
                f_min = model->fif_0;
                f_max = model->fif_0 + model->df * (model->nif - 1);
                n_mod = model->nif;
            } else {
                f_min = model->f[0];
                f_max = model->f[model->ngrid - 1];
                n_mod = model->ngrid;
            }
        }
        /*
         * Data lines which look like comments are ignored.  A
         * comment line starts with a '#', unless a
         * fit_data_format string explicitly starting with a '#'
         * has been defined.
         */
        if (strncmp(lbuf, "#", (size_t)1) == 0) {
            char *str = fit_data_format(NULL);
            if (str == NULL)
                continue;
            else if (strncmp(str, "#", (size_t)1) != 0)
                continue;
        }
        in_data_block = 1;
        if (fit_data_format(NULL) != NULL) { /* read data with sscanf() */
            double x[FIT_MAX_INPUT_FIELDS];
            int i, n, col;
            /*
             * Note: the args list for sscanf() must be edited if
             * FIT_MAX_INPUT_FIELDS changes.
             */
            n = sscanf(lbuf, fit_data_format(NULL), x, x+1, x+2, x+3);
            for (i = 0; i < n; ++i) {
                col = i + 1;
                if (col == fit_data->f_col)
                    f = x[i];
                if (col == fit_data->s_col)
                    s = x[i];
                if (col == fit_data->b_col)
                    b = x[i];
                if (col == fit_data->w_col)
                    w = x[i];
            }
        } else { /* read data directly by column number with strtok() */
            double x;
            char *tok, *endp;
            int col = 1;
            tok = strtok(lbuf, fit_data_delimiters(NULL));
            while (tok != NULL) {
                if ((col == fit_data->f_col) ||
                    (col == fit_data->s_col) ||
                    (col == fit_data->b_col) ||
                    (col == fit_data->w_col)) {
                    x = strtod(tok, &endp);
                    /* check for NaN or extra characters */
                    if (!(x == x) || (endp[0] != '\0')) {
                        errlog(143, file_linenum); 
                        return READ_FIT_DATA_ERROR;
                    }
                    if (col == fit_data->f_col)
                        f = x;
                    if (col == fit_data->s_col)
                        s = x;
                    if (col == fit_data->b_col)
                        b = x;
                    if (col == fit_data->w_col)
                        w = x;
                }
                ++col;
                tok = strtok(NULL, fit_data_delimiters(NULL));
            }
        }
        if (fabs(f) >= DBL_MAX || fabs(s) >= DBL_MAX)
            continue;
        if (b < 0.0) {
                errlog(42, 0);
                b = -b;
        }
        if (w < 0.0) {
                errlog(89, 0);
                w = -w;
        }
        /*
         * Convert the fit data to am native units.  The
         * weighting factor w = 1/sigma has units which are the
         * reciprocal of the spectral units.  No offset is
         * applied, because this is a differential quantity.
         */
        f += unit_tab[fit_data->f_unitnum].offset;
        f *= unit_tab[fit_data->f_unitnum].factor;
        s += unit_tab[fit_data->s_unitnum].offset;
        s *= unit_tab[fit_data->s_unitnum].factor;
        b += unit_tab[fit_data->f_unitnum].offset;
        b *= unit_tab[fit_data->f_unitnum].factor;
        w /= unit_tab[fit_data->s_unitnum].factor;
        /*
         * Skip to the next line if f is not within the model
         * frequency grid range.  For the special case of a
         * single-point model frequency grid, any data point
         * within 0.5 * model->df is accepted.
         */
        if (n_mod == 1) {
            if (f < (f_min - 0.5 * model->df) || f > (f_max + 0.5 * model->df))
                continue;
        } else {
            if (f < (f_min + 0.5 * b) || f > (f_max - 0.5 * b))
                continue;
        }
        /*
         * Increment the data point count, and enlarge the fit
         * data arrays if needed.
         */
        ++fit_data->npts;
        if (fit_data->nalloc < fit_data->npts) {
            if (grow_fit_data_arrays(fit_data)) {
                errlog(37, 0);
                close_fit_data_file(fit_data);
                return READ_FIT_DATA_ERROR;
            }
        }
        fit_data->f[fit_data->npts - 1] = f;
        fit_data->s[fit_data->npts - 1] = s;
        fit_data->b[fit_data->npts - 1] = b;
        fit_data->w[fit_data->npts - 1] = w;
    }
    if (status == READ_FIT_DATA_EOF) {
        close_fit_data_file(fit_data);
    }
    /*
     * If there are too few data points (and this wasn't just an
     * EOF on the line immediately after an EOB), signal this by
     * setting fit_data.npts to 0 and report an error.
     */
    if (fit_data->npts < (gridsize_t)simplex->n &&
        !(status == READ_FIT_DATA_EOF && !block_linenum && fit_data->blocknum)
        ) {
        fit_data->npts = 0;
        fprintf(stderr,
            "file \"%s\" data block %d: %ld data points < %d fit variables\n",
            fit_data->filename[fit_data->open_filenum],
            fit_data->blocknum,
            (long int)fit_data->npts,
            simplex->n);
        errlog(41, 0);
    }
    return status;
}   /* read_fit_data() */


/***********************************************************
* void report_fit_env_info(FILE *stream)
*
* Purpose:
*   Reports environment fit-related environment information
*   to a stream.
*
* Arguments:
*   FILE *stream - destination for output
************************************************************/

void report_fit_env_info(FILE *stream)
{
    fprintf(stream, "Fit I/O\n");
    fprintf(stream, "  %s = ", FIT_INPUT_PATH_ENV_STRING);
    if (!strcmp(fit_indir(), ""))
        fprintf(stream, "current directory (default setting)\n");
    else
        fprintf(stream, "%s\n", fit_indir());
    fprintf(stream, "  %s = ", FIT_OUTPUT_PATH_ENV_STRING);
    if (!strcmp(fit_outdir(), ""))
        fprintf(stream, "current directory (default setting)\n");
    else
        fprintf(stream, "%s\n", fit_outdir());
    return;
}   /* report_fit_env_info() */


/***********************************************************
* static void report_fit_status(
*   FILE *stream,
*   fit_data_t *fit_data,
*   simplex_t *simplex)
*
* Purpose:
*   Writes status information on the current fit iteration
*   to the specified stream.
*
* Arguments:
*   FILE *stream - stream to write status information to.
*   fit_data_t *fit_data - pointer to fit_data structure
*   simplex_t *simplex - pointer to simplex data structure
*
* Return:
*   none
************************************************************/

static void report_fit_status(
    FILE *stream,
    fit_data_t *fit_data,
    simplex_t *simplex)
{
    unsigned int i, j;

    fprintf(stream, "\n");
    fprintf(stream, "file %s : block %d\n",
        fit_data->filename[fit_data->open_filenum],
        fit_data->blocknum);
    fprintf(stream, "iter %d : ds %g (%g) : time %.3f s\n",
        simplex->iter,
        simplex_scaled_diameter(simplex),
        simplex->tol,
        fit_data->runtime);
    fprintf(stream, "vtx       E      ");
    for (j = 0; j < simplex->n; ++j)
        fprintf(stream, "     x%-3d    ", j);
    fprintf(stream, "\n");
    fprintf(stream, "--- ------------");
    for (j = 0; j < simplex->n; ++j)
        fprintf(stream, " ------------");
    fprintf(stream, "\n");
    for (i = 0; i <= simplex->n; ++i) {
        if (i == simplex->ihigh)
            fprintf(stream, "H");
        else if (i == simplex->ilow)
            fprintf(stream, "L");
        else
            fprintf(stream, " ");
        fprintf(stream, "%2d %12.4e", i, simplex->E[i]);
        for (j = 0; j < simplex->n; ++j) {
            print_with_unit(
                stream,
                " %12.4e",
                simplex->logarithmic ?
                    exp(simplex->vertex[i][j]) : simplex->vertex[i][j],
                simplex->unitnum[j]);
        }
        fprintf(stream, "\n");
    }
    fprintf(stream, "                ");
    for (j = 0; j < simplex->n; ++j)
        fprintf(stream, " ------------");
    fprintf(stream, "\n");
    fprintf(stream, "       max - min");
    for (j = 0; j < simplex->n; ++j) {
        print_differential_with_unit(
            stream,
            " %12.4e",
            simplex_variable_range(simplex, j),
            simplex->unitnum[j]);
    }
    fprintf(stream, "\n");
    fflush(stream);
    return;
}   /* report_fit_status() */


/***********************************************************
* static void compute_variance(
*   fit_data_t *fit_data,
*   simplex_t *simplex)
*
* Purpose:
*   If weights w=1/sigma have been specified for the fit
*   data points, this function computes the reduced
*   chi-squared statistic for the residuals.  For
*   unspecified weights, which default to unity, this will
*   instead be the mean variance of the residuals.
*
* Arguments:
*   fit_data_t *fit_data - pointer to fit_data structure
*   simplex_t *simplex - pointer to simplex data structure
*
* Return:
*   none
************************************************************/

static void compute_variance(
    fit_data_t *fit_data,
    simplex_t *simplex)
{
    gridsize_t i;

    fit_data->mean_var = 0.0;
    fit_data->mean_var_tracked = 0.0;
    for (i = 0; i < fit_data->npts; ++i) {
        double res, res_tracked;
        res = fit_data->s[i] - fit_data->s_mod[i];
        res_tracked = res - fit_data->res_est[i];
        fit_data->res[i] = res;
        res *= fit_data->w[i];
        res_tracked *= fit_data->w[i];
        fit_data->mean_var += res * res;
        fit_data->mean_var_tracked += res_tracked * res_tracked;
    }
    if (fit_data->npts <= (int)simplex->n) {
        fit_data->mean_var = 0.0;
        fit_data->mean_var_tracked = 0.0;
    } else {
        fit_data->mean_var /= fit_data->npts - simplex->n;
        fit_data->mean_var_tracked /= fit_data->npts - simplex->n;
    }
    return;
}   /* compute_variance() */


/***********************************************************
* static int update_simplex_vertex(
*   model_t *model,
*   model_t *lmodel,
*   fit_data_t *fit_data,
*   simplex_t *simplex,
*   double *vertex,
*   double *E)
*
* Purpose:
*
*   This function updates the model computation at the
*   specified simplex vertex, and computes a corresponding
*   goodness of fit estimator.
*
* Arguments:
*   model_t *model - pointer to model structure
*   fit_data_t *fit_data - pointer to fit_data structure
*   simplex_t *simplex - pointer to simplex data structure
*   double *vertex - pointer to vertex to be computed
*   double *E - pointer to corresponding goodness of fit estimator
*
* Return:
*   0 if OK, 1 otherwise
************************************************************/

static int update_simplex_vertex(
    model_t *model,
    model_t *lmodel,
    fit_data_t *fit_data,
    simplex_t *simplex,
    double *vertex,
    double *E)
{
    double *f_mod, *s_mod;
    gridsize_t i, n_mod;
    unsigned int iv;

    /*
     * Begin by substituting the vertex coordinates for the
     * corresponding model parameters.
     */
    for (iv = 0; iv < simplex->n; ++iv) {
        *simplex->varptr[iv] =
            simplex->logarithmic ? exp(vertex[iv]) : vertex[iv];
    }
    /*
     * Do the model computation.
     */
    if (compute_model(model, lmodel))
        return 1;
    /*
     * Interpolate the model output onto the fit data frequency
     * grid.  For the special case of a single-point grid, just
     * use the single grid point without interpolating.
     */
    f_mod = model->ifmode ? model->fif : model->f;
    s_mod = output[fit_data->data_type].spectrum;
    n_mod = model->ifmode ? model->nif : model->ngrid;
    for (i = 0; i < fit_data->npts; ++i) {
        if (n_mod == 1) {
            fit_data->s_mod[i] = s_mod[0];
        } else if (fit_data->b[i] < model->df * DBL_EPSILON) {
            /* just interpolate */
            double frac, fint;
            gridsize_t im;
            frac = modf((fit_data->f[i] - f_mod[0]) / model->df, &fint);
            im   = (gridsize_t)fint;
            /*
             * If the fit data point is right at the top of the
             * model frequency range, or slightly over due to
             * rounding, then back off by one model grid point
             * and do a slight extrapolation.  Note that the same
             * situation at the low end of the frequency range is
             * handled automatically, since modf() will just
             * return a small negative fractional part, and zero
             * integer part.
             */
            if (im == n_mod - 1) {
                --im;
                frac += 1.0;
            }
            fit_data->s_mod[i] = s_mod[im] + frac * (s_mod[im+1] - s_mod[im]);
        } else {
            /* do a trapezoidal integration over the channel bandwidth */
            double f1, f2;
            double frac1, frac2, fint1, fint2;
            double s1, s2;
            gridsize_t im1, im2;
            f1    = fit_data->f[i] - 0.5 * fit_data->b[i];
            f2    = f1 + fit_data->b[i];
            frac1 = modf((f1 - f_mod[0]) / model->df, &fint1);
            frac2 = modf((f2 - f_mod[0]) / model->df, &fint2);
            im1   = (gridsize_t)fint1;
            im2   = (gridsize_t)fint2;
            if (im2 == n_mod - 1) {
                --im2;
                frac2 += 1.0;
            }
            s1 = s_mod[im1] + frac1 * (s_mod[im1+1] - s_mod[im1]);
            s2 = s_mod[im2] + frac2 * (s_mod[im2+1] - s_mod[im2]);
            if (im1 == im2) {
                fit_data->s_mod[i] = 0.5 * (s1 + s2);
            } else {
                double norm;
                gridsize_t j;
                fit_data->s_mod[i]  = (1.0 - frac1) * (s1 + s_mod[im1+1]);
                fit_data->s_mod[i] += frac2 * (s_mod[im2] + s2);
                norm = 2.0 * (1.0 - frac1 + frac2);
                for (j = im1 + 1; j < im2; ++j) {
                    fit_data->s_mod[i] += s_mod[j] + s_mod[j+1];
                    norm += 2.0;
                }
                fit_data->s_mod[i] /= norm;
            }
        }
    }
    /*
     * Compute raw residuals.
     */
    for (i = 0; i < fit_data->npts; ++i)
        fit_data->res[i] = fit_data->s[i] - fit_data->s_mod[i];
    /*
     * Compute the fit estimator.  res_est[] are the estimated
     * residuals if residual tracking is active.  Otherwise, all
     * the res_est[] values are zero.
     */
    *E = 0.;
    switch (fit_data->estimator_type) {
    case FIT_ESTIMATOR_ABSOLUTE_RESIDUALS:
        for (i = 0; i < fit_data->npts; ++i) {
            double nres;
            nres  = fit_data->res[i] - fit_data->res_est[i];
            nres *= fit_data->w[i];
            *E += fabs(nres);
        }
        break;
    case FIT_ESTIMATOR_SQUARED_RESIDUALS:
        for (i = 0; i < fit_data->npts; ++i) {
            double nres;
            nres  = fit_data->res[i] - fit_data->res_est[i];
            nres *= fit_data->w[i];
            *E += nres * nres;
        }
        break;
    case FIT_ESTIMATOR_LOGARITHMIC:
        for (i = 0; i < fit_data->npts; ++i) {
            double nres;
            nres  = fit_data->res[i] - fit_data->res_est[i];
            nres *= fit_data->w[i];
            *E += log(1.0 + 0.5 * nres * nres);
        }
        break;
    default:
        errlog(48, 0);
        return 1;
    }
    return 0;
}   /* update_simplex_vertex() */


/***********************************************************
* static int write_fit_output_files(
*   model_t *model,
*   fit_data_t *fit_data,
*   simplex_t *simplex)
*
* Purpose:
*   Writes output files containing model configuration,
*   model spectra, and residual spectrum data.  The files to
*   be written are governed by bits of fit_data.output_mode.
*   The bit masks are
*
*     FIT_OUTPUT_CONFIG
*     FIT_OUTPUT_SPECTRUM
*     FIT_OUTPUT_RESIDUAL
*
*   defined in fit.h.
*   
*   The first time any of these files is written, a new file
*   is created, discarding any existing file of the same
*   name.  Subsequent writes, as may occur in the case of
*   fits to multiple blocks of data, are appended to the
*   file.  A blank line marks the boundary between blocks.
*
* Arguments:
*   model_t    *model    - pointer to model structure
*   fit_data_t *fit_data - pointer to fit_data structure
*   simplex_t  *simplex  - pointer to simplex data structure
*
* Return:
*   0 if OK, 1 otherwise
************************************************************/

static int write_fit_output_files(
    model_t *model,
    fit_data_t *fit_data,
    simplex_t *simplex)
{
    FILE *stream;
    char path[AM_MAX_PATH];
    static int config_first_write   = 1;
    static int spectrum_first_write = 1;
    static int residual_first_write = 1;

    /*
     * Since the output mode can change on-the-fly between data
     * blocks or data files, some extra logic is needed to ensure
     * that files are written, and old files overwritten, in a
     * consistent way.
     *
     * If this is data block 0 for this data file, the
     * first_write flags are set, indicating that the next write,
     * if any, to the corresponding output file will be the first
     * write for that file name.  On first write, any existing
     * file of the same name is overwritten, and the output data
     * block is written without a preceeding blank line.
     * Subsequent output blocks, if any, are appended, after a
     * blank line to separate blocks.
     */
    if (fit_data->blocknum == 0) {
        config_first_write   = 1;
        spectrum_first_write = 1;
        residual_first_write = 1;
    }
    /*
     * Write the model configuration file...
     */
    if (fit_data->output_mode & FIT_OUTPUT_CONFIG) {
        snprintf(path,
                sizeof(path),
                "%s%s%s",
                fit_outdir(),
                fit_data->filename[fit_data->open_filenum],
                MODEL_CONFIG_EXT);
        if ((stream = fopen(path, config_first_write ? "w" : "a+")) == NULL) {
            fprintf(stderr, "Could not open file \"%s\"\n", path);
            errlog(43, 0);
            return 1;
        }
        if (!config_first_write)
            fprintf(stream, "\n");
        config_first_write = 0;
        write_model_config_data(stream, model, fit_data, simplex);
        fclose(stream);
    }
    /*
     * the model spectrum file...
     */
    if (fit_data->output_mode & FIT_OUTPUT_SPECTRUM) {
        snprintf(path,
                sizeof(path),
                "%s%s%s",
                fit_outdir(),
                fit_data->filename[fit_data->open_filenum],
                MODEL_SPECTRUM_EXT);
        if ((stream = fopen(path, spectrum_first_write ? "w" : "a+")) == NULL) {
            fprintf(stderr, "Could not open file \"%s\"\n", path);
            errlog(43, 0);
            return 1;
        }
        if (!spectrum_first_write)
            fprintf(stream, "\n");
        spectrum_first_write = 0;
        write_model_spectra(stream, model, simplex);
        fclose(stream);
    }
    /*
     * and the residual spectrum file...
     */
    if (fit_data->output_mode & FIT_OUTPUT_RESIDUAL) {
        snprintf(path,
                sizeof(path),
                "%s%s%s",
                fit_outdir(),
                fit_data->filename[fit_data->open_filenum],
                RESIDUAL_SPECTRUM_EXT);
        if ((stream = fopen(path, residual_first_write ? "w" : "a+")) == NULL) {
            fprintf(stderr, "Could not open file \"%s\"\n", path);
            errlog(43, 0);
            return 1;
        }
        if (!residual_first_write)
            fprintf(stream, "\n");
        residual_first_write = 0;
        write_fit_residuals(stream, fit_data);
        fclose(stream);
    }
    return 0;
}   /* write_fit_output_files() */


/***********************************************************
* static void write_fit_summary_line(
*   fit_data_t *fit_data,
*   simplex_t *simplex)
*
* Purpose:
*   Writes to stdout a single line of the form:
*
*   filename block# #iterations stat par[1] range[1] ... par[n] range[n]
*
*   where stat is the reduced chi squared statistic, if
*   weights were given for the data points, or the standard
*   deviation of the residuals if not.  If residual tracking
*   is on, then stat is the reduced chi squared or standard
*   deviation for the corrected residuals.  Each par[i] is
*   the best fit value of the ith fit variable, and range[i]
*   is its range over the converged simplex.  If the fit
*   data source was the standard input, the filename is
*   omitted. 
*
* Arguments:
*   fit_data_t *fit_data - pointer to fit_data structure
*   simplex_t *simplex - pointer to simplex data structure
************************************************************/

static void write_fit_summary_line(
    fit_data_t *fit_data,
    simplex_t *simplex)
{
    unsigned int j;
    double stat;

    if ((strcmp(fit_data->filename[fit_data->open_filenum], "stdin") != 0) &&
        (strcmp(fit_data->filename[fit_data->open_filenum], "-") != 0)
        ) {
        printf("%s ", fit_data->filename[fit_data->open_filenum]);
    }
    if (fit_data->res_track_gain < 0.0)
        stat = fit_data->mean_var;
    else
        stat = fit_data->mean_var_tracked;
    if (fit_data->w_col == 0) /* no weights */
        stat = sqrt(stat);
    printf("%4d %4d %14.6e", fit_data->blocknum, simplex->iter, stat);
    for (j = 0; j < simplex->n; ++j) {
        print_with_unit(
            stdout,
            "  %#g",
            (simplex->logarithmic ?
                exp(simplex->vertex[0][j]) : simplex->vertex[0][j]),
            simplex->unitnum[j]);
        print_with_unit(
            stdout,
            " %#.3g",
            simplex_variable_range(simplex, j),
            simplex->unitnum[j]);
    }
    printf("\n");
    fflush(stdout);
    return;
}   /* write_fit_summary_line() */


/***********************************************************
* void reset_estimated_residuals(fit_data_t* fit_data)
*
* Purpose:
*   Resets the estimated residuals array to zero.
*
* Arguments:
*   fit_data_t *fit_data - pointer to fit_data structure
************************************************************/

void reset_estimated_residuals(fit_data_t* fit_data)
{
    gridsize_t i;

    for (i = 0; i < fit_data->nalloc; ++i)
        fit_data->res_est[i] = 0.0;
    return;
}   /* reset_estimated_residuals() */


/***********************************************************
* static void update_estimated_residuals(fit_data_t* fit_data)
*
* Purpose:
*   Compute the estimated residuals for the upcoming fit,
*   using a recursive filter based on the prior fit
*   residuals.
*
* Arguments:
*   fit_data_t *fit_data - pointer to fit_data structure
************************************************************/

static void update_estimated_residuals(fit_data_t* fit_data)
{
    gridsize_t i;

    if (fit_data->res_track_gain < 0.0)  /* tracking off */
        return;
    if (fit_data->res_track_gain == 0.0) /* nothing to do */
        return;
    for (i = 0; i < fit_data->npts; ++i) {
        fit_data->res_est[i] *= 1.0 - fit_data->res_track_gain;
        fit_data->res_est[i] += fit_data->res_track_gain * fit_data->res[i];
    }
    return;
}   /* update_estimated_residuals() */
