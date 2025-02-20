/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* output.c                   S. Paine rev. 2024 September 26
*
* Functions for program output.
************************************************************/

#include <ctype.h>
#ifdef _WIN32
    #include <fcntl.h>
#endif
#include <float.h>
#ifdef _WIN32
    #include <io.h>
#endif
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "abscoeff.h"
#include "am_sysdep.h"
#include "am_types.h"
#include "column.h"
#include "dcache.h"
#include "errlog.h"
#include "fit.h"
#include "ils.h"
#include "kcache.h"
#include "layer.h"
#include "linesum.h"
#include "mapping.h"
#include "model.h"
#include "nscale.h"
#include "output.h"
#include "phys_const.h"
#include "simplex.h"
#include "tags.h"
#include "units.h"
#include "version.h"

/*
 * Verbosity levels for reporting layer and column data
 */
enum {
    VERBOSITY_FULL,
    VERBOSITY_PARENT_LAYER,
    VERBOSITY_LAYER_NOT_IN_PATH
};

struct output_tabentry output[] = {
    {
        "",        /* text name of output array                             */
        NULL,      /* pointer to spectrum array                             */
        NULL,      /* pointer to Jacobian matrix                            */
        NULL,      /* pointer to Jacobian estimated rounding error matrix   */
        NULL,      /* pointer to Jacobian estimated truncation error matrix */
        NULL,      /* pointer to Jacobian sorted total error matrix         */
        UNIT_NONE, /* default unit number                                   */
        UNIT_NONE, /* user unit number                                      */
        0,         /* k_typenum for k output header                         */
        0,         /* flags                                                 */
        0          /* format                                                */
    },
    {
        "f",
        NULL,
        NULL,
        NULL,
        NULL,
        NULL,
        AM_UNIT_FREQUENCY,
        AM_UNIT_FREQUENCY,
        0,
        OUTPUT_USER|OUTPUT_ACTIVE,
        0
    },
    {
        "tau",
        NULL,
        NULL,
        NULL,
        NULL,
        NULL,
        AM_UNIT_OPACITY,
        AM_UNIT_OPACITY,
        0,
        JACOBIAN_ALLOWED,
        0
    },
    {
        "tx",
        NULL,
        NULL,
        NULL,
        NULL,
        NULL,
        AM_UNIT_NONE,
        AM_UNIT_NONE,
        0,
        OUTPUT_USER|OUTPUT_ACTIVE|ILS_ALLOWED|ILS_APPLIED|JACOBIAN_ALLOWED,
        0
    },
    {
        "I",
        NULL,
        NULL,
        NULL,
        NULL,
        NULL,
        AM_UNIT_RADIANCE,
        AM_UNIT_RADIANCE,
        0,
        ILS_ALLOWED|ILS_APPLIED|JACOBIAN_ALLOWED,
        0
    },
    {
        "I_diff",
        NULL,
        NULL,
        NULL,
        NULL,
        NULL,
        AM_UNIT_RADIANCE,
        AM_UNIT_RADIANCE,
        0,
        ILS_ALLOWED|ILS_APPLIED|JACOBIAN_ALLOWED,
        0
    },
    {
        "Tb",
        NULL,
        NULL,
        NULL,
        NULL,
        NULL,
        AM_UNIT_TEMPERATURE,
        AM_UNIT_TEMPERATURE,
        0,
        OUTPUT_USER|OUTPUT_ACTIVE|ILS_ALLOWED|ILS_APPLIED|JACOBIAN_ALLOWED,
        0
    },
    {
        "Trj",
        NULL,
        NULL,
        NULL,
        NULL,
        NULL,
        AM_UNIT_TEMPERATURE,
        AM_UNIT_TEMPERATURE,
        0,
        ILS_ALLOWED|ILS_APPLIED|JACOBIAN_ALLOWED,
        0
    },
    {
        "Tsys",
        NULL,
        NULL,
        NULL,
        NULL,
        NULL,
        AM_UNIT_TEMPERATURE,
        AM_UNIT_TEMPERATURE,
        0,
        ILS_ALLOWED|ILS_APPLIED|JACOBIAN_ALLOWED,
        0
    },
    {
        "Y",
        NULL,
        NULL,
        NULL,
        NULL,
        NULL,
        AM_UNIT_NONE,
        AM_UNIT_NONE,
        0,
        ILS_ALLOWED|ILS_APPLIED|JACOBIAN_ALLOWED,
        0
    },
    {
        "L",
        NULL,
        NULL,
        NULL,
        NULL,
        NULL,
        AM_UNIT_DELAY,
        AM_UNIT_DELAY,
        0,
        JACOBIAN_ALLOWED,
        0
    },
    {
        "tau_fsl",
        NULL,
        NULL,
        NULL,
        NULL,
        NULL,
        AM_UNIT_OPACITY,
        AM_UNIT_OPACITY,
        0,
        JACOBIAN_ALLOWED,
        0
    },
    {
        "k",
        NULL,
        NULL,
        NULL,
        NULL,
        NULL,
        UNIT_AUTO, /* molecular or binary coeff, determined at runtime */
        UNIT_AUTO,
        0,
        JACOBIAN_ALLOWED,
        0
    },
    {
        "END",
        NULL,
        NULL,
        NULL,
        NULL,
        NULL,
        UNIT_NONE,
        UNIT_NONE,
        0,
        0,
        OUTPUT_FORMAT_TEXT
    }
};

/*
 * The following array lists the user output columns in the order they
 * are to be written.  The initialization here determines the default
 * output spectra.
 */
int outcol[OUTPUT_END_OF_TABLE] = {
    OUTPUT_FREQUENCY,
    OUTPUT_TRANSMITTANCE,
    OUTPUT_TB_PLANCK
};

const char* output_format_name[] = {
    "text",
    "csv",
    "npy",
    "end"
};

/*
 * Constants related to output formats.
 */

enum {
    NPY_HEADER_ALIGNMENT = 64,   /* 16 in NEP 1, but 64 in NumPy docs */
    NPY_VERSION_MAJOR    = 1,
    NPY_VERSION_MINOR    = 0
};

static const char TXT_FIELD_SEPARATOR[]  = " ";
static const char TXT_RECORD_SEPARATOR[] = "\n";
static const char CSV_FIELD_SEPARATOR[]  = ",";
/*
 * The RFC 4180 (2005) specification for the CSV format calls for
 * "\r\n" as a record separator.  On Windows, "\n" will be
 * translated to "\r\n" when writing to a text-mode stream.
 */
#if defined _WIN32
    static const char CSV_RECORD_SEPARATOR[] = "\n";
#else
    static const char CSV_RECORD_SEPARATOR[] = "\r\n";
#endif
static const char NPY_MAGIC_STRING[]     = "\x93NUMPY";

static void write_column_data(
        FILE*,
        model_t*,
        simplex_t*,
        const int,
        const int,
        const int);
static void write_column_density_summary(FILE*, model_t*);
static void write_fit_data_delimiters(FILE *);
static void write_fit_variable_summary(FILE*, simplex_t*);
static void write_jacobian_differentiation_variables_summary(FILE*, simplex_t*);
static void write_jacobian_error_estimate_summary(FILE*, model_t*, simplex_t*);
static void write_layer_data(
        FILE*,
        const char*,
        model_t*,
        simplex_t*,
        const int,
        const int);
static void write_layer_lineshape_config(FILE*, const char*, layer_t*);
static void write_model_spectra_as_npy(
        FILE*,
        simplex_t*,
        gridsize_t,
        gridsize_t);
static void write_model_spectra_as_text(
        FILE*,
        model_t*,
        simplex_t*,
        gridsize_t,
        gridsize_t,
        const char*,
        const char*);
static void write_Nscale_config(FILE*, simplex_t*);
static void write_output_headers_as_text(
        FILE*,
        simplex_t*,
        const char*,
        const char*);
static void write_runtime_details(FILE*, model_t*);
static void write_variable_info(
        FILE*, const char*, double, double*, int, simplex_t*);


/***********************************************************
* int set_active_outputs(const int mode)
*
* Purpose:
*   Modifies the OUTPUT_ACTIVE flags in the global output[]
*   table according to mode.  mode is either the identifier
*   naming a single output array, or ALL_OUTPUTS.  In the
*   former case, only the computations needed to produce the
*   named array will be performed when compute_model() is
*   called.  This mode is used during fits to avoid
*   redundant computations.  If mode is ALL_OUTPUTS, then
*   all spectra which are program outputs or are being
*   fitted will be computed.
*   
*   The OUTPUT_ACTIVE bit in output[ALL_OUTPUTS].flags is
*   used to indicate whether or not all user-requested
*   outputs are currently active, which in turn indicates
*   whether a fit is taking place.
*
* Arguments:
*   int mode - either ALL_OUTPUTS, or one of the individual
*              output identifiers enumerated in output.h.
*
* Return:
*   0 if successful, 1 otherwise
************************************************************/

int set_active_outputs(const int mode)
{
    int i;

    if (mode == ALL_OUTPUTS) {
        for (i = 1; i < OUTPUT_END_OF_TABLE; ++i) {
            if (output[i].flags &
                    (OUTPUT_USER | OUTPUT_FITTED | OUTPUT_JACOBIAN))
                output[i].flags |= OUTPUT_ACTIVE;
        }
    } else if (mode > OUTPUT_NONE && mode < OUTPUT_END_OF_TABLE) {
        for (i = 0; i < OUTPUT_END_OF_TABLE; ++i)
            output[i].flags &= ~(OUTPUT_ACTIVE);
    } else {
        errlog(39, 0);
        return 1;
    }
    output[mode].flags |= OUTPUT_ACTIVE;
    return 0;
}   /* set_active_outputs() */


/***********************************************************
* void write_model_config_data(
*       FILE       *stream,
*       model_t    *model,
*       fit_data_t *fit_data,
*       simplex_t  *simplex)
*
* Purpose:
*   Writes the model configuration data to a stream, in
*   canonical config file format.
*
*   In general, in this function and those called by it,
*   user-supplied numbers are written with %.13g format to
*   recover all the digits supplied by the user.  (%.13g is
*   used rather than %.15g, because the latter can have base
*   conversion artifacts in the last one or two decimal
*   places.)  Computed numbers are written with an
*   appropriately-chosen truncation.
************************************************************/

void write_model_config_data(
        FILE *stream,
        model_t *model,
        fit_data_t *fit_data,
        simplex_t *simplex)
{
    unsigned int i, j;
    int lnum;
    int verbosity;

    /*
     * Print version, run timing, cache statistics, and fit
     * convergence information.
     */
    write_version_line(stream);
    if (output[ALL_OUTPUTS].flags & OUTPUT_FITTED) {
        fprintf(stream,
                "# fit run time %.3f s (%d iterations%s)\n",
                fit_data->runtime,
                simplex->iter,
                (simplex->iter == fit_data->max_iter) ?
                    ", not converged" : "");
        if (model->log_runtimes) {
            fprintf(stream,
                    "# last model computation run time %.3f s\n",
                    model->runtime);
            write_runtime_details(stream, model);
        }
        if (model->nkcache != 0)
            report_kcache_log_data(stream);
        report_dcache_log_data(stream);
        fprintf(stream,
                "# %ld data point%s, %d fit variable%s\n",
                (long int)fit_data->npts,
                fit_data->npts == 1 ? "" : "s",
                simplex->n,
                simplex->n == 1 ? "" : "s");
        write_fit_variable_summary(stream, simplex);
        if (fit_data->w_col == 0) {
            /*
             * No weights were defined, so report standard
             * deviation
             */
            const char *unitname = NULL;
            if (output[fit_data->data_type].default_unitnum != AM_UNIT_NONE) {
                unitname =
                    unit_tab[output[fit_data->data_type].default_unitnum].name; 
            }
            fprintf(stream,
                    "# standard deviation of residuals %.4g %s",
                    sqrt(fit_data->mean_var),
                    unitname == NULL ? "" : unitname);
            if (fit_data->res_track_gain >= 0.0) {
                fprintf(stream,
                        " (raw)  %.4g %s (tracked)\n",
                        sqrt(fit_data->mean_var_tracked),
                        unitname == NULL ? "" : unitname);
            } else {
                fprintf(stream, "\n");
            }
        } else {
            /*
             * Weights were defined, so report reduced
             * chi-squared
             */
            fprintf(stream,
                    "# reduced chi-squared %.4g",
                    fit_data->mean_var);
            if (fit_data->res_track_gain >= 0.0) {
                fprintf(stream,
                        " (raw)  %.4g (tracked)\n",
                        fit_data->mean_var_tracked);
            } else {
                fprintf(stream, "\n");
            }
        }
    } else {
        fprintf(stream, "# run time %.3f s\n", model->runtime);
        if (model->log_runtimes)
            write_runtime_details(stream, model);
        report_dcache_log_data(stream);
    }
    fprintf(stream, "\n");
    /*
     * If this was a fit, print fit configuration and kcache
     * information
     */
    if (output[ALL_OUTPUTS].flags & OUTPUT_FITTED) {
        fprintf(stream,
                "# fit %s %s\n",
                output[fit_data->data_type].name,
                fit_data->filename[fit_data->open_filenum]);
        fprintf(stream,
                "fit_data_columns %d %d %d %d\n",
                fit_data->f_col,
                fit_data->s_col,
                fit_data->b_col,
                fit_data->w_col);
        write_fit_data_delimiters(stream);
        if (fit_data_format(NULL) != NULL)
            fprintf(stream, "fit_data_format \"%s\"\n", fit_data_format(NULL));
        fprintf(stream,
                "fit_data_units %s %s\n",
                unit_tab[fit_data->f_unitnum].name,
                unit_tab[fit_data->s_unitnum].name);
        fprintf(stream,
                "fit_estimator %s (%.4g)\n",
                fit_estimator_type[fit_data->estimator_type].name,
                simplex->E[simplex->ilow]);
        fprintf(stream, "fit_iter %d\n", fit_data->max_iter);
        fprintf(stream, "fit_reinit %d\n", simplex->reinit);
        fprintf(stream,
                "fit_output config %d\n",
                fit_data->output_mode & FIT_OUTPUT_CONFIG ? 1 : 0);
        fprintf(stream,
                "fit_output spectrum %d\n",
                fit_data->output_mode & FIT_OUTPUT_SPECTRUM ? 1 : 0);
        fprintf(stream,
                "fit_output residual %d\n",
                fit_data->output_mode & FIT_OUTPUT_RESIDUAL ? 1 : 0);
        fprintf(stream,
                "fit_tol %.13g %d (%.4g)\n",
                simplex->tol,
                fit_data->max_restarts,
                simplex_scaled_diameter(simplex));
        if (fit_data->res_track_gain < 0.0) {
            fprintf(stream, "fit_track_residuals off\n");
        } else {
            fprintf(stream,
                    "fit_track_residuals %.13g\n",
                    fit_data->res_track_gain);
        }
        fprintf(stream,
                "fit_verbose %d\n",
                fit_data->output_mode & FIT_OUTPUT_VERBOSE ? 1 : 0);
        fprintf(stream, "kcache ");
        if (model->nkcache != 0) {
            print_with_unit(
                    stream,
                    "%.13g %s  ",
                    model->kcache_Tmin,
                    AM_UNIT_TEMPERATURE);
            print_with_unit(
                    stream,
                    "%.13g %s  ",
                    model->kcache_Tmax,
                    AM_UNIT_TEMPERATURE);
            print_with_unit(
                    stream,
                    "%.13g %s\n",
                    model->kcache_dT,
                    AM_UNIT_TEMPERATURE);
        } else {
            fprintf(stream, "off\n");
        }
        fprintf(stream, "simplex_log %d\n", simplex->logarithmic);
        fprintf(stream, "\n");
    }

    /*
     * Below, model configuration variables are printed using the
     * units given by the user where applicable.
     *
     * For fit variables, print the initial variable scale, and
     * the range of values for the variable over the vertices of
     * the converged simplex.  However, don't print these out
     * unless this was actually a fit (i.e.  if
     * output[ALL_OUTPUTS].flags & OUTPUT_FITTED is true).
     */

    /*
     * Frequency grid definition
     */
    print_with_unit(
            stream,
            "f %.13g %s",
            model->fmin,
            model->fmin_unitnum);
    print_with_unit(
            stream,
            "  %.13g %s",
            model->fmax,
            model->fmax_unitnum);
    print_with_unit(
            stream,
            "  %.13g %s\n",
            model->df,
            model->df_unitnum);
    /*
     * IF spectrum
     */
    if (model->ifmode) {
        if (model->ifmode & IFMODE_DSB) {
            fprintf(stream, "ifspec dsb");
        } else if (model->ifmode & IFMODE_USB) {
            fprintf(stream, "ifspec usb");
        } else if (model->ifmode & IFMODE_LSB) {
            fprintf(stream, "ifspec lsb");
        } else {
            fprintf(stream, "ifspec ???");
        }
        write_variable_info(
                stream,
                "",
                model->flo,
                &model->flo,
                model->flo_unitnum,
                simplex);
        fprintf(stream, "\n");
    }
    /*
     * restricted IF frequency range
     */
    if (model->fif_max >= 0.0) {
        print_with_unit(
                stream,
                "fif %.13g %s",
                model->fif_min,
                model->fif_min_unitnum);
        print_with_unit(
                stream,
                "  %.13g %s\n",
                model->fif_max,
                model->fif_max_unitnum);
    }
    /*
     * restricted output frequency range
     */
    if (model->fout_max >= 0.0) {
        print_with_unit(
                stream,
                "fout %.13g %s",
                model->fout_min,
                model->fout_min_unitnum);
        print_with_unit(
                stream,
                "  %.13g %s\n",
                model->fout_max,
                model->fout_max_unitnum);
    }
    /*
     * output definition
     */
    fprintf(stream,
            "output %s",
            output_format_name[output[ALL_OUTPUTS].format]);
    for (i = 0; outcol[i] != 0; ++i)
        fprintf(stream,
                "  %s %s",
                output[outcol[i]].name,
                unit_tab[output[outcol[i]].unitnum].name);
    fprintf(stream, "\n\n");
    /*
     * Jacobians
     */
    if (output[ALL_OUTPUTS].flags & OUTPUT_JACOBIAN) {
        fprintf(stream, "jacobian");
        for (i = 0; outcol[i] != 0; ++i) {
            if (output[outcol[i]].flags & OUTPUT_JACOBIAN)
                fprintf(stream, " %s", output[outcol[i]].name);
        }
        if (output[ALL_OUTPUTS].flags & OUTPUT_JACOBIAN_ERRS)
            fprintf(stream, " include_estimated_errors");
        fprintf(stream, "\n");
        /*
         * Print a simplex_log statement only if it hasn't been
         * printed already as part of a fit configuration.
         */
        if (!(output[ALL_OUTPUTS].flags & OUTPUT_FITTED))
            fprintf(stream, "simplex_log %d\n", simplex->logarithmic);
        fprintf(stream, "\n");
        fprintf(stream, "#\n");
        fprintf(stream, "# Summary data for Jacobians\n");
        fprintf(stream, "#\n");
        write_jacobian_differentiation_variables_summary(stream, simplex);
        write_jacobian_error_estimate_summary(stream, model, simplex);
        fprintf(stream, "\n");
    }
    /*
     * Instrumental line shape
     */
    if (model->ils != NULL) {
        fprintf(stream,
                "ILS %s",
                ils_type[model->ils_typenum].name);
        if (isvar(simplex, &(model->ils_fwhm))) {
            j = get_simplex_variable_index(simplex, &(model->ils_fwhm));
            print_with_unit(
                    stream,
                    " %g %s",
                    model->ils_fwhm,
                    model->ils_fwhm_unitnum);
            print_differential_with_unit(
                    stream,
                    "  %.13g %s",
                    simplex->scale[j],
                    model->ils_fwhm_unitnum);
        } else {
            print_with_unit(
                    stream,
                    " %.13g %s ",
                    model->ils_fwhm,
                    model->ils_fwhm_unitnum);
        }
        for (i = 0; outcol[i] != 0; ++i) {
            if (output[outcol[i]].flags & ILS_APPLIED)
                fprintf(stream, " %s", output[outcol[i]].name);
        }
        /*
         * Printing of the simplex range is deferred to after
         * the list of array names, so the syntax is preserved.
         */
        if (isvar(simplex, &(model->ils_fwhm))
                && output[ALL_OUTPUTS].flags & OUTPUT_FITTED) {
            j = get_simplex_variable_index(simplex, &(model->ils_fwhm));
            print_differential_with_unit(
                    stream,
                    "  (%.3g %s) ",
                    simplex_variable_range(simplex, j),
                    model->ils_fwhm_unitnum);
        }
        fprintf(stream, "\n");
        if (!(model->ilsmode & ILSMODE_NORMAL)) {
            if (model->ilsmode & ILSMODE_DSB) {
                fprintf(stream, "ilsmode dsb");
            } else if (model->ilsmode & ILSMODE_USB) {
                fprintf(stream, "ilsmode usb");
            } else if (model->ilsmode & ILSMODE_LSB) {
                fprintf(stream, "ilsmode lsb");
            } else {
                fprintf(stream, "ilsmode ???");
            }
            write_variable_info(
                    stream,
                    "",
                    model->ils_fif,
                    &model->ils_fif,
                    model->ils_fif_unitnum,
                    simplex);
            fprintf(stream, "\n");
        }
    }
    /*
     * Receiver noise temperature and gain correction factor
     */
    if (output[OUTPUT_TSYS].flags & OUTPUT_ACTIVE ||
            output[OUTPUT_Y].flags & OUTPUT_ACTIVE) {
        write_variable_info(
                stream,
                "Trx",
                model->Trx,
                &model->Trx,
                model->Trx_unitnum,
                simplex);
        fprintf(stream, "\n");
        write_variable_info(
                stream,
                "rx_gain_factor",
                model->rx_gain_factor,
                &model->rx_gain_factor,
                UNIT_NONE,
                simplex);
        fprintf(stream, "\n");
    }
    /*
     * Reference temperature Tref for Y-factor or I_diff
     */
    if (output[OUTPUT_Y].flags & OUTPUT_ACTIVE ||
            output[OUTPUT_RADIANCE_DIFF].flags & OUTPUT_ACTIVE) {
        write_variable_info(
                stream,
                "Tref",
                model->Tref,
                &model->Tref,
                model->Tref_unitnum,
                simplex);
        fprintf(stream, "\n");
    }
    /*
     * Sideband ratio for DSB ILS or DSB IF spectrum.
     */
    if (model->ilsmode & ILSMODE_DSB || model->ifmode & IFMODE_DSB) {
        write_variable_info(
                stream,
                "dsb_utol_ratio",
                model->dsb_utol_ratio,
                &model->dsb_utol_ratio,
                UNIT_NONE,
                simplex);
        fprintf(stream, "\n");
    }
    /*
     * Tolerances on fast linesum and self-broadening vmr.
     */
    fprintf(stream, "tol %.13g\n", model->tol);
    fprintf(stream, "selfbroad_vmr_tol %.13g\n", model->selfbroad_vmr_tol);
    /*
     * headers, runtime
     */
    if (output[ALL_OUTPUTS].flags & OUTPUT_HEADERS)
        fprintf(stream, "headers 1\n");
    if (model->log_runtimes != 0)
        fprintf(stream, "runtime %d\n", model->log_runtimes);
    fprintf(stream, "\n");
    /*
     * Geometry mode
     */
    if (model->geometry & GEOMETRY_PLANE_PARALLEL)
        fprintf(stream, "geometry plane-parallel\n");
    else if (model->geometry & GEOMETRY_SPHERICAL)
        fprintf(stream, "geometry spherical\n");
    else if (model->geometry & GEOMETRY_LIMB)
        fprintf(stream, "geometry limb\n");
    /*
     * Refraction model
     */
    if (model->geometry & GEOMETRY_REFRACT_NONE)
        fprintf(stream, "refract none\n");
    else if (model->geometry & GEOMETRY_REFRACT_RADIO)
        fprintf(stream, "refract radio\n");
    else if (model->geometry & GEOMETRY_REFRACT_OPTICAL)
        fprintf(stream, "refract optical\n");
    /*
     * zenith angle or secant(zenith angle).  If neither
     * has been user-defined, then the default zenith angle
     * is written, unless we're in limb mode.
     */
    if (model->geometry & GEOMETRY_SEC_ZA_USER_DEFINED) {
        write_variable_info(
                stream,
                "sec_za",
                model->sec_za,
                &model->sec_za,
                UNIT_NONE,
                simplex);
        fprintf(stream, "\n");
    } else if (!(model->geometry & GEOMETRY_LIMB)) {
        write_variable_info(
                stream,
                "za",
                model->za,
                &model->za,
                model->za_unitnum,
                simplex);
        fprintf(stream, "\n");
    }
    /*
     * reverse propagation flag
     */
    if (model->geometry & GEOMETRY_REVERSE)
        fprintf(stream, "reverse 1\n");
    /*
     * PTmode
     */
    write_PTmode(stream, model->PTmode);
    /*
     * Gravitational acceleration and gradient.
     */
    if (model->PTmode & PTMODE_HYDROSTATIC) {
        fprintf(stream, "\n");
        print_with_unit(
                stream,
                "g %.13g %s\n",
                model->g,
                model->g_unitnum);
        print_with_unit(
                stream,
                "dg_dz %.13g %s\n",
                model->dg_dz,
                model->dg_dz_unitnum);
    }
    /*
     * z0 and R0 get written for spherical and limb geometry
     * always, and for plane parallel geometry only if defined by
     * the user.
     */
    if (model->geometry &
            (GEOMETRY_SPHERICAL | GEOMETRY_LIMB |
             GEOMETRY_Z0_USER_DEFINED | GEOMETRY_R0_USER_DEFINED)) {
        fprintf(stream, "\n");
    }
    if (model->geometry &
            (GEOMETRY_SPHERICAL | GEOMETRY_LIMB | GEOMETRY_Z0_USER_DEFINED)) {
        write_variable_info(
                stream,
                "z0",
                model->z0,
                &model->z0,
                model->z0_unitnum,
                simplex);
        fprintf(stream, "\n");
    }
    if (model->geometry &
            (GEOMETRY_SPHERICAL | GEOMETRY_LIMB | GEOMETRY_R0_USER_DEFINED)) {
        write_variable_info(
                stream,
                "R0",
                model->R0,
                &model->R0,
                model->R0_unitnum,
                simplex);
        fprintf(stream, "\n");
    }
    /*
     * RH offset and scale.
     */
    if (model->RH_adj_flag) {
        fprintf(stream, "\n");
        if (isvar(simplex, &(model->RH_offset_exp))) {
            j = get_simplex_variable_index(simplex, &(model->RH_offset_exp));
            fprintf(stream,
                    "RH_offset %.4g%%  %.13g%%",
                    unmap_variable(model->RH_offset_exp, MAPPING_EXP),
                    unmap_differential(simplex->init[j], simplex->scale[j],
                        MAPPING_EXP)
                   );
            if (output[ALL_OUTPUTS].flags & OUTPUT_FITTED) {
                fprintf(stream,
                        "  (%.3g%%)",
                        unmap_differential(model->RH_offset_exp,
                            simplex_variable_range(simplex, j), MAPPING_EXP)
                       );
            }
        } else {
            fprintf(stream, "RH_offset %.13g%%",
                    unmap_variable(model->RH_offset_exp, MAPPING_EXP));
        }
        fprintf(stream, "\n");
        write_variable_info(
                stream,
                "RH_scale",
                model->RH_scale,
                &model->RH_scale,
                UNIT_NONE,
                simplex);
        fprintf(stream, "\n\n");
    }
    /*
     * Nscale
     */
    fprintf(stream, "\n");
    write_Nscale_config(stream, simplex);
    /*
     * Background temperature T0
     */
    write_variable_info(
            stream,
            "T0",
            model->T0,
            &model->T0,
            model->T0_unitnum,
            simplex);
    fprintf(stream, "\n");
    /*
     * User-defined interpolated levels
     */
    if (model->geometry & GEOMETRY_AT_LEAST_ONE_USER_DEFINED_LEVEL)
        fprintf(stream, "\n");
    if (model->geometry & GEOMETRY_POBS_USER_DEFINED) {
        write_variable_info(
                stream,
                "Pobs",
                model->Pobs,
                &model->Pobs,
                model->Pobs_unitnum,
                simplex);
        fprintf(stream, "\n");
    }
    if (model->geometry & GEOMETRY_ZOBS_USER_DEFINED) {
        write_variable_info(
                stream,
                "zobs",
                model->zobs,
                &model->zobs,
                model->zobs_unitnum,
                simplex);
        fprintf(stream, "\n");
    }
    if (model->geometry & GEOMETRY_PSOURCE_USER_DEFINED) {
        write_variable_info(
                stream,
                "Psource",
                model->Psource,
                &model->Psource,
                model->Psource_unitnum,
                simplex);
        if (model->geometry & GEOMETRY_SOURCE_NEAR)
            fprintf(stream, " near\n");
        else
            fprintf(stream, " far\n");
    }
    if (model->geometry & GEOMETRY_ZSOURCE_USER_DEFINED) {
        write_variable_info(
                stream,
                "zsource",
                model->zsource,
                &model->zsource,
                model->zsource_unitnum,
                simplex);
        if (model->geometry & GEOMETRY_SOURCE_NEAR)
            fprintf(stream, " near\n");
        else
            fprintf(stream, " far\n");
    }
    if (model->geometry & GEOMETRY_PTAN_USER_DEFINED) {
        write_variable_info(
                stream,
                "Ptan",
                model->Ptan,
                &model->Ptan,
                model->Ptan_unitnum,
                simplex);
        fprintf(stream, "\n");
    }
    if (model->geometry & GEOMETRY_ZTAN_USER_DEFINED) {
        write_variable_info(
                stream,
                "ztan",
                model->ztan,
                &model->ztan,
                model->ztan_unitnum,
                simplex);
        fprintf(stream, "\n");
    }
    /*
     * Layer data.  Layers above or below the propagation path
     * are printed in abbreviated form.
     */
    if (model->path_begin > 0 && model->path_end > 0)
        verbosity = VERBOSITY_LAYER_NOT_IN_PATH;
    else
        verbosity = VERBOSITY_FULL;
    for (lnum = 0; lnum < model->nlayers; ++lnum) {
        layer_t *layer = model->layer[lnum];
        fprintf(stream, "\n");
        if (verbosity == VERBOSITY_LAYER_NOT_IN_PATH &&
                (lnum == model->path_begin || lnum == model->path_end)) {
            fprintf(stream,
                    "#\n"
                    "# Layers above this line were not"
                        " in the propagation path.\n"
                    "#\n"
                    "\n");
            verbosity = VERBOSITY_FULL;
        }
        if (lnum == model->path_mid + 1) {
            fprintf(stream,
                    "#\n"
                    "# Layers below this line were not"
                        " in the propagation path.\n"
                    "#\n"
                    "\n");
            verbosity = VERBOSITY_LAYER_NOT_IN_PATH;
        }
        if (layer->type != LAYER_TYPE_DEF) {
            /*
             * Interpolated or extrapolated layer specifications are
             * written as comments.
             */
            fprintf(stream, "| layer");
            if (layer->tagnum)
                fprintf(stream, " %s", tag_string(layer->tagnum));
            write_layer_data(stream, "| ", model, simplex, lnum, verbosity); 
        } else if (
                lnum > 0 &&
                model->PTmode & PTMODE_HYDROSTATIC &&
                model->layer[lnum - 1]->type != LAYER_TYPE_DEF
                ) {
            /*
             * This is a parent layer for interpolated layer(s)
             * above.  Print its layer specification twice, first
             * as a commented layer specification with layer data
             * as modified by the interpolation, and second in
             * the unmodified state with no verbose data added.
             * This ensures that feeding the output config data
             * back into am is an indempotent operation.
             */
            fprintf(stream, "| layer");
            if (layer->tagnum)
                fprintf(stream, " %s", tag_string(layer->tagnum));
            write_layer_data(stream, "| ", model, simplex, lnum, verbosity);
            fprintf(stream, "\nlayer");
            if (layer->tagnum)
                fprintf(stream, " %s", tag_string(layer->tagnum));
            fprintf(stream, " | interpolated");
            if (lnum < model->nlayers - 1 &&
                    lnum == get_lnum_by_type(model, LAYER_TYPE_DEF))
                fprintf(stream, " and extrapolated");
            write_layer_data(
                    stream, "", model, simplex, lnum, VERBOSITY_PARENT_LAYER);
        } else if (
                lnum < model->nlayers - 1 &&
                model->PTmode & PTMODE_HYDROSTATIC &&
                lnum == get_lnum_by_type(model, LAYER_TYPE_DEF)
                ) {
            /*
             * This is a parent layer for extrapolated layer(s)
             * only.  This means it has not itself been modified,
             * so it only gets printed once.
             */
            fprintf(stream, "layer");
            if (layer->tagnum)
                fprintf(stream, " %s", tag_string(layer->tagnum));
            fprintf(stream, " | extrapolated");
            write_layer_data(stream, "", model, simplex, lnum, verbosity);
        } else {
            /*
             * Normal layer specification.
             */
            fprintf(stream, "layer");
            if (layer->tagnum)
                fprintf(stream, " %s", tag_string(layer->tagnum));
            write_layer_data(stream, "", model, simplex, lnum, verbosity);
        }
    }
    write_column_density_summary(stream, model);
    fflush(stream);
    return;
}   /* write_model_config_data() */


/***********************************************************
* static void write_column_data(
*     FILE       *stream,
*     model_t    *model,
*     simplex_t  *simplex,
*     const int  lnum,
*     const int  cnum
*     const int  verbosity)
*
* Purpose:
*   Writes a column configuration line, with a prepended
*   prefix string intended for commenting and indentation.
*   The verbosity flag controls printing of derived column
*   densities.
************************************************************/

static void write_column_data(
        FILE       *stream,
        model_t    *model,
        simplex_t  *simplex,
        const int  lnum,
        const int  cnum,
        const int  verbosity)
{
    layer_t  *layer  = model->layer[lnum];
    column_t *column = layer->column[cnum];

    fprintf(stream, "column %s ", col_type[column->col_typenum].name);
    if (column->N_mode & N_BY_VMR &&
            !(column->vmr_stat & VMR_BY_RH)) {
        /*
         * Column density by mixing ratio, but not by relative
         * humidity.
         */
        if (column->vmr_stat & VMR_USER_DEFINED) {
            if (column->N_mode & N_HYDROSTATIC) {
                fprintf(stream,
                        "hydrostatic %.13g",
                        unmap_variable(column->xvmr, MAPPING_VMR));
            } else if (column->N_mode & N_BY_DISTANCE) {
                fprintf(stream,
                        "vmr %.13g",
                        unmap_variable(column->xvmr, MAPPING_VMR));
            }
        } else if (verbosity == VERBOSITY_PARENT_LAYER) {
            /*
             * Automatically-computed mixing ratios are not
             * reported on a parent layer.
             */
            if (column->N_mode & N_HYDROSTATIC)
                fprintf(stream, "hydrostatic");
            else if (column->N_mode & N_BY_DISTANCE)
                fprintf(stream, "vmr");
        } else {
            if (column->N_mode & N_HYDROSTATIC) {
                fprintf(stream,
                        "hydrostatic (%g)",
                        unmap_variable(column->xvmr, MAPPING_VMR));
            } else if (column->N_mode & N_BY_DISTANCE) {
                fprintf(stream,
                        "vmr (%g)",
                        unmap_variable(column->xvmr, MAPPING_VMR));
            }
        }
        if (isvar(simplex, &(column->xvmr))) {
            int j = get_simplex_variable_index(simplex, &(column->xvmr));
            fprintf(stream,
                    "  %.13g",
                    unmap_differential(simplex->init[j], simplex->scale[j],
                        MAPPING_VMR)
                   );
            if (output[ALL_OUTPUTS].flags & OUTPUT_FITTED) {
                fprintf(stream,
                        "  (%.3g)",
                        unmap_differential(
                            simplex_variable_range(simplex, j),
                            column->xvmr, MAPPING_VMR)
                       );
            }
        }
        if (verbosity == VERBOSITY_FULL) {
            print_with_unit(
                    stream,
                    "  (%g %s)",
                    column->N,
                    column->N_unitnum);
        }
    } else if (column->N_mode & N_BY_VMR &&
            column->vmr_stat & VMR_BY_RH) {
        /*
         * Column density by RH
         */
        if (column->N_mode & N_INTERPOLATED &&
                verbosity != VERBOSITY_PARENT_LAYER) {
            /*
             * For an interpolated or extrapolated layer,
             * including a split parent layer, do not report RH.
             * This is because RH is converted to vmr at the
             * parent layer midpoint temperature before
             * interpolation or extrapolation, which is then done
             * holding composition constant.  (However, when the
             * original parent layer data are printed with
             * VERBOSITY_PARENT_LAYER, then we do want to see the
             * RH to vmr computation.)
             */
            fprintf(stream,
                    "vmr %.4g",
                    unmap_variable(column->xvmr, MAPPING_VMR));
            if (verbosity == VERBOSITY_FULL) {
                print_with_unit(
                        stream,
                        "  (%g %s)",
                        column->N,
                        column->N_unitnum);
            }
        } else {
            fprintf(stream,
                    "%s %g%%",
                    column->vmr_stat & VMR_BY_RH_LIQUID ? "RH" : "RHi",
                    column->RH);
            if (isvar(simplex, &(column->RH))) {
                int j = get_simplex_variable_index(simplex, &(column->RH));
                fprintf(stream,
                        "  %.13g%%",
                        simplex->scale[j]);
                if (output[ALL_OUTPUTS].flags & OUTPUT_FITTED) {
                    fprintf(stream,
                            "  (%.3g)",
                            simplex_variable_range(simplex, j));
                }
            }
            if (model->RH_adj_flag) {
                fprintf(stream,
                        "  (->%.4g%%  vmr %.4g)",
                        column->RH_adj,
                        unmap_variable(column->xvmr, MAPPING_VMR));
            } else {
                fprintf(stream,
                        "  (vmr %.4g)",
                        unmap_variable(column->xvmr, MAPPING_VMR));
            }
            if (verbosity == VERBOSITY_FULL) {
                print_with_unit(
                        stream,
                        "  (%g %s)",
                        column->N,
                        column->N_unitnum);
            }
        }
    } else {
        /*
         * Column density defined explicitly.
         */
        if (!(column->N_mode & N_INTERPOLATED) ||
                verbosity == VERBOSITY_PARENT_LAYER) {
            /*
             * If the layer was not split by an interpolated
             * level, or if we're reporting the defining parent
             * layer data, report the defining column density
             * N_def, which might be a fit or differentiation
             * variable.
             */
            if (isvar(simplex, &(column->N_def))) {
                int j = get_simplex_variable_index(simplex, &(column->N_def));
                print_with_unit(
                        stream,
                        "%g %s",
                        column->N_def,
                        column->N_unitnum);
                print_differential_with_unit(
                        stream,
                        "  %.13g %s",
                        simplex->scale[j],
                        column->N_unitnum);
                if (output[ALL_OUTPUTS].flags & OUTPUT_FITTED) {
                    print_differential_with_unit(
                            stream,
                            "  (%.3g %s)",
                            simplex_variable_range(simplex, j),
                            column->N_unitnum);
                }
            } else {
                print_with_unit(
                        stream,
                        "%.13g %s",
                        column->N_def,
                        column->N_unitnum);
            }
        } else {
            /*
             * Otherwise, this was an interpolated level.
             * Report N derived by splitting N_def across
             * levels.
             */
            print_with_unit(
                    stream,
                    "%.13g %s",
                    column->N,
                    column->N_unitnum);
        }
        /*
         * When a vmr was computed from an explicit column
         * density, provide it as a comment.
         */
        if (verbosity == VERBOSITY_FULL &&
                !(col_type[column->col_typenum].flags & COL_PARAMETRIC) &&
                (model->PTmode & PTMODE_HYDROSTATIC || layer->h >= 0.0)) {
            fprintf(stream,
                    "  (vmr %.4g)",
                    unmap_variable(column->xvmr, MAPPING_VMR));
        }
    }
    if (find_Nscale_list_entry(column->col_typenum, layer->tagnum) != NULL ||
            find_Nscale_list_entry(column->col_typenum, 0) != NULL) {
        fprintf(stream,
                "  (* %#g)",
                lookup_Nscale(column->col_typenum, layer->tagnum));
    }
    if (model->log_runtimes && verbosity == VERBOSITY_FULL)
        fprintf(stream, " ; %g us", 1.0e6 * column->runtime);
    fprintf(stream, "\n");
    if (column->unres_lines && verbosity == VERBOSITY_FULL) {
        fprintf(stream,
                "! Warning: Column included %d unresolved line%s.\n",
                column->unres_lines,
                column->unres_lines == 1 ? "" : "s");
    }
    /*
     * Report column mixing ratio computation problems.
     */
    if (column->vmr_stat & VMR_ERROR && verbosity != VERBOSITY_PARENT_LAYER) {
        if (column->vmr_stat & VMR_ZERO_DEFAULT) {
            fprintf(stream,
                    "! Warning: The volume mixing ratio for this column is"
                    " 0.0 by default.\n");
        } else if (column->vmr_stat & VMR_DEFAULT_ON_SCALED_COLUMN) {
            fprintf(stream,
                    "! Warning: The default volume mixing ratio is being"
                    " scaled for this column\n");
        }
    }
}   /* write_column_data() */


/***********************************************************
* static void write_column_density_summary(FILE *stream, model_t *model)
*
* Purpose:
*   Writes a summary of total column density by species
*   and layer tag, and for the entire model.  Also reports
*   airmass, nominal total airmass, total distance, and
*   total refraction when appropriate.
************************************************************/

static void write_column_density_summary(FILE *stream, model_t *model)
{
    int *tag_ord   = NULL;
    int *tag_hit   = NULL;
    int path_begin = model->path_begin;
    int path_mid   = model->path_mid;
    int path_end   = model->path_end;
    int path_min   = model->path_min;
    int lnum;
    int i, ntags;

    /*
     * Check that there is at least one defined column in the range of
     * layers covered by the propagation path.  If not, no column
     * density summary table will be printed.
     */
    for (lnum = path_min; lnum <= path_mid; ++lnum) {
        if (model->layer[lnum]->ncols)
            break;
    }
    if (lnum > path_mid)
        return;
    /*
     * The tag string table is ordered by first appearance of tag
     * strings anywhere in the configuration file, including in
     * Nscale statements and in layer statements.  However, the
     * column density summary table is easier to interpret if it
     * is written in model layer order.
     *
     * To do this, here we make a table of tag string index
     * numbers in model layer order.  Tag index 0, which stands
     * for all layers in the model, goes at the end of the tag
     * order table so totals for all layers will be printed last.
     */
    ntags = get_num_tag_strings();
    if ((tag_ord = (int*)malloc(ntags * sizeof(int))) == NULL ||
            (tag_hit = (int*)malloc(ntags * sizeof(int))) == NULL) {
        free(tag_ord);
        free(tag_hit);
        fprintf(stream,
                "# malloc() error generating column density summary table\n");
        return;
    }
    for (i = 0; i < ntags; ++i) {
        tag_ord[i] = 0;
        tag_hit[i] = 0;
    }
    /*
     * Re-count ntags as we go through the model layers, since we
     * don't want to include tags that occurred only in Nscale
     * statements.
     */
    ntags = 0;
    for (lnum = 0; lnum < model->nlayers; ++lnum) {
        int tagnum = model->layer[lnum]->tagnum;
        if (tagnum != 0 && !tag_hit[tagnum])
            tag_ord[ntags++] = tagnum;
        /*
         * tag_hit[0] counts all layers, tagged and untagged.
         */
        ++tag_hit[0];
        if (tagnum != 0)
            ++tag_hit[tagnum];
    }
    tag_ord[ntags++] = 0;
    /*
     * For each layer tag, add up the zenith and line-of-sight
     * column densities by column type.
     */
    for (i = 0; i < ntags; ++i) {
        double N[COL_TYPE_END_OF_TABLE];
        double N_los[COL_TYPE_END_OF_TABLE];
        double N_all, N_los_all;
        int los_col_def_flag[COL_TYPE_END_OF_TABLE];
        int zenith_col_def_flag[COL_TYPE_END_OF_TABLE];
        int ctype;
        int layer_count_down   = 0;
        int layer_count_up     = 0;
        int layer_count_zenith = 0;
        int tagnum = tag_ord[i];
        for (ctype = 0; ctype < COL_TYPE_END_OF_TABLE; ++ctype) {
            N[ctype]     = 0.0;
            N_los[ctype] = 0.0;
            los_col_def_flag[ctype]    = 0;
            zenith_col_def_flag[ctype] = 0;
        }
        N_all     = 0.0;
        N_los_all = 0.0;
        /*
         * Path airmass for a tagged set of layers is normalized
         * to zenith column density for that set of layers.  The
         * zenith totals include all layers between source and
         * observer.  If the path includes a tangent point, any
         * layers that are traversed twice are counted only once
         * in the zenith column density sum.
         */
        for (lnum = path_mid; lnum >= path_min; --lnum) {
            layer_t *layer = model->layer[lnum];
            int cnum;
            if (tagnum == 0 || tagnum == layer->tagnum) {
                ++layer_count_zenith;
                for (cnum = 0; cnum < layer->ncols; ++cnum) {
                    column_t *column = layer->column[cnum];
                    ctype = column->col_typenum;
                    ++zenith_col_def_flag[ctype];
                    if (!(col_type[ctype].flags & COL_PARAMETRIC)) {
                        N[ctype] += column->N_scaled;
                        N_all    += column->N_scaled;
                    }
                }
            }
        }
        for (lnum = path_begin; lnum <= path_mid; ++lnum) {
            layer_t *layer = model->layer[lnum];
            int cnum;
            if (tagnum == 0 || tagnum == layer->tagnum) {
                ++layer_count_down;
                for (cnum = 0; cnum < layer->ncols; ++cnum) {
                    column_t *column = layer->column[cnum];
                    ctype = column->col_typenum;
                    ++los_col_def_flag[ctype];
                    if (!(col_type[ctype].flags & COL_PARAMETRIC)) {
                        N_los[ctype] += column->N_scaled * layer->m;
                        N_los_all    += column->N_scaled * layer->m;
                    }
                }
            }
        }
        for (lnum = path_mid; lnum >= path_end; --lnum) {
            layer_t *layer = model->layer[lnum];
            int cnum;
            if (tagnum == 0 || tagnum == layer->tagnum) {
                ++layer_count_up;
                for (cnum = 0; cnum < layer->ncols; ++cnum) {
                    column_t *column = layer->column[cnum];
                    ctype = column->col_typenum;
                    ++los_col_def_flag[ctype];
                    if (!(col_type[ctype].flags & COL_PARAMETRIC)) {
                        N_los[ctype] += column->N_scaled * layer->m;
                        N_los_all    += column->N_scaled * layer->m;
                    }
                }
            }
        }
        /*
         * If this is the first tag in the table, print a header
         */
        if (i == 0) {
            fprintf(stream, "\n");
            fprintf(stream,
                    "#\n"
                    "# column densities [%s], including scale factors\n",
                    unit_tab[AM_UNIT_COL_DENSITY].name);
            if (model->geometry & GEOMETRY_DISPLAY_AIRMASS) {
                int obs_lnum, source_lnum, tan_lnum;
                if (model->geometry & GEOMETRY_REVERSE) {
                    source_lnum = get_obs_lnum(model);
                    obs_lnum    = get_source_lnum(model);
                } else {
                    obs_lnum    = get_obs_lnum(model);
                    source_lnum = get_source_lnum(model);
                }
                tan_lnum = get_tan_lnum(model);
                fprintf(stream,
                        "#\n"
                        "# airmass is line-of-sight / (zenith");
                if (path_min == 0) {
                    fprintf(stream, " above");
                } else if (path_min - 1 == obs_lnum) {
                    fprintf(stream, " between observing and");
                } else if (path_min - 1 == source_lnum) {
                    fprintf(stream, " between source and");
                } else {
                    fprintf(stream, " !error in %s at line %d",
                        __FILE__, __LINE__);
                }
                if (path_mid == obs_lnum) {
                    fprintf(stream, " observing level)\n");
                } else if (path_mid == source_lnum) {
                    fprintf(stream, " source level)\n");
                } else if (path_mid == tan_lnum) {
                    fprintf(stream, " tangent level)\n");
                } else {
                    fprintf(stream, " !error in %s at line %d)\n",
                        __FILE__, __LINE__);
                }
            }
            fprintf(stream,
                    "#\n"
                    "# %24s  %-18s %-18s",
                    "", "zenith", "line-of-sight");
            if (model->geometry & GEOMETRY_DISPLAY_AIRMASS)
                fprintf(stream, " %-8s", "airmass");
            fprintf(stream, "\n");
        }
        /*
         * Report the layer count for this tag
         */
        fprintf(stream,
                "# %s (%d layer%s",
                tagnum == 0 ? "total" : tag_string(tagnum),
                tag_hit[tagnum],
                tag_hit[tagnum] == 1 ? "" : "s");
        /*
         * If the propagation path does anything other than traverse
         * all the model layers once from end to end, also report
         * extra layer traversal information.
         */
        if (    (path_mid   != model->nlayers - 1)          ||
                (path_begin <= path_mid && path_begin != 0) ||
                (path_end   <= path_mid && path_end   != 0)
                ) {
            fprintf(stream,
                    ")\n#  (%d zenith, %d down, %d up",
                    layer_count_zenith,
                    layer_count_down,
                    layer_count_up);
        }
        fprintf(stream, "):\n");
        for (ctype = 1; ctype < COL_TYPE_END_OF_TABLE; ++ctype) {
            if (zenith_col_def_flag[ctype] == 0)
                continue;
            fprintf(stream, "# %24s", col_type[ctype].name);
            if (col_type[ctype].flags & COL_PARAMETRIC) {
                fprintf(stream, "  %-18s\n", "parametric");
            } else {
                fprintf(stream, "  %-18e", N[ctype]);
                if (los_col_def_flag[ctype]) {
                    fprintf(stream, " %-18e", N_los[ctype]);
                    if (model->geometry & GEOMETRY_DISPLAY_AIRMASS &&
                            N[ctype] > 0.0) {
                        double airmass = N_los[ctype] / N[ctype];
                        fprintf(stream, " %-8g", airmass);
                    }
                }
                fprintf(stream, "\n");
                if (col_type[ctype].common_unit != AM_UNIT_COL_DENSITY) {
                    char s[64];
                    snprint_with_unit(
                            s,
                            sizeof(s),
                            "(%g %.32s)",
                            N[ctype],
                            col_type[ctype].common_unit);
                    fprintf(stream, "# %24s   %-18s", "", s);
                    if (los_col_def_flag[ctype]) {
                        snprint_with_unit(
                                s,
                                sizeof(s),
                                "(%g %.32s)",
                                N_los[ctype],
                                col_type[ctype].common_unit);
                        fprintf(stream, " %-18s", s);
                    }
                    fprintf(stream, "\n");
                }
            }
        }
        if (model->geometry & GEOMETRY_DISPLAY_AIRMASS && tagnum == 0) {
            fprintf(stream, "#\n");
            fprintf(stream, "# %24s  %-18e %-18e",
                    "all", N_all, N_los_all);
            if (N_all > 0.0) {
                double airmass = N_los_all / N_all;
                fprintf(stream, " %-8g", airmass);
            }
            fprintf(stream, "\n");
        }
        fprintf(stream, "#\n");
    }
    /*
     * For finite-length paths without an endpoint at infinity,
     * report the refracted and geometric path distances.  Units
     * are picked up from model->R0;
     */
    if (    model->PTmode & PTMODE_HYDROSTATIC &&
            path_begin > 0 &&
            path_end > 0) {
        fprintf(stream, "# source to observer distance\n");
        print_with_unit(
                stream,
                "#  refracted : %g %s\n",
                source_to_obs_path_distance(model),
                model->R0_unitnum);
        print_with_unit(
                stream,
                "#  geometric : %g %s\n",
                source_to_obs_geometric_distance(model),
                model->R0_unitnum);
        print_with_unit(
                stream,
                "#  projected : %g %s\n",
                source_to_obs_projected_distance(model),
                model->R0_unitnum);
        fprintf(stream, "#\n");
    }
    if (model->geometry & GEOMETRY_DISPLAY_REFRACTION) {
        print_with_unit(
                stream,
                "# total refraction : %g %s\n",
                total_refraction(model),
                UNIT_ARCSEC);
        fprintf(stream, "#\n");
    }
    free(tag_ord);
    free(tag_hit);
    return;
}   /* write_column_density_summary() */


/***********************************************************
* static void write_fit_data_delimiters(FILE *stream)
*
* Purpose:
*   Writes a fit_data_delimiter config statement in the
*   form:
*
*     fit_data_delimiters "string" # "string_comment"
*
*   where string is the (quoted) delimiter string, excluding
*   the "\n\r" newline and carriage return characters which
*   are always appended internally.  If the delimiter string
*   includs a non-printing tab, a comment is included to
*   render the tab visible as \t.
************************************************************/

static void write_fit_data_delimiters(FILE *stream)
{
    int i;
    int hidden_char = 0;
    char *str = fit_data_delimiters(NULL);
    int len   = (int)strlen(str);

    fprintf(stream, "fit_data_delimiters \"");
    for (i = 0; i < len; ++i) {
        if (isprint(str[i])) {
            fprintf(stream, "%c", str[i]);
        } else if (str[i] == '\t') {
            hidden_char = 1;
            fprintf(stream, "%c", str[i]);
        }
    }
    if (hidden_char) {
        fprintf(stream, "\" # \"");
        for (i = 0; i < len; ++i) {
            if (isprint(str[i])) {
                fprintf(stream, "%c", str[i]);
            } else if (str[i] == '\t') {
                fprintf(stream, "\\t");
            }
        }
    }
    fprintf(stream, "\"\n");
    return;
}   /* write_fit_data_delimiters() */


/***********************************************************
* static void write_fit_variable_summary(FILE *stream, simplex_t *simplex)
*
* Purpose:
*   Writes a summary of fit variables by name and converged
*   value.
************************************************************/

static void write_fit_variable_summary(FILE *stream, simplex_t *simplex)
{
    unsigned int i;

    fprintf(stream, "#\n");
    fprintf(stream, "# %3s %15s %18s %8s  %s\n",
            "var",
            "init",
            "converged",
            "units",
            "description");
    for (i = 0; i < simplex->n; ++i) {
        fprintf(stream, "# %3d", i);
        print_with_unit(
                stream,
                " %15.13g%.0s",
                unmap_variable(simplex->init[i], simplex->mapping[i]),
                simplex->unitnum[i]);
        print_with_unit(
                stream,
                " %18g %8s",
                unmap_variable(*simplex->varptr[i], simplex->mapping[i]),
                simplex->unitnum[i]);
        fprintf(stream, "  %s\n", simplex->name[i]);
    }
    fprintf(stream, "#\n");
    return;
}   /* write_fit_variable_summary() */


/***********************************************************
* void write_fit_residuals(FILE *stream, fit_data_t *fit_data)
*
* Purpose:
*   Writes the fit residuals to a stream.  Frequency
*   is converted to same units used for model output
*   spectra.
************************************************************/

void write_fit_residuals(FILE *stream, fit_data_t *fit_data)
{
    gridsize_t i;
    int f_unitnum, s_unitnum;

    f_unitnum = output[OUTPUT_FREQUENCY].unitnum;
    s_unitnum = output[fit_data->data_type].unitnum;
    if (fit_data->res_track_gain < 0.0) { /* residual tracking off */
        for (i = 0; i < fit_data->npts; ++i) {
            fprintf(stream,
                    "%e % e %e %e %e\n",
                    fit_data->f[i] / unit_tab[f_unitnum].factor -
                        unit_tab[f_unitnum].offset,
                    fit_data->res[i] / unit_tab[s_unitnum].factor,
                    fit_data->s_mod[i] / unit_tab[s_unitnum].factor -
                        unit_tab[s_unitnum].offset,
                    fit_data->b[i] / unit_tab[f_unitnum].factor -
                        unit_tab[f_unitnum].offset,
                    fit_data->w[i] * unit_tab[s_unitnum].factor);
        }
    } else {
        for (i = 0; i < fit_data->npts; ++i) { /* tracking on */
            fprintf(stream,
                    "%e % e %e %e %e % e\n",
                    fit_data->f[i] / unit_tab[f_unitnum].factor -
                        unit_tab[f_unitnum].offset,
                    fit_data->res[i] / unit_tab[s_unitnum].factor,
                        fit_data->s_mod[i] / unit_tab[s_unitnum].factor -
                    unit_tab[s_unitnum].offset,
                        fit_data->b[i] / unit_tab[f_unitnum].factor -
                    unit_tab[f_unitnum].offset,
                        fit_data->w[i] * unit_tab[s_unitnum].factor,
                    fit_data->res_est[i] / unit_tab[s_unitnum].factor);
        }
    }
    fflush(stream);
    return;
}   /* write_fit_residuals() */


/***********************************************************
* static void write_jacobian_differentiation_variables_summary(
*   FILE *stream,
*   simplex_t *simplex)
*
* Purpose:
*   Writes a summary table of the user-defined differentiation
*   variables for Jacobians.
************************************************************/

static void write_jacobian_differentiation_variables_summary(
        FILE *stream,
        simplex_t *simplex)
{
    unsigned int i;

    fprintf(stream,
            "# Differentiation variable%s\n"
            "#\n" , simplex->n == 1 ? "" :  "s");
    if (simplex->n < 1) {
        fprintf(stream, "# none specified\n");
        fprintf(stream, "#\n");
        return;
    }
    fprintf(stream, "# %3s %15s %18s %8s  %s\n",
            "var",
            "value",
            "dvalue",
            "units",
            "description");
    for (i = 0; i < simplex->n; ++i) {
        fprintf(stream, "# %3d", i);
        print_with_unit(
                stream,
                " %15.13g%.0s",
                unmap_variable(*simplex->varptr[i], simplex->mapping[i]),
                simplex->unitnum[i]);
        print_differential_with_unit(
                stream,
                " %18g %8s",
                unmap_differential(*simplex->varptr[i], simplex->scale[i],
                    simplex->mapping[i]),
                simplex->unitnum[i]);
        fprintf(stream, "  %s\n", simplex->name[i]);
    }
    fprintf(stream, "#\n");
    return;
}   /* write_jacobian_differentiation_variables_summary() */


/***********************************************************
* static void write_jacobian_error_estimate_summary(
*   FILE *stream,
*   model_t *model,
*   simplex_t *simplex)
*   
* Purpose:
*   Writes a statistical summary table of numerical
*   derivative error estimates for Jacobians, by array
*   differentiation variable.  If no output spectra were
*   actually computed, as when am is invoked with the -a 
*   option, then this function does nothing.
************************************************************/

static void write_jacobian_error_estimate_summary(
        FILE *stream,
        model_t *model,
        simplex_t *simplex)
{
    unsigned int i, j, k;
    static const double percentile[] = {0, 5, 25, 50, 75, 95, 100};
    gridsize_t npts = (model->ifmode ? model->nif : model->ngrid);


    if (output[ALL_OUTPUTS].flags & OUTPUT_AM_ONLY)
        return;
    fprintf(stream,
            "#\n"
            "# Relative error estimates for numerical"
                " derivatives (percentiles)\n"
            "#\n"
            "# %6s %4s", "array", "var");
    for (i = 0; i < sizeof(percentile) / sizeof(double); ++i)
        fprintf(stream, "%8.0f", percentile[i]);
    fprintf(stream, "\n#\n");
    for (k = 0; outcol[k] != 0; ++k) {
        if (!(output[outcol[k]].flags & OUTPUT_JACOBIAN))
            continue;
        for (j = 0; j < simplex->n; ++j) {
            fprintf(stream, "# %6s %4d", j ? "" : output[outcol[k]].name, j);
            for (i = 0; i < sizeof(percentile) / sizeof(double); ++i) {
                double p = (double)npts * percentile[i] / 100.;
                double p_lo = floor(p);
                double p_hi = p_lo + 1.0;
                double q;
                gridsize_t i_lo = (gridsize_t)p_lo;
                gridsize_t i_hi = (gridsize_t)p_hi;
                if (p_hi < npts) {
                    q  = (p_hi - p) *
                     output[outcol[k]].jacobian_sorted_err[j][i_lo];
                    q += (p - p_lo) *
                     output[outcol[k]].jacobian_sorted_err[j][i_hi];
                } else {
                    q = output[outcol[k]].jacobian_sorted_err[j][npts - 1];
                }
                fprintf(stream, "%8.0e", q);
            }
            fprintf(stream, "\n");
        }
        fprintf(stream, "#\n");
    }
    return;
} /* write_jacobian_error_estimate_summary() */


/***********************************************************
* static void write_layer_data(
*     FILE       *stream,
*     const char *prefix,
*     model_t    *model,
*     simplex_t  *simplex,
*     const int  lnum
*     const int  verbosity)
*
* Purpose:
*   Writes layer configuration lines, with a prepended
*   prefix string intended for commenting and indentation.
*   The verbosity flag controls printing of derived
*   variables, so that an abbreviated layer specification
*   can be printed for layers that are not in the
*   propagation path.
*
*   This function assumes that the layer keyword, layer tag,
*   and any additional comments on the same line have
*   already been printed, without a newline.
************************************************************/

static void write_layer_data(
        FILE       *stream,
        const char *prefix,
        model_t    *model,
        simplex_t  *simplex,
        const int  lnum,
        const int  verbosity)
{
    layer_t *layer = model->layer[lnum];
    int cnum;
    int above_toa = model->PTmode & PTMODE_PBASE && layer->Pbase == 0.0;

    if (model->PTmode & PTMODE_HYDROSTATIC &&
            layer->h > 0.0)
        fprintf(stream, " (nonhydrostatic)");
    if (model->log_runtimes)
        fprintf(stream, " ; %g us", 1.0e6 * layer->runtime);
    fprintf(stream, "\n");
    if (layer->vmr_stat & VMR_SUM_EXCEEDS_UNITY)
        fprintf(stream,
                "! Warning: The sum of volume mixing ratios on this layer"
                " exceeds unity.\n");
    /*
     * Layer airmass is omitted for layers outside the
     * propagation path.
     */
    if (model->geometry & GEOMETRY_DISPLAY_AIRMASS &&
            verbosity == VERBOSITY_FULL) {
        fprintf(stream, "%s# airmass %g\n", prefix, layer->m); 
    }
    /* Base zenith angle is omitted for layers outside the
     * propagation path, unless the base level is a source or
     * observing level.
     */
    if (model->geometry & GEOMETRY_DISPLAY_ZA_BASE &&
            (verbosity   == VERBOSITY_FULL ||
             layer->type == LAYER_TYPE_OBS ||
             layer->type == LAYER_TYPE_SOURCE)) {
        fprintf(stream, "%s", prefix);
        print_with_unit(
                stream,
                "# za_base %g %s\n",
                layer->za_base,
                UNIT_DEG); 
    }
    /*
     * P
     */
    if (model->PTmode & PTMODE_P) {
        fprintf(stream, "%s", prefix);
        write_variable_info(
                stream,
                "P",
                layer->P,
                &layer->P,
                layer->P_unitnum,
                simplex);
        fprintf(stream, "\n");
    } else if (verbosity == VERBOSITY_FULL && !above_toa) {
        fprintf(stream, "%s", prefix);
        print_with_unit(
                stream,
                "# P %g %s\n",
                layer->P,
                layer->P_unitnum);
    }
    /*
     * T
     */
    if (model->PTmode & PTMODE_T) {
        fprintf(stream, "%s", prefix);
        write_variable_info(
                stream,
                "T",
                layer->T,
                &layer->T,
                layer->T_unitnum,
                simplex);
        fprintf(stream, "\n");
    } else if (verbosity == VERBOSITY_FULL && !above_toa) {
        fprintf(stream, "%s", prefix);
        print_with_unit(
                stream,
                "# T %g %s\n",
                layer->T,
                layer->T_unitnum);
    }
    /*
     * dP, gbase, and Pbase only get printed in hydrostatic PTmodes.
     */
    if (model->PTmode & PTMODE_HYDROSTATIC) {
        /*
         * Derived dP is omitted on parent layers and if Pbase == 0.
         */
        if (model->PTmode & PTMODE_DP) {
            fprintf(stream, "%s", prefix);
            if (layer->type == LAYER_TYPE_DEF) {
                write_variable_info(
                        stream,
                        "dP",
                        layer->dP_def,
                        &layer->dP_def,
                        layer->P_unitnum,
                        simplex);
            } else {
                print_with_unit(
                        stream,
                        "dP %g %s\n",
                        layer->dP,
                        layer->P_unitnum);
            }
            fprintf(stream, "\n");
        } else if (verbosity != VERBOSITY_PARENT_LAYER && !above_toa) {
            fprintf(stream, "%s", prefix);
            print_with_unit(
                    stream,
                    "# dP %g %s\n",
                    layer->dP,
                    layer->P_unitnum);
        }
        /*
         * gbase is printed only if dg_dz is nonzero, otherwise
         * it would be repetitious.  Also omit if Pbase <= 0.
         */
        if (    model->dg_dz != 0.0 &&
                verbosity != VERBOSITY_PARENT_LAYER &&
                layer->Pbase > 0.0) {
            fprintf(stream, "%s", prefix);
            print_with_unit(
                    stream,
                    "# gbase %g %s\n",
                    layer->gbase,
                    model->g_unitnum);
        }
        /*
         * Pbase is always printed for hydrostatic models, along
         * with z when appropriate.
         */
        fprintf(stream, "%s", prefix);
        if (model->PTmode & PTMODE_PBASE) {
            print_with_unit(
                    stream,
                    layer->type == LAYER_TYPE_DEF ?
                    "Pbase %.13g %s" : "Pbase %g %s",
                    layer->Pbase,
                    layer->P_unitnum);
        } else {
            print_with_unit(
                    stream,
                    "# Pbase %g %s",
                    layer->Pbase,
                    layer->P_unitnum);
        }
        /*
         * Except on parent layers, layer base height zbase is
         * always printed as a comment for spherical or limb
         * geometry.  For plane-parallel geometry, it is printed
         * only if the user has explicitly specified z0 or
         * defined a level by height.
         */
        if (verbosity != VERBOSITY_PARENT_LAYER &&
                (model->geometry & (GEOMETRY_SPHERICAL | GEOMETRY_LIMB) ||
                 model->geometry & GEOMETRY_AT_LEAST_ONE_USER_DEFINED_Z)) {
            if (layer->zbase < DBL_MAX) {
                print_with_unit(
                        stream,
                        "  # z = %.3f %s",
                        layer->zbase,
                        model->z0_unitnum);
            } else {
                fprintf(stream, "  # layer at infinity");
            }
        }
        /*
         * Label source, observing, and tangent levels with
         * comments.
         */
        if (!(model->geometry & GEOMETRY_LIMB)) {
            if (lnum == get_obs_lnum(model)) {
                if (model->geometry & GEOMETRY_REVERSE)
                    fprintf(stream, " (source level)");
                else
                    fprintf(stream, " (observing level)");
            }
            if (lnum == get_source_lnum(model)) {
                if (model->geometry & GEOMETRY_REVERSE)
                    fprintf(stream, " (observing level)");
                else
                    fprintf(stream, " (source level)");
            }
        }
        if (lnum == get_tan_lnum(model))
            fprintf(stream, " (tangent level)");
        fprintf(stream, "\n");
    }
    /*
     * Tbase only gets printed if Tbase is specified in the
     * current PTmode.  For interpolated layers, Tbase is a
     * derived quantity and is printed at default %g precision.
     * However, if the interpolated layer is at zero pressure,
     * as for layers generated by source, observing, or tangent
     * levels specified by height outside the atmosphere, Tbase
     * is omitted.
     */
    if (model->PTmode & PTMODE_TBASE) {
        if (layer->type == LAYER_TYPE_DEF) {
            fprintf(stream, "%s", prefix);
            write_variable_info(
                    stream,
                    "Tbase",
                    layer->Tbase,
                    &layer->Tbase,
                    layer->T_unitnum,
                    simplex);
            fprintf(stream, "\n");
        } else if (!above_toa) {
            fprintf(stream, "%s", prefix);
            print_with_unit(
                    stream,
                    "Tbase %g %s\n",
                    layer->Tbase,
                    layer->T_unitnum);
        }
    }
    /*
     * h only gets printed if it has been defined for this layer.
     */
    if (layer->h >= 0.0) {
        fprintf(stream, "%s", prefix);
        write_variable_info(
                stream,
                "h",
                layer->h,
                &layer->h,
                layer->h_unitnum,
                simplex);
        fprintf(stream, "\n");
    }   
    /*
     * Mair statement only gets printed for hydrostatic models,
     * and only on layers with dP > 0.0.  The derived value for
     * Mair is not printed for parent layers.
     */
    if (model->PTmode & PTMODE_HYDROSTATIC &&
            layer->dP > 0.0) {
        int flag = 0;
        int i;
        fprintf(stream, "%sMair %.13g", prefix, layer->M0);
        for (i = 1; i < COL_TYPE_END_OF_TABLE; ++i) {
            if (layer->Mair_flag[i]) {
                flag = 1;
                fprintf(stream, " %s", col_type[i].name);
            }
        }
        if (verbosity == VERBOSITY_PARENT_LAYER) {
            fprintf(stream, "\n");
        } else {
            if (flag)
                fprintf(stream, "  (%g)", layer->Mair);
            fprintf(stream, "\n");
            if (flag && layer->vmr_stat & VMR_SUM_EXCEEDS_UNITY) {
                fprintf(stream,
                        "! Warning: vmr sum exceeds unity.  Mair adjustment"
                        " may be inaccurate.\n");
            }
            if (flag && layer->vmr_stat & VMR_WEIGHTED_MEAN_MASS_UNDEFINED) {
                fprintf(stream,
                        "! Warning: Mair adjustment undefined.  Check for"
                        " related warnings.\n");
            }
        }
    }
    /*
     * Lineshape configuration data for this layer.
     */
    write_layer_lineshape_config(stream, prefix, layer);
    /*
     * Column data
     */
    for (cnum = 0; cnum < layer->ncols; ++cnum) {
        fprintf(stream, "%s", prefix);
        write_column_data(
                stream, model, simplex, lnum, cnum, verbosity);
    }
}   /* write_layer_data() */


/***********************************************************
* static void write_layer_lineshape_config(
*     FILE *stream, const char *prefix, layer_t *layer)
*
* Purpose:
*   Writes lineshape configuration data for a layer, with a
*   prepended prefix string intended for commenting and
*   indentation
************************************************************/

static void write_layer_lineshape_config(
        FILE *stream, const char *prefix, layer_t *layer)
{
    int cnum, knum;
    int lineshape, strict_selfbroad;
    int ktab[K_TYPE_END_OF_TABLE];

    for (knum = 0; knum < K_TYPE_END_OF_TABLE; ++knum)
        ktab[knum] = 0;
    for (cnum = 0; cnum < layer->ncols; ++cnum) {
        column_t *column;
        column = layer->column[cnum];
        for (knum = 0; knum < column->n_abscoeffs; ++knum) {
            int k_typenum;
            k_typenum = column->abscoeff[knum]->k_typenum;
            ktab[k_typenum] = k_type[k_typenum].dep_flags & DEP_ON_LSHAPE;
        }
    }
    for (lineshape = 0; lineshape < LINESHAPE_END_OF_TABLE; ++lineshape) {
        for (strict_selfbroad = 0;
                strict_selfbroad <= 1;
                ++strict_selfbroad) {
            int seen = 0;
            for (knum = 0; knum < K_TYPE_END_OF_TABLE; ++knum) {
                if (ktab[knum] &&
                        layer->lineshape[knum] == lineshape &&
                        layer->strict_selfbroad[knum] == strict_selfbroad) {
                    if (!seen) {
                        fprintf(stream, "%slineshape %s%s",
                                prefix,
                                lineshape_type[lineshape].name,
                                (strict_selfbroad ?
                                    " strict_self_broadening" : ""));
                        seen = 1;
                    }
                    fprintf(stream, " %s", k_type[knum].name);
                }
            }
            if (seen)
                fprintf(stream, "\n");
        }
    }
}   /* write_layer_lineshape_config() */


/***********************************************************
* void write_model_spectra(
*   FILE      *stream,
*   model_t   *model,
*   simplex_t *simplex)
*
* Purpose:
*   Writes the designated output spectra to a stream.
************************************************************/

void write_model_spectra(
        FILE      *stream,
        model_t   *model,
        simplex_t *simplex)
{
    double *f = output[OUTPUT_FREQUENCY].spectrum; /* model or IF grid */
    gridsize_t nout = (model->ifmode ? model->nif : model->ngrid);
    gridsize_t imin, imax;

    /*
     * Compute the range of output grid indices corresponding to
     * the output frequency range.  (If no output range has been
     * set, model->fout_min and model->fout_max are negative.)
     */
    if ((model->fout_min >= 0.0 && model->fout_min > f[nout - 1]) ||
            (model->fout_max >= 0.0 && model->fout_max < f[0])) {
        errlog(87, 0); /* no overlap */
        return;
    }
    if (model->fout_min < 0.0 || model->fout_min <= f[0]) {
        imin = 0;
    } else {
        imin = (gridsize_t)ceil(
                ((model->fout_min - f[0]) * (1.0 - DBL_EPSILON)) / model->df);
    }
    if (model->fout_max < 0.0 || model->fout_max >= f[nout - 1]) {
        imax = nout - 1;
    } else {
        imax = (gridsize_t)floor(
                ((model->fout_max - f[0]) * (1.0 + DBL_EPSILON)) / model->df);
    }
    switch (output[ALL_OUTPUTS].format) {
    case OUTPUT_FORMAT_TEXT:
        write_model_spectra_as_text(
                stream, model, simplex, imin, imax,
                TXT_FIELD_SEPARATOR, TXT_RECORD_SEPARATOR);
        break;
    case OUTPUT_FORMAT_CSV:
        write_model_spectra_as_text(
                stream, model, simplex, imin, imax,
                CSV_FIELD_SEPARATOR, CSV_RECORD_SEPARATOR);
        break;
    case OUTPUT_FORMAT_NPY:
        write_model_spectra_as_npy(stream, simplex, imin, imax);
        break;
    default:
        errlog(204, output[ALL_OUTPUTS].format);
        break;
    }
    fflush(stream);
    return;
}   /* write_model_spectra() */


/***********************************************************
* static void write_model_spectra_as_npy(
*   FILE       *stream,
*   simplex_t  *simplex
*   gridsize_t imin,
*   gridsize_t imax)
*
* Purpose:
*   Writes the designated output spectra to a stream in npy
*   binary format, version 1.0.  The specification of the
*   format is available at 
*
*   https://numpy.org/devdocs/reference/generated/numpy.lib.format.html#
************************************************************/

static void write_model_spectra_as_npy(
        FILE       *stream,
        simplex_t  *simplex,
        gridsize_t imin,
        gridsize_t imax)
{
    int header_cnt, header_len, header_data_len;
    int ncols, nrows;
    gridsize_t i;
    unsigned int j;
    double x;
    char dbl_dtype_str[8];

#ifdef _WIN32
    /*
     * On Windows, set the output stream to binary mode to
     * prevent translation of '\n' to "\r\n".
     */
    if(_setmode(_fileno(stream), _O_BINARY) == -1)
        errlog(210,0);
#endif

    /*
     * Count the number of data columns
     */
    ncols = 0;
    for (j = 0; outcol[j] != 0; ++j) {
        ncols += 1;
        if (output[outcol[j]].flags & OUTPUT_JACOBIAN) {
            unsigned int k;
            if (output[ALL_OUTPUTS].flags & OUTPUT_JACOBIAN_ERRS) {
                for (k = 0; k < simplex->n; ++k)
                    ncols += 3;
            } else {
                for (k = 0; k < simplex->n; ++k)
                    ncols += 1;
            }
        }
    }
    nrows = (imax < imin ? 0 : imax - imin + 1);
    /*
     * The npy header data must be written in the same byte order
     * on any machine.  In contrast, the array data can be
     * written in any byte order, so long as it is documented in
     * the dtype string in the header.  We will write the array
     * data (which is all 8-byte IEEE-754 doubles) in the native
     * byte order of the machine, and record the following dtype
     * in the header.
     */
    strncpy(dbl_dtype_str,
            little_endian_machine() ? "<f8" : ">f8",
            sizeof(dbl_dtype_str) - 1);
    dbl_dtype_str[sizeof(dbl_dtype_str) - 1] = '\0';
    if (output[ALL_OUTPUTS].flags & OUTPUT_HEADERS) {
        /*
         * If output headers are on, write a header for a
         * structured array with nrows records.  Field names are
         * output array names with units.
         */
        const char HEADER_FMT_STR_BEGIN[] =
            "%s%c%c%c%c{'descr': [";
        const char STRUCT_ARR_FIELD_FMT_K[] =
            "('%s(%s)[%s]', '%s'), ";
        const char STRUCT_ARR_FIELD_FMT[] =
            "('%s[%s]', '%s'), ";
        const char STRUCT_ARR_FIELD_FMT_K_JACOBIAN[] =
            "('d%s(%s)/dvar[%d][%s]/[%s]', '%s'), ";
        const char STRUCT_ARR_FIELD_FMT_JACOBIAN[] =
            "('d%s/dvar[%d][%s]/[%s]', '%s'), ";
        const char STRUCT_ARR_FIELD_FMT_JACOBIAN_ERRS[] =
            "('round_err[%d][%d]', '%s'), ('trunc_err[%d][%d]', '%s'), ";
        const char HEADER_FMT_STR_END[] =
            "], 'fortran_order': False, 'shape': (%d,), }";
        header_data_len = 0;
        header_cnt = snprintf(
                NULL,
                (size_t)0, 
                HEADER_FMT_STR_BEGIN,
                NPY_MAGIC_STRING,
                (char)NPY_VERSION_MAJOR,
                (char)NPY_VERSION_MINOR,
                header_data_len & 0xff,
                (header_data_len & 0xff00) >> 8);
        for (j = 0; outcol[j] != 0; ++j) {
            if (outcol[j] == OUTPUT_K) {
                header_cnt += snprintf(
                        NULL,
                        (size_t)0,
                        STRUCT_ARR_FIELD_FMT_K,
                        output[outcol[j]].name,
                        k_type[output[outcol[j]].k_typenum].name,
                        unit_tab[output[outcol[j]].unitnum].name,
                        dbl_dtype_str);
            } else {
                header_cnt += snprintf(
                        NULL,
                        (size_t)0,
                        STRUCT_ARR_FIELD_FMT,
                        output[outcol[j]].name,
                        unit_tab[output[outcol[j]].unitnum].name,
                        dbl_dtype_str);
            }
            if (output[outcol[j]].flags & OUTPUT_JACOBIAN) {
                unsigned int k;
                for (k = 0; k < simplex->n; ++k) {
                    if (outcol[j] == OUTPUT_K) {
                        header_cnt += snprintf(
                                NULL,
                                (size_t)0,
                                STRUCT_ARR_FIELD_FMT_K_JACOBIAN,
                                output[outcol[j]].name,
                                k_type[output[outcol[j]].k_typenum].name,
                                k,
                                unit_tab[output[outcol[j]].unitnum].name,
                                unit_tab[simplex->unitnum[k]].name,
                                dbl_dtype_str);
                    } else {
                        header_cnt += snprintf(
                                NULL,
                                (size_t)0,
                                STRUCT_ARR_FIELD_FMT_JACOBIAN,
                                output[outcol[j]].name,
                                k,
                                unit_tab[output[outcol[j]].unitnum].name,
                                unit_tab[simplex->unitnum[k]].name,
                                dbl_dtype_str);
                    }
                    if (output[ALL_OUTPUTS].flags & OUTPUT_JACOBIAN_ERRS)
                        header_cnt += snprintf(
                                NULL,
                                (size_t)0,
                                STRUCT_ARR_FIELD_FMT_JACOBIAN_ERRS,
                                j, k, dbl_dtype_str, j, k, dbl_dtype_str);
                }
            }
        }
        header_cnt += snprintf(NULL, (size_t)0, HEADER_FMT_STR_END, nrows); 
        header_len = header_cnt +
            NPY_HEADER_ALIGNMENT - (header_cnt % NPY_HEADER_ALIGNMENT);
        header_data_len = header_len - (int)strlen(NPY_MAGIC_STRING) - 4;
        fprintf(stream,
                HEADER_FMT_STR_BEGIN,
                NPY_MAGIC_STRING,
                (char)NPY_VERSION_MAJOR,
                (char)NPY_VERSION_MINOR,
                header_data_len & 0xff,
                (header_data_len & 0xff00) >> 8);
        for (j = 0; outcol[j] != 0; ++j) {
            if (outcol[j] == OUTPUT_K) {
                fprintf(stream,
                        STRUCT_ARR_FIELD_FMT_K,
                        output[outcol[j]].name,
                        k_type[output[outcol[j]].k_typenum].name,
                        unit_tab[output[outcol[j]].unitnum].name,
                        dbl_dtype_str);
            } else {
                fprintf(stream,
                        STRUCT_ARR_FIELD_FMT,
                        output[outcol[j]].name,
                        unit_tab[output[outcol[j]].unitnum].name,
                        dbl_dtype_str);
            }
            if (output[outcol[j]].flags & OUTPUT_JACOBIAN) {
                unsigned int k;
                for (k = 0; k < simplex->n; ++k) {
                    if (outcol[j] == OUTPUT_K) {
                        fprintf(stream,
                                STRUCT_ARR_FIELD_FMT_K_JACOBIAN,
                                output[outcol[j]].name,
                                k_type[output[outcol[j]].k_typenum].name,
                                k,
                                unit_tab[output[outcol[j]].unitnum].name,
                                unit_tab[simplex->unitnum[k]].name,
                                dbl_dtype_str);
                    } else {
                        fprintf(stream,
                                STRUCT_ARR_FIELD_FMT_JACOBIAN,
                                output[outcol[j]].name,
                                k,
                                unit_tab[output[outcol[j]].unitnum].name,
                                unit_tab[simplex->unitnum[k]].name,
                                dbl_dtype_str);
                    }
                    if (output[ALL_OUTPUTS].flags & OUTPUT_JACOBIAN_ERRS)
                        fprintf(stream,
                                STRUCT_ARR_FIELD_FMT_JACOBIAN_ERRS,
                                j, k, dbl_dtype_str, j, k, dbl_dtype_str);
                }
            }
        }
        fprintf(stream, HEADER_FMT_STR_END, nrows); 
    } else {
        /*
         * If output headers are off, write a header for a simple
         * array with dimensions nrows x ncols.
         */
        const char HEADER_FMT_STR[] = 
                "%s%c%c%c%c{'descr': '%s', 'fortran_order': False, "
                "'shape': (%d, %d), }";
        header_data_len = 0;
        header_cnt = snprintf(
                NULL,
                (size_t)0, 
                HEADER_FMT_STR,
                NPY_MAGIC_STRING,
                (char)NPY_VERSION_MAJOR,
                (char)NPY_VERSION_MINOR,
                header_data_len & 0xff,
                (header_data_len & 0xff00) >> 8,
                dbl_dtype_str,
                nrows,
                ncols);
        header_len = header_cnt +
            NPY_HEADER_ALIGNMENT - (header_cnt % NPY_HEADER_ALIGNMENT);
        header_data_len = header_len - (int)strlen(NPY_MAGIC_STRING) - 4;
        fprintf(stream, 
                HEADER_FMT_STR,
                NPY_MAGIC_STRING,
                (char)NPY_VERSION_MAJOR,
                (char)NPY_VERSION_MINOR,
                header_data_len & 0xff,
                (header_data_len & 0xff00) >> 8,
                dbl_dtype_str,
                nrows,
                ncols);
    }
    /*
     * Write the header alignment padding and terminating newline.
     */
    while (++header_cnt < header_len)
        fputc(' ', stream);
    fputc('\n', stream);
    /*
     * Write the spectral data.  The layout is the same for
     * simple and structured arrays.
     */
    for (i = imin; i <= imax; ++i) {
        for (j = 0; outcol[j] != 0; ++j) {
            x = output[outcol[j]].spectrum[i]
                / unit_tab[output[outcol[j]].unitnum].factor;
            x -= unit_tab[output[outcol[j]].unitnum].offset;
            fwrite(&x, sizeof(double), (size_t)1, stream);
            if (output[outcol[j]].flags & OUTPUT_JACOBIAN) {
                unsigned int k;
                if (output[ALL_OUTPUTS].flags & OUTPUT_JACOBIAN_ERRS) {
                    for (k = 0; k < simplex->n; ++k) {
                        x = output[outcol[j]].jacobian[k][i] *
                            unit_tab[simplex->unitnum[k]].factor /
                            unit_tab[output[outcol[j]].unitnum].factor;
                        fwrite(&x, sizeof(double), (size_t)1, stream);
                        x = output[outcol[j]].jacobian_round_err[k][i] *
                            unit_tab[simplex->unitnum[k]].factor /
                            unit_tab[output[outcol[j]].unitnum].factor;
                        fwrite(&x, sizeof(double), (size_t)1, stream);
                        x = output[outcol[j]].jacobian_trunc_err[k][i] *
                            unit_tab[simplex->unitnum[k]].factor /
                            unit_tab[output[outcol[j]].unitnum].factor;
                        fwrite(&x, sizeof(double), (size_t)1, stream);
                    }
                } else {
                    for (k = 0; k < simplex->n; ++k) {
                        x = output[outcol[j]].jacobian[k][i] *
                            unit_tab[simplex->unitnum[k]].factor /
                            unit_tab[output[outcol[j]].unitnum].factor;
                        fwrite(&x, sizeof(double), (size_t)1, stream);
                    }
                }
            }
        }
    }
    return;
}   /* write_model_spectra_as_npy() */


/***********************************************************
* static void write_model_spectra_as_text(
*   FILE       *stream,
*   model_t    *model,
*   simplex_t  *simplex
*   gridsize_t imin,
*   gridsize_t imax,
*   const char *fs,
*   const char *rs)
*
* Purpose:
*   Writes the designated output spectra to a stream in text
*   format, with field separator string *fs and record
*   separator string *rs.
************************************************************/

static void write_model_spectra_as_text(
        FILE       *stream,
        model_t    *model,
        simplex_t  *simplex,
        gridsize_t imin,
        gridsize_t imax,
        const char *fs,
        const char *rs)
{
    double *f = output[OUTPUT_FREQUENCY].spectrum;
    int freq_dig;
    double x;
    gridsize_t i;

    /*
     * Compute the precision which will be used for printing
     * frequencies.  The goal is to represent df to at least 1%
     * accuracy when printing the highest value of f, subject to
     * a minimum precision of FLT_DIG, and a maximum precision of
     * DBL_DIG.
     */
    freq_dig = 2 + (int)ceil(log10(f[imax] / model->df));
    freq_dig = freq_dig > FLT_DIG ? freq_dig : FLT_DIG;
    freq_dig = freq_dig > DBL_DIG ? DBL_DIG : freq_dig;
    if (output[ALL_OUTPUTS].flags & OUTPUT_HEADERS)
        write_output_headers_as_text(stream, simplex, fs, rs);
    /*
     * Write the spectra.
     */
    for (i = imin; i <= imax; ++i) {
        unsigned int j;
        for (j = 0; outcol[j] != 0; ++j) {
            x = output[outcol[j]].spectrum[i]
                / unit_tab[output[outcol[j]].unitnum].factor;
            x -= unit_tab[output[outcol[j]].unitnum].offset;
            fprintf(stream,
                    "%s%.*e",
                    (j == 0) ? "" : fs,
                    outcol[j] == OUTPUT_FREQUENCY ? freq_dig: FLT_DIG,
                    x);
            if (output[outcol[j]].flags & OUTPUT_JACOBIAN) {
                unsigned int k;
                if (output[ALL_OUTPUTS].flags & OUTPUT_JACOBIAN_ERRS) {
                    for (k = 0; k < simplex->n; ++k) {
                        fprintf(stream,
                                "%s%e%s%e%s%e",
                                fs,
                                output[outcol[j]].jacobian[k][i] *
                                unit_tab[simplex->unitnum[k]].factor / 
                                unit_tab[output[outcol[j]].unitnum].factor,
                                fs,
                                output[outcol[j]].jacobian_round_err[k][i] *
                                unit_tab[simplex->unitnum[k]].factor / 
                                unit_tab[output[outcol[j]].unitnum].factor,
                                fs,
                                output[outcol[j]].jacobian_trunc_err[k][i] *
                                unit_tab[simplex->unitnum[k]].factor / 
                                unit_tab[output[outcol[j]].unitnum].factor);
                    }
                } else {
                    for (k = 0; k < simplex->n; ++k) {
                        fprintf(stream,
                                "%s%e",
                                fs,
                                output[outcol[j]].jacobian[k][i] *
                                unit_tab[simplex->unitnum[k]].factor / 
                                unit_tab[output[outcol[j]].unitnum].factor);
                    }
                }
            }
        }
        fprintf(stream, "%s", rs);
    }
    return;
}   /* write_model_spectra_as_text() */


/***********************************************************
* static void write_Nscale_config(
*   FILE *stream,
*   simplex_t *simplex)
*
* Purpose:
*   Walks through the Nscale list, writing out corresponding
*   configuration statements.
************************************************************/

static void write_Nscale_config(
        FILE *stream,
        simplex_t *simplex)
{
    Nscale_list_t *ptr;
    unsigned int j;

    for (ptr = Nscale_list_head(); ptr != NULL; ptr = ptr->next) {
        if (ptr->tagnum != 0)
            fprintf(stream, "Nscale %s", tag_string(ptr->tagnum));
        else
            fprintf(stream, "Nscale");
        if (isvar(simplex, &(ptr->Nscale))) {
            j = get_simplex_variable_index(simplex, &(ptr->Nscale));
            fprintf(stream,
                    " %s %#g  %.13g",
                    col_type[ptr->col_typenum].name,
                    ptr->Nscale,
                    simplex->scale[j]);
            if (output[ALL_OUTPUTS].flags & OUTPUT_FITTED) {
                fprintf( stream,
                        "  (%.3g)",
                        simplex_variable_range(simplex, j));
            }
            fprintf(stream, "\n");
        } else {
            fprintf(stream,
                    " %s %.13g\n",
                    col_type[ptr->col_typenum].name,
                    ptr->Nscale);
        }
        if (ptr->next == NULL)
            fprintf(stream, "\n");
    }
    return;
}   /* write_Nscale_config() */


/***********************************************************
* static void write_output_headers_as_text(
*   FILE *stream,
*   simplex_t *simplex,
*   const char *fs,
*   const char *rs)
*
* Purpose:
*   Writes column headers for all designated output spectra
*   to a stream in text format, with field separator string
*   *fs and record separator string *rs.
************************************************************/

static void write_output_headers_as_text(
        FILE *stream,
        simplex_t *simplex,
        const char *fs,
        const char *rs)
{
    int j;
    for (j = 0; outcol[j] != 0; ++j) {
        if (outcol[j] == OUTPUT_K) {
            fprintf(stream,
                    "%s%s(%s)[%s]",
                    (j == 0) ? "" : fs,
                    output[outcol[j]].name,
                    k_type[output[outcol[j]].k_typenum].name,
                    unit_tab[output[outcol[j]].unitnum].name);
        } else {
            fprintf(stream,
                    "%s%s[%s]",
                    (j == 0) ? "" : fs,
                    output[outcol[j]].name,
                    unit_tab[output[outcol[j]].unitnum].name);
        }
        if (output[outcol[j]].flags & OUTPUT_JACOBIAN) {
            unsigned int i;
            for (i = 0; i < simplex->n; ++i) {
                if (outcol[j] == OUTPUT_K) {
                    fprintf(stream,
                            "%sd%s(%s)/dvar[%d][%s]/[%s]",
                            fs,
                            output[outcol[j]].name,
                            k_type[output[outcol[j]].k_typenum].name,
                            i,
                            unit_tab[output[outcol[j]].unitnum].name,
                            unit_tab[simplex->unitnum[i]].name
                           );
                } else {
                    fprintf(stream,
                            "%sd%s/dvar[%d][%s]/[%s]",
                            fs,
                            output[outcol[j]].name,
                            i,
                            unit_tab[output[outcol[j]].unitnum].name,
                            unit_tab[simplex->unitnum[i]].name
                           );
                }
                if (output[ALL_OUTPUTS].flags & OUTPUT_JACOBIAN_ERRS)
                    fprintf(stream, "%sround_err%strunc_err", fs, fs);
            }
        }
    }
    fprintf(stream, "%s", rs);
    return;
}   /* write_output_headers_as_text() */


/***********************************************************
* void write_PTmode(FILE *stream, int PTmode)
*
* Purpose:
*   Writes a newline-terminated PTmode statement to a
*   stream, specifying the PTmode.  In addition, if the
*   user has specified the base extrapolation mode, a
*   corresponding second PTmode statement is written.
*
* Arguments:
*   FILE *stream - destination for output
*   int PTmode - PTmode, encoded as defined in model.h
************************************************************/

void write_PTmode(FILE *stream, int PTmode)
{
    fprintf(stream, "PTmode ");
    if ((PTmode & PTMODE_P) && (PTmode & PTMODE_T)) {
        fprintf(stream, "midpoint\n");
    } else {
        if (PTmode & PTMODE_DP)
            fprintf(stream, "dP ");
        else if (PTmode & PTMODE_P)
            fprintf(stream, "P ");
        else if (PTmode & PTMODE_PBASE)
            fprintf(stream, "Pbase ");
        else
            fprintf(stream, "? ");
        if (PTmode & PTMODE_T)
            fprintf(stream, "T\n");
        else if (PTmode & PTMODE_TBASE)
            fprintf(stream, "Tbase\n");
        else
            fprintf(stream, "?\n");
    }
    if (PTmode & PTMODE_EXTEND_USER_DEFINED) {
        fprintf(stream, "PTmode extend ");
        if (PTmode & PTMODE_EXTEND_ISOTHERMAL)
            fprintf(stream, "isothermal\n");
        else if (PTmode & PTMODE_EXTEND_DRY_ADIABATIC)
            fprintf(stream, "dry_adiabatic\n");
        else if (PTmode & PTMODE_EXTEND_MOIST_ADIABATIC)
            fprintf(stream, "moist_adiabatic\n");
        else if (PTmode & PTMODE_EXTEND_ENVIRONMENTAL)
            fprintf(stream, "environmental\n");
        else
            fprintf(stream, "?\n");
    }
    fflush(stream);
    return;
}   /* write_PTmode() */


/***********************************************************
* static void write_runtime_details(FILE *stream, model_t *model)
*
* Purpose:
*   Writes logged model performance timings.
************************************************************/

static void write_runtime_details(FILE *stream, model_t *model)
{
    fprintf(stream, "#   atmospheric model  %9.0f us\n",
            1.0e6 * model->am_runtime);
    fprintf(stream, "#   optical depths     %9.0f us\n",
            1.0e6 * model->od_runtime);
    fprintf(stream, "#   radiative transfer %9.0f us\n",
            1.0e6 * model->rt_runtime);
    fprintf(stream, "#   derived spectra    %9.0f us\n",
            1.0e6 * model->spec_runtime);
    return;
}   /* write_runtime_details() */


/***********************************************************
* void write_variable_info(
*       char *var_name,
*       double *var,
*       int var_unitnum)
*
* Purpose:
*   Prints a line in a standard format for variables,
*   including any units.  If the variable is a
*   differentiation variable, the characteristic scale is
*   printed next to the value.  If the variable is a fit
*   variable, the initial simplex scale and range over the
*   converged simplex are printed.  User-supplied numbers
*   are written with %.13g format to recover all digits
*   supplied by the user.  Computed numbers are written
*   with an appropriate truncation.
*
* Arguments:
*   FILE *stream       - destination for output
*   char *var_name     - variable name string
*   double var         - variable value (possibly unmapped)
*   double var_ptr     - pointer to variable
*   int var_unitnum    - index number for units
*                        UNIT_NONE if the variable is unitless
*   simplex_t *simplex - pointer to simplex structure
************************************************************/

void write_variable_info(
        FILE *stream,
        const char *var_name,
        double var,
        double *var_ptr,
        int var_unitnum,
        simplex_t *simplex)
{
    fprintf(stream, "%s", var_name);
    if (var_unitnum == UNIT_NONE) {
        if (isvar(simplex, var_ptr)) {
            int j = get_simplex_variable_index(simplex, var_ptr);
            fprintf(stream,
                    " %g  %.13g",
                    var,
                    simplex->scale[j]);
            if (output[ALL_OUTPUTS].flags & OUTPUT_FITTED) {
                fprintf(stream,
                        "  (%.3g)",
                        simplex_variable_range(simplex, j));
            }
        } else {
            fprintf(stream, " %.13g", var);
        }
    } else {
        if (isvar(simplex, var_ptr)) {
            int j = get_simplex_variable_index(simplex, var_ptr);
            print_with_unit(
                    stream,
                    " %g %s",
                    var,
                    var_unitnum);
            print_differential_with_unit(
                    stream,
                    "  %.13g %s",
                    simplex->scale[j],
                    var_unitnum);
            if (output[ALL_OUTPUTS].flags & OUTPUT_FITTED) {
                print_differential_with_unit(
                        stream,
                        "  (%.3g %s)",
                        simplex_variable_range(simplex, j),
                        var_unitnum);
            }
        } else {
            print_with_unit(
                    stream,
                    " %.13g %s",
                    var,
                    var_unitnum);
        }
    }
    return;
}   /* write_variable_info() */
