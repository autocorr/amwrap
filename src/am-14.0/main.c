/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* main.c                         S. Paine rev. 2024 April 17
*
* Main file for the am atmospheric model
************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef _OPENMP
    #include <omp.h>
#endif

#include "am_alloc.h"
#include "am_sysdep.h"
#include "am_types.h"
#include "config.h"
#include "dcache.h"
#include "doc.h"
#include "errlog.h"
#include "fit.h"
#include "jacobian.h"
#include "kcache.h"
#include "model.h"
#include "nscale.h"
#include "output.h"
#include "planck.h"
#include "simplex.h"
#include "tags.h"
#include "transform.h"
#include "units.h"
#include "version.h"

static void benchmarks(void);
static int  report_environment(FILE*);

int main(int argc, char *argv[])
{
    model_t        model = MODEL_INIT;
    model_t       lmodel = MODEL_INIT;
    fit_data_t  fit_data = FIT_DATA_INIT;
    simplex_t    simplex = SIMPLEX_INIT;
    double        tstart = 0.0;
    int        benchmark = 0;
    int          errstat = 0;
    char       *progname = argv[0];

    if (argc <= 1) {
        version(stderr);
        usage(stderr, progname);
        return 0;
    }
    while (argc > 1) {
        if (
                !strcmp(argv[1], "-a") ||
                !strcmp(argv[1], "--atmosphere")
                ) {
            output[ALL_OUTPUTS].flags |= OUTPUT_AM_ONLY;
            if (--argc == 1) {
                errlog(170, 0);
                return print_errlog();
            }
            ++argv;
        } else if (
                !strcmp(argv[1], "-b") ||
                !strcmp(argv[1], "--benchmark")
                ) {
            if (--argc == 1) {
                benchmarks();
                return 0;
            }
            benchmark = 1;
            ++argv;
        } else if (
                !strcmp(argv[1], "-e") ||
                !strcmp(argv[1], "--environment")
                ) {
            return report_environment(stderr);
        } else if (
                !strcmp(argv[1], "-h") ||
                !strcmp(argv[1], "-?") ||
                !strcmp(argv[1], "--help")
                ) {
            usage(stdout, progname);
            return 0;
        } else if (
                !strcmp(argv[1], "-l") ||
                !strcmp(argv[1], "--legal")
                ) {
            legal();
            return 0;
        } else if (
                !strcmp(argv[1], "-r") ||
                !strcmp(argv[1], "--references")
                ) {
            references();
            return 0;
        } else if (
                !strcmp(argv[1], "-v") ||
                !strcmp(argv[1], "--version")
                ) {
            version(stdout);
            return 0;
        } else {
            /*
             * All command-line switches have been consumed.
             * argv[0] now points to the first command-line
             * parameter for substitution into the config file,
             * and argc is the number of parameters.
             */
            break;
        }
    }
    errstat = parse_config_file(argc, argv, &model, &fit_data, &simplex);
    /*
     * For Jacobians and fits, the structure lmodel ("last
     * model") holds copies of the scalar variables in the
     * structure model from the last round of spectral
     * computations performed.  When these scalar variables are
     * changed for a subsequent round of spectral computations
     * (as when computing Jacobians with respect to a particular
     * scalar variable), these copies support dependency checking
     * to avoid redundant spectral computations.
     */
    if (!errstat &&
            output[ALL_OUTPUTS].flags & (OUTPUT_FITTED | OUTPUT_JACOBIAN)) {
        /*
         * Here, the layer and column dimensions of model are
         * replicated in lmodel.  But first, to establish the
         * initial shape of model, any interpolated levels need
         * to be added by calling setup_atmospheric_model().
         * These interpolated levels include user-defined source
         * or observing levels, and auto-generated or
         * user-defined path tangent levels.  Later changes in
         * shape associated with level reorderings (which might
         * occur, for example, when differentiating with respect
         * to a level position) are synchronized between model
         * and lmodel.
         */
        errstat = setup_atmospheric_model(&model, NULL);
        if (!errstat)
            copy_model_dimensions(&model, &lmodel);
    }
    /*
     * Setup is complete.  Initialize the random number
     * generator, start the performance clock, and run the
     * requested computation.
     */
    srand((unsigned int)time(NULL));
    tstart = am_timer(0.0);
    if (errstat) {
        /*
         * Errors occurred above; do nothing.
         */
    } else if (output[ALL_OUTPUTS].flags & OUTPUT_AM_ONLY) {
        /*
         * Atmospheric model only, no radiative transfer
         */
        if (!setup_atmospheric_model(&model, NULL)) {
            model.am_runtime = model.runtime = am_timer(tstart);
            if (benchmark) {
                printf("%.0f us \n", 1.0e6 * model.runtime);
            } else {
                write_model_config_data(stderr, &model, &fit_data, &simplex);
                print_with_unit(stdout, "%-6.13g", model.za, UNIT_DEGREE);
                printf(" %#8g", total_airmass(&model));
                print_with_unit(stdout, " %8.2f",
                        total_refraction(&model), UNIT_ARCSEC);
                if (
                        model.PTmode & PTMODE_HYDROSTATIC &&
                        model.path_begin > 0 &&
                        model.path_end   > 0
                        ) {
                    print_with_unit(stdout, " %#8g",
                            source_to_obs_path_distance(&model),
                            model.R0_unitnum);
                    print_with_unit(stdout, " %#8g",
                            source_to_obs_geometric_distance(&model),
                            model.R0_unitnum);
                    print_with_unit(stdout, " %#8g\n",
                            source_to_obs_projected_distance(&model),
                            model.R0_unitnum);
                } else {
                    printf("\n");
                }
            }
        }
    } else if (output[ALL_OUTPUTS].flags & OUTPUT_FITTED) {
        /*
         * Model fit to spectral data
         */
        fit(&model, &lmodel, &fit_data, &simplex);
        if (benchmark)
            printf("%.0f us \n", 1.0e6 * am_timer(tstart));
    } else if (output[ALL_OUTPUTS].flags & OUTPUT_JACOBIAN) {
        /*
         * Jacobians.  The Jacobian spectra are computed first,
         * after which the model spectra are recomputed at the
         * reference state.
         */
        if (    !compute_jacobians(&model, &lmodel, &simplex) &&
                !compute_model(&model, &lmodel)) {
            model.runtime = am_timer(tstart);
            if (benchmark) {
                printf("%.0f us \n", 1.0e6 * model.runtime);
            } else {
                write_model_config_data(stderr, &model, &fit_data, &simplex);
                write_model_spectra(stdout, &model, &simplex);
            }
        }
    } else {
        /*
         * Single-pass radiative transfer
         */
        if (!compute_model(&model, NULL)) {
            model.runtime = am_timer(tstart);
            if (benchmark) {
                printf("%.0f us \n", 1.0e6 * model.runtime);
            } else {
                write_model_config_data(stderr, &model, &fit_data, &simplex);
                write_model_spectra(stdout, &model, &simplex);
            }
        }
    }
    free_fit_data_entities(&fit_data);
    free_model_entities(&model);
    free_model_entities(&lmodel);
    free_jacobians(&simplex);
    free_simplex_entities(&simplex);
    kcache_free_all();
    free_Nscale_list();
    free_tag_string_table();
    return print_errlog();
}   /* main() */


/***********************************************************
* static void benchmarks(void)
*
* Purpose:
*   Runs internal performance benchmarks on time-critical
*   functions.  This facility is useful for algorithm
*   development and for measuring the performance effect of
*   OpenMP directives and compiler options.
************************************************************/

static void benchmarks(void)
{
    printf(
            "Running internal benchmarks.\n\n"
            "Environment\n"
            "-----------\n"
            );
    report_environment(stdout);
    printf(
            "Timings\n"
            "-------\n"
            );
    planck_benchmarks();
    printf("\n");
    ft_benchmarks();
    printf("\n");
    return;
}   /* benchmarks() */


/***********************************************************
* static int report_environment(FILE *stream)
*
* Purpose:
*   Prints information about the runtime environment to a
*   stream.
*
* Arguments:
*   FILE *stream - destination for output
*
* Return:
*   0 if OK
*   1 if any errors were logged.
************************************************************/

static int report_environment(FILE *stream)
{
#ifdef _OPENMP
#if   (_OPENMP < 200805)
    /*
     * In OpenMP prior to version 3.0 (200805), nested
     * parallelism was enabled or disabled via the OMP_NESTED
     * environment variable, with no control on nesting depth.
     * This was a dangerous design.
     */
    fprintf(stream,
            "OpenMP\n"
            "  Available processors = %d\n"
            "  OMP_NUM_THREADS = %d\n"
            "  OMP_NESTED = %s\n",
            omp_get_num_procs(),
            omp_get_max_threads(),
            omp_get_nested() ? "true" : "false"
            );
#elif (_OPENMP <= 201511)
    /*
     * OpenMP version 3.0 (200805) introduced the environment
     * variable OMP_MAX_ACTIVE_LEVELS, which sets a limit on the
     * maximum depth of nested parallel regions.  This was an
     * improvement, but made OMP_NESTED redundant, since setting
     * OMP_MAX_ACTIVE_LEVELS = 1 is sufficient to disable nested
     * parallelism.  Moreover, in some implementations setting
     * OMP_MAX_ACTIVE_LEVELS > 1 will override setting OMP_NESTED
     * to false.  With that in mind, here we report the actual
     * run-time state based on OpenMP's internal control
     * variables.
     *
     * This was the state of the OpenMP specification through
     * OpenMP version 4.5 (201511).
     */
    fprintf(stream,
            "OpenMP\n"
            "  Available processors = %d\n"
            "  OMP_NUM_THREADS = %d\n"
            "  OMP_NESTED = %s\n"
            "  OMP_MAX_ACTIVE_LEVELS = %d\n",
            omp_get_num_procs(),
            omp_get_max_threads(),
            omp_get_nested() ? "true" : "false",
            omp_get_max_active_levels()
            );
#else
    /*
     * In OpenMP 5.0 (201811), the redundant OMP_NESTED
     * environment variable and related API calls were
     * deprecated.  Partial implementations (e.g. icc) that
     * define nonstandard version macros between 201511 and
     * 201811 are also included here.
     *
     * If the user has set the environment variable OMP_NESTED to
     * any value, a deprecation message is printed.  Here we
     * report the actual value read from the environment, rather
     * than the effective state from omp_get_nested(), since the
     * deprecated omp_get_nested() call will be eliminated
     * entirely in OpenMP 6.0 per Technical Report 12 (Nov 2023).
     */
    if (getenv("OMP_NESTED") == NULL) {
        fprintf(stream,
                "OpenMP\n"
                "  Available processors = %d\n"
                "  OMP_NUM_THREADS = %d\n"
                "  OMP_MAX_ACTIVE_LEVELS = %d\n",
                omp_get_num_procs(),
                omp_get_max_threads(),
                omp_get_max_active_levels()
                );
    } else {
        fprintf(stream,
                "OpenMP\n"
                "  Available processors = %d\n"
                "  OMP_NUM_THREADS = %d\n"
                "  OMP_NESTED = %s (deprecated, use OMP_MAX_ACTIVE_LEVELS to\n"
                "    control nested parallelism)\n"
                "  OMP_MAX_ACTIVE_LEVELS = %d\n",
                omp_get_num_procs(),
                omp_get_max_threads(),
                getenv("OMP_NESTED"),
                omp_get_max_active_levels()
                );
    }
#endif
#endif
    report_dcache_env_info(stream);
    report_kcache_env_info(stream);
    report_fit_env_info(stream);
    fprintf(stream, "\n");
    return print_errlog();
}   /* report_environment() */
