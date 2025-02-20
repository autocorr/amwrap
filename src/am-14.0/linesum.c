/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* linesum.c                       S. Paine rev. 2024 July 25
*
* Computation of line-by-line absorption coefficient.
************************************************************/

#include <float.h>
#include <math.h>
#include <stdlib.h>
#ifdef _OPENMP
    #include <omp.h>
#endif

#include "am_sysdep.h"
#include "am_types.h"
#include "errlog.h"
#include "linesum.h"
#include "math_const.h"
#include "phys_const.h"
#include "specfunc.h"

/*
 * LINESUM_MIN_THD_BLOCKSIZE is the smallest number of frequency grid
 * points that will be given to a thread to do a block of a
 * line-by-line computation.  The rationale for setting a lower limit
 * is that blocks any smaller than this will not efficiently amortize
 * the fixed overhead of creating a thread.
 */
#ifndef LINESUM_MIN_THD_BLOCKSIZE
    #define LINESUM_MIN_THD_BLOCKSIZE 8
#endif

/*
 * Each block of the line-by-line computation, except for the last, is
 * a multiple of LINESUM_BLOCK_ALIGN points long.  This can be used to
 * control the alignment of the subarrays within each block relative
 * to the alignment of the full arrays.  For example, for arrays of
 * 8-byte doubles initially allocated with 16-byte alignment, this
 * alignment will be preserved on subarrays if LINESUM_BLOCK_ALIGN is
 * set to 2.
 */
#ifndef LINESUM_BLOCK_ALIGN
    #define LINESUM_BLOCK_ALIGN 2
#endif

/*
 * Additions to the lineshape type table between program version
 * number changes should be made at the end, since the enum value in
 * linesum.h associated with each lineshape type is part of the cache
 * header data. 
 */
struct lineshape_type_tabentry lineshape_type[] = {
    {"none",                0.0},
    {"Doppler",             0.0},
    {"Full_Lorentz",        0.0},
    {"Gross",               0.0},
    {"Lorentz",             0.0},
    {"VVW",                 0.0},
    {"VVW_coupled",         0.0},
    {"VVH",                 0.0},
    {"VVH_750",           750.0},
    {"Voigt-Kielkopf",      0.0},
    {"END",                 0.0},
};

/*
 * An array of the following structures contains precomputed line data
 * for the current P, T, and vmr.
 */
typedef struct line_data_t line_data_t;

struct line_data_t {
    double f0;          /* frequency, adjusted for pressure shift   */
    double S;           /* line strength at current temperature     */
    double alpha;       /* Doppler width                            */
    double gamma;       /* collisional width                        */
    double Voigt_hwhm_reciprocal; /* 1 / (Voigt HWHM)               */
    double Voigt_eta;   /* Lorentzian coeff. in Voigt approximation */
    double Y;           /* line mixing parameter                    */
    double G;           /* line strength mixing parameter           */
};

static void bracket_in_band_lines(
        const double*,
        const gridsize_t,
        const line_data_t*,
        const int,
        int*,
        int*);

static int linesum_block(
        double*,
        double*,
        double*,
        double*,
        const double*,
        const double*,
        const double,
        const gridsize_t,
        const line_data_t*,
        const int,
        const int,
        const double,
        const double);

static int Qtab_interp(
        const double,
        const double*,
        const int,
        const int,
        double*);


/***********************************************************
* static void bracket_in_band_lines(
*         const double *f,
*         const gridsize_t ngrid,
*         const line_data_t* line_data,
*         const int nlines,
*         int *p,
*         int *q)
*
* Purpose:
*   Finds line catalog indices p, q such that that p and q
*   index the first and last lines, respectively, in the
*   closed interval [f[0], f[ngrid-1]].  If there are no
*   lines within this interval, this function exits with
*   q = p - 1.
*
* Arguments:
*   const double *f              - frequency grid [GHz]
*   const gridsize_t ngrid       - number of frequency grid
*                                  points
*   const line_data_t *line_data - precomputed line data
*                                  table
*   const int nlines             - number of catalog lines
*   int *p                       - index of lowest line in
*                                  frequency grid range
*   int *q                       - index of highest line in
*                                  frequency grid range
************************************************************/

static void bracket_in_band_lines(
        const double *f,
        const gridsize_t ngrid,
        const line_data_t *line_data,
        const int nlines,
        int *p,
        int *q)
{
    if (nlines <= 0) { /* empty line catalog */
        *p =  0;
        *q = -1;
    } else if (f[0] > line_data[nlines-1].f0) {
        /* All lines below f[0]. */
        *p = nlines;
        *q = nlines - 1;
    } else if (f[ngrid-1] < line_data[0].f0) {
        /* All lines above f[ngrid-1]. */
        *p =  0;
        *q = -1;
    } else {
        if (f[0] < line_data[0].f0) {
            /* no lines below f[0] */
            *p = 0;
        } else {
            /* f[0] bracketed by catalog.  Search for p. */
            int i, mid;
            i = 0;
            *p = nlines - 1;
            while (*p - i > 1) {
                mid = (*p + i) / 2;
                if (f[0] > line_data[mid].f0)
                    i = mid;
                else
                    *p = mid;
            }
        }
        if (f[ngrid-1] > line_data[nlines-1].f0) {
            /* No lines above f[ngrid-1] */
            *q = nlines - 1;
        } else {
            /* f[ngrid-1] bracketed by catalog.  Search for q. */
            int i, mid;
            *q = *p - 1;
            i = nlines - 1;
            while (i - *q > 1) {
                mid = (*q + i) / 2;
                if (f[ngrid-1] > line_data[mid].f0)
                    *q = mid;
                else
                    i = mid;
            }
        }
    }
    return;
}   /* bracket_in_band_lines() */


/***********************************************************
* int linesum(
*         double *k,
*         const double *f,
*         const double *f2,
*         const double df,
*         const gridsize_t ngrid,
*         const cat_entry_t *cat,
*         const line_coupling_table_entry_t *line_coupling_table,
*         const int nlines,
*         const int lineshape,
*         const int iso,
*         const double P,
*         const double vmr,
*         const double T,
*         const double Tref,
*         const double *Qrat,
*         const double *abundance_tab,
*         const double *mass_tab,
*         const double tol,
*         int *unresolved_line_count)
*
* Purpose:
*   Performs a line-by-line sum to compute a spectral
*   absorption coefficient k[] over the frequency grid f[].
*   Line parameters are taken from the line catalog cat[],
*   and line widths, strengths, and shifts are computed from
*   the catalog parameters based on pressure P, volume
*   mixing ratio vmr, temperature T, and partition function
*   ratio Qrat[].  The latter is an array, indexed by HITRAN
*   isotopologue number, giving for each isotopologue the
*   ratio Q(Tref) / Q(T), where Tref is the catalog
*   reference temperature.
*
*   If the parameter tol is non-zero, the computation is
*   accelerated by dropping lines from the summation, while
*   attempting to maintain a maximum fractional error
*   dk/k < tol at each grid point, as decribed in
*   linesum_block(), below.
*
*   The pointer unresolved_line_count references an integer
*   which accumulates a count of all in-band lines for which
*   the linewidth is less than the frequency grid spacing
*   df.
*
* Arguments:
*   double *k              - absorption coefficient [cm^2]
*   const double *f        - frequency grid [GHz]
*   const double *f2       - squared frequency grid [GHz^2]
*   const double df        - frequency grid interval [GHz]
*   const gridsize_t ngrid - number of frequency grid points
*   const cat_entry_t *cat - line catalog
*   const line_coupling_table_entry_t *line_coupling_table
*                          - line coupling coefficients
*   const int nlines       - number of catalog lines
*   const int lineshape    - index number of lineshape
*   const int iso          - HITRAN isotopologue index
*                            (0 means all isotopologues)
*   const double P         - air pressure [mbar]
*   const double vmr       - volume mixing ratio for self-
*                            broadening
*   const double T         - temperature [K]
*   const double Tref      - catalog reference temperature
*   const double *Qrat     - table of Q(Tref) / Q(T) for all
*                            isotopologues
*   const double *abundance_tab - isotopologue abundances
*   const double *mass_tab      - isotopologue masses
*   const double tol            - error tolerance for fast
*                                 line sum
*   int *unresolved_line_count  - cumulative count of
*                                 unresolved lines
*
* Return:
*   0 if OK, 1 if an error occurred
************************************************************/

int linesum(
    double *k,
    const double *f,
    const double *f2,
    const double df,
    const gridsize_t ngrid,
    const cat_entry_t *cat,
    const line_coupling_table_entry_t *line_coupling_table,
    const int nlines,
    const int lineshape,
    const int iso,
    const double P,
    const double vmr,
    const double T,
    const double Tref,
    const double *Qrat,
    const double *abundance_tab,
    const double *mass_tab,
    const double tol,
    int *unresolved_line_count)
{
    double dktol;
    line_data_t *line_data = NULL;
    double *dk  = NULL;
    double *ta1 = NULL;
    double *ta2 = NULL;
    gridsize_t blocksize, offset;
    int n_iso_lines;
    int i, p, q, retval;

    /*
     * Allocate memory for precomputed line data to be derived from
     * catalog data.
     */
    if ((line_data = (line_data_t*)malloc(nlines * sizeof(line_data_t)))
            == NULL) {
        errlog(59, 0);
        return 1;
    }
    /*
     * Compute line strengths adjusted for temperature, and line
     * frequencies adjusted for linear pressure shift.
     *
     * For the single-isotopologue case, the HITRAN abundance
     * weighting of the line strength is removed, and the line
     * strength for all other isotopologues is set to zero.
     *
     * Note that the units of Elo are [K], and Elo < 0 is a flag
     * indicating no temperature adjustment of line strength is to be
     * made.
     */
    if (iso == 0) {
        n_iso_lines = nlines;
        #ifdef _OPENMP
        #pragma omp parallel for schedule(guided, 4096)
        #endif
        for (i = 0; i < nlines; ++i) {
            int flag = cat[i].Elo < 0.0;
            double S = cat[i].S;
            S *= Qrat[cat[i].iso];
            S *= exp((cat[i].Elo * (T - Tref)) / (T * Tref));
            S *= 1. - exp(-(H_ON_KB * cat[i].freq) / T);
            S /= 1. - exp(-(H_ON_KB * cat[i].freq) / Tref);
            line_data[i].S =  !flag * S +  flag * cat[i].S;
            line_data[i].f0 = cat[i].freq + P * cat[i].delta_air;
        }
    } else {
        n_iso_lines = 0;
        #ifdef _OPENMP
        #pragma omp parallel for schedule(guided, 4096)\
            reduction(+:n_iso_lines)
        #endif
        for (i = 0; i < nlines; ++i) {
            int flag = cat[i].Elo < 0.0;
            int flag_iso = iso == cat[i].iso;
            double S = cat[i].S;
            S *= Qrat[cat[i].iso];
            S *= exp((cat[i].Elo * (T - Tref)) / (T * Tref));
            S *= 1. - exp(-(H_ON_KB * cat[i].freq) / T);
            S /= 1. - exp(-(H_ON_KB * cat[i].freq) / Tref);
            S =  !flag * S +  flag * cat[i].S;
            S *= flag_iso;
            line_data[i].S  = S / abundance_tab[iso];
            line_data[i].f0 = cat[i].freq + P * cat[i].delta_air;
            n_iso_lines += flag_iso;
        }
    }
    /*
     * Compute Doppler widths.
     */
    if (
            lineshape == LINESHAPE_DOPPLER        ||
            lineshape == LINESHAPE_VOIGT_KIELKOPF
            ) {
        if (T < DBL_MIN) {
            errlog(81, 0);
            free(line_data);
            return 1;
        }
        #ifdef _OPENMP
        #pragma omp parallel for schedule(guided, 4096)
        #endif
        for (i = 0; i < nlines; ++i) {
            line_data[i].alpha = line_data[i].f0 *
                sqrt((TWOKB_ON_AMU_CSQUARED * T) / mass_tab[cat[i].iso]);
        }
    }
    /*
     * Compute collisional widths.  Following Eq. A12 of the 1996
     * HITRAN manual, assume that the temperature dependence of
     * gam_self is the same as that for gam_air.
     */
    if (
            lineshape == LINESHAPE_FULL_LORENTZ   ||
            lineshape == LINESHAPE_GROSS          ||
            lineshape == LINESHAPE_LORENTZ        ||
            lineshape == LINESHAPE_VVW            ||
            lineshape == LINESHAPE_VVW_COUPLED    ||
            lineshape == LINESHAPE_VVH            ||
            lineshape == LINESHAPE_VVH_750        ||
            lineshape == LINESHAPE_VOIGT_KIELKOPF
            ) {
        if (P < DBL_MIN) {
            errlog(82, 0);
            free(line_data);
            return 1;
        }
        if (T < DBL_MIN) {
            errlog(125, 0);
            free(line_data);
            return 1;
        }
        #ifdef _OPENMP
        #pragma omp parallel for schedule(guided, 4096)
        #endif
        for (i = 0; i < nlines; ++i) {
            /*
             * In the O2 catalog data, "air" means N2, so as to
             * maintain the validity of this formula.
             */
            line_data[i].gamma = 
                ((1.0 - vmr) * cat[i].gam_air + vmr * cat[i].gam_self) *
                P * pow(Tref / T, cat[i].nair);
        }
    }
    /*
     * Compute Voigt parameters.
     */
    if (lineshape == LINESHAPE_VOIGT_KIELKOPF) {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(guided, 4096)
        #endif
        for (i = 0; i < nlines; ++i) {
            const double ln2 = 0.6931471805599453;
            const double eps = 0.0990;
            const double r1 = eps * ln2;
            const double r2 = 1.0 + r1;
            const double r3 = (1.0 - r1) * (1.0 - r1);
            const double alpha = line_data[i].alpha;
            const double gamma = line_data[i].gamma;
            const double a = gamma / alpha;
            const double el = 2.0 / (r2 + sqrt(r3 + 4.0 * ln2 / (a * a)));
            const double g2 = (1.0 / ln2) * (1.0 - r2 * el + r1 * el * el);
            double r4;
            if (a <= 1.5) {
                const double b0 =  0.47047;
                const double b1 =  0.61686;
                const double b2 = -0.16994;
                const double b3 =  1.32554;
                const double t  =  1.0 / (1.0 + b0 * a);
                r4 = t * (b1 + t * (b2 + t * b3));
            } else {
                r4 = 1.0 /
                    (a + 0.5 /
                        (a + 1.0 /
                            (a + 1.5 /
                                (a + 2.0 /
                                    (a + 2.5 / 
                                        (a + 3.0 / a))))));
            }
            line_data[i].S *= a * r4 / (PI * gamma);
            line_data[i].Voigt_hwhm_reciprocal = el / gamma;
            line_data[i].Voigt_eta = el / (el + g2);
        }
    }
    /*
     * Compute line mixing parameters.
     */
    if (lineshape == LINESHAPE_VVW_COUPLED) {
        if (line_coupling_table == NULL ) {
            #ifdef _OPENMP
            #pragma omp parallel for schedule(guided, 4096)
            #endif
            for (i = 0; i < nlines; ++i) {
                line_data[i].Y = 0.0;
                line_data[i].G = 0.0;
            }
        } else {
            #ifdef _OPENMP
            #pragma omp parallel for schedule(guided, 4096)
            #endif
            for (i = 0; i < nlines; ++i) {
                double Y, G, dnu, u; 
                u  = line_coupling_table[i].Tref / T;
                Y  = line_coupling_table[i].y0;
                Y += line_coupling_table[i].y1 * (u - 1.0);
                Y *= P * pow(u, line_coupling_table[i].yx);
                line_data[i].Y = Y;
                G  = line_coupling_table[i].g0;
                G += line_coupling_table[i].g1 * (u - 1.0);
                G *= P * P * pow(u, line_coupling_table[i].gx);
                line_data[i].G = G;
                dnu  = line_coupling_table[i].d0;
                dnu += line_coupling_table[i].d1 * (u - 1.0);
                dnu *= P * P * pow(u, line_coupling_table[i].dx);
                line_data[i].f0 += dnu;
            }   
        }
    }
    /*
     * Count unresolved in-band lines
     */
    bracket_in_band_lines(f, ngrid, line_data, nlines, &p, &q);
    if (lineshape == LINESHAPE_DOPPLER) {
        int unres_lines = 0;
        #ifdef _OPENMP
        #pragma omp parallel for schedule(guided, 4096)\
            reduction(+: unres_lines)
        #endif
        for (i = p; i <= q; ++i) {
            if (line_data[i].S == 0.0)
                continue;
            if (line_data[i].alpha < df)
                ++unres_lines;
        }
        *unresolved_line_count = unres_lines;
    }
    else if (
            lineshape == LINESHAPE_FULL_LORENTZ ||
            lineshape == LINESHAPE_GROSS        ||
            lineshape == LINESHAPE_LORENTZ      ||
            lineshape == LINESHAPE_VVW          ||
            lineshape == LINESHAPE_VVW_COUPLED  ||
            lineshape == LINESHAPE_VVH          ||
            lineshape == LINESHAPE_VVH_750        
            ) {
        int unres_lines = 0;
        #ifdef _OPENMP
        #pragma omp parallel for schedule(guided, 4096)\
            reduction(+: unres_lines)
        #endif
        for (i = p; i <= q; ++i) {
            if (line_data[i].S == 0.0)
                continue;
            if (line_data[i].gamma < df)
                ++unres_lines;
        }
        *unresolved_line_count = unres_lines;
    }
    else if (lineshape == LINESHAPE_VOIGT_KIELKOPF) {
        int unres_lines = 0;
        #ifdef _OPENMP
        #pragma omp parallel for schedule(guided, 4096)\
            reduction(+: unres_lines)
        #endif
        for (i = p; i <= q; ++i) {
            if (line_data[i].S == 0.0)
                continue;
            if (line_data[i].Voigt_hwhm_reciprocal * df >= 1.0)
                ++unres_lines;
        }
        *unresolved_line_count = unres_lines;
    }
    /*
     * Allocate scratch arrays for use within linesum_block().
     * dk[] is a buffer for single-line contributions to an
     * accumulating absorption coefficient.  ta1[] and ta2[] are
     * temporary arrays used for loop splitting to improve
     * vectorizability.
     */
    if ((dktol = tol / n_iso_lines) > DBL_EPSILON &&
        lineshape != LINESHAPE_VVW_COUPLED) {
        if ((dk  = (double*)malloc(ngrid * sizeof(double))) == NULL) {
            errlog(160, 0);
            free(line_data);
            return 1;
        }
    }
    if (
            lineshape == LINESHAPE_DOPPLER        ||
            lineshape == LINESHAPE_VOIGT_KIELKOPF
            ) {
        if ((ta1 = (double*)malloc(ngrid * sizeof(double))) == NULL) {
            errlog(160, 0);
            free(dk);
            free(line_data);
            return 1;
        }
    }
    if (lineshape == LINESHAPE_VOIGT_KIELKOPF) {
        if ((ta2 = (double*)malloc(ngrid * sizeof(double))) == NULL) {
            errlog(160, 0);
            free(ta1);
            free(dk);
            free(line_data);
            return 1;
        }
    }
    /*
     * Break the computation over the frequency grid into blocks that
     * will fit in L1 cache.  Under OpenMP, smaller blocks are used as
     * needed to share the work among threads, subject to a minimum
     * block size per thread of LINESUM_MIN_THD_BLOCKSIZE.
     */
    blocksize = ((L1_CACHE_BYTES) == 0) ? ngrid :
        (gridsize_t)((L1_CACHE_BYTES) / ((L1_CACHE_WAYS) * sizeof(double)));
    #ifdef _OPENMP
    {
        gridsize_t thd_blocksize;
        int nthreads;
        int nblocks;
        nthreads = omp_get_max_threads();
        nblocks = (int)(ngrid / LINESUM_MIN_THD_BLOCKSIZE);
        nblocks = nblocks == 0 ? 1 : nblocks;
        nblocks = nthreads < nblocks ? nthreads : nblocks;
        /*
         * Avoid a small remainder block by rounding up on division. 
         */
        thd_blocksize = (ngrid + (nblocks - 1)) / nblocks;
        blocksize = (blocksize < thd_blocksize) ? blocksize : thd_blocksize;
    }
    #endif
    if (blocksize > LINESUM_BLOCK_ALIGN) {
        volatile gridsize_t aligned_blocksize = blocksize;
        aligned_blocksize /= LINESUM_BLOCK_ALIGN;
        aligned_blocksize *= LINESUM_BLOCK_ALIGN;
        blocksize = aligned_blocksize;
    }
    retval = 0;
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 1) reduction(+: retval)
    #endif
    for (offset = 0; offset < ngrid; offset += blocksize) {
        double *kb, *dkb, *ta1b, *ta2b, *fb, *f2b;
        gridsize_t nblock;
        nblock = ngrid - offset;
        if (nblock > blocksize)
            nblock = blocksize;
        kb   = k   + offset;
        dkb  = (dk  != NULL ? dk  + offset : NULL);
        ta1b = (ta1 != NULL ? ta1 + offset : NULL);
        ta2b = (ta2 != NULL ? ta2 + offset : NULL);
        fb  = (double*)f  + offset; 
        f2b = (double*)f2 + offset; 
        retval += linesum_block(
                kb,
                dkb,
                ta1b,
                ta2b,
                fb,
                f2b,
                df,
                nblock,
                line_data,
                nlines,
                lineshape,
                T,
                dktol);
    }
    free(line_data);
    free(dk);
    free(ta1);
    free(ta2);
    return (retval > 0);
}   /* linesum() */


/***********************************************************
* static int linesum_block(
*         double *k,
*         double *dk,
*         double *ta1,
*         double *ta2,
*         const double *f,
*         const double *f2,
*         const double df,
*         const gridsize_t ngrid,
*         const line_data_t *line_data,
*         const int nlines,
*         const int lineshape,
*         const double T,
*         const double dktol)
*
* Purpose:
*   Used by the function linesum() to compute a block of a
*   line-by-line absorption coefficient.  The line by line
*   sum is broken into three passes as follows:
*
*   In the first pass, all lines falling within the
*   frequency subgrid for this block are added in,
*   regardless of line strength.  
*
*   In the second pass, lines centered below the minimum
*   subgrid frequency are added in, in reverse frequency
*   order (i.e.  increasing detuning).  Initially, the
*   minimum line strength threshold is set to zero.  As
*   lines are added in, if a line is found to produce a
*   relative change smaller than dktol for for every k[i],
*   the strength of this line becomes the new (increased)
*   cutoff threshold.  Subsequent lines that are weaker than
*   the threshold are ignored.
*
*   For the third (last) pass, the threshold is reset to
*   zero, and lines centered above the highest subgrid
*   frequency are added in, in ascending frequency order.
*   Again, the threshold is adjusted upwards if possible as
*   the sum proceeds.
*
*   This tolerance checking is bypassed if dk == NULL, in
*   which case every catalog line is accumulated directly
*   into k[].  Note that for coupled lines, in particular
*   for lineshape == VVW_COUPLED, this function should
*   always be called with dk == NULL.
*
*   If lineshape == LINESHAPE_DOPPLER, the temporary array
*   pointer ta1 must point to allocated space for ngrid
*   doubles.
*
*   If lineshape == LINESHAPE_VOIGT_KIELKOPF, then both ta1
*   and ta2 must point to allocated space for ngrid doubles.
*
*   For other lineshapes, ta1 and ta2 can be NULL.
*
* Arguments:
*   double *k              - absorption coefficient [cm^2]
*   double *dk             - single-line abs. coeff. buffer
*   double *ta1            - temporary array
*   double *ta2            - temporary array
*   const double *f        - frequency grid [GHz]
*   const double *f2       - squared frequency grid [GHz^2]
*   const double df        - frequency grid interval [GHz]
*   const gridsize_t ngrid - number of grid points
*   const line_data_t *line_data
*                          - precomputed line data table
*   const int nlines       - number of catalog lines
*   const int lineshape    - index number of lineshape
*   const double T         - temperature [K]
*   const double dktol     - tolerance / num lines in sum
*
* Return:
*   0 if OK, 1 if an error occurred
************************************************************/

static int linesum_block(
        double *k,
        double *dk,
        double *ta1,
        double *ta2,
        const double *f,
        const double *f2,
        const double df,
        const gridsize_t ngrid,
        const line_data_t *line_data,
        const int nlines,
        const int lineshape,
        const double T,
        const double dktol)
{
    int p, q;       /* catalog index range for in-band lines */
    int pass;       /* line sum pass number (0-2, see below) */
    gridsize_t j;

    /*
     * Initialize the absorption coefficient array.
     */
    for (j = 0; j < ngrid; ++j)
        k[j] = 0.0;
    /*
     * Get the catalog index range p -> q covering in-band lines for
     * this block.
     */
    bracket_in_band_lines(f, ngrid, line_data, nlines, &p, &q);
    /*
     * LOOP OVER CATALOG SUBSETS
     * 
     * The catalog lines are divided into three subsets relative to
     * this block of the frequency grid:
     *
     * pass  line center freq.  cat. index range  (empty subset p,q)
     * ---------------------------------------------------------------
     *  0   | in block         | p   -> q        | q = p-1
     *  1   | below block      | p-1 -> 0        | p = 0
     *  2   | above block      | q+1 -> nlines-1 | q = nlines-1
     *
     * During pass 0, the spectral absorption coefficient for each
     * line is always summed directly into the total spectral
     * absorption coefficient k[].  During passes 1 and 2, lines are
     * first buffered into the array dk[] for tolerance checking and
     * possible line strength threshold adjustment, unless dk == NULL,
     * indicating that tolerance checking will not take place
     */
    for (pass = 0; pass <= 2; ++pass) {
        double Smin = 0.0; /* adaptive line strength threshold */
        int i, start, stop, step;
        if (pass == 0) {
            start = p    ; stop = q + 1 ; step =  1;
        } else if (pass == 1) {
            start = p - 1; stop = -1    ; step = -1;
        } else {
            start = q + 1; stop = nlines; step =  1;
        }
        /*
         * LINE CATALOG LOOP
         */
        for (i = start; i != stop; i += step) {
            double f0 = line_data[i].f0;
            double  S = line_data[i].S;
            /*
             * [jmin, jmax) is the range of frequency channels that
             * get touched by a line.  The initial values here cover
             * the entire frequency block.  Truncated lines (like
             * VVH_750) will modify these for each line.
             */
            gridsize_t jmin = 0;
            gridsize_t jmax = ngrid;
            if (S <= Smin)
                continue;
            switch (lineshape) {
            case LINESHAPE_DOPPLER:
                {
                    double r0, r1;
                    r0 = 1.0 / line_data[i].alpha;
                    r1 = S * r0 * ONE_ON_SQRT_PI;
                    for (j = 0; j < ngrid; ++j) {
                        double r2 = (f[j] - f0) * r0;
                        ta1[j] = r2 * r2;
                    }
                    if (dk == NULL || pass == 0)
                        for (j = 0; j < ngrid; ++j)
                            k[j] += r1 * am_exp(-ta1[j]);
                    else
                        for (j = 0; j < ngrid; ++j)
                            dk[j] = r1 * am_exp(-ta1[j]);
                }
                break;
            case LINESHAPE_FULL_LORENTZ:
                {
                    double gamma = line_data[i].gamma;
                    double r0, r1, r2;
                    r0 = S * FOUR_ON_PI * gamma;
                    r1 = gamma * gamma;
                    r2 = 2. * f0;
                    r1 += f0 * f0;
                    if (dk == NULL || pass == 0) {
                        for (j = 0; j < ngrid; ++j) {
                            double rff, r3, r4, r5, r6;
                            rff = f[j];
                            r3 = rff * r2;
                            rff *= rff;
                            r4 = rff + r1;
                            r5 = r4 - r3;
                            r6 = r4 + r3;
                            k[j] += r0 / (r5 * r6);
                        }
                    } else {
                        for (j = 0; j < ngrid; ++j) {
                            double rff, r3, r4, r5, r6;
                            rff = f[j];
                            r3 = rff * r2;
                            rff *= rff;
                            r4 = rff + r1;
                            r5 = r4 - r3;
                            r6 = r4 + r3;
                            dk[j] = r0 / (r5 * r6);
                        }
                    }
                }
                break;
            case LINESHAPE_GROSS:
                {
                    double gamma = line_data[i].gamma;
                    double r0, r1, r2;
                    r0 = S * FOUR_ON_PI * gamma;
                    r1 = 4. * gamma * gamma;
                    r2 = f0 * f0;
                    if (dk == NULL || pass == 0) {
                        for (j = 0; j < ngrid; ++j) {
                            double r3, r4;
                            r3 = r1 * f2[j];
                            r4 = f2[j] - r2;
                            r4 *= r4;
                            k[j] += r0 / (r4 + r3);
                        }
                    } else {
                        for (j = 0; j < ngrid; ++j) {
                            double r3, r4;
                            r3 = r1 * f2[j];
                            r4 = f2[j] - r2;
                            r4 *= r4;
                            dk[j] = r0 / (r4 + r3);
                        }
                    }
                }
                break;
            case LINESHAPE_LORENTZ:
                {
                    double gamma = line_data[i].gamma;
                    double r0, r1;
                    r0 = S * ONE_ON_PI * gamma;
                    r1 = gamma * gamma;
                    if (dk == NULL || pass == 0) {
                        for (j = 0; j < ngrid; ++j) {
                            double r2;
                            r2 = f[j] - f0;
                            r2 *= r2;
                            k[j] += r0 / (r1 + r2);
                        }
                    } else {
                        for (j = 0; j < ngrid; ++j) {
                            double r2;
                            r2 = f[j] - f0;
                            r2 *= r2;
                            dk[j] = r0 / (r1 + r2);
                        }
                    }
                }
                break;
            case LINESHAPE_VVW:
                {
                    double gamma = line_data[i].gamma;
                    double r0, r1, r2;
                    r0 = f0 * f0;
                    r1 = S * TWO_ON_PI * gamma / r0;
                    r0 += gamma * gamma;
                    r2 = 2. * f0;
                    if (dk == NULL || pass == 0) {
                        for (j = 0; j < ngrid; ++j) {
                            double rff, r3, r4, r5, r6;
                            rff = f[j];
                            r3 = rff * r2;
                            rff *= rff;
                            r4 = rff + r0;
                            r5 = r4 - r3;
                            r6 = r4 + r3;
                            k[j] += r1 * r4 / (r5 * r6);
                        }
                    } else {
                        for (j = 0; j < ngrid; ++j) {
                            double rff, r3, r4, r5, r6;
                            rff = f[j];
                            r3 = rff * r2;
                            rff *= rff;
                            r4 = rff + r0;
                            r5 = r4 - r3;
                            r6 = r4 + r3;
                            dk[j] = r1 * r4 / (r5 * r6);
                        }
                    }
                }
                break;
            case LINESHAPE_VVW_COUPLED:
                {
                    double gamma = line_data[i].gamma;
                    double r0, r1, r2;
                    double Y = line_data[i].Y;
                    r0 = S / (PI * f0 * f0);
                    r1 = gamma * gamma;
                    r2 = gamma * (1.0 + line_data[i].G);
                    for (j = 0; j < ngrid; ++j) {
                        double rf, r3, r4, r5, r6;
                        rf = f[j];
                        r3 = rf - f0;
                        r4 = rf + f0;
                        r5 = r2 + Y * r3;
                        r3 *= r3;
                        r3 += r1;
                        r6 = r2 - Y * r4;
                        r4 *= r4;
                        r4 += r1;
                        k[j] += r0 * (r5 / r3 + r6 / r4);
                    }
                }
                break;
            case LINESHAPE_VVH:
                {
                    double gamma = line_data[i].gamma;
                    double r0, r1, r2;
                    r0 = f0 * f0 + gamma * gamma;
                    r1 = S * TWO_ON_PI * gamma;
                    r1 /= f0 * tanh(0.5 * H_ON_KB * f0 / T);
                    r2 = 2. * f0;
                    if (dk == NULL || pass == 0) {
                        for (j = 0; j < ngrid; ++j) {
                            double rf, r3, r4, r5, r6;
                            rf = f[j];
                            r3 = rf * r2;
                            r4 = rf * rf + r0;
                            r5 = r4 - r3;
                            r6 = r4 + r3;
                            k[j] += r1 * r4 / (r5 * r6);
                        }
                    } else {
                        for (j = 0; j < ngrid; ++j) {
                            double rf, r3, r4, r5, r6;
                            rf = f[j];
                            r3 = rf * r2;
                            r4 = rf * rf + r0;
                            r5 = r4 - r3;
                            r6 = r4 + r3;
                            dk[j] = r1 * r4 / (r5 * r6);
                        }
                    }
                }
                break;
            case LINESHAPE_VVH_750:
                {
                    double gamma = line_data[i].gamma;
                    double r0, r1, r2;
                    double fc, f_min, f_max;
                    gridsize_t jmax_neg_line;
                    fc = lineshape_type[LINESHAPE_VVH_750].f_cutoff;
                    r0 = S * ONE_ON_PI * gamma;
                    r0 /= f0 * tanh(0.5 * H_ON_KB * f0 / T);
                    r1 = gamma * gamma;
                    r2 = r0 / (r1 + fc * fc);
                    f_min = f0 - fc;
                    f_max = f0 + fc;
                    /*
                     * If the positive frequency line is outside the
                     * grid, so are the negative frequency line and
                     * all remaining lines in this pass.  Break out
                     * with jmax < 0 to signal the end of this catalog
                     * pass.
                     */
                    if (f_min > f[ngrid-1] || f_max < f[0]) {
                        jmin = 0;
                        jmax = -1;
                        break;
                    }
                    jmin = f_min < f[0] ?
                        0 : (gridsize_t)ceil((f_min - f[0]) / df);
                    jmax = f_max > f[ngrid-1] ?
                        ngrid :
                        ngrid - (gridsize_t)floor((f[ngrid-1] - f_max) / df);
                    if (dk == NULL || pass == 0) {
                        for (j = jmin; j < jmax; ++j) {
                            double r3;
                            r3 = f[j] - f0;
                            r3 *= r3;
                            k[j] += r0 / (r1 + r3) - r2;
                        }
                    } else {
                        for (j = jmin; j < jmax; ++j) {
                            double r3;
                            r3 = f[j] - f0;
                            r3 *= r3;
                            dk[j] = r0 / (r1 + r3) - r2;
                        }
                    }
                    f_max = -f0 + fc;
                    /*
                     * If jmin > 0, then the truncated positive
                     * frequency line lies entirely at f > 0.
                     * This means the truncated negative line
                     * lies entirely at f < 0 and makes no
                     * contribution, so break out here.
                     */
                    if (jmin > 0)
                        break;
                    /*
                     * If there is a negative line contribution,
                     * it will always start at f[0], but the
                     * upper truncation point may be below
                     * f[ngrid-1].  Here, find the maximum grid
                     * index.
                     */
                    jmax_neg_line = f_max > f[ngrid-1] ?
                        ngrid :
                        ngrid - (gridsize_t)floor((f[ngrid-1] - f_max) / df);
                    if (dk == NULL || pass == 0) {
                        for (j = 0; j < jmax_neg_line; ++j) {
                            double r3;
                            r3 = f[j] + f0;
                            r3 *= r3;
                            k[j] += r0 / (r1 + r3) - r2;
                        }
                    } else {
                        /* Sum into dk[] in line overlap area */
                        for (j = 0; j < jmax_neg_line; ++j) {
                            double r3;
                            r3 = f[j] + f0;
                            r3 *= r3;
                            dk[j] += r0 / (r1 + r3) - r2;
                        }
                    }
                }
                break;
            case LINESHAPE_VOIGT_KIELKOPF:
                {
                    double eta = line_data[i].Voigt_eta;
                    double r1 = 1.0 - eta;
                    double r2 = line_data[i].Voigt_hwhm_reciprocal;
                    for (j = 0; j < ngrid; ++j) {
                        double x = r2 * (f[j] - f0);
                        double x2 = x * x;
                        ta1[j] = x2;
                    }
                    for (j = 0; j < ngrid; ++j) {
                        const double ln2 = 0.6931471805599453;
                        ta2[j] = am_exp(-ln2 * ta1[j]);
                    }
                    if (dk == NULL || pass == 0) {
                        for (j = 0; j < ngrid; ++j) {
                            double x2 = ta1[j];
                            double G  = ta2[j];
                            double L  = 1.0 / (1.0 + x2);
                            double E  = (0.8029 - 0.4207 * x2) /
                                (1.0 + x2 * (0.2030 + 0.07335 * x2));
                            double U  = r1 * G + eta * L +
                                eta * r1 * E * (G - L);
                            k[j] += U * S;
                        }
                    } else {
                        for (j = 0; j < ngrid; ++j) {
                            double x2 = ta1[j];
                            double G  = ta2[j];
                            double L  = 1.0 / (1.0 + x2);
                            double E  = (0.8029 - 0.4207 * x2) /
                                (1.0 + x2 * (0.2030 + 0.07335 * x2));
                            double U  = r1 * G + eta * L +
                                eta * r1 * E * (G - L);
                            dk[j] = U * S;
                        }
                    }
                }
                break;
            default:
                errlog(60, lineshape);
                return 1;
            }
            /*
             * For truncated lines (like VVH_750), jmax < 0 indicates
             * that for this line and all remaining lines in this
             * catalog subset pass, the truncated line wings won't
             * touch this frequency interval. If so, the pass is
             * finished.
             */
            if (jmax < 0)
                break;
            /*
             * If tolerance is being monitored on this pass, a line is
             * accumulated into the total spectral absorption
             * coefficient k[] here.
             *
             * If the fractional change dk/k < dktol for every
             * frequency channel in this frequency block, the strength
             * of the current line becomes the new threshold strength
             * for ignoring lines at greater detuning from the block.
             */
            if (dk != NULL && pass != 0) {
                for (j = jmin; j < jmax; ++j)
                    k[j] += dk[j];
                for (j = jmin; j < jmax; ++j)
                    if (dk[j] > dktol * k[j])
                        break;
                if (j >= jmax)
                    Smin = S * (1. + DBL_EPSILON);
            }
        } /* END OF LINE CATALOG LOOP */
    } /* END OF CATALOG SUBSET LOOP */
    /*
     * Common factors that depend on f only.
     */
    switch (lineshape) {
    case LINESHAPE_DOPPLER:
    case LINESHAPE_LORENTZ:
    case LINESHAPE_VOIGT_KIELKOPF:
        break;
    case LINESHAPE_FULL_LORENTZ:
    case LINESHAPE_GROSS:
    case LINESHAPE_VVW:
    case LINESHAPE_VVW_COUPLED:
        for (j = 0; j < ngrid; ++j)
            k[j] *= f2[j];
        break;
    case LINESHAPE_VVH:
    case LINESHAPE_VVH_750:
        {
            double r0 = 0.5 * H_ON_KB / T;
            for (j = 0; j < ngrid; ++j)
                k[j] *= f[j] * tanh(r0 * f[j]);
        }
        break;
    default:
        errlog(60, lineshape);
        return 1;
    }
    return 0;
}   /* linesum_block() */


/***********************************************************
* int Qratio(
*         const double T,
*         const double Tref,
*         const double *Qtab,
*         const int nrows,
*         const int ncols,
*         double *Qrat)
*
* Purpose:
*   For a molecule with partition sum table *Qtab, computes
*   the ratio Q(Tref) / Q(T) for each isotopologue.  See
*   the function Qtab_interp(), below, for a description of
*   the layout of *Qtab.
*
* Arguments:
*   const double T     - temperature [K]
*   const double Tref  - line catalog reference
*                        temperature [K]
*   const double *Qtab - partition sum table
*   const int nrows    - number of table rows
*   const int ncols    - number of table columns
*   double *Qrat       - pointer to array to hold computed
*                        Q ratio
*
* Return:
*   0 if OK, 1 if an error occurred
************************************************************/

int Qratio(
        const double T,
        const double Tref,
        const double *Qtab,
        const int nrows,
        const int ncols,
        double *Qrat
    )
{
    double Q[MAX_ISO];
    int iso;

    if (Qtab_interp(Tref, Qtab, nrows, ncols, Qrat))
        return 1;
    if (Qtab_interp(T, Qtab, nrows, ncols, Q))
        return 1;
    for (iso = 1; iso < ncols; ++iso)
        Qrat[iso] /= Q[iso];
    return 0;
}   /* Qratio() */


/***********************************************************
* static int Qtab_interp(
*         const double T,
*         const double *Qtab,
*         const int nrows,
*         const int ncols,
*         double *Q)
*
* Purpose:
*   Linearly interpolates the partition sum table for a
*   given molecule to a temperature T.  The interpolated
*   values are stored in the array *Q, indexed by HITRAN
*   isotopologue number.
*
*   The partition sum table *Qtab is an array of doubles,
*   conceptually arranged in rows of the form:
*
*       T0 Q1(T0) Q2(T0) ... Q_ncols-1(T0)
*       T1 Q1(T1) Q2(T1) ... Q_ncols-1(T1)
*       .     .      .          .
*
*   and so on.  The temperatures T0 .. T_nrows-1 are in
*   ascending order, but not necessarily equally spaced.
*   Note that the number of tabulated isotopologues is
*   ncols-1, and the column indices 1 .. ncols-1 are the
*   HITRAN isotopologue number.
*
*   If T falls outside the tabulated range, the return
*   value is set to 1.
*
* Arguments:
*   const double T     - temperature [K]
*   const double *Qtab - partition sum table
*   const int nrows    - number of table rows
*   const int ncols    - number of table columns
*   double *Q          - pointer to array to receive
*                        interpolated Q(T) values,
*                        indexed by isotopologue.
*
* Return:
*   0 if OK, 1 if an error occurred
************************************************************/

static int Qtab_interp(
        const double T,
        const double *Qtab,
        const int nrows,
        const int ncols,
        double *Q
    )
{
    int ilow, ihigh;
    int iso;

    /*
     * Check that T is in the tabulated range
     */
    if (T < Qtab[0] || T > Qtab[(nrows - 1) * ncols])
        return 1;
    /*
     * Binary search for table rows bounding T.  During
     * the search, ilow and ihigh are row numbers.
     */
    ilow  = 0;
    ihigh = nrows - 1;
    while (ihigh - ilow > 1) {
        int imid = (ilow + ihigh) / 2;
        if (T > Qtab[imid * ncols])
            ilow  = imid;
        else
            ihigh = imid;
    }
    /*
     * Scale ilow and ihigh by ncols to point directly to
     * the first column in their respective rows.
     */
    ilow  *= ncols;
    ihigh *= ncols;
    /*
     * Interpolate the table for each isotopologue.
     */
    for (iso = 1; iso < ncols; ++iso) {
        Q[iso]  = Qtab[ihigh + iso] - Qtab[ilow + iso];
        /*
         * Qtab is a static table constructed such that this
         * denominator is always non-zero.
         */
        Q[iso] /= Qtab[ihigh] - Qtab[ilow];
        Q[iso] *= T - Qtab[ilow];
        Q[iso] += Qtab[ilow + iso];
    }
    return 0;
}   /* Qtab_interp() */
