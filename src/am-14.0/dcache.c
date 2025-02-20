/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* dcache.c                        S. Paine rev. 2024 July 23
*
* Disk caching of absorption coefficients.
************************************************************/

#include <float.h>
#include <math.h>
#include <signal.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "abscoeff.h"
#include "am_types.h"
#include "dcache.h"
#include "errlog.h"
#include "fileops.h"
#include "phys_const.h"
#include "version.h"

enum {
    CACHE_FILES_PER_HASHVAL = 4,
    HASH_MODULUS_DEFAULT    = 1021,
    HASH_MODULUS_MAX        = 16384,
    HASH_MULTIPLIER         = 37
};

static const char DCACHE_DIRPATH_ENV_STRING[] = "AM_CACHE_PATH";
static const char DCACHE_HASHMOD_ENV_STRING[] = "AM_CACHE_HASH_MODULUS";
static const char DCACHE_FILE_ID_STRING[]     = "am cache file\x04\x1A";

struct cache_header {
    unsigned int version;
    unsigned int hashval;
    int          k_typenum;
    int          lineshapenum;
    int          unres_lines;
    gridsize_t   ngrid;
    double       fmin;
    double       fmax;
    double       df;
    double       tol;
    double       P;
    double       vmr;
    double       T;
};

/*
 * Log data
 */
static unsigned long dcache_hit_count     = 0;
static unsigned long dcache_miss_count    = 0;
static unsigned long dcache_discard_count = 0;


static char         *dcache_path(void);
static int          dcache_path_valid(void);
static unsigned int dcache_version(void);
static unsigned int hash(
        double, double, double, double, int, int, unsigned int);
static unsigned int hash_modulus(void);
static void         init_cache_hdr_struct(
        struct cache_header*, const model_t*, int, int, int, double);
static int          insert_as_mru(
        struct cache_header*, char*new_mrufile, char*);
static int          is_cache_hit(struct cache_header*, struct cache_header*);
static int          make_cache_file_pathstr(
        char*, const char*, unsigned int, int);
static int          promote_to_mru(char*, unsigned int, int);
static int          read_cache_file_data(
        model_t*, double*, int*, struct cache_header*, FILE*);
static int          read_cache_file_hdr(struct cache_header*, FILE*);
static int          write_cache_file(struct cache_header*, double*, FILE*);


/***********************************************************
* void dcache_log(int event)
*
* Purpose:
*   Keeps a count of dcache events, including hits, misses
*   and discards.
*
* Arguments:
*   int event - event type (defined in dcache.h)
************************************************************/

void dcache_log(int event)
{
    switch (event) {
        case DCACHE_HIT:
            ++dcache_hit_count;
            break;
        case DCACHE_MISS:
            ++dcache_miss_count;
            break;
        case DCACHE_DISCARD:
            ++dcache_discard_count;
            break;
        case DCACHE_CLEAR_LOG:
            dcache_hit_count     = 0;
            dcache_miss_count    = 0;
            dcache_discard_count = 0;
            break;
        default:
            break;
    }
    return;
}   /* dcache_log() */


/***********************************************************
* int dcache_lookup_absorption_coefficient(
*         model_t *model,
*         int lnum,
*         int cnum,
*         int knum,
*         double *k,
*         double T)
*
* Purpose:
*   Looks in the disk cache for a file containing the
*   currently needed absorption coefficient knum, of column
*   cnum, of layer lnum.  If a matching file is found, the
*   absorption coefficient is read out of the file, and the
*   file is promoted to the MRU slot of its hash bucket.
*
*   The address k of the absorption coefficient array, and
*   the temperature T, are passed as separate parameters to
*   accommodate computation of absorption coefficients on the
*   fixed T grid used for the memory cache (kcache).
*
* Arguments:
*   model_t *model - pointer to model structure
*   int lnum       - layer number
*   int cnum       - column number
*   int knum       - absorption coefficient number
*   double *k      - pointer to receive absorption
*                    coefficient array
*   double T       - temperature [K]
*
* Return:
*   0 if the absorption coefficient was recovered from the
*     cache
*   1 otherwise
************************************************************/

int dcache_lookup_absorption_coefficient(
        model_t *model,
        int lnum,
        int cnum,
        int knum,
        double *k,
        double T)
{
    char       *cache_dirpath;
    struct     cache_header requested_hdr, cached_hdr;
    char       cache_filepath[AM_MAX_PATH];
    FILE       *cache_fileptr;
    int        slot;
    abscoeff_t *abscoeff = model->layer[lnum]->column[cnum]->abscoeff[knum];

    /*
     * If there is no cache directory, or if the hash modulus is
     * invalid, do nothing.
     */
    if (((cache_dirpath = dcache_path()) == NULL) || (hash_modulus() == 0))
        return 1;
    /*
     * If this is not a cacheable absorption coefficient type, do
     * nothing.
     */
    if (!(k_type[abscoeff->k_typenum].comp_flags & CACHEABLE))
        return 1;
    /*
     * Set up a cache_header structure with the parameters for
     * the requested absorption coefficient, including their hash
     * value.
     */
    init_cache_hdr_struct(&requested_hdr, model, lnum, cnum, knum, T);
    /*
     * Look for a cache hit in the appropriate hash bucket.
     * Conceptually, each hash bucket has slots for files
     * numbered from 1 up to a maximum of CACHE_FILES_PER_HASHVAL
     * in order of time of use, with slot 1 being the most
     * recently used (MRU) file.  (Slot 0 holds a temporary lock
     * file, used to prevent more than one am process reordering
     * a hash bucket simultaneously.)
     */
    for (slot = 1; slot <= CACHE_FILES_PER_HASHVAL; ++slot) {
        if (make_cache_file_pathstr(
                cache_filepath,
                cache_dirpath,
                requested_hdr.hashval,
                slot)) {
            return 1;
        }
        cache_fileptr = fopen(cache_filepath, "rb");
        /*
         * If there is no file in this slot, continue on to the
         * next.  That is, don't assume there are no filled slots
         * beyond this one-- there might be a hole left by a
         * deleted file.
         */
        if (cache_fileptr == NULL)
            continue;
        /*
         * Read the file header.  If this fails, the file might
         * be corrupt or binary-incompatible.  Continue on to the
         * next file.  Note that we don't want to try deleting a
         * corrupt file now, because the association between file
         * name and file pointer is not protected.  Instead, am
         * relies on any corrupt file (which will always be a
         * cache miss) eventually dropping to the LRU slot and
         * out of the cache.
         */
        if (read_cache_file_hdr(&cached_hdr, cache_fileptr)) {
            fclose(cache_fileptr);
            continue;
        }
        /*
         * Got a good file-- see if it has the requested data.
         * If so, read the data and promote this file to the MRU
         * slot.  Otherwise, close this file and continue to the
         * file in the next slot.
         */
        if (is_cache_hit(&requested_hdr, &cached_hdr)) {
            void (*saved_handler)(int);
            if (read_cache_file_data(
                        model,
                        k,
                        &(abscoeff->unres_lines),
                        &cached_hdr,
                        cache_fileptr)) {
                fclose(cache_fileptr);
                continue;
            }
            fclose(cache_fileptr);
            /*
             * The call to promote_to_mru() is protected from
             * SIGINT (Ctrl-C) to guard against leaving an
             * abandoned lock file.
             */
            saved_handler = signal(SIGINT, SIG_IGN);
            promote_to_mru(cache_dirpath, requested_hdr.hashval, slot);
            signal(SIGINT, saved_handler);
            dcache_log(DCACHE_HIT);
            return 0;
        } else {
            fclose(cache_fileptr);
        }
    }
    dcache_log(DCACHE_MISS);
    return 1;
}   /* dcache_lookup_absorption_coefficient() */


/***********************************************************
* static char *dcache_path(void)
*
* Purpose:
*   On the first call, looks for an environment string
*   naming a path to a cache directory.  If a path is found,
*   a pointer to a static string containing the path is
*   returned.  Subsequent calls to dcache_path() will return
*   the same pointer, without checking the environment.
*
*   If a cache directory is not named in the environment, or
*   if the environment string is empty, NULL is returned.
*
* Return:
*   Pointer to cache path string, or NULL if no cache path
*   is named in the environment.
************************************************************/

static char *dcache_path(void)
{
    static char *pathptr;
    static char path[AM_MAX_DIRPATH];
    static int  initialized = 0;

    if (initialized)
        return pathptr;

    pathptr = getenv(DCACHE_DIRPATH_ENV_STRING);
    if (pathptr != NULL) {
        if (pathptr[0] == '\0') { /* empty path string */
            pathptr = NULL;
        } else {
            if ((size_t)snprintf(
                        path,
                        sizeof(path),
                        "%s",
                        pathptr
                        ) > sizeof(path)) {
                errlog(70, 0); /* environment path string too long */
                pathptr = NULL;
            } else if (check_for_dir_separator(path)) {
                errlog(70, 0); /* path too long to append separator */
                pathptr = NULL;
            } else {
                pathptr = (char *)path;
            }
        }
    }
    initialized = 1;
    return pathptr;
}   /* dcache_path() */


/***********************************************************
* static int dcache_path_valid(void)
*
* Purpose:
*   Tests the validity of the dcache path, if any, named in
*   the environment by attempting to open, close, and delete
*   a temporary file in the named directory.
*
* Return:
*   1 if the test succeeds
*   0 if the test fails, or if no dcache path is named in
*     the environment.
************************************************************/

static int dcache_path_valid(void)
{
    FILE *fp;
    char tmpfpath[AM_MAX_PATH];

    if (dcache_path() == NULL)
        return 0;
    if ((fp = am_tmpfile(dcache_path(), tmpfpath)) == NULL) {
        return 0;
    }
    fclose(fp);
    if (remove_with_retry(tmpfpath)) {
        errlog(75, __LINE__);
        return 0;
    }
    return 1;
}   /* dcache_path_valid() */


/***********************************************************
* int dcache_save_absorption_coefficient(
*         model_t *model,
*         int lnum,
*         int cnum,
*         int knum,
*         double *k,
*         double T)
*
* Purpose:
*   Writes a temporary cache file containing the header
*   information and data for absorption coefficient knum, of
*   column cnum, of layer lnum of the model, then inserts
*   this file into the MRU slot of the appropriate hash
*   bucket.
*
*   The address k of the absorption coefficient array, and
*   the temperature T, are passed as separate parameters to
*   accommodate computation of cached absorption coefficients
*   on a fixed T grid.
*
* Arguments:
*   model_t *model - pointer to model structure
*   int lnum       - layer number
*   int cnum       - column number
*   int knum       - absorption coefficient number
*   double *k      - pointer to absorption coefficient array
*                    to be saved
*   double T       - temperature [K]
*
* Return:
*   0 if successful
*   1 otherwise
************************************************************/

int dcache_save_absorption_coefficient(
        model_t *model,
        int lnum,
        int cnum,
        int knum,
        double *k,
        double T)
{
    struct cache_header hdr;
    void (*saved_handler)(int);
    FILE *fp;
    char tmpfpath[AM_MAX_PATH];
    char *cache_dirpath;
    int  status;
    abscoeff_t *abscoeff = model->layer[lnum]->column[cnum]->abscoeff[knum];

    /*
     * If there is no cache directory, or the hash modulus is
     * invalid, do nothing.
     */
    if (((cache_dirpath = dcache_path()) == NULL) || (hash_modulus() == 0))
        return 1;
    /*
     * If this is not a cacheable absorption coefficient type, do
     * nothing.
     */
    if (!(k_type[abscoeff->k_typenum].comp_flags & CACHEABLE))
        return 1;
    /*
     * Write a temporary cache file.  If it cannot be created, or
     * if an error occurs on write, log an error and return.
     */
    if ((fp = am_tmpfile(cache_dirpath, tmpfpath)) == NULL) {
        errlog(23, 0);
        return 1;
    }
    init_cache_hdr_struct(&hdr, model, lnum, cnum, knum, T);
    if (write_cache_file(&hdr, k, fp)) {
        errlog(22, 0);
        fclose(fp);
        if (remove_with_retry(tmpfpath)) {
            errlog(75, __LINE__);
        }
        return 1;
    }
    fclose(fp);
    /*
     * Insert the temporary file into the MRU slot of the
     * appropriate cache bucket.  The call to insert_as_mru() is
     * protected from SIGINT (Ctrl-C) to guard against leaving an
     * abandoned lock file.
     */
    saved_handler = signal(SIGINT, SIG_IGN);
    status = insert_as_mru(&hdr, tmpfpath, cache_dirpath);
    if (status) {
        if (remove_with_retry(tmpfpath)) {
            errlog(75, __LINE__);
        }
    }
    signal(SIGINT, saved_handler);
    return status;
}   /* dcache_save_absorption_coefficient() */


/***********************************************************
* static unsigned int dcache_version(void)
*
* Purpose:
*   Returns a version number for dcache files.  This
*   function is a wrapper around the program version
*   numbers, and can be used to implement different dcache
*   file versioning schemes as needed.
*
*   In particular, in early versons of am, the major version
*   number was used to track changes in spectrosopic data,
*   so cache files were tagged by major version only.  This
*   is no longer the case, and dcache files are now tagged
*   with an integer composed of both the major and minor
*   version numbers.
*
* Return:
*   dcache version number, as an unsigned int.
************************************************************/

static unsigned int dcache_version(void)
{
    unsigned int version_num;

    version_num = (unsigned int)AM_VERSION_MAJOR;
    version_num *= 0x10000;
    version_num += (unsigned int) AM_VERSION_MINOR;
    return version_num;
}   /* dcache_version() */


/***********************************************************
* static unsigned int hash(
*         double df,
*         double P,
*         double vmr,
*         double T,
*         int lineshapenum,
*         int k_typenum,
*         int mod)
*
* Purpose:
*   Compute a hash value for indexing the disk cache of
*   computed absorption coefficients.  The hash value, in
*   the range 0 to mod-1, is a function of the frequency
*   grid spacing, pressure, mixing ratio, temperature,
*   column type, and lineshape type associated with the
*   absorption coefficient.  Frequency, pressure, and
*   temperature are in am native units.
*
*   The hash value is computed from the individual bytes of
*   these numbers, without regard for the byte order, or
*   integer size of the machine, so this function will
*   produce different hash values on different machine
*   architectures.
*
* Arguments:
*   double df        - frequency grid spacing
*   double P         - pressure
*   double vmr       - volume mixing ratio
*   double T         - temperature
*   int lineshapenum - lineshape type number
*   int k_typenum    - absorption coefficient type number
*   unsigned int mod - hash modulus
*
* Return:
*   hash value in the range [0, mod-1], as an unsigned int
************************************************************/

static unsigned int hash(
        double df,
        double P,
        double vmr,
        double T,
        int lineshapenum,
        int k_typenum,
        unsigned int mod)
{
    unsigned int i;
    unsigned int h = 0;
    unsigned char *byte;

    byte = (unsigned char*)(&df);
    for (i = 0; i < sizeof(df); ++i) {
        h = HASH_MULTIPLIER * h + byte[i];
    }
    byte = (unsigned char*)(&P);
    for (i = 0; i < sizeof(P); ++i) {
        h = HASH_MULTIPLIER * h + byte[i];
    }
    byte = (unsigned char*)(&vmr);
    for (i = 0; i < sizeof(vmr); ++i) {
        h = HASH_MULTIPLIER * h + byte[i];
    }
    byte = (unsigned char*)(&T);
    for (i = 0; i < sizeof(T); ++i) {
        h = HASH_MULTIPLIER * h + byte[i];
    }
    byte = (unsigned char*)(&lineshapenum);
    for (i = 0; i < sizeof(lineshapenum); ++i) {
        h = HASH_MULTIPLIER * h + byte[i];
    }
    byte = (unsigned char*)(&k_typenum);
    for (i = 0; i < sizeof(k_typenum); ++i) {
        h = HASH_MULTIPLIER * h + byte[i];
    }
    return h % mod;
}   /* hash() */


/***********************************************************
* static unsigned int hash_modulus(void)
*
* Purpose:
*   On the first call, checks whether the hash modulus has
*   been specified in the environment.  If so, the value is
*   stored in an static variable and returned.  If not, a
*   default value is assigned.  In either case, subsequent
*   calls return the same value without checking the
*   environment again.
*
* Return:
*   hash modulus, as an unsigned int
************************************************************/

static unsigned int hash_modulus(void)
{
    static unsigned int hashmod = HASH_MODULUS_DEFAULT;
    static int initialized = 0;
    char *str;
    char *endp;

    if (initialized)
        return hashmod;

    if ((str = getenv(DCACHE_HASHMOD_ENV_STRING)) != NULL) {
        hashmod = (unsigned int)strtoul(str, &endp, 0);
        if (endp[0] != '\0' || hashmod < 1) {
            errlog(73, 0);
            hashmod = 0;
        } else if (hashmod > HASH_MODULUS_MAX) {
            errlog(74, (int)HASH_MODULUS_MAX);
            hashmod = 0;
        }
    }
    initialized = 1;
    return (unsigned int)hashmod;
}   /* hash_modulus() */


/***********************************************************
* static void init_cache_hdr_struct(
*         struct cache_header *hdr,
*         model_t *model,
*         int lnum,
*         int cnum,
*         int knum,
*         double T)
*
* Purpose:
*   Fills in a cache_hdr_data structure with model data for
*   column cnum of layer lnum, program version information,
*   and the hash value generated from the model data.  The
*   temperature T is passed separately to accommodate
*   computation on the fixed temperature grid points of the
*   kcache.
*
* Arguments:
*   struct cache_hdr *hdr - pointer to cache_hdr structure
*   model_t *model        - pointer to model structure
*   int lnum              - layer number
*   int cnum              - column number
*   int knum              - absorption coefficient number
*   double T              - Temperature [K]
************************************************************/

static void init_cache_hdr_struct(
        struct cache_header *hdr,
        const model_t *model,
        int lnum,
        int cnum,
        int knum,
        double T)
{
    int        k_typenum;
    int        dep_flags;
    layer_t    *layer;
    column_t   *column;
    abscoeff_t *abscoeff;

    layer     = model->layer[lnum];
    column    = layer->column[cnum];
    abscoeff  = column->abscoeff[knum];
    k_typenum = abscoeff->k_typenum;
    dep_flags = k_type[k_typenum].dep_flags;
    /*
     * Header fields upon which the absorption coefficient does
     * not depend are defaulted to zero.
     */
    hdr->version      = dcache_version();
    hdr->k_typenum    = k_typenum;
    hdr->lineshapenum =
        dep_flags & DEP_ON_LSHAPE ? layer->lineshape[k_typenum] : 0;
    hdr->unres_lines  = abscoeff->unres_lines;
    hdr->ngrid        = model->ngrid;
    hdr->fmin         = model->fmin;
    hdr->fmax         = model->fmax;
    hdr->df           = model->df;
    hdr->tol          = dep_flags & DEP_ON_TOL ? model->tol : 0.0;
    hdr->P            = dep_flags & DEP_ON_P ? layer->P : 0.0;
    hdr->vmr          =
        dep_flags & DEP_ON_VMR_SELFBROAD ? abscoeff->vmr_selfbroad : 0.0;
    hdr->T            = dep_flags & DEP_ON_T ? T : 0.0;
    hdr->hashval      = hash(
                            hdr->df,
                            hdr->P,
                            hdr->vmr,
                            hdr->T,
                            hdr->lineshapenum,
                            hdr->k_typenum,
                            hash_modulus());
    return;
}   /* init_cache_hdr_struct() */


/***********************************************************
* static int insert_as_mru(
*         struct cache_header *new_hdr,
*         char *new_file,
*         char *cache_dirpath)
*
* Purpose:
*   Inserts the file named by new_file into the MRU slot of
*   its proper hash bucket in the cache directory named by
*   cache_dirpath.  The hash value is contained in the
*   cache_header structure new_hdr.  If the hash bucket is
*   full, the file in the LRU slot is evicted.
*
*   On success, either (1) new_file will have been renamed
*   to the MRU slot name, or (2) if the data in new_file had
*   already been placed in the cache, new_file, having first
*   been renamed to the lock file slot for the hash bucket,
*   will be removed to clear the lock.
*
*   On failure, new_file will retain its original name.  If
*   it was a temporary file, it will need to be removed by
*   the calling function.
*
* Arguments:
*   struct cache_header *new_hdr - header data struct for
*                         file to be inserted into MRU slot
*   char *new_file      - full path of file to be inserted
*   char *cache_dirpath - path string for cache directory
*
* Return:
*   0 if the new file was inserted into the MRU slot, or if
*     the data was already placed in the cache by another
*     process.
*   1 otherwise.
************************************************************/

static int insert_as_mru(
        struct cache_header *new_hdr,
        char *new_file,
        char *cache_dirpath)
{
    FILE *fp;
    struct cache_header cached_hdr;
    char lockpath[AM_MAX_PATH];
    char slotpath[AM_MAX_PATH];
    char old_slotpath[AM_MAX_PATH];
    int slot;

    /*
     * Check for a lock file in slot 0.  If one exists, another
     * process is working with this hash bucket.  Give up on
     * inserting the new file.
     */
    if (make_cache_file_pathstr(
                lockpath,
                cache_dirpath,
                new_hdr->hashval,
                0)) {
        return 1;
    }
    if ((fp = fopen(lockpath, "rb")) != NULL) {
        fclose(fp);
        return 1;
    }
    /*
     * Got past the lock file check, so try putting the new file
     * into the lock slot.  If the rename() fails, another
     * process could have won a race for the lock file.  If so,
     * give up on inserting the new file.
     */
    if (rename(new_file, lockpath))
        return 1;
    /*
     * Now that we have the lock on this hash bucket, make sure
     * that another process hasn't already inserted a cache entry
     * which would make this one redundant.  If so, clear the
     * lock by deleting the file, and return without modifying
     * the cache bucket.
     */
    for (slot = 1; slot <= CACHE_FILES_PER_HASHVAL; ++slot) {
        if (make_cache_file_pathstr(
                    slotpath,
                    cache_dirpath,
                    new_hdr->hashval,
                    slot)) {
            return 1;
        }
        if ((fp = fopen(slotpath, "rb")) == NULL)
            continue;
        if (read_cache_file_hdr(&cached_hdr, fp)) {
            fclose(fp);
            continue;
        }
        fclose(fp);
        if (is_cache_hit(new_hdr, &cached_hdr)) {
            if (remove_with_retry(lockpath)) {
                errlog(75, __LINE__);
            }
            return 0;
        }
    }
    /*
     * All checks are satisfied; insert the file.  Start by
     * removing the file (if any) in the LRU slot of this hash
     * bucket.  Then, increment the slot number of all other
     * existing files.  The last one to be incremented is the new
     * file in slot 0 (the lock slot), which gets put into the
     * MRU slot (slot 1), clearing the lock.
     */
    if (make_cache_file_pathstr(
            slotpath,
            cache_dirpath,
            new_hdr->hashval,
            CACHE_FILES_PER_HASHVAL)) {
        return 1;
    }
    if ((fp = fopen(slotpath, "rb")) != NULL) {
        fclose(fp);
        if (remove_with_retry(slotpath)) {
            errlog(75, __LINE__);
        }
        dcache_log(DCACHE_DISCARD);
    }
    for (slot = CACHE_FILES_PER_HASHVAL; slot > 0; --slot) {
        if (make_cache_file_pathstr(
                old_slotpath,
                cache_dirpath,
                new_hdr->hashval,
                slot - 1)) {
            return 1;
        }
        if (make_cache_file_pathstr(
                slotpath,
                cache_dirpath,
                new_hdr->hashval,
                slot)) {
            return 1;
        }
        if ((fp = fopen(old_slotpath, "rb")) != NULL) {
            fclose(fp);
            if (rename_with_retry(old_slotpath, slotpath)) {
                /*
                 * If rename fails, log a warning, clear the
                 * lock, and return.
                 */
                errlog(56, 0);
                if (remove_with_retry(lockpath)) {
                    errlog(75, __LINE__); /* couldn't clear lock */
                }
                return 1;
            }
        }
    }
    return 0;
}   /* insert_as_mru() */


/***********************************************************
* static int is_cache_hit(
*         struct cache_header *requested,
*         struct cache_header *cached)
*
* Purpose:
*   Compares two cache_header structures, corresponding to
*   the requested data and the data in a cache file, to
*   determine if the file data constitutes a cache hit.
*
* Arguments:
*   struct cache_header *requested - structure containing
*                   parameters for requested data
*   struct cache_header *cached    - structure containing
*                   parameters for the data in a cache file
*
* Return:
*   0 if not a cache hit
*   1 if a cache hit
************************************************************/

static int is_cache_hit(
        struct cache_header *requested,
        struct cache_header *cached)
{
    double requested_imin, cached_imin;
    double requested_imax, cached_imax;
    double df = requested->df;

    if (
            (requested->version != cached->version) ||
            (requested->k_typenum != cached->k_typenum) ||
            (requested->lineshapenum != cached->lineshapenum) ||
            (fabs(requested->P - cached->P) >
             fabs(requested->P * DBL_EPSILON)) ||
            (fabs(requested->vmr - cached->vmr) >
             fabs(requested->vmr * DBL_EPSILON)) ||
            (fabs(requested->T - cached->T) >
             fabs(requested->T * DBL_EPSILON)) ||
            (fabs(df - cached->df) > fabs(df * DBL_EPSILON)) ||
            ((fabs(requested->tol) * (1. + DBL_EPSILON)) <
             fabs(cached->tol))
            ) {
        return 0;
    }
    /*
     * Compute the aligned frequency grid positions corresponding
     * to the requested and cached fmin and fmax.  If there were
     * unresolved lines when the cached absorption coefficient
     * was computed, then only a request for the full array
     * counts as a cache hit, since this is the only case for
     * which the unresolved line count can be guaranteed to be
     * correct.  (It is also the most common case in practice.)
     * If there were no unresolved lines, then any subarray
     * counts as a cache hit.
     */
    requested_imin = ceil( (requested->fmin * (1.0 - DBL_EPSILON)) / df);
    requested_imax = floor((requested->fmax * (1.0 + DBL_EPSILON)) / df);
    cached_imin    = ceil( (cached->fmin    * (1.0 - DBL_EPSILON)) / df);
    cached_imax    = floor((cached->fmax    * (1.0 + DBL_EPSILON)) / df);
    if (cached->unres_lines != 0) {
        if((requested_imin != cached_imin) ||
                (requested_imax != cached_imax)) {
            return 0;
        }
    } else {
        if ((requested_imin < cached_imin) ||
                (requested_imax > cached_imax)) {
            return 0;
        }
    }
    return 1;
}   /* is_cache_hit() */


/***********************************************************
* static int make_cache_file_pathstr(
*         char *cache_filepath,
*         const char *cache_dirpath,
*         const unsigned int hashval,
*         const int slot)
*
* Purpose:
*   Writes the full path to a cache file, to the buffer
*   named by cache_filepath.  The path is composed of the
*   cache directory path, and an appended file name
*   containing integer fields for the hash value and slot
*   number within the corresponding hash bucket.
*
* Arguments:
*   char *cache_filepath      - buffer to receive string
*   const char *cache_dirpath - path to cache directory
*   unsigned int hashval      - hash value
*   int slot                  - slot number in hash bucket
*
* Return:
*   0 on success, 1 on failure
************************************************************/

static int make_cache_file_pathstr(
        char *cache_filepath,
        const char *cache_dirpath,
        unsigned int hashval,
        int slot)
{
    unsigned int hmax;
    int hdigits;
    /*
     * In the hash file name, hash values are padded, with
     * leading zeros, to the same number of hex digits as the
     * maximum hash value.
     */
    hmax = hash_modulus() - 1;
    hdigits = 1;
    while (hmax >>= 4)
        ++hdigits;
    if (snprintf(cache_filepath,
            AM_MAX_PATH,
            "%sam_%0*x_%x",
            cache_dirpath,
            hdigits,
            hashval,
            slot) >= AM_MAX_PATH) {
        return 1;
    }
    return 0;
}   /* make_cache_file_pathstr() */


/***********************************************************
* static int promote_to_mru(
*         char *cache_dirpath,
*         unsigned int hashval,
*         int old_slot)
*
* Purpose:
*   Promotes the file in slot old_slot, in the cache bucket
*   associated with hashval, to the MRU slot.  To free up
*   the MRU slot, the files having slot numbers numbers from
*   1..old_slot-1 have their slot numbers incremented.
*
* Arguments:
*   char *cache_dirpath  - path string for cache directory
*   unsigned int hashval - hash value for the hash bucket
*                          being reordered
*   int old_slot         - slot number of file to be
*                          promoted to MRU
*
* Return:
*   0 if successful
*   1 otherwise (The cache bucket may have been locked.)
************************************************************/

static int promote_to_mru(
        char *cache_dirpath,
        unsigned int hashval,
        int old_slot)
{
    FILE *fp;
    char lockpath[AM_MAX_PATH];
    char slotpath[AM_MAX_PATH];
    char old_slotpath[AM_MAX_PATH];
    int slot;

    /*
     * If the file to be promoted is already in the MRU slot,
     * just return.
     */
    if (old_slot == 1)
        return 0;
    /*
     * Check for a lock file in slot 0.  If one exists, another
     * process is working with this hash bucket.  Give up on
     * promoting the file.
     */
    if (make_cache_file_pathstr(
            lockpath,
            cache_dirpath,
            hashval,
            0)) {
        return 1;
    }
    if ((fp = fopen(lockpath, "rb")) != NULL) {
        fclose(fp);
        return 1;
    }
    /*
     * Got past the lock file check, so try moving the file from
     * old_slot to the lock slot (slot 0). If the rename() fails,
     * another process could have won a race for the lock file.
     * If so, give up on promoting the file.
     */
    if (make_cache_file_pathstr(
                old_slotpath,
                cache_dirpath,
                hashval,
                old_slot)) {
        return 1;
    }
    if (rename(old_slotpath, lockpath))
        return 1;
    /*
     * Next, increment the slot number of the files in slots
     * old_slot - 1 to 0.  The last one to be incremented is the
     * promoted file, in slot 0, which gets put into slot 1,
     * clearing the lock.  To accommodate the possibility that a
     * file has been deleted, test for existence before calling
     * rename_with_retry().
     */
    for (slot = old_slot - 1; slot >= 0; --slot) {
        if (make_cache_file_pathstr(
                    old_slotpath,
                    cache_dirpath,
                    hashval,
                    slot)) {
            return 1;
        }
        if (make_cache_file_pathstr(
                    slotpath,
                    cache_dirpath,
                    hashval,
                    slot + 1)) {
            return 1;
        }
        if ((fp = fopen(old_slotpath, "rb")) != NULL) {
            fclose(fp);
            if (rename_with_retry(old_slotpath, slotpath)) {
                /*
                 * If rename fails, log a warning, clear the
                 * lock, and return.
                 */
                errlog(55, 0);
                if (remove_with_retry(lockpath)) {
                    errlog(75, __LINE__); /* couldn't clear lock */
                }
                return 1;
            }
        }
    }
    return 0;
}   /* promote_to_mru() */


/***********************************************************
* static int read_cache_file_data(
*         model_t *model,
*         double *k,
*         int *unres_lines,
*         struct cache_header *hdr,
*         FILE *fp)
*
* Purpose:
*   Reads the requested data from the cache file containing
*   it into the absorption coefficient array pointed to by
*   k.  This function must follow a call to
*   read_cache_file_hdr(), which will move the current file
*   position to the end of the header area.
*
*   The unresolved line count associated with the absorption
*   coefficent computation is restored along with the
*   absorption coefficient.  Note that cache entries with
*   nonzero unresolved line counts are only counted as cache
*   hits if the entire k array is requested, since it is not
*   possible to restore a correct unresolved line count on a
*   subarray.
*
* Arguments:
*   model_t model            - model data structure
*   double *k                - pointer to array to receive
*                              absorption coefficient data
*   int *unres_lines         - unresolved line count for
*                              cached array
*   struct cache_header *hdr - structure containing cache
*                              file header data
*   FILE *fp                 - file pointer for cache file
*
* Return:
*   0 if successful
*   1 otherwise
************************************************************/

static int read_cache_file_data(
        model_t *model,
        double *k,
        int *unres_lines,
        struct cache_header *hdr,
        FILE *fp)
{
    gridsize_t nskipped;

    /*
     * The requested data array can be a subarray of the cached
     * data array (i.e. same df, but narrower range from fmin to
     * fmax).  Compute the number of points at the beginning of
     * the cached array which will be skipped.  Frequency grids
     * in am are always aligned to the origin.
     */
    nskipped =
        (gridsize_t)(ceil(model->fmin * (1.0 - DBL_EPSILON) / model->df)
                - ceil(hdr->fmin * (1.0 - DBL_EPSILON) / hdr->df));
    if (fseek(fp, nskipped * (long int)sizeof(double), SEEK_CUR))
        return 1;
    *unres_lines = hdr->unres_lines;
    if (fread((void*)k, sizeof(double), (size_t)model->ngrid, fp)
            != (size_t)model->ngrid) {
        return 1;
    }
    return 0;
}   /* read_cache_file_data() */


/***********************************************************
* static int read_cache_file_hdr(
*         struct cache_header *hdr
*         FILE *fp)
*
* Purpose:
*   Reads the header data from a cache file into a
*   cache_header structure.  This function leaves the file
*   position at the next byte after the header data area.
*
* Arguments:
*   struct cache_header *hdr - structure to receive header
*                              data
*   FILE *fp                 - file pointer for cache file
*
* Return:
*   0 if successful
*   1 otherwise (bad fp, or the cache file was corrupt or
*     binary-incompatable.)
************************************************************/

static int read_cache_file_hdr(
        struct cache_header *hdr,
        FILE *fp)
{
    unsigned int hashval;
    char id_string[sizeof(DCACHE_FILE_ID_STRING)];

    rewind(fp);
    if (fread(
            (void*)id_string,
            sizeof(char),
            sizeof(DCACHE_FILE_ID_STRING),
            fp) != sizeof(DCACHE_FILE_ID_STRING)) {
        return 1;
    }
    if (strncmp(
        id_string, DCACHE_FILE_ID_STRING,
        sizeof(DCACHE_FILE_ID_STRING)) != 0) {
        errlog(69, 0);
        return 1;
    }
    /*
     * Check the version first, in case the header data has
     * changed between versions.
     */
    if(fread((void*)&hdr->version, sizeof(unsigned int), (size_t)1, fp) != 1) {
        errlog(69, 0);
        return 1;
    }
    if (hdr->version != dcache_version()) {
        /*
         * Fail quietly.  Old files will eventually be evicted
         * from the cache as new ones are added.
         */
        return 1;
    }
    if (
            (fread((void*)&hdr->hashval, sizeof(unsigned int),
                   (size_t)1, fp) != 1) ||
            (fread((void*)&hdr->k_typenum,    sizeof(int),
                   (size_t)1, fp) != 1) ||
            (fread((void*)&hdr->lineshapenum, sizeof(int),
                   (size_t)1, fp) != 1) ||
            (fread((void*)&hdr->unres_lines,  sizeof(int),
                   (size_t)1, fp) != 1) ||
            (fread((void*)&hdr->ngrid, sizeof(gridsize_t),
                   (size_t)1, fp) != 1) ||
            (fread((void*)&hdr->fmin, sizeof(double),
                   (size_t)1, fp) != 1) ||
            (fread((void*)&hdr->fmax, sizeof(double),
                   (size_t)1, fp) != 1) ||
            (fread((void*)&hdr->df,   sizeof(double),
                   (size_t)1, fp) != 1) ||
            (fread((void*)&hdr->tol,  sizeof(double),
                   (size_t)1, fp) != 1) ||
            (fread((void*)&hdr->P,    sizeof(double),
                   (size_t)1, fp) != 1) ||
            (fread((void*)&hdr->vmr,  sizeof(double),
                   (size_t)1, fp) != 1) ||
            (fread((void*)&hdr->T,    sizeof(double),
                   (size_t)1, fp) != 1)
        ) {
        errlog(69, 0);
        return 1;
    }
    /*
     * To guard against corrupt or binary-incompatable cache
     * files, recompute the hash value and check against the
     * stored value.
     */
    hashval = hash(
            hdr->df,
            hdr->P,
            hdr->vmr,
            hdr->T,
            hdr->lineshapenum,
            hdr->k_typenum,
            hash_modulus());
    if (hashval != hdr->hashval) {
        errlog(69, 0);
        return 1;
    }
    return 0;
}   /* read_cache_file_hdr() */


/***********************************************************
* void report_dcache_env_info(FILE *stream)
*
* Purpose:
*   Reports environment information pertaining to the disk
*   cache to a stream.
*
* Arguments:
*   FILE *stream - destination for output
************************************************************/

void report_dcache_env_info(FILE *stream)
{
    fprintf(stream, "Disk cache\n");
    if (dcache_path() == NULL) {
        fprintf(stream,
            "  %s is not set, so am is running without a disk cache.\n",
            DCACHE_DIRPATH_ENV_STRING);
    } else {
        fprintf(stream, "  %s = %s", DCACHE_DIRPATH_ENV_STRING, dcache_path());
        if (dcache_path_valid())
            fprintf(stream, "\n");
        else
            fprintf(stream, " ! Warning: not a valid directory.\n");
    }
    fprintf(stream, "  %s = %d", DCACHE_HASHMOD_ENV_STRING, hash_modulus());
    if (getenv(DCACHE_HASHMOD_ENV_STRING) == NULL)
        fprintf(stream, " (default setting)\n");
    else
        fprintf(stream, "\n");
    return; 
}   /* report_dcache_env_info() */


/**********************************************************
* void report_dcache_log_data(FILE *stream)
*
* Purpose:
*   Reports dcache log data to a stream.
*
* Arguments:
*   FILE *stream - destination for output
************************************************************/

void report_dcache_log_data(FILE *stream)
{
    if (dcache_path() != NULL)
        fprintf(stream,
            "# dcache hit: %lu  miss: %lu  discard: %lu\n",
            dcache_hit_count, dcache_miss_count, dcache_discard_count);
    return;
}   /* report_dcache_log_data() */


/***********************************************************
* static int write_cache_file(
*         struct cache_header *hdr,
*         double *k,
*         FILE *fp)
*
* Purpose:
*   Writes header data and absorption coefficient data for
*   column cnum of layer lnum of the model to the stream fp.
*
* Arguments:
*   struct cache_header hdr - file header data
*   double *k               - pointer to array containing
*                             the spectral absorption
*                             coefficient to be written
*   FILE *fp                - stream to write onto
*
* Return:
*   0 if successful
*   1 otherwise
************************************************************/

static int write_cache_file(
        struct cache_header *hdr,
        double *k,
        FILE *fp)
{
    /*
     * File ID string.
     */
    if (fwrite(
                (void*)DCACHE_FILE_ID_STRING,
                sizeof(char),
                sizeof(DCACHE_FILE_ID_STRING),
                fp) != sizeof(DCACHE_FILE_ID_STRING)) {
        return 1;
    }
    /*
     * Header items.  These get written out member by member to
     * keep the struct alignment from affecting the file layout.
     * Cache files should be binary compatable within a given
     * machine architecture and integer size.
     */
    if (
            (fwrite((void*)&hdr->version, sizeof(unsigned int),
                    (size_t)1, fp) != 1) ||
            (fwrite((void*)&hdr->hashval, sizeof(unsigned int),
                    (size_t)1, fp) != 1) ||
            (fwrite((void*)&hdr->k_typenum,    sizeof(int),
                    (size_t)1, fp) != 1) ||
            (fwrite((void*)&hdr->lineshapenum, sizeof(int),
                    (size_t)1, fp) != 1) ||
            (fwrite((void*)&hdr->unres_lines,  sizeof(int),
                    (size_t)1, fp) != 1) ||
            (fwrite((void*)&hdr->ngrid, sizeof(gridsize_t),
                    (size_t)1, fp) != 1) ||
            (fwrite((void*)&hdr->fmin, sizeof(double),
                    (size_t)1, fp) != 1) ||
            (fwrite((void*)&hdr->fmax, sizeof(double),
                    (size_t)1, fp) != 1) ||
            (fwrite((void*)&hdr->df,   sizeof(double),
                    (size_t)1, fp) != 1) ||
            (fwrite((void*)&hdr->tol,  sizeof(double),
                    (size_t)1, fp) != 1) ||
            (fwrite((void*)&hdr->P,    sizeof(double),
                    (size_t)1, fp) != 1) ||
            (fwrite((void*)&hdr->vmr,  sizeof(double),
                    (size_t)1, fp) != 1) ||
            (fwrite((void*)&hdr->T,    sizeof(double),
                    (size_t)1, fp) != 1)
       ) {
        return 1;
    }
    /*
     * The absorption coefficient array.
     */
    if (fwrite(
                (void*)k,
                sizeof(double),
                (size_t)hdr->ngrid,
                fp) != (size_t)hdr->ngrid) {
        return 1;
    }
    fflush(fp);
    return 0;
}   /* write_cache_file() */
