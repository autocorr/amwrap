/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* kcache.c                     S. Paine rev. 2021 January 26
*
* Memory management functions for the kcache.
************************************************************/

#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#include "am_types.h"
#include "errlog.h"
#include "kcache.h"

#ifndef SIZE_MAX
    #define SIZE_MAX ((size_t) (-1))
#endif
/*
 * KCACHE_INIT_BLOCK_SIZE is the number of kcache entries in the
 * first block of memory allocated for the kcache.
 *
 * KCACHE_BLOCKSIZE_GROW is the factor by which successively
 * allocated blocks of kcache memory grow in size.
 *
 * KCACHE_MAX_NUM_BLOCKS is the maximum number of blocks which
 * will ever be allocated, and is the size of the static block
 * table.
 *
 * 1 / KCACHE_DISCARD_FRAC is the fraction of kcache entries
 * which will be deleted from the kcache to make room for new
 * entries if memory runs out, or if a user-defined memory limit
 * defined in the environment is reached.
 */

enum {
    KCACHE_INIT_BLOCK_SIZE =  8,
    KCACHE_BLOCKSIZE_GROW  =  2,
    KCACHE_MAX_NUM_BLOCKS  = 32,
    KCACHE_DISCARD_FRAC    =  4
};

static const double KCACHE_ENTRY_FREE   = 0.0;
static const double KCACHE_ENTRY_IN_USE = 1.0;
static const double KCACHE_ENTRY_LOCKED = 2.0;
static const char KCACHE_MEM_LIMIT_ENV_STRING[] = "AM_KCACHE_MEM_LIMIT";

/*
 * kcache_block_table[] is a table of pointers to blocks obtained
 * with malloc().  **kcache_pointer_table is a table of pointers
 * into these blocks, one for each kcache entry.  kcache entries
 * are arrays, of length equal to the model frequency grid,
 * containing an absorption coefficient.  The pointers in the
 * table point to the zeroth element of each kcache array.  An
 * extra element in front of the array is allocated to track
 * whether an array is free, in use, or locked to prevent it
 * being discarded.  
 */
static double      *kcache_block_table[KCACHE_MAX_NUM_BLOCKS];
static double     **kcache_pointer_table    = NULL;
static size_t       kcache_pointer_tabsize  = 0;
static size_t       kcache_num_free_entries = 0;
static unsigned int kcache_num_blocks       = 0;
static size_t       kcache_num_bytes        = 0;

/*
 * Log data
 */
static unsigned long kcache_hit_count     = 0;
static unsigned long kcache_miss_count    = 0;
static unsigned long kcache_discard_count = 0;

static size_t grow_kcache(model_t*);
static size_t discard_kcache_entries(model_t*, size_t);
static size_t kcache_mem_limit(void);


/***********************************************************
* size_t discard_kcache_entries(model_t* model, size_t last_alloc)
*
* Purpose:
*   This function is called when no more memory is available
*   for kcache entries, either because a malloc() request
*   for an additional memory block has failed, or because
*   the user-defined limit set by the environment variable
*   KCACHE_MEM_LIMIT has been reached.
*
*   To make space available for new kcache entries, this
*   function discards a fraction 1/KCACHE_DISCARD_FRAC of
*   the existing entries.  The discarded entries are taken
*   sequentially from the kcache pointer table, starting
*   from the position just after the last
*   successfully-allocated entry and wrapping around the end
*   of the table if needed.  After marking the entries as
*   free, the corresponding pointers in the model are found
*   and set to NULL.
*
* Arguments:
*   model_t *model - pointer to model structure
*   last_alloc     - index in the pointer table of the
*                    most recently allocated kcache entry
*
* Return:
*   number of discarded entries
************************************************************/

size_t discard_kcache_entries(model_t* model, size_t last_alloc)
{
    size_t i, ndiscard_try, ndiscard;
    int lnum;
    /*
     * Mark the entries to be freed.  Skip an entry if it is
     * locked.
     */
    ndiscard_try = kcache_pointer_tabsize / KCACHE_DISCARD_FRAC;
    ndiscard = 0;
    for (i = 1; i <= ndiscard_try; ++i) {
        size_t imod = (last_alloc + i) % kcache_pointer_tabsize;
        if (kcache_pointer_table[imod][-1] != KCACHE_ENTRY_LOCKED) {
            kcache_pointer_table[imod][-1] = KCACHE_ENTRY_FREE;
            kcache_log(KCACHE_DISCARD);
            ++ndiscard;
        }
    }
    /*
     * Iterate through the model, finding the freed entries and
     * setting their pointers to NULL.
     */
    for (lnum = 0; lnum < model->nlayers; ++lnum) {
        layer_t *layer = model->layer[lnum];
        int cnum;
        for (cnum = 0; cnum < layer->ncols; ++cnum) {
            column_t *column = layer->column[cnum];
            int knum;
            for (knum = 0; knum < column->n_abscoeffs; ++knum) {
                abscoeff_t *abscoeff = column->abscoeff[knum];
                int j;
                if (abscoeff->kcache == NULL)
                    continue;
                for (j = 0; j < model->nkcache; ++j) {
                    if ((abscoeff->kcache[j] != NULL) &&
                        (abscoeff->kcache[j][-1] == KCACHE_ENTRY_FREE))
                            abscoeff->kcache[j] = NULL;
                }
            }
        }
    }
    kcache_num_free_entries += ndiscard;
    return ndiscard;
}   /* discard_kcache_entries() */


/***********************************************************
* static size_t grow_kcache(model_t *model)
*
* Purpose:
*   Grow the kcache by requesting a block of memory with
*   malloc, larger than the last block by a factor
*   KCACHE_BLOCKSIZE_GROW, and enlarging the kcache pointer
*   table with entries pointing into the new block.
*
* Arguments:
*   model_t *model - pointer to model structure
*
* Return:
*   Number of entries in new block.
************************************************************/

static size_t grow_kcache(model_t *model)
{
    static size_t last_blocksize = 0;
    size_t new_blocksize, max_blocksize;
    size_t new_block_bytes, max_block_bytes;
    double **tptr;
    size_t i, last_tabsize;

    /*
     * Make sure there is room for another block address in the
     * block table
     */
    if (kcache_num_blocks >= KCACHE_MAX_NUM_BLOCKS)
        return 0;
    /*
     * Compute the maximum size of a new block.
     */
    if ((max_block_bytes = kcache_mem_limit()) == 0)
        max_block_bytes = SIZE_MAX;
    max_block_bytes -= kcache_num_bytes;
    max_blocksize = max_block_bytes / ((1 + model->ngrid) * sizeof(double));
    if (max_blocksize == 0)
        return 0;
    /*
     * If this is the initial block, the size is
     * KCACHE_INIT_BLOCK_SIZE.  Otherwise, try to grow the block
     * size by KCACHE_BLOCKSIZE_GROW, subject to the maximum new
     * block size computed above.
     */
    if (kcache_num_blocks == 0) {
        if (max_blocksize >= KCACHE_INIT_BLOCK_SIZE)
            new_blocksize = KCACHE_INIT_BLOCK_SIZE;
        else
            return 0;
    } else {
        if (last_blocksize > (max_blocksize / KCACHE_BLOCKSIZE_GROW))
            new_blocksize = max_blocksize;
        else
            new_blocksize = last_blocksize * KCACHE_BLOCKSIZE_GROW;
    }
    /*
     * malloc() the new block.  If this fails, return 0 as the
     * number of new entries.
     */
    new_block_bytes = new_blocksize * (1 + model->ngrid) * sizeof(double);
    if ((kcache_block_table[kcache_num_blocks] =
        (double*)malloc(new_block_bytes)) == NULL) {
        return 0;
    }
    /*
     * realloc() the pointer table to hold the pointers into the
     * new block.  If this fails, free the just-allocated block
     * and return 0.
     */
    if ((tptr = (double**)realloc(kcache_pointer_table,
        (kcache_pointer_tabsize + new_blocksize) * sizeof(double*))) == NULL) {
        free(kcache_block_table[kcache_num_blocks]);
        return 0;
    } else {
        kcache_pointer_table = tptr;
    }
    last_tabsize = kcache_pointer_tabsize;
    kcache_pointer_tabsize  += new_blocksize;
    kcache_num_free_entries += new_blocksize;
    /*
     * Put the addresses of the new arrays into the kcache table,
     * and mark the arrays as free.
     */
    kcache_pointer_table[last_tabsize] =
        1 + kcache_block_table[kcache_num_blocks];
    for (i = last_tabsize + 1; i < kcache_pointer_tabsize; ++i)
        kcache_pointer_table[i] = kcache_pointer_table[i-1] + model->ngrid + 1;
    for (i = last_tabsize; i < kcache_pointer_tabsize; ++i)
        kcache_pointer_table[i][-1] = KCACHE_ENTRY_FREE;
    /*
     * Update the counts of blocks and bytes, and return the new
     * block size.
     */
    ++kcache_num_blocks;
    kcache_num_bytes += new_block_bytes;
    last_blocksize    = new_blocksize;
    return new_blocksize;
}   /* grow_kcache() */


/***********************************************************
* double *kcache_alloc(model_t* model)
*
* Purpose:
*   Returns a pointer to space for a kcache entry.  If
*   sufficient space is not available in memory, it is
*   obtained by discarding a fraction of the existing kcache
*   entries, and setting their pointers in the model
*   structure to NULL.  Consequently, kcache pointers should
*   always be checked before use after any calls to
*   kcache_alloc().  This is a decidedly single-threaded
*   design; because of this, the kcache cannot be accessed
*   in a parallel block.
*
* Arguments:
*   model_t *model - pointer to model structure
*
* Return:
*   Pointer to free space, or NULL on failure
************************************************************/

double *kcache_alloc(model_t* model)
{
    static size_t last_alloc = KCACHE_INIT_BLOCK_SIZE - 1;
    /*
     * If there are no free entries, try to get a new block.
     */
    if (kcache_num_free_entries == 0) {
        if (grow_kcache(model) == 0) {
            if (kcache_num_blocks == 0) {
                /*
                 * Failed on first block.  Log an error and
                 * return NULL.
                 */
                errlog(114, 0);
                return NULL;
            } else {
                /*
                 * Ran out of memory.  Free space by discarding
                 * older kcache entries.
                 */
                discard_kcache_entries(model, last_alloc);
            }
        }
    }
    /*
     * The next free space for a kcache entry is found by a
     * cyclical linear search of the pointer table starting with
     * the entry just beyond the last to have been allocated.
     * Most likely, the first location checked will be free,
     * unless the only remaining free space consists of random
     * holes left by calls to kcache_free().
     */
    do {
        ++last_alloc;
        last_alloc = last_alloc % kcache_pointer_tabsize;
    } while (kcache_pointer_table[last_alloc][-1] != KCACHE_ENTRY_FREE);
    /*
     * Once a free space has been found, mark it in use, and
     * return the array address from the pointer table.
     */
    kcache_pointer_table[last_alloc][-1] = KCACHE_ENTRY_IN_USE;
    --kcache_num_free_entries;
    return kcache_pointer_table[last_alloc];
}   /* kcache_alloc() */


/***********************************************************
* void kcache_free(double *ptr)
*
* Purpose:
*   Frees an individual entry in the kcache.  The entry will
*   be freed even if it is marked as locked.
*
* Arguments:
*   double *ptr - pointer to a kcache entry.
************************************************************/

void kcache_free(double *ptr)
{
    if (ptr != NULL)
        ptr[-1] = KCACHE_ENTRY_FREE;
    return;
}   /* kcache_free() */


/***********************************************************
* void kcache_free_all(void)
*
* Purpose:
*   Frees all blocks which were allocated to the kcache, and
*   frees the pointer table.
************************************************************/

void kcache_free_all(void)
{
    while (kcache_num_blocks > 0)
        free(kcache_block_table[--kcache_num_blocks]);
    free(kcache_pointer_table);
    kcache_pointer_tabsize = 0;
    kcache_num_free_entries = 0;
    kcache_num_blocks = 0;
    kcache_num_bytes = 0;
    return;
}   /* kcache_free_all() */


/***********************************************************
* void kcache_lock_entry(double *ptr)
*
* Purpose:
*   Locks a kcache entry, to prevent it being discarded if
*   discard_kcache_entries() is called to free memory for
*   new entries.  The entry will only be locked if it is
*   in use.  If the entry is not in use, or if ptr == NULL,
*   this function does nothing.
*
* Arguments:
*   double *ptr - pointer to a kcache entry.
************************************************************/

void kcache_lock_entry(double *ptr)
{
    if ((ptr != NULL) && (ptr[-1] == KCACHE_ENTRY_IN_USE))
        ptr[-1] = KCACHE_ENTRY_LOCKED;
    return;
}   /* kcache_lock_entry() */


/***********************************************************
* void kcache_log(int event)
*
* Purpose:
*   Keeps a count of kcache events, including hits, misses
*   and discards.
*
* Arguments:
*   int event - event type (defined in kcache.h)
************************************************************/

void kcache_log(int event)
{
    switch (event) {
    case KCACHE_HIT:
        ++kcache_hit_count;
        break;
    case KCACHE_MISS:
        ++kcache_miss_count;
        break;
    case KCACHE_DISCARD:
        ++kcache_discard_count;
        break;
    case KCACHE_CLEAR_LOG:
        kcache_hit_count = 0;
        kcache_miss_count = 0;
        kcache_discard_count = 0;
        break;
    default:
        break;
    }
    return;
}   /* kcache_log() */


/***********************************************************
* static size_t kcache_mem_limit(void)
*
* Purpose:
*   On the first call, checks whether a memory size limit
*   for the kcache has been specified in the environment.
*   If so, the value is stored in an internal static
*   variable and returned.  If not, 0 is assigned, meaning
*   there is no defined limit.  In either case, subsequent
*   calls return the same value without checking the
*   environment again.
*
* Arguments:
*   None
*
* Return:
*   memory limit in bytes, as a size_t
************************************************************/

static size_t kcache_mem_limit(void)
{
    static unsigned long mem_limit = 0;
    static int initialized = 0;
    char *str;
    char *endp;

    if (initialized)
        return (size_t)mem_limit;
    if ((str = getenv(KCACHE_MEM_LIMIT_ENV_STRING)) != NULL) {
        errno = 0;
        mem_limit = strtoul(str, &endp, 0);
        if ((endp[0] != '\0') ||
            ((mem_limit == ULONG_MAX) && (errno == ERANGE))) {
            errlog(113,0);
            mem_limit = 0;
        }
    }
    if (mem_limit > SIZE_MAX)
        mem_limit = SIZE_MAX;
    initialized = 1;
    return (size_t)mem_limit;
}   /* kcache_mem_limit() */
    

/***********************************************************
* void kcache_unlock_entry(double *ptr)
*
* Purpose:
*   Unlocks a kcache entry, and marks it in use.  (Only
*   entries which were in use could have been locked.)  If
*   the entry was not locked, or if ptr == NULL, this
*   function does nothing.
*
* Arguments:
*   double *ptr - pointer to a kcache entry.
************************************************************/

void kcache_unlock_entry(double *ptr)
{
    if ((ptr != NULL) && (ptr[-1] == KCACHE_ENTRY_LOCKED))
        ptr[-1] = KCACHE_ENTRY_IN_USE;
    return;
}   /* kcache_unlock_entry() */


/**********************************************************
* void report_kcache_env_info(FILE *stream)
*
* Purpose:
*   Reports environment information pertaining to the kcache
*   to a stream.
*
* Arguments:
*   FILE *stream - destination for output
************************************************************/

void report_kcache_env_info(FILE *stream)
{
    size_t mem_limit = kcache_mem_limit();
    fprintf(stream, "kcache\n");
    if (mem_limit == 0) {
        fprintf(stream, "  %s = 0 (use maximum available memory)\n",
            KCACHE_MEM_LIMIT_ENV_STRING);
    } else {
        fprintf(stream, "  %s = %lu bytes",
            KCACHE_MEM_LIMIT_ENV_STRING, (unsigned long)mem_limit);
        if ((mem_limit >> 10) == 0)
            fprintf(stream, "\n");
        else if ((mem_limit >> 20) == 0)
            fprintf(stream, " (%.1f KB)\n", (double)mem_limit / (1 << 10));
        else if ((mem_limit >> 30) == 0)
            fprintf(stream, " (%.1f MB)\n", (double)mem_limit / (1 << 20));
        else
            fprintf(stream, " (%.1f GB)\n", (double)mem_limit / (1 << 30));
    }
    return;
}   /* report_kcache_env_info() */


/**********************************************************
* void report_kcache_log_data(FILE *stream)
*
* Purpose:
*   Reports kcache log data to a stream.
*
* Arguments:
*   FILE *stream - destination for output
************************************************************/

void report_kcache_log_data(FILE *stream)
{
    fprintf(stream,
        "# kcache hit: %lu  miss: %lu  discard: %lu\n",
        kcache_hit_count, kcache_miss_count, kcache_discard_count);
    return;
}   /* report_kcache_log_data() */
