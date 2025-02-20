/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* fileops.c                       S. Paine rev. 2024 July 21
*
* Miscellaneous file operations.
************************************************************/

#include <errno.h>
#include <signal.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "am_sysdep.h"
#include "fileops.h"

/*
 * Format string for am temporary file names
 */
static const char TMPNAME_FMT[] = "%sam_tmp_%04x";

/*
 * Initial retry delays for failed file operations.  Retry delays
 * double on each attempt, up to the maximum number of retries.
 */
static const double REMOVE_WAIT_INIT = 0.125;   /* seconds */
static const double RENAME_WAIT_INIT = 0.125;   /* seconds */

/*
 * Maximum number of retries for failed file operations.
 */
enum {
    TMPFILE_NUM_TRIES  = 8,
    REMOVE_NUM_RETRIES = 8,
    RENAME_NUM_RETRIES = 8
};

/*
 * recognized directory separators:
 * '/' - Unix, POSIX, OS X.  Also works in Windows standard C library.
 * '\\' - Windows native
 * ':' - Older MacOS
 * ']' - VMS
 */
static const char DIR_SEPARATORS[] = "/\\:]";
static const char DEFAULT_DIR_SEPARATOR[] = "/";


/***********************************************************
* FILE *am_tmpfile(const char *dirpath, char *fpath)
*
* Purpose:
*   Creates a temporary file, opened in binary update mode,
*   and returns a file pointer.  The file is created in the
*   directory named by the string dirpath.  If dirpath is
*   NULL, or points to an empty string, the file is created
*   in the current directory.  The full path to the
*   temporary file is written to the buffer fpath, which
*   should be at least AM_MAX_PATH characters long.
*
*   If the file cannot be created, NULL is returned, and an
*   empty string is written to fpath.
*
*   The last character of dirpath must be a directory
*   separator appropriate for the current environment.
*
* Arguments:
*   const char *path - path to directory
*   char *fpath      - pointer to buffer to receive path to
*                      file.
*
* Return:
*   file pointer, or NULL if the file could not be created.
************************************************************/

FILE *am_tmpfile(const char *dirpath, char *fpath)
{
    int i;
    FILE *stream;

    for (i = 0; i < TMPFILE_NUM_TRIES; ++i) {
        /*
         * Generate a random file name.
         */
        snprintf(fpath, AM_MAX_PATH,
                TMPNAME_FMT, ((dirpath == NULL) ? "" : dirpath), rand());
        if ((stream = fopen(fpath, "rb")) != NULL) {
            /*
             * File already exists.  Close and try another name.
             */
            fclose(stream);
        } else {
            /*
             * Didn't exist already, so try opening in binary
             * update mode.  If, by chance, this process lost a
             * race for this file name, the outcome depends on
             * the implementation of fopen().  If "wb+" opens
             * with exclusive access, then this fopen() fails and
             * all is well.  However, another possibility is that
             * this fopen() causes the other process to lose its
             * file name to file pointer association.  In am, a
             * possible consequence of this would be a file
             * getting copied to the wrong hash bucket of the
             * disk cache. This situation resolves itself--the
             * misplaced file will always be a cache miss based
             * on its header data, and it will eventually be
             * evicted from the cache.
             */
            if ((stream = fopen(fpath, "wb+")) != NULL)
                return stream;
        }
    }
    /*
     * Too many tries.  Give up and return error status.
     */
    fpath[0] = '\0';
    return NULL;
}   /* am_tmpfile() */


/***********************************************************
* int check_for_dir_separator(char *dirpath)
*
* Purpose:
*   Checks that dirpath ends in a directory separator
*   character.  If not, a default one is appended, unless
*   the resulting string would be longer than
*   AM_MAX_DIRPATH.
*
* Arguments:
*   char *dirpath - path to directory
*
* Return:
*   0 on success
*   1 otherwise
************************************************************/

int check_for_dir_separator(char *dirpath)
{
    size_t length;

    if (dirpath == NULL)
        return 1; /* NULL path doesn't exist */
    if ((length = strlen(dirpath)) == 0)
        return 0; /* empty string means use the current working dir */
    /*
     * Look for a valid directory separator character at the
     * end of dirpath.  If one isn't found, append a default.
     */
    if (strchr(DIR_SEPARATORS, dirpath[length - 1]) == NULL) {
        if (length + sizeof(DEFAULT_DIR_SEPARATOR) < AM_MAX_DIRPATH) {
            strncat(dirpath,
                    DEFAULT_DIR_SEPARATOR,
                    AM_MAX_DIRPATH - length);
        } else {
            return 1;
        }
    }
    return 0;
}   /* check_for_dir_separator() */


/***********************************************************
* int remove_with_retry(const char *filename)
*
* Purpose:
*   This function is a wrapper around the standard C library
*   function remove().  It calls remove(), and if the call
*   fails keeps trying with exponential back-off for up
*   to REMOVE_NUM_RETRIES attempts.
*
* Arguments:
*   const char *filename - file name or path string
*
* Return:
*   0 if the remove() call eventually succeeds
*   1 if the timeout expires
************************************************************/

int remove_with_retry(const char *filename)
{
    double t_wait;
    int i;

    if (remove(filename) == 0)
        return 0;
    t_wait = REMOVE_WAIT_INIT;
    for (i = 0; i < REMOVE_NUM_RETRIES; ++i) {
        am_sleep((unsigned int)t_wait + 1);
        errno = 0;
        if (remove(filename) == 0)
            return 0;
        /*
         * See if the file was already removed.  Note that
         * ENOENT is defined by POSIX, not by ANSI C.
         */
        if (errno == ENOENT)
            return 0;
        t_wait *= 2.;
    }
    return 1;
}   /* remove_with_retry() */


/***********************************************************
* int rename_with_retry(const char *oldname, const char *newname)
*
* Purpose:
*   This function is a wrapper around the standard C library
*   function rename().  It calls rename(), and, if the call
*   fails, keeps trying with exponential back-off for up
*   to RENAME_NUM_RETRIES attempts.
*
* Arguments:
*   const char *oldname - old file name or path string
*   const char *newname - new file name or path string
*
* Return:
*   0 if the rename() call eventually succeeds
*   1 if the timeout expires
************************************************************/

int rename_with_retry(const char *oldname, const char *newname)
{
    double t_wait;
    int i;

    if (rename(oldname, newname) == 0)
        return 0;
    t_wait = RENAME_WAIT_INIT;
    for (i = 0; i < RENAME_NUM_RETRIES; ++i) {
        am_sleep((unsigned int)t_wait + 1);
        if (rename(oldname, newname) == 0)
            return 0;
        t_wait *= 2.;
    }
    return 1;
}   /* rename_with_retry() */
