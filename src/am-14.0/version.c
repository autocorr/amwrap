/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* version.c                   S. Paine rev. 2021 December 23
*
* Functions reporting program version and build date.
*
* In the Makefile, version.o depends upon this file and all
* other object files, ensuring that the build date and time
* are always updated, even for partial builds.
************************************************************/

#include <limits.h>
#include <stdio.h>

#include "version.h"


/***********************************************************
* void version(FILE* stream)
*
* Purpose:
*   Prints version information to a stream.
************************************************************/

void version(FILE* stream)
{
    fprintf(stream, "am version %d.%d\n", AM_VERSION_MAJOR, AM_VERSION_MINOR);
    fprintf(stream, "build date %s %s\n", __DATE__, __TIME__);
#if _WIN64
    fprintf(stream, "64-bit ");
#elif (ULONG_MAX == 0xFFFFFFFFFFFFFFFF)
    fprintf(stream, "64-bit ");
#elif (ULONG_MAX == 0xFFFFFFFF)
    fprintf(stream, "32-bit ");
#endif
#ifdef _OPENMP
    fprintf(stream, "multi-threaded (OpenMP version %d)\n\n", _OPENMP);
#else
    fprintf(stream, "single-threaded\n\n");
#endif
    return;
}   /* version() */


/***********************************************************
* void write_version_line(FILE *stream)
*
* Purpose:
*   Writes a single line with the program version and build date
*   to a stream.
************************************************************/

void write_version_line(FILE *stream)
{
    fprintf(stream,
        "# am version %d.%d (build date %s %s)\n"
        "\n",
        AM_VERSION_MAJOR,
        AM_VERSION_MINOR,
        __DATE__,
        __TIME__
    );
    return;
}   /* write_version_line() */
