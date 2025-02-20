/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* mapping.c                  S. Paine rev. 2016 September 26
*
* Certain model parameter values and their differentials
* are mapped from their physical domain to the positive
* real domain.  This facilitates their representation in
* a log simplex space during fits.  These functions do the
* forward and reverse mappings.
************************************************************/

#include <float.h>
#include <stdio.h>
#include <math.h>

#include "errlog.h"
#include "mapping.h"

/*
 * In the functions below, x is an unmapped variable in the
 * physical domain, and y is the corresponding mapped value.
 * The possible mappings are:
 * 
 * MAPPING_NONE
 *    y = x
 *   dy = dx
 *
 * MAPPING_VMR
 *   maps volume mixing ratios x in [0, 1] to y in [0, 1 / FLT_EPSILON]
 *    y = x / [(1 + FLT_EPSILON) - x)
 *   dy = [dx * (1 + FLT_EPSILON)] / [(1 + FLT_EPSILON) - x]^2
 *
 * MAPPING_EXP
 *   maps a (typically small) signed x in [-inf,inf] to y in [0,inf]
 *    y = exp(x)
 *   dy = exp(x) * dx
 */


/***********************************************************
* double map_differential(double x, double dx, int mapping)
*
* Purpose:
*   Given x and dx, returns dy
************************************************************/

double map_differential(double x, double dx, int mapping)
{
    double dy;
    volatile double r1;

    switch (mapping) {
    case MAPPING_NONE:
        dy = dx;
        break;
    case MAPPING_VMR:
        x = x > 1.0 ? 1.0 : (x < 0.0 ? 0.0 : x);
        dx = dx > 1.0 ? 1.0 : (dx < 0.0 ? 0.0 : dx);
        r1 = 1.0 - x;
        r1 += FLT_EPSILON;
        r1 *= r1;
        dy = dx * (1.0 + FLT_EPSILON) / r1;
        break;
    case MAPPING_EXP:
        dy = dx * exp(x);
        break;
    default:
        errlog(25, mapping);
        dy = 0.0;
        break;
    }
    return dy;
}   /* map_differential() */


/***********************************************************
* double map_variable(double x, int mapping)
*
* Purpose:
*   Given x, returns y
************************************************************/

double map_variable(double x, int mapping)
{
    double y;
    volatile double r1;

    switch (mapping) {
    case MAPPING_NONE:
        y = x;
        break;
    case MAPPING_VMR:
        x = x > 1.0 ? 1.0 : (x < 0.0 ? 0.0 : x);
        r1 = 1.0 - x;
        r1 += FLT_EPSILON;
        y = x / r1;
        break;
    case MAPPING_EXP:
        y = exp(x);
        break;
    default:
        errlog(123, mapping);
        y = 0.0;
        break;
    }
    return y;
}   /* map_variable() */


/***********************************************************
* double unmap_differential(double y, double dy, int mapping)
*
* Purpose:
*   Given y and dy, returns dx
************************************************************/

double unmap_differential(double y, double dy, int mapping)
{
    double dx;
    volatile double r1;

    switch (mapping) {
    case MAPPING_NONE:
        dx = dy;
        break;
    case MAPPING_VMR:
        r1 = 1.0 + y;
        r1 *= r1;
        dx = dy * ((1.0 + FLT_EPSILON) / r1);
        break;
    case MAPPING_EXP:
        dx = dy / y;
        break;
    default:
        errlog(65, mapping);
        dx = 0.0;
        break;
    }
    return dx;
}   /* unmap_differential() */


/***********************************************************
* double unmap_variable(double y, int mapping)
*
* Purpose:
*   Given y, returns x
************************************************************/

double unmap_variable(double y, int mapping)
{
    double x;

    switch (mapping) {
    case MAPPING_NONE:
        x = y;
        break;
    case MAPPING_VMR:
        x = (y / (1.0 + y)) * (1.0 + FLT_EPSILON);
        break;
    case MAPPING_EXP:
        x = log(y);
        break;
    default:
        errlog(124, mapping);
        x = 0.0;
        break;
    }
    return x;
}   /* unmap_variable() */

#ifdef UNIT_TEST

/***********************************************************
* int main(int argc, char **argv)
*
* Purpose:
*   Tests the functions in this file.
************************************************************/

int main(int argc, char **argv)
{
    struct testvec {
        double x;
        double dx;
        double y;
        double dy;
    };
    struct testvec test_none[] = {
        {-1.0, 0.125, -1.0, 0.125},
        { 0.0, 0.125,  0.0, 0.125},
        { 1.0, 0.125,  1.0, 0.125},
        { 0.0, 0.0, 0.0, 0.0}
    };
    struct testvec test_vmr[] = {
        { 0.0, 0.125, 0.0, 0.125 / (1.0 + FLT_EPSILON)},
        { FLT_EPSILON, 0.125, FLT_EPSILON, (1.0 + FLT_EPSILON) * 0.125},
        { 0.5, 0.125, 1.0 / (1.0 + 2.0 * FLT_EPSILON),
            0.5 * (1.0 + FLT_EPSILON) / ((1.0 + 2.0 * FLT_EPSILON) *
            (1.0 + 2.0 * FLT_EPSILON))},
        { 1.0, 0.125, 1.0 / FLT_EPSILON,
            0.125 * ((1.0 + FLT_EPSILON) / (FLT_EPSILON * FLT_EPSILON))},
        { 1.0, 0.125, 8796102459392.0, 0.0},
        { 1.0, 0.125, DBL_MAX, 0.0},
        { 0.0, 0.0, 0.0, 0.0}
    };
    struct testvec test_exp[] = {
        { -10.0, 0.125, 4.5399929762484854e-5,  0.125 * 4.5399929762484854e-6},
        {   0.0, 0.125,         1.0,       0.125},
        {  10.0, 0.125,  2.2026465794806718e4,  0.125 * 2.2026465794806718e3},
        { 0.0, 0.0, 0.0, 0.0}
    };
    int mode[] =
        { MAPPING_NONE,   MAPPING_VMR,   MAPPING_EXP};
    char mode_name[sizeof(mode) / sizeof(int)][16] =
        {"MAPPING_NONE", "MAPPING_VMR", "MAPPING_EXP"};
    struct testvec *test[] =
        {    test_none,      test_vmr,      test_exp};
    int i, j;
    for (i = 0; i < sizeof(mode) / sizeof(int); ++i) {
        printf("\n");
        for (j = 0; test[i][j].dx != 0.0; ++j) {
            printf("      map_variable(%.*e, %16s) = %.*e\n",
                DBL_DIG + 1,
                test[i][j].x,
                mode_name[i],
                DBL_DIG + 1,
                map_variable(test[i][j].x, mode[i]));
            printf("    unmap_variable(%.*e, %16s) = %.*e\n",
                DBL_DIG + 1,
                test[i][j].y,
                mode_name[i],
                DBL_DIG + 1,
                unmap_variable(test[i][j].y, mode[i]));
            printf("  map_differential(%.*e, %.*e, %16s) = %.*e\n",
                DBL_DIG + 1,
                test[i][j].x,
                DBL_DIG + 1,
                test[i][j].dx,
                mode_name[i],
                DBL_DIG + 1,
                map_differential(test[i][j].x, test[i][j].dx, mode[i]));
            printf("unmap_differential(%.*e, %.*e, %16s) = %.*e\n",
                DBL_DIG + 1,
                test[i][j].y,
                DBL_DIG + 1,
                test[i][j].dy,
                mode_name[i],
                DBL_DIG + 1,
                unmap_differential(test[i][j].y, test[i][j].dy, mode[i]));
            printf("\n");
        }
    }
    return 0;
}   /* main() */

/***********************************************************
* int errlog(const int errnum, const int data)
*
* Purpose:
*   dummy errlog function for unit test
************************************************************/

int errlog(const int errnum, const int data)
{
    printf("errlog(%d, %d)\n", errnum, data);
    return 0;
}   /* errlog() */

#endif /* UNIT_TEST */
