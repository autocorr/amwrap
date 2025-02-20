/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* simplex.c                       S. Paine rev. 2024 June 27
*
* A simplex in a subspace of model parameter space is used
* for carrying out model fits by downhill simplex
* minimization.  For computing Jacobians, parts of the same
* simplex structure are used to hold differentiation
* variables and associated data.  This file contains
* functions for simplex operations.
************************************************************/

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "am_types.h"
#include "errlog.h"
#include "simplex.h"


static int increment_simplex_dimension(simplex_t*);

/***********************************************************
* int add_var_to_simplex(
*   simplex_t *simplex,
*   const char *name,
*   double *varptr,
*   double init,
*   double scale,
*   int unitnum,
*   int mapping)
*
* Purpose:
*   Adds a new variable to the variable space of a simplex.
*   Memory is realloc()'ed for the added dimension
*   associated with the new variable, with a call to
*   increment_simplex_dimension().  The simplex variable
*   tables and initial vertex vectors are then updated.
*
* Arguments:
*   simplex_t *simplex - pointer to simplex structure.
*   char *name         - descriptive variable name
*   double *varptr     - pointer to new variable
*   double init        - initial value of new variable
*   double scale       - characteristic scale of new
*                        variable
*   int unitnum        - unit number for new variable
*   int mapping        - mapping from physical domain to
*                        simplex coordinates
*
* Return:
*   0 if successful, 1 otherwise
************************************************************/

int add_var_to_simplex(
    simplex_t *simplex,
    const char *name,
    double *varptr,
    double init,
    double scale,
    int unitnum,
    int mapping)
{
    unsigned int i, m, n;

    /*
     * The initial scale of the parameter is assumed positive, and
     * must be non-zero.
     */
    scale = fabs(scale);
    if (scale < DBL_MIN) {
        errlog(61, 0);
        return 1;
    }
    /*
     * If log mapping of simplex coordinates is in effect, the
     * initial parameter value must be positive.
     */
    if (simplex->logarithmic) {
        if (init < DBL_MIN) {
            errlog(46, 0);
            return 1;
        }
    }
    /*
     * Check for a redundant variable.
     */
    if (isvar(simplex, varptr)) {
        errlog(106, 0);
        return 1;
    }
    if (increment_simplex_dimension(simplex))
        return 1;
    m = simplex->n - 1; /* dimensions run from 0 to m */
    n = simplex->n;     /* vertices run from 0 to n   */

    if ((simplex->name[m] =
                (char*)malloc((strlen(name) + 1) * sizeof(char))) == NULL) {
        errlog(121, 0);
    }
    strcpy(simplex->name[m], name);
    simplex->init[m]    = init;
    simplex->scale[m]   = scale;
    simplex->unitnum[m] = unitnum;
    simplex->mapping[m] = mapping;
    simplex->varptr[m]  = varptr;

    /*
     * Fill in the elements of the new vertex vector
     */
    for (i = 0; i < m; ++i) {
        simplex->vertex[n][i] =
            simplex->logarithmic ? log(simplex->init[i]) : simplex->init[i];
    }
    simplex->vertex[n][m] =
        simplex->logarithmic ? log(init + scale) : init + scale;
    /*
     * Fill in the new last element of the old vertex vectors
     */
    for (i = 0; i < n; ++i) {
        simplex->vertex[i][m] =
            simplex->logarithmic ? log(init) : init;
    }
    return 0;
}   /* add_var_to_simplex() */


/***********************************************************
* int create_null_simplex(simplex_t *simplex)
*
* Purpose:
*   Create a zero-dimensional simplex, with no fit variables
*   assigned.  The only non-NULL entities are those with
*   n+1 dimensions-- the vertex table has a single entry
*   pointing to a NULL vertex vector, and the table of
*   goodness of fit estimators E has a single entry for
*   the fit estimator at this vertex.  For a null fit,
*   this will just be the fit estimator computed for the
*   model as given in the configuration file.
*
* Arguments:
*   simplex_t *simplex - pointer to simplex structure
*
* Return:
*   0 if successful
*   1 otherwise
************************************************************/

int create_null_simplex(simplex_t *simplex)
{
    if (simplex->vertex != NULL || simplex->E != NULL) {
        errlog(95, 0);
        return 1;
    }
    if ((simplex->E      =  (double*)malloc(sizeof(double)))  == NULL ||
        (simplex->vertex = (double**)malloc(sizeof(double*))) == NULL) {
        errlog(18, 0);
        return 1;
    }
    simplex->vertex[0] = NULL;
    return 0;
}   /* create_null_simplex() */


/***********************************************************
* unsigned int get_simplex_variable_index(simplex_t *simplex, double *x)
*
* Purpose:
*   Given a pointer to a simplex variable, returns the index number
*   of the variable.  If the pointer does not correspond to a
*   simplex variable, simplex->n is returned.
*
* Arguments:
*   simplex_t *simplex - pointer to simplex structure.
*   double *x          - pointer to variable
*
* Return:
*   index number of variable, as an unsigned int
************************************************************/

unsigned int get_simplex_variable_index(simplex_t *simplex, double *x)
{
    unsigned int j;

    for (j = 0; j < simplex->n; ++j) {
        if (x == simplex->varptr[j])
            break;
    }
    return j;
}   /* get_simplex_variable_index() */


/***********************************************************
* static int increment_simplex_dimension(simplex_t *simplex)
*
* Purpose:
*   Increase the dimensionality of an existing simplex.  All
*   arrays in the simplex structure get realloc()'ed to
*   accommodate the increased dimensions.
*
* Arguments:
*   simplex_t *simplex - pointer to simplex structure
*
* Return:
*   0 if successful
*   1 otherwise
************************************************************/

static int increment_simplex_dimension(simplex_t *simplex)
{
    int i, n;
    void *tptr;

    if (simplex->n == 0) {
        if (create_null_simplex(simplex))
            return 1;
    }
    n = simplex->n + 1;
    if ((tptr = realloc(simplex->name, n * sizeof(char*))) == NULL) {
        errlog(18, n);
        return 1;
    }
    simplex->name = (char**)tptr;
    simplex->name[n-1] = NULL;
    if ((tptr = realloc(simplex->init, n * sizeof(double))) == NULL) {
        errlog(18, n);
        return 1;
    }
    simplex->init = (double*)tptr;
    if ((tptr = realloc(simplex->scale, n * sizeof(double))) == NULL) {
        errlog(18, n);
        return 1;
    }
    simplex->scale = (double*)tptr;
    if ((tptr = realloc(simplex->unitnum, n * sizeof(int))) == NULL) {
        errlog(18, n);
        return 1;
    }
    simplex->unitnum = (int*)tptr;
    if ((tptr = realloc(simplex->mapping, n * sizeof(int))) == NULL) {
        errlog(18, n);
        return 1;
    }
    simplex->mapping = (int*)tptr;
    if ((tptr = realloc(simplex->varptr, n * sizeof(double*))) == NULL) {
        errlog(18, n);
        return 1;
    }
    simplex->varptr = (double**)tptr;
    if ((tptr = realloc(simplex->vertex, (n + 1) * sizeof(double*))) == NULL) {
        errlog(18, n);
        return 1;
    }
    simplex->vertex = (double**)tptr;
    /*
     * set the new vertex vector pointer to NULL before calling
     * realloc() to increase the dimensions of all the vertices.
     */
    simplex->vertex[n] = NULL;
    for (i = 0; i <= n; ++i) {
        if ((tptr = realloc(simplex->vertex[i], n * sizeof(double))) == NULL) {
            errlog(18, n);
            return 1;
        }
        simplex->vertex[i] = (double*)tptr;
    }
    if ((tptr = realloc(simplex->E, (n + 1) * sizeof(double))) == NULL) {
        errlog(18, n);
        return 1;
    }
    simplex->E = (double*)tptr;
    if ((tptr = realloc(simplex->pbar, n * sizeof(double))) == NULL) {
        errlog(18, n);
        return 1;
    }
    simplex->pbar = (double*)tptr;
    if ((tptr = realloc(simplex->p1, n * sizeof(double))) == NULL) {
        errlog(18, n);
        return 1;
    }
    simplex->p1 = (double*)tptr;
    if ((tptr = realloc(simplex->p2, n * sizeof(double))) == NULL) {
        errlog(18, n);
        return 1;
    }
    simplex->p2 = (double*)tptr;
    if ((tptr = realloc(simplex->pc, n * sizeof(double))) == NULL) {
        errlog(18, n);
        return 1;
    }
    simplex->pc = (double*)tptr;
    simplex->n = n;
    return 0;
}   /* increment_simplex_dimension() */


/***********************************************************
* int isvar(simplex_t *simplex, double *x)
*
* Purpose:
*   Given a pointer x, this function looks for a matching
*   pointer in the array of pointers simplex.varptr.  If a
*   match is found, then x points to a simplex variable.
*
* Arguments:
*   simplex_t *simplex - pointer to simplex structure.
*   double *x          - pointer to variable
*
* Return:
*   1 if x points to a simplex variable, 0 otherwise
************************************************************/

int isvar(simplex_t *simplex, double *x)
{
    unsigned int j;

    for (j = 0; j < simplex->n; ++j) {
        if (x == simplex->varptr[j])
            return 1;
    }
    return 0;
}   /* isvar() */


/***********************************************************
* void free_simplex_entities(simplex_t *simplex)
*
* Purpose:
*   Frees the memory allocated to all the entities within
*   a simplex data structure and NULLs the pointers.
*
* Arguments:
*   simplex_t *simplex - pointer to simplex structure.
*
* Return:
*   none
************************************************************/

void free_simplex_entities(simplex_t *simplex)
{
    free(simplex->pc);
    simplex->pc   = NULL;
    free(simplex->p2);
    simplex->p2   = NULL;
    free(simplex->p1);
    simplex->p1   = NULL;
    free(simplex->pbar);
    simplex->pbar = NULL;
    free(simplex->E);
    simplex->E    = NULL;
    if (simplex->n) {
        unsigned int i;
        for (i = 0; i < simplex->n; ++i) {
            free(simplex->name[i]);
            simplex->name[i] = NULL;
        }
        for (i = 0; i <= simplex->n; ++i) {
            free(simplex->vertex[i]);
            simplex->vertex[i] = NULL;
        }
    }
    free(simplex->name);
    simplex->name    = NULL;
    free(simplex->vertex);
    simplex->vertex  = NULL;
    free(simplex->mapping);
    simplex->mapping = NULL;
    free(simplex->varptr);
    simplex->varptr  = NULL;
    free(simplex->unitnum);
    simplex->unitnum = NULL;
    free(simplex->scale);
    simplex->scale   = NULL;
    free(simplex->init);
    simplex->init = NULL;
    simplex->n = 0;
    return;
}   /* free_simplex_entities() */


/***********************************************************
* void reset_simplex_vertices(simplex_t *simplex)
*
* Purpose:
*   Resets a simplex to its initial dimensions, with edges
*   parallel to the coordinate axes.  If this is the start
*   of a fit, and simplex->reinit is set, then vertex[0]
*   is reset to the initial trial parameters.  Otherwise,
*   (as for restarts, or at the start of a new fit when
*   simplex->reinit is not set) the simplex is reset about
*   the current vertex[0].
*
* Arguments:
*   simplex_t *simplex - pointer to simplex structure.
*
* Return:
*   0 if OK, 1 otherwise
************************************************************/

void reset_simplex_vertices(simplex_t *simplex)
{
    unsigned int i, j;

    /*
     * If simplex.reinit is set, reset the original initial
     * values from the config file.  Otherwise, leave vertex[0]
     * where it was an the end of the last fit (but undo the log
     * mapping, if log mapping is in effect).
     */
    if ((simplex->restart_count == 0) && simplex->reinit) {
        for (j = 0; j < simplex->n; ++j)
            simplex->vertex[0][j] = simplex->init[j];
    } else if (simplex->logarithmic) {
        for (j = 0; j < simplex->n; ++j)
            simplex->vertex[0][j] = exp(simplex->vertex[0][j]);
    }
    /*
     * If any of the coordinates of vertex[0] are less than the
     * scale for the corresponding variable, set them equal to
     * the scale.  When log mapping is in effect, this keeps a
     * variable from getting stuck at zero from one fit to the
     * next.  If log mapping is turned off, this will restore  a
     * variable which converged to a negative value on a previous
     * fit to positive territory.  (Most am fit variables
     * correspond to positive physical quantities, but a
     * parameter near enough to zero could end up converging to a
     * negative value if log mapping is turned off.)
     */
    for (j = 0; j < simplex->n; ++j) {
        if (simplex->vertex[0][j] < simplex->scale[j])
            simplex->vertex[0][j] = simplex->scale[j];
    }
    /*
     * For the other vertices, vertex[i] = vertex[0] + e[i-1] *
     * scale[i-1], where e[i] is the unit vector corresponding to
     * the ith fit variable.
     */
    for (i = 1; i <= simplex->n; ++i) {
        for (j = 0; j < simplex->n; ++j)
            simplex->vertex[i][j] = simplex->vertex[0][j];
    }
    for (i = 1; i <= simplex->n; ++i)
        simplex->vertex[i][i-1] += simplex->scale[i-1];
    /*
     * Now return to logarithmic coordinates, if needed.
     */
    if (simplex->logarithmic) {
        for (i = 0; i <= simplex->n; ++i)
            for (j = 0; j < simplex->n; ++j)
                simplex->vertex[i][j] = log(simplex->vertex[i][j]);
    }
    return;
}   /* reset_simplex_vertices() */


/***********************************************************
* double simplex_scaled_diameter(simplex_t *simplex)
*
* Purpose:
*   Reports the scaled diameter of a simplex, defined as the
*   square root of the sum of the squared dimensionless
*   distances (max - min) / scale over all the simplex
*   variables.  This quantity is used as a convergence
*   metric during fits.
*
* Arguments:
*   simplex_t *simplex - pointer to simplex structure.
*   
* Return:
*   scaled diameter of the simplex, as a double.
************************************************************/

double simplex_scaled_diameter(simplex_t *simplex)
{
    double r1 = 0.0;
    unsigned int i;

    for (i = 0; i < simplex->n; ++i) {
        double r2 = simplex_variable_range(simplex, i) / simplex->scale[i];
        r1 += r2 * r2;
    }
    return sqrt(r1);
}   /* simplex_scaled_diameter() */


/***********************************************************
* double simplex_scaled_distance(simplex_t *simplex, double *v1, double *v2)
*
* Purpose:
*   Reports the scaled distance between two vectors in
*   simplex space, defined as the square root of the sum of
*   the squared differences between the dimensionless
*   distances (v1[i] - v2[i]) / scale[i] over all the
*   simplex coordinates.  This quantity is used to check for
*   consistent convergence after a fit restart.
*
* Arguments:
*   simplex_t *simplex - pointer to simplex structure.
*   double *v1         - first vector
*   double *v2         - second vector
*   
* Return:
*   scaled distance from v1 to v2, as a double.
************************************************************/

double simplex_scaled_distance(simplex_t *simplex, double *v1, double *v2)
{
    double r1 = 0.0;
    unsigned int i;

    for (i = 0; i < simplex->n; ++i) {
        double r2 = (v1[i] - v2[i]) / simplex->scale[i];
        r1 += r2 * r2;
    }
    return sqrt(r1);
}   /* simplex_scaled_distance() */


/***********************************************************
* double simplex_variable_range(simplex_t *simplex, unsigned int nvar)
*
* Purpose:
*   Reports the (max - min) range of a specified simplex
*   variable.
*
* Arguments:
*   simplex_t *simplex - pointer to simplex structure.
*   nvar               - index of variable
*
* Return:
*   max - min range of variable nvar, as a double.  If nvar
*   is invalid, logs an error and returns 0.0.
************************************************************/

double simplex_variable_range(simplex_t *simplex, unsigned int nvar)
{
    unsigned int i;
    double min, max;

    if (nvar >= simplex->n) {
        errlog(26, 0);
        return 0.0;
    }
    min = max = simplex->vertex[0][nvar];
    for (i = 1; i <= simplex->n; ++i) {
        if (min > simplex->vertex[i][nvar])
            min = simplex->vertex[i][nvar];
        if (max < simplex->vertex[i][nvar])
            max = simplex->vertex[i][nvar];
    }
    if (simplex->logarithmic) {
        min = exp(min);
        max = exp(max);
    }
    return max - min;
}   /* simplex_variable_range() */
