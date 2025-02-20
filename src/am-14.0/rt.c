/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* rt.c                           S. Paine rev. 2022 April 21
*
* Radiative transfer and Planck functions.
************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "am_sysdep.h"
#include "am_types.h"
#include "column.h"
#include "errlog.h"
#include "layer.h"
#include "model.h"
#include "phys_const.h"
#include "planck.h"
#include "specfunc.h"
#include "spectra.h"
#include "rt.h"

static int  compute_column_zenith_opacities(model_t*, model_t*);
static int  compute_layer_line_of_sight_opacities(model_t*);
static void compute_layer_Planck_functions(model_t*, model_t*);
static void layer_radiative_transfer(model_t*, const int, const int);


/***********************************************************
* static int compute_column_zenith_opacities(
*   model_t *model, model_t *lmodel)
*
* Purpose:
*   Computes (or updates when lmodel != NULL) the zenith
*   optical depth for all columns on all layers of the
*   model.
*
*   For an initial model computation, layer->updateflag will
*   be set on every layer.  For an update computation,
*   layer->updateflag will be set only on layers on which at
*   least one column optical depth has changed.
*
*
* Arguments:
*   model_t *model  - pointer to model data structure
*   model_t *lmodel - pointer to a model data structure
*       containing scalar data from a prior computation
*
* Return:
*   0 on success, 1 on error.
************************************************************/

static int compute_column_zenith_opacities(
    model_t *model, model_t *lmodel)
{
    int lnum;

    /*
     * Column opacity computations, and the underlying spectral
     * absorption coefficient computations, are done in sequence.
     * They are relatively expensive per frequency grid point, so they
     * are parallelized at a lower level, over frequency grid
     * iterations.  This allows model computations (or recomputations)
     * involving as little as one spectral absorption coefficient to
     * take advantage of multiple processors.
     */
    for (lnum = model->path_min; lnum <= model->path_mid; ++lnum) {
        layer_t *layer = model->layer[lnum];
        double layer_tstart = 0.0;
        int cnum;
        if (layer->tau == NULL) /* empty layer */
            continue;
        if (model->log_runtimes)
            layer_tstart = am_timer(0.0);
        /*
         * layer->updateflag will indicate if the layer opacity needs
         * to be updated.  It is set on an initial computation of the
         * model (lmodel == NULL), and on subsequent passes if the
         * zenith angle or airmass changes, or if the zenith opacity
         * of any column on the layer is changed in column_opacity().
         */
        layer->updateflag = (lmodel == NULL ||
            fabs(layer->m - lmodel->layer[lnum]->m)
                > (DBL_EPSILON * layer->m) ||
            fabs(model->sec_za - lmodel->sec_za)
                > (DBL_EPSILON * model->sec_za));
        for (cnum = 0; cnum < layer->ncols; ++cnum) {
            column_t *column = layer->column[cnum];
            double col_tstart = 0.0;
            if (model->log_runtimes)
                col_tstart = am_timer(0.0);
            /*
             * If Nscale is zero and nonvariable for this column type,
             * it is empty, and alloc_model_arrays() will not have
             * allocated memory for its opacity.  Skip it.
             */
            if (column->ztau == NULL)
                continue;
            if (column_opacity(model, lmodel, lnum, cnum, &(layer->updateflag)))
                return 1;
            if (model->log_runtimes)
                column->runtime += am_timer(col_tstart);
        }
        if (model->log_runtimes)
            layer->runtime += am_timer(layer_tstart);
    }
    return 0;
}   /* compute_column_zenith_opacities() */


/***********************************************************
* static int compute_layer_line_of_sight_opacities(model_t *model)
*
* Purpose:
*   This function is called to accumulate the column
*   opacities in each layer and compute the layer's
*   line-of-sight opacity.
*
* Arguments:
*   model_t *model  - pointer to model data structure
*
* Return:
*   0 on success, 1 on error.
************************************************************/

static int compute_layer_line_of_sight_opacities(model_t *model)
{
    int lnum;

    /*
     * Here, the zenith spectral opacities for the columns on each
     * layer are adjusted by their respective airmass or zenith angle
     * dependencies, and accumulated to produce the total
     * line-of-sight opacities for the layer.
     *
     * The tight loops over the frequency grid do not parallelize
     * efficiently, so instead these computations are parallelized
     * across layers, enabling multi-layer models to take advantage of
     * up to as many processors as there are layers.
     */
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 1) \
            if (model->isub[0].max - model->isub[0].min >= 16384)
    #endif
    for (lnum = model->path_min; lnum <= model->path_mid; ++lnum) {
        layer_t *layer = model->layer[lnum];
        double layer_tstart = 0.0;
        int i, cnum;
        if (layer->tau == NULL) /* empty layer */
            continue;
        if (model->log_runtimes)
            layer_tstart = am_timer(0.0);
        if (layer->updateflag) {
            for (i = 0; i < model->ngrid; ++i)
                layer->tau[i] = 0.0;
            for (cnum = 0; cnum < layer->ncols; ++cnum) {
                column_t *column = layer->column[cnum];
                int col_typenum = column->col_typenum;
                if (column->ztau == NULL) /* skip empty column */
                    continue;
                switch (col_type[col_typenum].zadep) {
                    double x;
                case ZADEP_NONE:
                    for (i = model->isub[0].min; i <= model->isub[0].max; ++i)
                        layer->tau[i] += column->ztau[i];
                    for (i = model->isub[1].min; i <= model->isub[1].max; ++i)
                        layer->tau[i] += column->ztau[i];
                    break;
                case ZADEP_AIRMASS:
                    for (i = model->isub[0].min; i <= model->isub[0].max; ++i)
                        layer->tau[i] += layer->m * column->ztau[i];
                    for (i = model->isub[1].min; i <= model->isub[1].max; ++i)
                        layer->tau[i] += layer->m * column->ztau[i];
                    break;
                case ZADEP_SEC_ZA:
                    for (i = model->isub[0].min; i <= model->isub[0].max; ++i)
                        layer->tau[i] += model->sec_za * column->ztau[i];
                    for (i = model->isub[1].min; i <= model->isub[1].max; ++i)
                        layer->tau[i] += model->sec_za * column->ztau[i];
                    break;
                case ZADEP_SIN_SQUARED_ZA:
                    x = model->sec_za * model->sec_za;
                    x = (x - 1.0) / x;
                    for (i = model->isub[0].min; i <= model->isub[0].max; ++i)
                        layer->tau[i] += x * column->ztau[i];
                    for (i = model->isub[1].min; i <= model->isub[1].max; ++i)
                        layer->tau[i] += x * column->ztau[i];
                    break;
                case ZADEP_SIN_SQUARED_2ZA:
                    x = model->sec_za * model->sec_za;
                    x = 4.0 * (x - 1.0) / (x * x);
                    for (i = model->isub[0].min; i <= model->isub[0].max; ++i)
                        layer->tau[i] += x * column->ztau[i];
                    for (i = model->isub[1].min; i <= model->isub[1].max; ++i)
                        layer->tau[i] += x * column->ztau[i];
                    break;
                default:
                    #ifdef _OPENMP
                    #pragma omp critical
                    #endif
                    {
                    errlog(96, col_type[col_typenum].zadep);
                    }
                    break;
                }
            }
            for (i = model->isub[0].min; i <= model->isub[0].max; ++i)
                layer->tx[i] = am_exp(-layer->tau[i]);
            for (i = model->isub[1].min; i <= model->isub[1].max; ++i)
                layer->tx[i] = am_exp(-layer->tau[i]);
        }
        if (model->log_runtimes)
            layer->runtime += am_timer(layer_tstart);
    }
    if (errtest(96)) /* an error occurred in the parallel block */
        return 1;
    return 0;
}   /* compute_layer_line_of_sight_opacities() */


/***********************************************************
* static void compute_layer_Planck_functions(
*   model_t *model, model_t *lmodel)
*
* Purpose:
*   Computes spectral Planck functions B(f) on all the model
*   layers.  These are computed for all layers on the first
*   pass through the model, and on subsequent passes only if
*   the layer midpoint temperature, base temperature, or
*   grid range has changed.
*
*   For models defined by midpoint temperature, B(f) is
*   evaluated at the midpoint temperature.  For models
*   defined by base temperatures, B(f) is evaluated at the
*   base temperature, and the layer radiances are computed
*   using the linear-in-tau Planck function approximation.
*
* Arguments:
*   model_t *model  - pointer to model data structure
*   model_t *lmodel - pointer to a model data structure
*       containing scalar data from a prior computation.
************************************************************/

static void compute_layer_Planck_functions(
    model_t *model, model_t *lmodel)
{
    int lnum, lnum_min;

    /*
     * In PTMODE_TBASE, to support the linear-in-tau radiative
     * transfer approximation, we need the Planck function on the base
     * of the layer above the highest layer in the propagation path,
     * unless that layer is the top model layer, which is always
     * isothermal.
     */
    if (model->path_min && (model->PTmode & PTMODE_TBASE)) {
        lnum_min = model->path_min - 1;
    } else {
        lnum_min = 0;
    }
    /*
     * The function B_Planck() is parallelized over frequency grid
     * iterations, so no need to parallelize over layers here.
     */
    for (lnum = lnum_min; lnum <= model->path_mid; ++lnum) {
        layer_t *layer = model->layer[lnum];
        double layer_tstart = 0.0;
        /*
         * Planck functions on empty layers aren't needed in midpoint
         * T mode.
         */
        if ((layer->tau == NULL) && (model->PTmode & PTMODE_T))
            continue;
        if (model->log_runtimes)
            layer_tstart = am_timer(0.0);
        if (lmodel == NULL ||
            fabs(layer->T - lmodel->layer[lnum]->T)
                > (DBL_EPSILON * layer->T) ||
            fabs(layer->Tbase - lmodel->layer[lnum]->Tbase)
                > (DBL_EPSILON * layer->Tbase) ||
            compare_spectral_subgrid_ranges(model, lmodel)) {
            double T;
            if (model->PTmode & PTMODE_TBASE)
                T = layer->Tbase;
            else
                T = layer->T;
            B_Planck(
                &layer->B[model->isub[0].min],
                &model->f[model->isub[0].min],
                model->df,
                model->isub[0].max - model->isub[0].min + 1,
                T);
            B_Planck(
                &layer->B[model->isub[1].min],
                &model->f[model->isub[1].min],
                model->df,
                model->isub[1].max - model->isub[1].min + 1,
                T);
        }
        if (model->log_runtimes)
            layer->runtime += am_timer(layer_tstart);
    }
    return;
}   /* compute_layer_Planck_functions() */


/***********************************************************
* void compute_spectral_radiance(model_t *model, model_t *lmodel)
*
* Purpose:
*   Performs an iterative radiative transfer computation
*   through the model layers in the propagation path to
*   compute the emerging spectral radiance.
*
* Arguments:
*   model_t *model  - pointer to model data structure
*   model_t *lmodel - pointer to a model data structure containing
*                     scalar data from a prior computation
************************************************************/

void compute_spectral_radiance(model_t *model, model_t *lmodel)
{
    gridsize_t i;
    int lnum;

    /*
     * Initialize the background radiance I0.
     */
    if (lmodel == NULL ||
        fabs(model->T0 - lmodel->T0) > (DBL_EPSILON * model->T0) ||
        compare_spectral_subgrid_ranges(model, lmodel)) {
        B_Planck(
            &model->I0[model->isub[0].min],
            &model->f[model->isub[0].min],
            model->df,
            model->isub[0].max - model->isub[0].min + 1,
            model->T0);
        B_Planck(
            &model->I0[model->isub[1].min],
            &model->f[model->isub[1].min],
            model->df,
            model->isub[1].max - model->isub[1].min + 1,
            model->T0);
    }
    for (i = model->isub[0].min; i <= model->isub[0].max; ++i)
        model->I[i] = model->I0[i];
    for (i = model->isub[1].min; i <= model->isub[1].max; ++i)
        model->I[i] = model->I0[i];
    /*
     * Initialize layer Planck functions
     */
    compute_layer_Planck_functions(model, lmodel);
    /*
     * Iterate through the layers.  Note that the path indices may be
     * such that one or both of these loops is skipped.
     */
    for (lnum = model->path_begin; lnum <= model->path_mid; ++lnum)
        layer_radiative_transfer(model, lnum, RT_DOWN);
    for (lnum = model->path_mid;   lnum >= model->path_end; --lnum)
        layer_radiative_transfer(model, lnum, RT_UP);
    return;
}   /* compute_spectral_radiance() */


/***********************************************************
* int compute_opacity_spectrum(model_t *model, model_t *lmodel)
*
* Purpose:
*   Computes the total optical depth for the model propagation
*   path.
*
* Arguments:
*   model_t *model  - pointer to model structure
*   model_t *lmodel - pointer to a model data structure containing
*                     scalar data from a prior computation
************************************************************/

int compute_opacity_spectrum(model_t *model, model_t *lmodel)
{
    int lnum;
    gridsize_t i;

    /*
     * Compute zenith and line-of-sight optical depths for all
     * the model layers
     */
    if (compute_column_zenith_opacities(model, lmodel))
        return 1;
    if (compute_layer_line_of_sight_opacities(model))
        return 1;
    /*
     * Compute the total line-of-sight optical depth through the
     * whole propagation path.
     */
    for (i = model->isub[0].min; i <= model->isub[0].max; ++i)
        model->tau[i] = 0.0;
    for (i = model->isub[1].min; i <= model->isub[1].max; ++i)
        model->tau[i] = 0.0;
    for (lnum = model->path_begin; lnum <= model->path_mid; ++lnum) {
        layer_t *layer = model->layer[lnum];
        if (layer->tau == NULL) /* empty layer */
            continue;
        for (i = model->isub[0].min; i <= model->isub[0].max; ++i)
            model->tau[i] += layer->tau[i];
        for (i = model->isub[1].min; i <= model->isub[1].max; ++i)
            model->tau[i] += layer->tau[i];
    }
    for (lnum = model->path_mid; lnum >= model->path_end; --lnum) {
        layer_t *layer = model->layer[lnum];
        if (layer->tau == NULL) /* empty layer */
            continue;
        for (i = model->isub[0].min; i <= model->isub[0].max; ++i)
            model->tau[i] += layer->tau[i];
        for (i = model->isub[1].min; i <= model->isub[1].max; ++i)
            model->tau[i] += layer->tau[i];
    }
    return 0;
}   /* compute_opacity_spectrum() */


/***********************************************************
* static void layer_radiative_transfer(
*     model_t *model,
*     const int lnum,
*     const int direction)
*
* Purpose:
*   Computes radiative transfer across a single layer,
*   updating the radiance array model->I.
*
* Arguments:
*   model_t *model - pointer to model data structure
*   const int lnum - index of layer.
*   const int direction - direction of propagation, either
*     RT_UP or RT_DOWN
************************************************************/

static void layer_radiative_transfer(
    model_t *model,
    const int lnum,
    const int direction)
{
    int subgrid;
    gridsize_t i;
    layer_t *layer = model->layer[lnum];

    if (layer->tau == NULL) /* empty layer */
        return;
    /*
     * Models defined by base temperatures use the linear-in-tau
     * Planck function approximation to compute the layer
     * radiance for all layers except layer 0, which is defined
     * to be isothermal in am.
     *
     * Otherwise, the layer radiance is computed using the Planck
     * function evaluated at the layer midpoint temperature.
     */
    if (lnum && (model->PTmode & PTMODE_TBASE)) {
        double *Bin, *Bout;
        if (direction == RT_DOWN) {
            Bin = model->layer[lnum-1]->B;
            Bout = layer->B;
        } else {
            Bin = layer->B;
            Bout = model->layer[lnum-1]->B;
        }
        for (subgrid = 0; subgrid <= 1; ++subgrid) {
            for (i  = model->isub[subgrid].min;
                 i <= model->isub[subgrid].max;
                 ++i) {
                double u, v;
                double tau = layer->tau[i];
                double tx = layer->tx[i];
                if (tau > 2.0e-3) { /* threshold for Taylor series */
                    v = (1.0 - tx) / tau;
                    u = 1.0 - v;
                    v -= tx;
                } else {
                    u = tau * (0.5 - tau * ((1. / 6.) - (1. / 24.) *  tau));
                    v = tau * (0.5 - tau * ((1. / 3.) - 0.125 * tau));
                }
                model->I[i] *= tx;
                model->I[i] += Bout[i] * u + Bin[i] * v;
            }
        }
    } else {
        for (subgrid = 0; subgrid <= 1; ++subgrid) {
            for (i  = model->isub[subgrid].min;
                 i <= model->isub[subgrid].max;
                 ++i) {
                double tau = layer->tau[i];
                double tx = layer->tx[i];
                model->I[i] *= tx;
                model->I[i] += tau > DBL_EPSILON ?
                    layer->B[i] * (1.0 - tx) :
                    layer->B[i] * tau;
            }
        }
    }
    return;
} /* layer_radiative_transfer() */
