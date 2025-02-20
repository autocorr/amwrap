/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* layer.c                         S. Paine rev. 2024 July 19
*
* Operations on model layers and their scalar variables.
************************************************************/

#include <float.h>
#include <math.h>
#include <stdlib.h>

#include "abscoeff.h"
#include "am_alloc.h"
#include "am_types.h"
#include "column.h"
#include "errlog.h"
#include "h2o_psat.h"
#include "interp.h"
#include "layer.h"
#include "mapping.h"
#include "model.h"
#include "molecules.h"
#include "nscale.h"
#include "phys_const.h"
#include "units.h"

static double adjusted_selfbroad_vmr_func(
                  const double, const double, const double);
static void   compute_hydrostatic_layer_column_densities(layer_t *layer);
static void   compute_nonhydrostatic_layer_column_densities(layer_t *layer);
static int    get_parent_lnum(const model_t*, const int);
static double layer_base_dry_lapse_rate(const model_t*, const int);
static double layer_base_moist_lapse_rate(const model_t*, const int);
static double layer_environmental_lapse_rate(const model_t*, int);
static double layer_h2o_vmr(const model_t*, const int);


/***********************************************************
* static double adjusted_selfbroad_vmr_func(
*   const double vmr,
*   const double vmr_def,
*   const double tol);
*
* Purpose:
*   Computes an adjusted volume mixing ratio that is used
*   when computing collisional widths.
*
*   To optimize re-use of previously computed or cached
*   absorption coefficient spectra, it is desirable to
*   ignore negligible changes in self broadening.  For
*   example, if a hydrostatic model contains an adjustable
*   scale factor on a water vapor profile, a change in this
*   scale factor will cause all other mixing ratios in the
*   model to change.  These changes will result in dcache
*   misses, or will trigger recomputation of spectra already
*   in memory.
*
*   The performance optimization adopted here is to use an
*   adjusted mixing ratio value when computing self
*   broadening.  For small deviations from the (possibly
*   zero) default vmr for a given species, the adjusted
*   value is held constant at the default value.  The
*   user-modifiable parameter selfbroad_vmr_tol sets the
*   maximum deviation for which the adjusted mixing ratio is
*   held constant.  For larger deviations, the adjusted
*   value varies linearly with the exact value, catching up
*   with the exact value at the ends of the interval [0,1].
*
* Arguments:
*   const double vmr     - exact mixing ratio
*   const double vmr_def - default mixing ratio
*   const double tol     - tolerance parameter
*
* Return:
*   Adjusted volume mixing ratio
************************************************************/

static double adjusted_selfbroad_vmr_func(
    const double vmr,
    const double vmr_def,
    const double tol)
{
    double vmr_adj;

    if (tol == 0.0) {
        vmr_adj = vmr;
    } else if (vmr < vmr_def - tol - DBL_EPSILON) {
        vmr_adj = vmr * vmr_def / (vmr_def - tol);
    } else if (vmr <= vmr_def + tol + DBL_EPSILON) {
        vmr_adj = vmr_def;
    } else {
        double u = 1.0 - vmr;
        double u_def = 1.0 - vmr_def;
        vmr_adj = 1.0 - u * u_def / (u_def - tol);
    }
    return vmr_adj;
}   /* adjusted_selfbroad_vmr_func() */


/***********************************************************
* void compute_adjusted_selfbroad_vmrs(model_t *model)
*
* Purpose:
*   On all model layers, computes the adjusted mixing ratios
*   that are used for computing self broadening.  The
*   purpose of the adjustment is to optimize re-use of
*   spectral absorption coefficients that are already in
*   memory or in the dcache.  For details, see the comment
*   in the function adjusted_selfbroad_vmr_func(), which is
*   called by this function.
*
*   In a model configuration file, this adjustment can be
*   disabled globally by setting selfbroad_vmr_tol to zero.
*   It can also be disabled for a given lineshape and
*   species on a given model layer using the
*   strict_selfbroad keyword in a lineshape statement.
*
* Arguments:
*   model_t *model - pointer to model data structure
************************************************************/

void compute_adjusted_selfbroad_vmrs(model_t *model)
{
    int lnum;

    for (lnum = 0; lnum < model->nlayers; ++lnum) {
        layer_t *layer = model->layer[lnum];
        int cnum;
        for (cnum = 0; cnum < layer->ncols; ++cnum) {
            column_t *column = layer->column[cnum];
            int col_typenum = column->col_typenum;
            int knum;
            /*
             * If the strict self broadening flag is set for this
             * absorption coefficient type on this layer, the scaled
             * vmr computed in set_def_layer_heights_and_col_densities()
             * is used unmodified. The VMR_ERROR bit in vmr_stat is set
             * here if there were any other error flag bits set during
             * that computation.
             *
             * If the strict self broadening flag is not set, modify
             * the mixing ratio for self broadening according to the
             * value of selfbroad_vmr_tol.
             */
            for (knum = 0; knum < column->n_abscoeffs; ++knum) {
                abscoeff_t *abscoeff = column->abscoeff[knum];
                if (layer->strict_selfbroad[abscoeff->k_typenum]) {
                    abscoeff->vmr_selfbroad =
                        column->vmr_scaled * col_type[col_typenum].pmr[knum];
                    if ((column->N_scaled > 0.0) &&
                        (column->vmr_stat & VMR_ALL_ERRS)) {
                        column->vmr_stat |= VMR_ERROR;
                        errlog(84, 0);
                    }
                } else {
                    abscoeff->vmr_selfbroad = adjusted_selfbroad_vmr_func(
                            column->vmr_scaled *
                                col_type[col_typenum].pmr[knum],
                            col_type[col_typenum].default_vmr *
                                col_type[col_typenum].pmr[knum],
                            model->selfbroad_vmr_tol);
                }
            }
        }
    }
    return;
}   /* compute_adjusted_selfbroad_vmrs() */


/***********************************************************
* static void compute_hydrostatic_layer_column_densities(layer_t *layer)
*
* Purpose:
*   For hydrostatic models, computes mean molecular mass
*   and column densities on a layer.  Mixing ratios are
*   normalized by scaling default mixing ratios; a warning
*   is set if this is not possible.
*
* Arguments:
*   layer_t *layer - pointer to layer structure
************************************************************/

static void compute_hydrostatic_layer_column_densities(layer_t *layer)
{
    double N0;
    double N;
    double Nbal;
    double vmr_bal;
    double r1 = 0.0;
    double r2 = 0.0;
    double r3 = 0.0;
    double r4 = 0.0;
    double r5;
    int cnum;

    /*
     * If a layer has negligible dP, set vmr-based column
     * densities to zero, and log a vmr normalization warning if
     * there are any non-zero explicity-defined column densities.
     */
    if (layer->dP < DBL_EPSILON * layer->Pbase) {
        for (cnum = 0; cnum < layer->ncols; ++cnum) {
            column_t *column = layer->column[cnum];
            int col_typenum = column->col_typenum;
            if (col_type[col_typenum].flags & COL_PARAMETRIC) {
                continue;
            } else if (column->N_mode & N_EXPLICIT) {
                /* N_base sets the scale for negligible N. */
                double N_base =
                    layer->Pbase / (layer->M0 * PCONST_AMU * layer->gbase);
                if (column->N_scaled > N_base * DBL_EPSILON) {
                    if (!(layer->vmr_stat & VMR_SUM_EXCEEDS_UNITY))
                        errlog(101, 0);
                    layer->vmr_stat |= VMR_SUM_EXCEEDS_UNITY;
                }
            } else {
                column->N = 0.0;
                column->N_scaled = 0.0;
            }
        }
        /*
         * For a zero-dP layer, Mair just defaults to M0.
         */
        layer->Mair = layer->M0;
        return;
    }
    /*
     * N0 is the total hydrostatic column density assuming the
     * mean relative molecular mass is equal to the default value
     * M0.
     */
    N0  = layer->dP / (layer->M0 * PCONST_AMU * layer->gbase);
    N0 *= unit_tab[UNIT_ONE_ON_MSQUARED].factor;
    N0 *= unit_tab[UNIT_M].factor;
    N0 /= unit_tab[UNIT_PA].factor;
    /*
     * Compute layer->Mair, the mean relative molecular mass
     * after adjusting for all columns that have been flagged to
     * be included in the average mass computation.  N is the
     * total hydrostatic column density corresponding to Mair.
     */
    for (cnum = 0; cnum < layer->ncols; ++cnum) {
        double r6;
        column_t *column = layer->column[cnum];
        int col_typenum = column->col_typenum;
        if (col_type[col_typenum].flags & COL_PARAMETRIC)
            continue;
        if (layer->Mair_flag[col_typenum])
            r6 = 1.0 - col_type[col_typenum].Mr / layer->M0;
        else
            r6 = 0.0;
        if (column->N_mode & N_EXPLICIT) {
            r1 += column->N_scaled;
            r2 += r6 * column->N_scaled / N0;
        } else if (column->vmr_stat & VMR_USER_DEFINED) {
            r3 += column->vmr_scaled;
            r4 += r6 * column->vmr_scaled;
        }
    }
    /*
     * If the user has specified a set of column densities such
     * that the weighted mean molecular mass is undefined,
     * negative, or zero, log a warning and don't attempt to
     * compute a mass adjustment.
     *
     * Note that significantly super-hydrostatic combinations of
     * vmr's and explicit column densities won't necessarily
     * trigger this warning, but will trigger the vmr sum warning
     * below.  In some cases, it is desirable to have a little
     * numeric headroom above perfect hydrostatic balance, as
     * when computing Jacobians with no automatically-scaled vmrs
     * on a layer to obtain true partial derivatives.
     */
    r2 = 1.0 + r2;
    r4 = 1.0 - r4;
    if ((r2 <= 0.0) || (r4 <= 0.0)) {
        errlog(157, 0);
        layer->vmr_stat |= VMR_WEIGHTED_MEAN_MASS_UNDEFINED;
        r5 = 1.0;
    } else {
        r5 = r2 / r4;
    }
    N = N0 * r5;
    layer->Mair = layer->M0 / r5;
    /*
     * Nbal is the balance remaining after subtracting
     * user-defined column densities from the total hydrostatic
     * column density.  Default vmrs get scaled by a factor
     * (Nbal / N), in addition to any applied Nscale factor.
     *
     * If Nbal is negative, then the combination of user-defined
     * column densities and mixing ratios is super-hydrostatic.
     * If so, log a warning and clamp Nbal at zero.
     */
    Nbal = N * (1.0 - r3) - r1;
    if (Nbal < 0.0) {
        double vmr_sum_err = -(Nbal / N);
        if (vmr_sum_err > VMR_SUM_WARNING_TOL) {
            if (!(layer->vmr_stat & VMR_SUM_EXCEEDS_UNITY))
                errlog(101, 0);
            layer->vmr_stat |= VMR_SUM_EXCEEDS_UNITY;
        }
        Nbal = 0.0;
    }
    vmr_bal = Nbal / N;
    /*
     * Now finish up column density and mixing ratio
     * computations.  For user-defined column density, compute
     * the corresponding mixing ratio; for user-defined mixing
     * ratio, compute the column density.  Scale default mixing
     * ratios by Nbal / N (in addition to any Nscale factor that
     * may have been applied) and compute the corresponding
     * column density.
     */
    for (cnum = 0; cnum < layer->ncols; ++cnum) {
        column_t *column = layer->column[cnum];
        int col_typenum = column->col_typenum;
        if (col_type[col_typenum].flags & COL_PARAMETRIC) {
            continue;
        } else if (column->N_mode & N_EXPLICIT) {
            double vmr;
            if (column->N <= N * DBL_EPSILON)
                vmr = 0.0;
            else if (N < column->N * DBL_EPSILON)
                vmr = 1.0 / DBL_EPSILON;
            else
                vmr = column->N / N;
            column->xvmr = map_variable(vmr, MAPPING_VMR);
            if (column->N_scaled <= N * DBL_EPSILON)
                column->vmr_scaled = 0.0;
            else if (N < column->N_scaled * DBL_EPSILON)
                column->vmr_scaled = 1.0 / DBL_EPSILON;
            else
                column->vmr_scaled = column->N_scaled / N;
        } else if (column->vmr_stat & VMR_USER_DEFINED) {
            column->N = N * unmap_variable(column->xvmr, MAPPING_VMR);
            column->N_scaled = column->vmr_scaled * N;
        } else {
            double vmr = vmr_bal * unmap_variable(column->xvmr, MAPPING_VMR);
            column->xvmr = map_variable(vmr, MAPPING_VMR);
            column->vmr_scaled *= vmr_bal;
            column->N = vmr * N;
            column->N_scaled = column->vmr_scaled * N;
        }
    }
    return;
}   /* compute_hydrostatic_layer_column_densities() */


/***********************************************************
* static void compute_nonhydrostatic_layer_column_densities(layer_t *layer)
*
* Purpose:
*   For nonhydrostatic models, compute column densities on a
*   layer.  Mixing ratios are normalized by scaling default
*   mixing ratios; a warning is set if this is not possible.
*
* Arguments:
*   layer_t *layer - pointer to layer structure
************************************************************/

static void compute_nonhydrostatic_layer_column_densities(layer_t *layer)
{
    int cnum;

    /*
     * If h is not defined on this layer, then only explicit
     * column density definitions are possible.  Mixing ratios
     * for self broadening are assumed to be the default dry air
     * values.  Set an advisory flag if the default mixing ratio
     * is zero, or if the default ratio is being assumed on a
     * scaled column.
     */
    if (layer->h < 0.0) {
        for (cnum = 0; cnum < layer->ncols; ++cnum) {
            column_t *column = layer->column[cnum];
            int col_typenum = column->col_typenum;
            if (col_type[col_typenum].flags & COL_PARAMETRIC) {
                continue;
            } else if (col_type[col_typenum].default_vmr == 0.0) {
                column->vmr_scaled = 0.0;
                column->vmr_stat |= VMR_ZERO_DEFAULT;
            } else {
                column->vmr_scaled = col_type[col_typenum].default_vmr;
                if ((find_Nscale_list_entry(col_typenum, layer->tagnum)
                    != NULL) ||
                    (find_Nscale_list_entry(col_typenum, 0) != NULL)) {
                    column->vmr_stat |= VMR_DEFAULT_ON_SCALED_COLUMN;
                }
            }
        }
    } else { /* h is defined on this layer */
        double N;
        double vmr_bal, user_vmr_sum = 0.0;
        N = (layer->P / P_STP) * (T_STP / layer-> T) * N_STP * layer->h;
        for (cnum = 0; cnum < layer->ncols; ++cnum) {
            column_t *column = layer->column[cnum];
            int col_typenum = column->col_typenum;
            if (col_type[col_typenum].flags & COL_PARAMETRIC) {
                continue;
            } else if (column->N_mode & N_EXPLICIT) {
                double vmr;
                if (column->N <= N * DBL_EPSILON)
                    vmr = 0.0;
                else if (N < column->N * DBL_EPSILON)
                    vmr = 1.0 / DBL_EPSILON;
                else
                    vmr = column->N / N;
                column->xvmr = map_variable(vmr, MAPPING_VMR);
                if (column->N_scaled <= N * DBL_EPSILON)
                    column->vmr_scaled = 0.0;
                else if (N < column->N_scaled * DBL_EPSILON)
                    column->vmr_scaled = 1.0 / DBL_EPSILON;
                else
                    column->vmr_scaled = column->N_scaled / N;
                user_vmr_sum += column->vmr_scaled;
            } else if (column->vmr_stat & VMR_USER_DEFINED) {
                user_vmr_sum += column->vmr_scaled;
                column->N_scaled = N * column->vmr_scaled;
                column->N = N * unmap_variable(column->xvmr, MAPPING_VMR);
            } else {
                continue;
            }
        }
        if (user_vmr_sum > 1.0) {
            double vmr_sum_err = user_vmr_sum - 1.0;
            if (vmr_sum_err > VMR_SUM_WARNING_TOL) {
                if (!(layer->vmr_stat & VMR_SUM_EXCEEDS_UNITY))
                    errlog(101, 0);
                layer->vmr_stat |= VMR_SUM_EXCEEDS_UNITY;
            }
            vmr_bal = 0.0;
        } else {
            vmr_bal = 1.0 - user_vmr_sum;
        }
        for (cnum = 0; cnum < layer->ncols; ++cnum) {
            column_t *column = layer->column[cnum];
            int col_typenum = column->col_typenum;
            if (col_type[col_typenum].flags & COL_PARAMETRIC) {
                continue;
            } else if (column->N_mode & N_EXPLICIT) {
                continue;
            } else if (column->vmr_stat & VMR_USER_DEFINED) {
                continue;
            } else {
                double vmr =
                    vmr_bal * unmap_variable(column->xvmr, MAPPING_VMR);
                column->xvmr = map_variable(vmr, MAPPING_VMR);
                column->vmr_scaled *= vmr_bal;
                column->N = vmr * N;
                column->N_scaled = column->vmr_scaled * N;
            }
        }
    }
    return;
} /* compute_nonhydrostatic_layer_column_densities() */


/***********************************************************
* static int get_parent_lnum(const model_t *model, const int clnum)
*
* Purpose:
*   Finds the parent model definition layer associated with
*   layer clnum.  For a layer associated with an
*   interpolated level this is first model definition layer
*   with lnum >= clnum.  (Note that this means that a model
*   definition layer is its own parent layer.)  For layer
*   associated with a level extrapolated below the model
*   definition layers, the parent layer is the lowest model
*   definition layer.
*
* Arguments:
*   const model_t *model - pointer to model data structure
*   const int clnum       - child layer number
*
* Return:
*   Parent layer number
*   If there are no model layers, or if clnum is outside
*   the range of model layers, returns -1
************************************************************/

static int get_parent_lnum(const model_t *model, const int clnum)
{
    int plnum;
    if (clnum < 0 || clnum >= model->nlayers)
        return -1;
    for (plnum = clnum; plnum < model->nlayers; ++plnum) {
        if (model->layer[plnum]->type == LAYER_TYPE_DEF)
            return plnum;
    }
    for (plnum = clnum; plnum >= 0; --plnum) {
        if (model->layer[plnum]->type == LAYER_TYPE_DEF)
            return plnum;
    }
    return plnum;
}   /* get_parent_lnum() */


/***********************************************************
* int get_lnum_by_type(const model_t *model, const int type)
*
* Purpose:
*   Seek for the first model layer of a given type, starting
*   from the base of the layer stack (highest index).
*
* Arguments:
*   const model_t *model - pointer to model data structure
*   const int type       - layer type, defined in am_types.h
*
* Return:
*   Returns the index number of the highest-index model
*   layer having the specified type.  If no such layer
*   exists, returns -1.
************************************************************/

int get_lnum_by_type(const model_t *model, const int type)
{
    int lnum;
    for (lnum = model->nlayers - 1; lnum >= 0; --lnum) {
        if (model->layer[lnum]->type == type)
            break;
    }
    return lnum;
}   /* get_lnum_by_type() */


/***********************************************************
* int get_obs_lnum(model_t *model)
*
* Purpose:
*   Returns the layer number of the observing layer (the
*   layer bounded below by the observing level) in a
*   hydrostatic model.  This is either a user-defined
*   interpolated (or extrapolated) layer, or the default
*   (last model layer).  Returns -1 if there are no layers
*   at all.
*
* Arguments:
*   model_t *model - pointer to model data structure
************************************************************/

int get_obs_lnum(model_t *model)
{
    int lnum;
    if ((lnum = get_lnum_by_type(model, LAYER_TYPE_OBS)) >= 0)
        return lnum;
    else
        return get_lnum_by_type(model, LAYER_TYPE_DEF);
}   /* get_obs_lnum() */


/***********************************************************
* int get_source_lnum(model_t *model)
*
* Purpose:
*   Returns the layer number of the source layer (the layer
*   bounded below by the source level) in a hydrostatic
*   model, if one has been defined. Returns -1 if no source
*   level has been defined, meaning that the source is at
*   infinity, at Psource = 0.
*
* Arguments:
*   model_t *model - pointer to model data structure
************************************************************/

int get_source_lnum(model_t *model)
{
    int lnum;
    if ((lnum = get_lnum_by_type(model, LAYER_TYPE_SOURCE)) >= 0)
        return lnum;
    else
        return -1;
}   /* get_source_lnum() */


/***********************************************************
* int get_tan_lnum(model_t *model)
*
* Purpose:
*   Returns the layer number of the tangent layer in a
*   hydrostatic model if the optical path includes a tangent
*   point, otherwise returns -1.
*
* Arguments:
*   model_t *model - pointer to model data structure
************************************************************/

int get_tan_lnum(model_t *model)
{
    int lnum;
    if ((lnum = get_lnum_by_type(model, LAYER_TYPE_TAN)) >= 0)
        return lnum;
    else if (model->geometry & GEOMETRY_LIMB)
        /*
         * In limb mode, the tangent level is the lowest model
         * definition level.  get_lnum_by_type() searches from
         * the bottom up.
         */
        return get_lnum_by_type(model, LAYER_TYPE_DEF);
    else
        return -1;
}   /* get_tan_lnum() */


/***********************************************************
* int insert_interpolated_level(
*         model_t *model,
*         model_t *lmodel,
*         const int type,
*         const double zbase,
*         const double Pbase,
*         const double nbase)
*
* Purpose:
*   Inserts or updates an interpolated level, by inserting
*   or updating a layer bounded below by that level.  The
*   new layer is positioned in the layer table at the lowest
*   index consistent with height order, and the base
*   pressure Pbase and refractive index nbase are set on the
*   layer.  (Pbase may be zero for layers outside the
*   atmosphere.)  This function assumes that zbase is up to
*   date on all model definition layers.
*
*   The layer type argument is one of:
*
*     LAYER_TYPE_OBS    - layer bounded below by an
*                         interpolated observing level
*
*     LAYER_TYPE_SOURCE - layer bounded below by an
*                         interpolated source level
*
*     LAYER_TYPE_TAN    - layer bounded below by an
*                         interpolated tangent level
*
*   No more than one inserted layer each of these types is
*   allowed.  If an existing layer of the same type exists,
*   and is already in the correct position in the layer
*   stack, then zbase and Pbase are simply updated (on
*   model, but not lmodel).  If an existing layer is found
*   but it is out of position, it is deleted; a new layer is
*   then inserted with its base level at the correct
*   position (on both model and lmodel).
*
*   The model definition layer on which the user-defined
*   level is interpolated, or from which it is extrapolated,
*   is called the parent layer.  A newly-inserted layer
*   inherits the column structure of the parent layer,
*   unless its base pressure Pbase == 0.0.
*
*   Setting (on *model, but not *lmodel) of all other layer
*   and column variables derived by copying, splitting,
*   interpolating, or extrapolating parent layer variables
*   is deferred until all interpolated or extrapolated
*   levels and their associated layers have been inserted or
*   updated.
*
* Arguments:
*   model_t *model     - pointer to model data structure
*   model_t *lmodel    - pointer to a model data structure
*                        containing scalar data from a
*                        prior computation.
*   const int type     - layer type
*   const double zbase - base height of new layer
*   const double Pbase - base pressure of new layer
*   const double nbase - base refractive index of new layer
*
* Return:
*   1 if there were errors, 0 otherwise.
************************************************************/

int insert_interpolated_level(
        model_t *model,
        model_t *lmodel,
        const int type,
        const double zbase,
        const double Pbase,
        const double nbase)
{
    int lnum;
    layer_t *layer  = NULL;
    layer_t *llayer = NULL;

    /*
     * Validate the layer type.
     */
    if (type <= LAYER_TYPE_NONE || type >= LAYER_TYPE_END ||
            type == LAYER_TYPE_DEF) {
        errlog(186, type);
        return 1;
    }
    /*
     * Look for an existing layer of the same type.
     */
    lnum = get_lnum_by_type(model, type);
    if (lnum > 0) {
        /* A layer of the same type already exists. */
        if (
                (lnum == 0 ||
                    model->layer[lnum - 1]->zbase > zbase) &&
                (lnum == model->nlayers - 1 ||
                    model->layer[lnum + 1]->zbase <= zbase)
                ) {
            /*
             * The existing layer is already in correct
             * order in the layer table.  Update zbase and
             * Pbase for the layer in *model (but not
             * *lmodel), and return.
             */
            layer        = model->layer[lnum];
            layer->zbase = zbase;
            layer->Pbase = Pbase;
            return 0;
        } else {        
            /*
             * The existing layer is no longer in correct
             * position.  Delete it from both model and
             * lmodel.  It will be re-inserted as a new
             * layer below.  Return an error if layer deletion
             * fails.
             */
            if (delete_layer(model, lnum) ||
                    (lmodel != NULL && delete_layer(lmodel, lnum))) {
                errlog(187, lnum);
                return 1;
            }
        }
    }
    /*
     * Insert a new layer into model and lmodel.  Initially,
     * add_layer() will place the layer at the end of the table
     * of layer pointers.  Make a copy of the pointer, then shift
     * table entries to re-insert it in correct order.
     */
    if (add_layer(model, type))
        return 1;
    layer = model->layer[model->nlayers - 1];
    if (lmodel != NULL) {
        if (add_layer(lmodel, type))
            return 1;
        llayer = lmodel->layer[lmodel->nlayers - 1];
    }
    for (lnum = model->nlayers - 1; lnum > 0; --lnum) {
        if (model->layer[lnum - 1]->zbase <= zbase) {
            model->layer[lnum] = model->layer[lnum - 1];
            if (lmodel != NULL)
                lmodel->layer[lnum] = lmodel->layer[lnum - 1];
        } else {
            break;
        }
    }
    model->layer[lnum] = layer;
    if (lmodel != NULL)
        lmodel->layer[lnum] = llayer;
    /*
     * On the new layer, set zbase, Pbase, and nbase on *model,
     * but not on *lmodel.
     */
    layer->zbase = zbase;
    layer->Pbase = Pbase;
    layer->nbase = nbase;
    /*
     * Set the layer dimensions and array allocations to match
     * the parent layer, unless this is a zero-pressure layer.
     * Zero-pressure layers get no added columns and a basic
     * array allocation.
     */
    if (layer->Pbase > 0.0) {
        int plnum = get_parent_lnum(model, lnum);
        if (plnum >= 0) {
            if (copy_layer_dimensions(model->layer[plnum], layer) ||
                (lmodel != NULL &&
                 copy_layer_dimensions(lmodel->layer[plnum], llayer))) {
                return 1;
            }
            if (copy_layer_allocations(model, model->layer[plnum], layer))
                return 1;
        }
    } else {
        alloc_layer_arrays(model, lnum);
    }
    return 0;
}   /* insert_interpolated_level() */


/***********************************************************
* static double layer_base_dry_lapse_rate(
*   const model_t *model, const int lnum)
*
* Purpose:
*   Computes the dry adiabatic lapse rate -dT/dz in K/cm at
*   the base of a hydrostatic layer.  In midpoint T mode, or
*   for a nonhydrostatic layer, a lapse rate of zero is
*   returned
*
*   This function assumes that Mair and gbase have already
*   been set on the layer by a call to
*   set_def_layer_heights_and_col_densities().
*
* Arguments:
*   model_t *model - pointer to model data structure
*   const int lnum - index of layer.
*
* Return:
*   dry adiabatic lapse rate [K / cm]
************************************************************/

static double layer_base_dry_lapse_rate(
        const model_t *model, const int lnum)
{
    double R, cp;
    layer_t *layer = model->layer[lnum];

    if (   !(model->PTmode & PTMODE_HYDROSTATIC) ||
            (model->PTmode & PTMODE_T)) {
        return 0.0;
    }
    R = PCONST_KB / (layer->M0 * PCONST_AMU);
    /*
     * Here we assume the atmospheric composition is dominated by
     * diatomic molecules.
     */
    cp = (7./2.) * R;
    /*
     * Convert cp  from [m^2 s^-2 K^-1] to [cm^2 s^2 K^-1] to
     * compute g / cp in native am units [K cm^-1].
     */
    cp *= unit_tab[UNIT_M].factor * unit_tab[UNIT_M].factor;
    return layer->gbase / cp;
}   /* layer_base_dry_lapse_rate() */


/***********************************************************
* double layer_base_lapse_rate(
*   const model_t *model, const int lnum)
*
* Purpose:
*   Returns the lapse rate -dT/dz at the base of layer lnum,
*   to be used for computing the position of extrapolated
*   model levels.  The lapse rate is computed according to
*   the model base extrapolation mode setting.
*
*   This function assumes that Mair and gbase have already
*   been set on the layer by a call to
*   set_def_layer_heights_and_col_densities().
*
* Arguments:
*   model_t *model - pointer to model data structure
*   const int lnum - index of layer.
*
* Return:
*   lapse rate -dT/dz [K/cm]
************************************************************/

double layer_base_lapse_rate(
        const model_t *model, const int lnum)
{
    double gamma = 0.0;

    if (model->PTmode & PTMODE_T)
        gamma = 0.0;
    else if (model->PTmode & PTMODE_EXTEND_ISOTHERMAL)
        gamma = 0.0;
    else if (model->PTmode & PTMODE_EXTEND_DRY_ADIABATIC)
        gamma = layer_base_dry_lapse_rate(model, lnum);
    else if (model->PTmode & PTMODE_EXTEND_MOIST_ADIABATIC)
        gamma = layer_base_moist_lapse_rate(model, lnum);
    else if (model->PTmode & PTMODE_EXTEND_ENVIRONMENTAL)
        gamma = layer_environmental_lapse_rate(model, lnum);
    else
        errlog(206, model->PTmode);
    return gamma;
}   /* layer_base_lapse_rate() */


/***********************************************************
* static double layer_base_moist_lapse_rate(
*   const model_t *model, const int lnum)
*
* Purpose:
*   Computes the moist adiabatic lapse rate -dT/dz in K/cm at
*   at the base of a hydrostatic layer.  In midpoint T mode,
*   or for a nonhydrostatic layer, a lapse rate of zero is
*   returned
*
*   This function assumes that Mair and gbase have already
*   been set on the layer by a call to
*   set_def_layer_heights_and_col_densities().
*
* Arguments:
*   model_t *model - pointer to model data structure
*   const int lnum - index of layer.
*
* Return:
*   moist adiabatic lapse rate [K / cm]
************************************************************/

static double layer_base_moist_lapse_rate(
        const model_t *model, const int lnum)
{
    double cp, L, mu_s, R, R_v, T, x;
    layer_t *layer = model->layer[lnum];

    if (   !(model->PTmode & PTMODE_HYDROSTATIC) ||
            (model->PTmode & PTMODE_T)) {
        return 0.0;
    }
    T = layer->Tbase;
    /*
     * Here we ignore the temperature dependence of the latent
     * heat of vaporization is ignored, and approximate it by its
     * value at 0 C.
     */
    L = H2O_LATENT_HEAT_VAP_0C;
    /*
     * mu_s is the saturated mass mixing ratio.
     */
    mu_s  = H2O_liquid_Psat(T, COL_TYPE_H2O) / layer->Pbase;
    mu_s *= MASS_H2O / layer->Mair;
    R     = PCONST_KB / (layer->M0 * PCONST_AMU);
    R_v   = PCONST_KB / (MASS_H2O  * PCONST_AMU);
    /*
     * Here we assume the atmospheric composition is dominated by
     * diatomic molecules.
     */
    cp = (7./2.) * R;
    x  = 1.0 + (L * mu_s) / (R * T);
    x /= 1.0 + (L * L * mu_s) / (cp * R_v * T * T);
    /*
     * x is now the dimensionless ratio of the moist adiabatic
     * lapse rate to the dry adiabatic lapse rate g / cp.
     * Convert cp  from [m^2 s^-2 K^-1] to [cm^2 s^2 K^-1] to
     * compute g / cp in native am units [K cm^-1].
     */
    cp *= unit_tab[UNIT_M].factor * unit_tab[UNIT_M].factor;
    return x * layer->gbase / cp;
}   /* layer_base_moist_lapse_rate() */


/***********************************************************
* static double layer_environmental_lapse_rate(
*   const model_t *model, int lnum)
*
* Purpose:
*   In Tbase mode, computes the environmental lapse rate
*   -dT/dZ across a layer in am native units [K /cm].
*
*   On the top model layer, which is always isothermal, or
*   in midpoint T mode, this function returns zero.
*
*   This function assumes that zbase has been set on model
*   definition layers by a call to 
*   set_def_layer_heights_and_col_densities().
*
* Arguments:
*   model_t *model - pointer to model data structure
*   const int lnum - index of layer.
*
* Return:
*   environmental lapse rate [K / cm]
************************************************************/

static double layer_environmental_lapse_rate(
        const model_t *model, int lnum)
{
    layer_t *layer = model->layer[lnum];
    double T_lo = layer->Tbase;
    double z_lo = layer->zbase;

    if (   !(model->PTmode & PTMODE_HYDROSTATIC) ||
            (model->PTmode & PTMODE_T)) {
        return 0.0;
    }
    /*
     * Seek upwards for another model definition level with
     * nonzero height difference from this layer's base level.
     * If found, compute the lapse rate between levels and return
     * the result.
     */
    while (--lnum >= 0) {
        layer = model->layer[lnum];
        if (layer->type == LAYER_TYPE_DEF) {
            double T_hi = layer->Tbase;
            double dz   = layer->zbase - z_lo;
            if (dz > DBL_EPSILON)
                return (T_lo - T_hi) / dz;
        }
    }
    /*
     * Not found, so this was the highest model definition
     * layer, which is isothermal.
     */
    return 0.0;
}   /* layer_environmental_lapse_rate() */


/***********************************************************
* static double layer_h2o_vmr(const model_t *model, const int lnum)
*
* Purpose:
*   Computes the total h2o mixing ratio on a layer, for use
*   in computing the refractive index.  This function
*   assumes that scaled mixing ratios have already been
*   computed on all layers.
*
* Arguments:
*   model_t *model - pointer to model data structure
*   const int lnum - index of layer.
*
* Return:
*   total h2o mixing ratio
************************************************************/

static double layer_h2o_vmr(const model_t *model, const int lnum)
{
    double h2o_vmr = 0.0;
    layer_t *layer = model->layer[lnum];
    int cnum;

    /*
     * The total h2o vmr is taken to be the sum of vmr's for all
     * column types that include H2_16O.
     */
    for (cnum = 0; cnum < layer->ncols; ++cnum) {
        column_t *column = layer->column[cnum];
        switch (column->col_typenum) {
        case COL_TYPE_H2O:
        case COL_TYPE_H2_16O_LINES_PLUS_CONTINUUM:
        case COL_TYPE_H2O_OPTICAL_REFRACTIVITY:
            h2o_vmr += column->vmr_scaled;
            break;
        default:
            break;
        }
    }
    return h2o_vmr;
}   /* layer_h2o_vmr() */


/***********************************************************
* void reset_layer_performance_timers(model_t *model)
*
* Purpose:
*   Resets layer and column runtimes to zero.
************************************************************/

void reset_layer_performance_timers(model_t *model)
{
    int lnum;
    for (lnum = 0; lnum < model->nlayers; ++lnum) {
        layer_t *layer = model->layer[lnum];
        int cnum;
        layer->runtime = 0.0;
        for (cnum = 0; cnum < layer->ncols; ++cnum)
            layer->column[cnum]->runtime = 0.0;
    }
    return;
}   /* reset_layer_performance_timers() */


/***********************************************************
* void restore_def_layer_explicit_col_densities(model_t *model)
*
* Purpose:
*   
*   Restores explicity-defined column densities N on
*   model definition layers.
*
*   The defining value of an explicit column density may be
*   split to accomodate interpolated levels and their
*   associated layers.  This function restores the defining
*   value on subsequent model computation passes.  This is
*   needed because the defining value may be a fit or
*   differentiation variable subject to change between model
*   computations.
*
* Arguments:
*   model_t *model - pointer to model structure
************************************************************/

void restore_def_layer_explicit_col_densities(model_t *model)
{
    int lnum;
    for (lnum = 0; lnum < model->nlayers; ++lnum) {
        layer_t *layer = model->layer[lnum];
        if (layer->type == LAYER_TYPE_DEF) {
            int cnum;
            for (cnum = 0; cnum < layer->ncols; ++cnum) {
                column_t *column = layer->column[cnum];
                    column->N = column->N_def;
            }
        }
    }
    return;
}   /* restore_def_layer_explicit_col_densities() */


/***********************************************************
* void set_def_layer_heights_and_col_densities(model_t *model)
*
* Purpose:
*  Computes zenith column densities for all layers.  For
*  hydrostatic models, layer heights and gravitational
*  acceleration at layer base levels are also computed.
*  This function assumes that layer pressure and temperature
*  variables are already up to date.
*
* Arguments:
*   model_t *model - pointer to model structure
************************************************************/

void set_def_layer_heights_and_col_densities(model_t *model)
{
    if (model->PTmode & PTMODE_HYDROSTATIC) {
        int lnum;
        layer_t* layer_ref;
        /*
         * In hydrostatic models, the initial height reference is
         * the base level of the lowest user-defined model level.
         * Set the inital zbase and gbase for this layer, and
         * compute its column densities and mean molecular mass.
         */
        if ((lnum = get_lnum_by_type(model, LAYER_TYPE_DEF)) < 0)
            return; /* no layers */
        layer_ref = model->layer[lnum];
        layer_ref->zbase = model->z0;
        /*
         * Gravity adjustment is always applied in hydrostatic
         * models, though the default gradient dg_dz is zero.
         * Our conventions are that model->g is the gravitational
         * acceleration at height z = 0, and that a user-defined
         * gradient dg_dz will normally have negative sign.
         */
        layer_ref->gbase  = model->g;
        layer_ref->gbase += model->dg_dz * layer_ref->zbase;
        if (layer_ref->gbase < 0.0)
            layer_ref->gbase = 0.0;
        if (layer_ref->h >= 0.0) /* plume or instrument layer, see below */
            compute_nonhydrostatic_layer_column_densities(layer_ref);
        else
            compute_hydrostatic_layer_column_densities(layer_ref);
        /*
         * Work upwards through the model definition layers,
         * computing layer height increments and gravity.
         */
        while (--lnum >= 0) {
            layer_t *layer = model->layer[lnum];
            double dz;
            if (layer->type == LAYER_TYPE_DEF) {
                if (layer->Pbase <= 0.0) {
                    /*
                     * Model definition layers with Pbase == 0 are
                     * placed at infinity.
                     */
                    layer->zbase = DBL_MAX;
                    layer->gbase = 0.0;
                    layer->Mair  = 0.0;
                    continue;
                } else if (layer->h >= 0.0) {
                    /*
                     * In hydrostatic models, a non-hydrostatic layer
                     * can be embedded in the model by assigning a
                     * thickness h to a layer with dP == 0.  This
                     * facility is used for modeling things like
                     * plumes or instrument layers.  Such a layer is
                     * always assumed to have negligible height dz
                     * irrespective of h.
                     */
                    layer->gbase = layer_ref->gbase;
                    layer->zbase = layer_ref->zbase;
                    layer->Mair  = layer->M0;
                    dz = 0.0;
                    compute_nonhydrostatic_layer_column_densities(layer);
                } else {
                    /*
                     * For a normal hydrostatic layer, the base level
                     * height is computed relative to the current
                     * reference layer.  When computing the
                     * hydrostatic layer thickness, the gravitational
                     * acceleration gbase at the layer base level is
                     * used rather than a midpoint value.  This is
                     * normally a very accurate approximation, and
                     * avoids an iterative or extrapolated solution.
                     */
                    double R, Tbar;
                    double P_ref = layer_ref->Pbase;
                    double P     = layer->Pbase;
                    if (model->PTmode & PTMODE_TBASE) {
                        /*
                         * Assuming T is linear in log(P), Tbar is the
                         * mean of T_ref and T.  See the discussion on
                         * layer thickness in Chapter 1 of the manual.
                         */
                        double T_ref = layer_ref->Tbase;
                        double T     = layer->Tbase;
                        Tbar = 0.5 * (T_ref + T);
                    } else {
                        /*
                         * In midpoint T mode, model layers are
                         * isothermal.
                         */
                        Tbar = layer_ref->Tbase;
                    }
                    R   = PCONST_KB / (layer_ref->Mair * PCONST_AMU);
                    R  *= unit_tab[UNIT_M].factor * unit_tab[UNIT_M].factor;
                    dz  = R * Tbar / layer_ref->gbase;
                    dz *= log(P_ref / P);
                    layer->zbase  = layer_ref->zbase + dz;
                    layer->gbase  = model->g;
                    layer->gbase += model->dg_dz * layer->zbase;
                    if (layer->gbase < 0.0)
                        layer->gbase = 0.0;
                    compute_hydrostatic_layer_column_densities(layer);
                }
                layer_ref = layer;
            }
        }
    } else { /* Nonhydrostatic models */
        int lnum;
        for (lnum = 0; lnum < model->nlayers; ++lnum)
            compute_nonhydrostatic_layer_column_densities(model->layer[lnum]);
    }
    return;
}   /* set_def_layer_heights_and_col_densities() */


/***********************************************************
* void set_def_layer_pressures(model_t *model)
*
* Purpose:
*   Sets additional pressure variables derived from the
*   user-specified ones on all model definition layers.
*
* Arguments:
*   model_t *model - pointer to model structure
************************************************************/

void set_def_layer_pressures(model_t *model)
{
    /*
     * In each pair of loops below, the first loop sets P = Pbase on
     * the first model definition layer.
     */
    if (model->PTmode & PTMODE_DP) {
        int lnum;
        layer_t *llayer = NULL;
        for (lnum = 0; lnum < model->nlayers; ++lnum) {
            llayer = model->layer[lnum];
            if (llayer->type == LAYER_TYPE_DEF) {
                llayer->Pbase = llayer->dP_def;
                llayer->P     = llayer->dP_def;
                break;
            }
        }
        for (++lnum; lnum < model->nlayers; ++lnum) {
            layer_t *layer = model->layer[lnum];
            if (layer->type == LAYER_TYPE_DEF) {
                layer->Pbase = llayer->Pbase + layer->dP_def;
                layer->P     = 0.5 * (llayer->Pbase + layer->Pbase);
                llayer = layer;
            }
        }
    } else if (model->PTmode & PTMODE_PBASE) {
        /*
         * In PTMODE_PBASE, set dP_def on model definition
         * layers as the defining value of dP, for later use
         * if the layer is split by an interpolated level.
         */
        int lnum;
        layer_t *llayer = NULL;
        for (lnum = 0; lnum < model->nlayers; ++lnum) {
            llayer = model->layer[lnum];
            if (llayer->type == LAYER_TYPE_DEF) {
                llayer->dP     = llayer->Pbase;
                llayer->dP_def = llayer->dP;
                llayer->P      = llayer->Pbase;
                break;
            }
        }
        for (++lnum; lnum < model->nlayers; ++lnum) {
            layer_t *layer = model->layer[lnum];
            if (layer->type == LAYER_TYPE_DEF) {
                layer->dP     = layer->Pbase - llayer->Pbase;
                layer->dP_def = layer->dP;
                layer->P      = 0.5 * (llayer->Pbase + layer->Pbase);
                if (fabs(layer->dP) < DBL_EPSILON * layer->P)
                    layer->dP = 0.0;
                llayer = layer;
            }
        }
    }
    return;
} /* set_def_layer_pressures() */


/***********************************************************
* void set_def_layer_refractive_indices(model_t *model)
*
* Purpose:
*   On all model definition layers, set the layer base
*   refractive index.  The refractive index profile is
*   anchored to the base levels of the model definition
*   layers, and interpolated linearly in density (P/T)
*   elsewhere.
*
*   This function assumes that pressures, temperatures, and
*   scaled mixing ratios have already been set on the model
*   definition layers.
*
* Arguments:
*   model_t *model - pointer to model data structure
************************************************************/

void set_def_layer_refractive_indices(model_t *model)
{
    int lnum;

    for (lnum = 0; lnum < model->nlayers; ++lnum) {
        layer_t *layer = model->layer[lnum];
        if (layer->type == LAYER_TYPE_DEF) {
            double P, Pdry, Ph2o;
            double T;
            double h2o_vmr;
            double n;
            if (model->PTmode & PTMODE_HYDROSTATIC)
                P = layer->Pbase;
            else
                P = layer->P;
            if (model->PTmode & PTMODE_TBASE)
                T = layer->Tbase;
            else
                T = layer->T;
            h2o_vmr = layer_h2o_vmr(model, lnum);
            Pdry    = P * (1.0 - h2o_vmr);
            Ph2o    = P * h2o_vmr;
            n       = 1.0;
            if (model->geometry & GEOMETRY_REFRACT_OPTICAL) {
                n +=    (  
                        OPTICAL_REFRACTIVITY_K1 * Pdry +
                        OPTICAL_REFRACTIVITY_K2 * Ph2o
                        ) / T;
            }
            if (model->geometry & GEOMETRY_REFRACT_RADIO) {
                n +=    (
                        RADIO_REFRACTIVITY_K1 * Pdry +
                        RADIO_REFRACTIVITY_K2 * Ph2o +
                        RADIO_REFRACTIVITY_K3 * Ph2o / T
                        ) / T;
            }
            layer->nbase = n;
        }
    }
    return;
}   /* set_def_layer_refractive_indices() */


/***********************************************************
* void set_def_layer_temperatures(model_t *model)
*
* Purpose:
*   When model definition layer temperatures are set by base
*   temperature, this function sets interpolated midpoint
*   temperatures.  Interpolation is linear in log(P) for
*   hydrostatic models, and simply linear otherwise.
*
*   When model definition layer temperatures are set by
*   midpoint temperature, layers are assumed isothermal, and
*   base temperatures are set to match the midpoint
*   temperature.  This convention avoids the need to check
*   the PTmode in some cases.
*
*   This function should be preceded by a call to
*   set_def_layer_pressures() to ensure that base pressures
*   are up to date.
*
* Arguments:
*   model_t *model - pointer to model structure
************************************************************/

void set_def_layer_temperatures(model_t *model)
{
    if (model->PTmode & PTMODE_TBASE) {
        int lnum;
        layer_t *llayer = NULL;
        /*
         * On the top model definition layer, the midpoint
         * temperature is set equal to the base temperature.
         */
        for (lnum = 0; lnum < model->nlayers; ++lnum) {
            llayer = model->layer[lnum];
            if (llayer->type == LAYER_TYPE_DEF) {
                llayer->T = llayer->Tbase;
                break;
            }
        }
        if (model->PTmode & PTMODE_HYDROSTATIC) {
            while (++lnum < model->nlayers) {
                layer_t *layer = model->layer[lnum];
                if (layer->type == LAYER_TYPE_DEF) {
                    layer->T = log_x_interp(
                            llayer->Pbase,
                            llayer->Tbase,
                            layer->Pbase,
                            layer->Tbase,
                            layer->P);
                    llayer = layer;
                }
            }
        } else {
            while (++lnum < model->nlayers) {
                layer_t *layer = model->layer[lnum];
                if (layer->type == LAYER_TYPE_DEF) {
                    layer->T = 0.5 * (llayer->Tbase + layer->Tbase);
                    llayer = layer;
                }
            }
        }
    } else {
        int lnum;
        for (lnum = 0; lnum < model->nlayers; ++lnum) {
            layer_t *layer = model->layer[lnum];
            if (layer->type == LAYER_TYPE_DEF)
                layer->Tbase = layer->T;
        }
    }
    return;
} /* set_def_layer_temperatures() */


/***********************************************************
* void set_def_layer_vmrs_and_apply_Nscale_factors(
*   model_t *model)
*
* Purpose:
*   On all model definition layers, this function
*   initializes user-defined and default volume mixing
*   ratios, including user-defined mixing ratios defined by
*   relative humidity.  Any applicable Nscale factors are
*   then applied to these mixing ratios and to user-defined
*   explicit column densities.
*
*   Note that the default mixing ratios applied in this
*   function are subject to subsequent scaling to normalize
*   the volume mixing ratios to unity in calls to either
*
*       compute_hydrostatic_layer_column_densities()
*   or
*       compute_nonhydrostatic_layer_column_densities()
*
* Arguments:
*   model_t *model - pointer to model structure
************************************************************/

void set_def_layer_vmrs_and_apply_Nscale_factors(
    model_t *model)
{
    int lnum;

    for (lnum = 0; lnum < model->nlayers; ++lnum) {
        layer_t *layer = model->layer[lnum];
        if (layer->type == LAYER_TYPE_DEF) {
            int cnum;
            for (cnum = 0; cnum < layer->ncols; ++cnum) {
                double Nscale;
                column_t *column = layer->column[cnum];
                int col_typenum = column->col_typenum;
                /*
                 * Clear vmr error flags
                 */
                layer->column[cnum]->vmr_stat &= ~(VMR_ALL_ERRS);
                /*
                 * Compute any vmrs defined by relative humidity
                 */
                if (column->vmr_stat & VMR_BY_RH) {
                    double Psat, vmr;
                    if (column->vmr_stat & VMR_BY_RH_LIQUID)
                        Psat = H2O_liquid_Psat(layer->T, col_typenum);
                    else
                        Psat = H2O_ice_Psat(layer->T, col_typenum);
                    column->RH_adj  = column->RH -
                        unmap_variable(model->RH_offset_exp, MAPPING_EXP);
                    column->RH_adj *= model->RH_scale;
                    if (column->RH_adj < 0.0) {
                        column->RH_adj = 0.0;
                        errlog(118, 0);
                    }
                    vmr = 0.01 * column->RH_adj * (Psat / layer->P);
                    if (vmr > 1.0) {
                        errlog(122, 0);
                        vmr = 1.0;
                    }
                    column->xvmr = map_variable(vmr, MAPPING_VMR);
                }
                /*
                 * Apply any Nscale factors to mixing ratios and column
                 * densities.
                 */
                Nscale = lookup_Nscale(col_typenum, layer->tagnum);
                if (column->N_mode & N_EXPLICIT) {
                    /*
                     * User-defined column density and parametric column
                     * types
                     */
                    column->N_scaled = Nscale * column->N;
                } else if (column->vmr_stat & VMR_USER_DEFINED) {
                    /*
                     * User-defined mixing ratio
                     */
                    column->vmr_scaled  = unmap_variable(
                            column->xvmr, MAPPING_VMR);
                    column->vmr_scaled *= Nscale;
                } else {
                    /*
                     * Default mixing ratio.
                     *
                     * For columns defined with default mixing ratios,
                     * the default mixing ratio is subject to subsequent
                     * re-normalization to account for user-defined
                     * column densities or mixing ratios.  Here, restore
                     * the default.  Log a warning if Nscale is applied
                     * to a default mixing ratio, since the user may not
                     * have expected this.
                     */
                    double vmr;
                    vmr = col_type[col_typenum].default_vmr;
                    column->xvmr = map_variable(vmr, MAPPING_VMR);
                    column->vmr_scaled = vmr * Nscale;
                    if (
                            (find_Nscale_list_entry(col_typenum, 0) != NULL) ||
                            (find_Nscale_list_entry(col_typenum, layer->tagnum)
                                != NULL )
                            ) {
                        errlog(112,0);
                    }
                }
            }
        }
    }
    return;
}   /* set_def_layer_vmrs_and_apply_Nscale_factors() */


/***********************************************************
* void set_interpolated_layer_variables(model_t *model)
*
* Purpose:
*   In a hydrostatic model, after interpolated levels have
*   been placed among the model definition layers, this
*   function sets all variables on the corresponding layers
*   and the columns contained therein.  It also sets
*   adjusted variables on the associated parent layers that
*   were split to create these added layers.
*
* Arguments:
*   model_t *model - pointer to model structure
*
* Return:
*   1 if there were errors, 0 otherwise.
************************************************************/

int set_interpolated_layer_variables(model_t *model)
{
    int lnum, def_lnum_top, def_lnum_base;
    layer_t *player;

    if (model->nlayers <= 0) {
        /*
         * no layers
         */
        return 0;
    } else if (!(model->PTmode & PTMODE_HYDROSTATIC)) {
        /*
         * not a hydrostatic model
         */
        return 0;
    } else {
        /*
         * Set P and dP on all layers to reflect any changes from
         * added levels.  Here, we take care to follow the rule
         * that the highest layer with non-zero base pressure is
         * isobaric with P = Pbase, unless the layer above is a
         * dummy model definition layer with Pbase = 0, in which
         * case P = Pbase / 2.
         *
         * Because of the check above, we know there is at least
         * one layer.  That's hard for the compiler to figure
         * out, so initialize llayer to NULL.
         */
        layer_t *llayer = NULL;
        for (lnum = 0; lnum < model->nlayers; ++lnum) {
            llayer = model->layer[lnum];
            llayer->dP = llayer->Pbase;
            llayer->P  = llayer->Pbase;
            if (llayer->Pbase > 0.0 || llayer->type == LAYER_TYPE_DEF)
                break;
        }
        for (++lnum; lnum < model->nlayers; ++lnum) {
            layer_t *layer = model->layer[lnum];
            layer->dP = layer->Pbase - llayer->Pbase;
            layer->P  = 0.5 * (layer->Pbase + llayer->Pbase);
            llayer = layer;
        }
    }
    /*
     * Get indices for the model definition layers at the top and
     * the base of the layer stack.  These are used in the
     * temperature setting loops below to set loop ranges
     * covering interpolated layers above, within, and below the
     * range spanned by model definition levels.
     */
    def_lnum_top  = get_parent_lnum(model, 0);
    def_lnum_base = get_parent_lnum(model, model->nlayers - 1);
    /*
     * Set interpolated layer temperatures.
     */
    if (model->PTmode & PTMODE_TBASE) {
        layer_t *pplayer;
        /*
         * In PTMODE_TBASE, temperatures are extrapolated
         * isothermally at the top of the atmosphere,
         * interpolated in log(P) between model definition levels
         * within the atmosphere, and extrapolated according the
         * the base extrapolation mode setting at the bottom of
         * the atmosphere.
         */
        player = model->layer[def_lnum_top];
        for (lnum = 0; lnum < def_lnum_top; ++lnum) {
            layer_t *layer = model->layer[lnum];
            layer->T       = player->Tbase;
            layer->Tbase   = player->Tbase;
        }
        /*
         * pplayer tracks the previous parent layer, for
         * interpolating between model definition levels.
         */
        pplayer = player;
        while (++lnum < def_lnum_base) {
            layer_t *layer = model->layer[lnum];
            if (layer->type == LAYER_TYPE_DEF) {
                pplayer = layer; 
            } else {
                int plnum = get_parent_lnum(model, lnum);
                player    = model->layer[plnum];
                layer->Tbase = log_x_interp(
                        pplayer->Pbase,
                        pplayer->Tbase,
                        player->Pbase,
                        player->Tbase,
                        layer->Pbase);
                layer->T = log_x_interp(
                        pplayer->Pbase,
                        pplayer->Tbase,
                        player->Pbase,
                        player->Tbase,
                        layer->P);
                if (plnum - lnum == 1) {
                    /*
                     * If this is the last interpolated level
                     * within this parent layer, adjust the
                     * midpoint temperature of the parent layer.
                     */
                    player->T = log_x_interp(
                            pplayer->Pbase,
                            pplayer->Tbase,
                            player->Pbase,
                            player->Tbase,
                            player->P);
                }
            }
        }
        player = model->layer[def_lnum_base];
        /*
         * Assign extrapolated layer base temperatures according
         * the the base extrapolation mode setting, then set the
         * midpoint temperature by interpolating in log(P)
         * between the current layer and the base of the layer
         * above (which could be another extrapolated layer).
         */
        while (++lnum < model->nlayers) {
            double gamma;
            layer_t *layer = model->layer[lnum];
            gamma = layer_base_lapse_rate(model, def_lnum_base);
            layer->Tbase = player->Tbase +
                gamma * (player->zbase - layer->zbase);
            layer->T = log_x_interp(
                    model->layer[lnum-1]->Pbase,
                    model->layer[lnum-1]->Tbase,
                    layer->Pbase,
                    layer->Tbase,
                    layer->P);
        }
    } else {
        /*
         * In PTMODE_T, temperatures are simply picked up from
         * the parent model definition layers.
         */
        player = model->layer[def_lnum_top];
        for (lnum = 0; lnum < def_lnum_top; ++lnum) {
            layer_t *layer = model->layer[lnum];
            layer->T       = player->T;
            layer->Tbase   = player->T;
        }
        while (++lnum < def_lnum_base) {
            layer_t *layer = model->layer[lnum];
            if (layer->type != LAYER_TYPE_DEF) {
                int plnum    = get_parent_lnum(model, lnum);
                player       = model->layer[plnum];
                layer->T     = player->T;
                layer->Tbase = player->T;
            }
        }
        player = model->layer[def_lnum_base];
        while (++lnum < model->nlayers) {
            layer_t *layer = model->layer[lnum];
            layer->T       = player->T;
            layer->Tbase   = player->T;
        }
    }
    /*
     * Column densities on model definition layers have all been
     * computed previously, both those defined explicitly and
     * those defined by mixing ratios.  Here, split the column
     * densities on parent layers containing interpolated levels.
     * This includes parametric column types, which behave as
     * absorption distributed uniformly across a layer.
     */
    for (lnum = 0; lnum < def_lnum_base; ++lnum) {
        layer_t *layer = model->layer[lnum];
        if (layer->type != LAYER_TYPE_DEF) {
            int    cnum;
            int    plnum = get_parent_lnum(model, lnum);
            double P_frac;
            player = model->layer[plnum];
            P_frac = layer->dP / player->dP_def;
            for (cnum = 0; cnum < layer->ncols; ++cnum) {
                column_t *column  = layer->column[cnum];
                column_t *pcolumn = player->column[cnum];
                column->N         = P_frac * pcolumn->N;
                column->N_scaled  = P_frac * pcolumn->N_scaled;
                column->N_mode    = pcolumn->N_mode | N_INTERPOLATED;
                column->N_unitnum = pcolumn->N_unitnum;
            }
            if (plnum - lnum == 1) {
                /*
                 * If this is the last interpolated level within
                 * this parent layer, adjust the column densities
                 * on the parent layer.
                 */
                P_frac = player->dP / player->dP_def;
                for (cnum = 0; cnum < player->ncols; ++cnum) {
                    column_t *pcolumn  = player->column[cnum];
                    pcolumn->N         = P_frac * pcolumn->N;
                    pcolumn->N_scaled  = P_frac * pcolumn->N_scaled;
                    pcolumn->N_mode    = pcolumn->N_mode | N_INTERPOLATED;
                }
            }
        }
    }
    /*
     * Next, compute column densities on extrapolated layers.
     * Note that the parent layer column density may have been
     * split if it contained an interpolated level, so here we
     * compute P_frac with dP, not dP_def.
     */
    player = model->layer[def_lnum_base];
    while (++lnum < model->nlayers) {
        int     cnum;
        layer_t *layer = model->layer[lnum];
        double  P_frac = layer->dP / player->dP;
        for (cnum = 0; cnum < layer->ncols; ++cnum) {
            column_t *column  = layer->column[cnum];
            column_t *pcolumn = player->column[cnum];
            column->N         = P_frac * pcolumn->N;
            column->N_scaled  = P_frac * pcolumn->N_scaled;
            column->N_mode    = pcolumn->N_mode | N_INTERPOLATED;
            column->N_unitnum = pcolumn->N_unitnum;
        }
    }
    /*
     * Copy intensive parent layer variables and set base level
     * gravity.
     */
    for (lnum = 0; lnum < model->nlayers; ++lnum) {
        layer_t *layer = model->layer[lnum];
        if (layer->type != LAYER_TYPE_DEF) {
            int i, cnum;
            int plnum = get_parent_lnum(model, lnum);
            player = model->layer[plnum];
            layer->Mair = player->Mair;
            layer->M0   = player->M0;
            for (i = 0; i < K_TYPE_END_OF_TABLE; ++i)
                layer->lineshape[i] = player->lineshape[i];
            for (i = 0; i < K_TYPE_END_OF_TABLE; ++i)
                layer->strict_selfbroad[i] = player->strict_selfbroad[i];
            for (i = 0; i < COL_TYPE_END_OF_TABLE; ++i)
                layer->Mair_flag[i] = player->Mair_flag[i];
            layer->tagnum    = player->tagnum;
            layer->P_unitnum = player->P_unitnum;
            layer->T_unitnum = player->T_unitnum;
            layer->h_unitnum = player->h_unitnum;
            layer->vmr_stat  = player->vmr_stat;
            for (cnum = 0; cnum < layer->ncols; ++cnum) {
                column_t *column   = layer->column[cnum];
                column_t *pcolumn  = player->column[cnum];
                column->xvmr       = pcolumn->xvmr;
                column->vmr_scaled = pcolumn->vmr_scaled;
                column->RH         = pcolumn->RH;
                column->RH_adj     = pcolumn->RH_adj;
                column->vmr_stat   = pcolumn->vmr_stat;
            }
            layer->gbase  = player->gbase;
            layer->gbase += (layer->zbase - player->zbase) * model->dg_dz; 
        }
    }
    return 0;
}   /* set_interpolated_layer_variables() */
