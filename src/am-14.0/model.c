/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* model.c                         S. Paine rev. 2024 July 16
*
* Atmospheric model and radiative transfer computations.
************************************************************/

#include <float.h>
#include <math.h>

#include "am_sysdep.h"
#include "am_types.h"
#include "column.h"
#include "errlog.h"
#include "ils.h"
#include "interp.h"
#include "layer.h"
#include "math_const.h"
#include "model.h"
#include "output.h"
#include "phys_const.h"
#include "rt.h"
#include "spectra.h"
#include "units.h"

/*
 * The following two constants control iterative solution
 * for the tangent height of paths below the astronomical
 * horizon.  Convergence is deemed to have occurred when the
 * absolute value of the fractional change in ztan from one
 * iteration to the next is less than ZTAN_TOL.
 * ZTAN_ITER_MAX is the maximum number of iterations
 * allowed, after which convergence is deemed to have
 * failed.
 */
#define ZTAN_TOL      1e-10
#define ZTAN_ITER_MAX 100

/*
 * If the sine of the zenith angle computed from the
 * refractive invariant at the base of a layer exceeds
 * unity, then we're in a duct on the layer below.  However,
 * especially near a tangent point, there may be a small
 * numerical overshoot.  If sin(za) exceeds unity by less
 * than SIN_ZA_DUCT_TOL, it is simply clamped to unity.
 *
 * This tolerance can be somewhat smaller than ZTAN_TOL,
 * above, since ZTAN_TOL enters the refractive invariant via
 * the refractivity n-1.
 */
#define SIN_ZA_DUCT_TOL 1e-12

/*
 * For thin layers, TRIG_TOL sets the point at which path
 * trigonometry computations cannot be done to useful
 * accuracy at double precision, and a thin layer limit
 * is used instead.
 */
#define TRIG_TOL 1e-12

/*
 * Constants for set_hydrostatic_pressure_limits()
 *
 * The product of TOA_EPSILON and the base pressure of the
 * highest-altitude model definition layer is considered the
 * "top of atmosphere" pressure level.  A user-defined
 * source, observing, or tangent level placed at a lower
 * pressure will generate a warning, but will be
 * accommodated.  A user-defined level set at a height that
 * would correspond to a lower pressure (such as the orbital
 * height zobs of a satellite) is placed at P = 0.  When
 * such a level is inserted into the model, it is sorted
 * into height order with any other such levels.
 *
 * The product of P_OVERSHOOT_LIMIT and the base pressure of
 * the lowest model definition layer is the maximum allowed
 * pressure for any extrapolated level.
 */
#ifndef TOA_EPSILON
    #define TOA_EPSILON       0.1353352832366127
#endif

#ifndef P_OVERSHOOT_LIMIT
    #define P_OVERSHOOT_LIMIT 2.718281828459045
#endif

static int  compute_airmass_and_refraction(model_t*);
static int  compute_refractive_invariant(model_t*);
static int  copy_model_scalars(model_t*, model_t*);
static int  find_n_from_z(double*, const model_t*, const double);
static int  find_P_from_z(double*, const model_t*, const double);
static int  find_z_from_P(double*, const model_t*, const double);
static void set_hydrostatic_pressure_limits(model_t*);
static int  set_user_defined_levels(model_t*, model_t*);
static int  set_propagation_path_indices(model_t*);
static int  solve_for_tan_level(model_t*, model_t*);


/***********************************************************
* static int compute_airmass_and_refraction(model_t *model)
*
* Purpose:
*   Compute airmass and refraction for the propagation path.
*
*   An excellent and highly readable general reference on
*   this topic is the review:
*
*     A. T. Young 2006, "Understanding astronomical
*     refraction." The Observatory 126:82-115.
*
*   For spherical ray tracing, we use the approximations
*   from:
*
*     S. N. Kivalov 2007, "Improved ray tracing air mass
*     numbers model."  Applied Optics 46:7091-7098.
*
*   The Kivalov formulation has the very useful property
*   that the path segment in a layer is computed using the
*   mean of the radii and zenith angles at the upper and
*   lower layer boundaries.  This enables a midpoint
*   integration of the path as a function of height,
*   avoiding singularities at the tangent level.
*
*   For plane-parallel geometry, a similar approximation is
*   used for the airmass, which is computed as the secant of
*   the mean of the zenith angles at the upper and lower
*   layer boundaries.
*
* Arguments:
*   model_t *model     - pointer to model data structure
*
* Return:
*   1 if there were errors, 0 otherwise.
************************************************************/

static int compute_airmass_and_refraction(model_t *model)
{
    if (model->geometry & GEOMETRY_PLANE_PARALLEL) {
        double sin_za, za_hi;
        int lnum;
        sin_za = model->p;
        if (sin_za > 1.0) { /* duct */
            errlog(172, 0);
            return 1;
        }
        za_hi = asin(sin_za);
        for (lnum = 0; lnum <= model->path_mid; ++lnum) {
            layer_t *layer = model->layer[lnum];
            double za_lo, za_bar;
            sin_za = model->p / layer->nbase;
            if (sin_za > 1.0 + SIN_ZA_DUCT_TOL) { /* duct */
                errlog(172, lnum);
                return 1;
            }
            sin_za = sin_za > 1.0 ? 1.0 : sin_za;
            za_lo  = layer->za_base = asin(sin_za);
            za_bar = 0.5 * (za_hi + za_lo);
            za_hi  = za_lo;
            /*
             * An instrument or plume layer on a hydrostatic
             * model is assigned a fixed airmass factor of unity.
             */
            if ((model->PTmode & PTMODE_HYDROSTATIC) && (layer->h > 0.0))
                layer->m = 1.0;
            else
                layer->m = 1.0 / cos(za_bar);
        }
    } else if (model->geometry & (GEOMETRY_SPHERICAL | GEOMETRY_LIMB)) {
        double z_hi  = 0.0; /* inits to suppress compiler warning */
        double r_hi  = 0.0;
        double za_hi = 0.0;
        int lnum;
        /*
         * Here, scan down through the layer stack to find the
         * highest layer with its base at finite height and upper
         * boundary at infinity.  The airmass of this layer is
         * asymptotically equal to unity.  The central angle dphi
         * is asymptotically equal to the zenith angle at the
         * layer base; in the total refraction computation these
         * will cancel, and there will be no net refractive
         * bending in the top layer of spherical models.
         */
        for (lnum = 0; lnum <= model->path_mid; ++lnum) {
            layer_t *layer = model->layer[lnum];
            if (layer->Pbase > 0.0 || layer->type != LAYER_TYPE_DEF) {
                double sin_za;
                z_hi   = layer->zbase;
                r_hi   = model->R0 + z_hi;
                sin_za = model->p / (layer->nbase * r_hi);
                if (sin_za > 1.0 + SIN_ZA_DUCT_TOL) { /* duct */
                    errlog(172, lnum);
                    return 1;
                }
                sin_za         = sin_za > 1.0 ? 1.0 : sin_za;
                za_hi          = asin(sin_za);
                layer->m       = 1.0;
                layer->dphi    = za_hi;
                layer->za_base = za_hi;
                break;
            }
            /*
             * Dummy layer at infinity.
             */
            layer->m       = 1.0;
            layer->dphi    = 0.0;
            layer->za_base = 0.0;
        }
        /*
         * Log a warning if the entire path lies within a single
         * layer with upper bound at infinity (P == 0), since the
         * user might not have expected this.
         */
        if (lnum == model->path_mid)
            errlog(169, 0);
        /*
         * Now work down through the remaining layers.
         */
        for (++lnum; lnum <= model->path_mid; ++lnum) {
            layer_t *layer = model->layer[lnum];
            double z_lo, dz;
            double r_lo;
            double sin_za;
            double za_lo, za_bar;
            /*
             * A layer with h > 0 in a hydrostatic model
             * (spherical models are always hydrostatic) is an
             * embedded non-hydrostatic instrument or plume
             * layer. It is treated as a zero-thickness layer
             * having no effect on the path geometry.  The
             * airmass factor is set to unity.
             */
            if (layer->h > 0.0) {
                layer->m       = 1.0;
                layer->dphi    = 0.0;
                layer->za_base = model->layer[lnum - 1]->za_base;
                continue;
            }
            z_lo   = layer->zbase;
            dz     = z_hi - z_lo;
            r_lo   = model->R0 + z_lo;
            sin_za = model->p / (layer->nbase * r_lo);
            if (sin_za > 1.0 + SIN_ZA_DUCT_TOL) { /* duct */
                errlog(172, lnum);
                return 1;
            }
            sin_za = sin_za > 1.0 ? 1.0 : sin_za;
            za_lo  = asin(sin_za);
            za_bar = 0.5 * (za_lo + za_hi);
            /*
             * Thin layer limits are used when the computation
             * cannot be done at double precision.  The airmass
             * is undefined on a zero-thickness layer, and is
             * simply set to zero.
             */
            if ((PI_ON_TWO - za_bar < TRIG_TOL) ||
                    (dz < r_lo * TRIG_TOL)) {
                layer->m       = 0.0;
                layer->dphi    = 0.0;
                layer->za_base = model->layer[lnum - 1]->za_base;
                continue;
            }
            /*
             * For a normal layer, compute dphi and m using
             * Kivalov Eq.(14) and Eq.(16).
             */
            layer->dphi    = 2.0 * atan(dz * tan(za_bar) / (r_lo + r_hi));
            layer->m       = sqrt(r_lo * r_lo + r_hi * r_hi -
                    2.0 * r_lo * r_hi * cos(layer->dphi)) / dz;
            layer->za_base = za_lo;
            z_hi  = z_lo;
            r_hi  = r_lo;
            za_hi = za_lo;
        }
    } else {
        return 1; /* model geometry flags were invalid */
    }
    return 0;
}   /* compute_airmass_and_refraction() */


/***********************************************************
* int compute_model(model_t *model, model_t *lmodel_in)
*
* Purpose:
*   This is the top-level function for setting up the
*   atmospheric model and running radiative transfer
*   computations.
*
*   If lmodel_in != NULL, and lmodel_in->ngrid is nonzero,
*   then this is an update of a prior model computation, and
*   *lmodel_in holds the values of scalar model variables at
*   the time of that prior computation.  In this case,
*   expensive spectral computations are only redone if they
*   depend on changed scalar variables.  (Scalar variables
*   are always recomputed, because these computations take
*   relatively negligible time.)
*
*   After all computations are finished, scalar variables in
*   *lmodel_in are updated with the corresponding values
*   from *model to save the most recently computed model
*   state.
*
* Arguments:
*   model_t *model     - pointer to model data structure
*   model_t *lmodel_in - pointer to a model data structure
*                        holding scalar data from a prior
*                        computation, or NULL for a single-
*                        pass computation.
*
* Return:
*   1 if there were errors, 0 otherwise.
************************************************************/

int compute_model(model_t *model, model_t *lmodel_in)
{
    model_t *lmodel;
    double tmark = 0.0;

    if (model == NULL) {
        errlog(76, 0);
        return 1;
    }
    /*
     * *lmodel_in is ignored for scalar variable dependency
     * checking if it does not yet reflect a completed model
     * computation (i.e., if ngrid == 0 for non-NULL lmodel_in).
     * This is done by setting lmodel to NULL.  Otherwise,
     * lmodel is just a copy of lmodel_in.
     */
    lmodel = (lmodel_in != NULL && lmodel_in->ngrid == 0) ? NULL : lmodel_in;
    /*
     * Start the model run performance timer.
     */
    if (model->log_runtimes)
        model->runtime = tmark = am_timer(0.0);
    /*
     * Set up the atmospheric model.  This step may insert or
     * reorder model layers to support interpolated observing,
     * source, or tangent levels.  Because of this, we operate on
     * lmodel_in also, even though it is not being used here for
     * dependency checking.
     */
    if (setup_atmospheric_model(model, lmodel_in))
        return 1;
    /*
     * Log the model setup timing.  Start the layer performance
     * timers now that all layers are in place.
     */
    if (model->log_runtimes) {
        tmark += (model->am_runtime = am_timer(tmark));
        reset_layer_performance_timers(model);
    }
    /*
     * In some cases (IF spectra, particularly), we don't need to
     * compute all spectral channels.  Here, the sub-grids to be
     * computed are defined (or updated).
     */
    if (set_spectral_subgrid_ranges(model, lmodel))
        return 1;
    /*
     * Compute the optical depth (opacity) spectrum for the
     * entire propagation path.
     */
    if (compute_opacity_spectrum(model, lmodel))
        return 1;
    if (model->log_runtimes)
        tmark += (model->od_runtime = am_timer(tmark));
    /*
     * Do the radiative transfer computation along the
     * propagation path to compute the spectral radiance.
     */
    compute_spectral_radiance(model, lmodel);
    if (model->log_runtimes)
        tmark += (model->rt_runtime = am_timer(tmark));
    /*
     * Compute other spectra derived from radiance and opacity as
     * needed.
     */
    if (output[OUTPUT_TRANSMITTANCE].flags & OUTPUT_ACTIVE)
        compute_transmittance(model);
    if (output[OUTPUT_TB_PLANCK].flags & OUTPUT_ACTIVE)
        compute_Tb(model);
    if (output[OUTPUT_TB_RAYLEIGH_JEANS].flags & OUTPUT_ACTIVE ||
            output[OUTPUT_Y].flags & OUTPUT_ACTIVE ||
            output[OUTPUT_TSYS].flags & OUTPUT_ACTIVE)
        compute_Trj(model);
    if (output[OUTPUT_TSYS].flags & OUTPUT_ACTIVE)
        compute_spectral_Tsys(model);
    if (output[OUTPUT_Y].flags & OUTPUT_ACTIVE)
        compute_spectral_Y_factor(model);
    if (output[OUTPUT_RADIANCE_DIFF].flags & OUTPUT_ACTIVE)
        compute_radiance_difference_spectrum(model, lmodel);
    if (output[OUTPUT_DELAY].flags & OUTPUT_ACTIVE) {
        if (compute_delay_spectrum(model))
            return 1;
    }
    /*
     * Compute other output spectra as needed.
     */
    if (output[OUTPUT_FREE_SPACE_LOSS].flags & OUTPUT_ACTIVE)
        compute_free_space_loss(model);
    if (output[OUTPUT_K].flags & OUTPUT_ACTIVE)
        compute_k_out(model);
    /*
     * Apply modifications to spectra that model the response of
     * spectrometers or receivers.
     */
    if (model->ils != NULL) {
        if (apply_instrumental_line_shape(model, lmodel))
            return 1;
    }
    if (model->ifmode) {
        compute_if_spectra(model);
    }
    if (model->log_runtimes)
        tmark += (model->spec_runtime = am_timer(tmark));
    /*
     * Finished.  Log the total model compute time and update
     * *lmodel_in.
     */
    if (model->log_runtimes)
        model->runtime = am_timer(model->runtime);
    if (lmodel_in != NULL)
        return copy_model_scalars(model, lmodel_in);
    else
        return 0;
}   /* compute_model() */


/***********************************************************
* static int compute_refractive_invariant(model_t *model)
*
* Purpose:
*   Computes the refractive invariant, which is n * sin(za)
*   in plane-parallel mode, and n * r * sin(za) in spherical
*   mode.  This function assumes that layer heights have
*   already been computed.
*
* Arguments:
*   model_t *model - pointer to model structure
*
* Return:
*   0 if successful, 1 otherwise
************************************************************/

static int compute_refractive_invariant(model_t *model)
{
    int lnum;

    /*
     * Find the index number of the observing level for
     * plane-parallel or spherical models, or of the tangent
     * level for limb models.  The default observing or tangent
     * level if not otherwise defined is the base of the lowest
     * defined model layer.
     */
    if (model->geometry & GEOMETRY_LIMB)
        lnum = get_tan_lnum(model);
    else
        lnum = get_obs_lnum(model);
    /*
     * Compute the refractive invariant.
     */
    if (model->geometry & GEOMETRY_PLANE_PARALLEL) {
        model->p = model->layer[lnum]->nbase * sin(model->za);
    } else if (model->geometry & GEOMETRY_SPHERICAL) {
        double r = model->R0 + model->layer[lnum]->zbase;
        double sin_za = sin(model->za);
        /*
         * Here, clean up numerical rounding error for the
         * common nadir case, za == PI.
         */
        if (fabs(sin_za) < DBL_EPSILON)
            sin_za = 0.0;
        model->p = model->layer[lnum]->nbase * r * sin_za;
    } else if (model->geometry & GEOMETRY_LIMB) {
        double r = model->R0 + model->layer[lnum]->zbase;
        model->p = model->layer[lnum]->nbase * r;
    } else {
        return 1;
    }
    return 0;
}   /* compute_refractive_invariant() */


/***********************************************************
* static int copy_model_scalars(model_t *model, model_t *lmodel)
*
* Purpose:
*   Copies scalar variables from *model to *lmodel, where
*   *lmodel is a model structure recording the state of
*   a prior model computation.  The variables that need to
*   be copied are those that are compared from one model
*   computation to the next when determining whether
*   expensive spectral computations need to be redone, and
*   those that are involved in insertion and removal of
*   interpolated layers.
*
* Arguments:
*   model_t *model, *lmodel - pointers to model structures
*
* Return:
*   0 if successful, 1 otherwise
************************************************************/

static int copy_model_scalars(model_t *model, model_t *lmodel)
{
    int lnum, cnum, knum;

    if (model->nlayers != lmodel->nlayers) {
        errlog(35, 0);
        return 1;
    }
    lmodel->flo            = model->flo;
    lmodel->fif_min        = model->fif_min;
    lmodel->fif_max        = model->fif_max;
    lmodel->ils_fwhm       = model->ils_fwhm;
    lmodel->ils_fif        = model->ils_fif;
    lmodel->dsb_utol_ratio = model->dsb_utol_ratio;
    lmodel->T0             = model->T0;
    lmodel->Tref           = model->Tref;
    lmodel->sec_za         = model->sec_za;
    lmodel->isub[0]        = model->isub[0];
    lmodel->isub[1]        = model->isub[1];
    lmodel->ifmode         = model->ifmode;
    lmodel->ilsmode        = model->ilsmode;
    lmodel->ils_typenum    = model->ils_typenum;
    /*
     * Layer variables are copied only from layers that lie
     * within the range covered by the propagation path indices.
     * Outside this range they are not guaranteed to be updated.
     */
    for (lnum = model->path_min; lnum <= model->path_mid; ++lnum) {
        layer_t *layer  = model->layer[lnum];
        layer_t *llayer = lmodel->layer[lnum];
        if (layer->ncols != llayer->ncols) {
            errlog(35, 0);
            return 1;
        }
        llayer->P       = layer->P;
        llayer->T       = layer->T;
        llayer->dP      = layer->dP;
        llayer->dP_def  = layer->dP_def;
        llayer->Pbase   = layer->Pbase;
        llayer->Tbase   = layer->Tbase;
        llayer->h       = layer->h;
        llayer->m       = layer->m;
        llayer->type    = layer->type;
        for (cnum = 0; cnum < layer->ncols; ++cnum) {
            column_t *column  = layer->column[cnum];
            column_t *lcolumn = llayer->column[cnum];
            if (column->n_abscoeffs != lcolumn->n_abscoeffs) {
                errlog(35, 0);
                return 1;
            }
            lcolumn->N          = column->N;
            lcolumn->N_def      = column->N_def;
            lcolumn->N_scaled   = column->N_scaled;
            lcolumn->xvmr       = column->xvmr;
            lcolumn->vmr_scaled = column->vmr_scaled;
            lcolumn->RH = column->RH;
            for (knum = 0; knum < column->n_abscoeffs; ++knum) {
                abscoeff_t *abscoeff  = column->abscoeff[knum];
                abscoeff_t *labscoeff = lcolumn->abscoeff[knum];
                labscoeff->vmr_selfbroad = abscoeff->vmr_selfbroad;
            }
        }
    }
    /*
     * ngrid is initialized to 0, but is always non-zero in a
     * valid model.  It is copied here to indicate that valid
     * model data have been copied to *lmodel.
     */
    lmodel->ngrid = model->ngrid;
    return 0;
}   /* copy_model_scalars() */


/***********************************************************
* int find_n_from_z(
*         double *n, const model_t *model, const double z)
*
* Purpose:
*   Given a height z, finds the refractive index at that
*   height by interpolating in density (P/T) between model
*   definition levels.  If z lies outside the model
*   definition levels, n-1 is scaled by P/T.
*
*   This interpolation scheme ensures stability of path
*   geometry computations.  It avoids inadvertent creation
*   of ducts caused by composition changes above thin
*   interpolated layers, and maintains continuity of the
*   refractive index profile when solving for the tangent
*   level of paths below the astronomical horizon.
*
*   This function assumes that layer base heights are up to
*   date on model definition layers, and that model pressure
*   limits have been set.
*
* Arguments:
*   double *n      - pointer to receive refractive index
*   model_t *model - pointer to model data structure
*   const double z - height
*
* Return:
*   0 on success, 1 on failure
************************************************************/

int find_n_from_z(
        double *n, const model_t *model, const double z)
{
    int lnum_lo, lnum_hi;
    layer_t *layer_lo = NULL;
    layer_t *layer_hi = NULL;
    double P;

    if (model->Pmax == 0.0) {
        /*
         * There are no model definition layers at non-zero
         * pressure.
         */
        *n = 1.0;
        return 0;
    }
    /*
     * find_P_from_z() clamps P at a maximum extrapolated value
     * of model->Pmax.  Below that level, we assume constant P
     * for the purpose of computing n.  This situation arises
     * when the tangent point lies below the planetary surface.
     * This is normal for nadir sounding paths, and will not
     * cause trouble since in such cases the source level is on
     * the near side of the tangent point and no tangent level
     * will be inserted.  For this reason, we ignore the error
     * return unless P == 0.0.  (Note that a non-error return
     * with P == 0.0 is also possible.)
     */
    if (find_P_from_z(&P, model, z) && P == 0.0)
        return 1;
    /*
     * Scan for two model definition layers with nonzero base
     * pressures that have zbase values that bracket z.  Here,
     * lnum_hi and lnum_lo refer to indices of model definition
     * levels above and below the interpolated level,
     * respectively.
     */
    lnum_hi = -1;
    for (lnum_lo = 0; lnum_lo < model->nlayers; ++lnum_lo) {
        layer_lo = model->layer[lnum_lo];
        if (layer_lo->type != LAYER_TYPE_DEF) {
            /* skip if not a model definition layer */
            continue;
        } else if (layer_lo->Pbase == 0.0) {
            /* skip zero-pressure model definition layer */
            continue;
        } else if (layer_lo->dP == 0.0 && layer_lo->h >= 0.0) {
            /* skip plume layer */
            continue;
        } else if (layer_lo->zbase > z) {
            /* update upper z boundary */
            lnum_hi  = lnum_lo;
            layer_hi = model->layer[lnum_hi];
            continue;
        } else {
            /* found lower z boundary */
            break;
        }
    }
    if (lnum_hi < 0 && lnum_lo == model->nlayers) {
        /*
         * No model definition layers above or below z with
         * nonzero base pressure.  Not a normal case.
         */
        *n = 1.0;
    } else if (lnum_hi < 0) {
        /*
         * z is in the highest model definition layer with
         * nonzero base pressure.  Extrapolate isothermally from
         * the base of this layer upwards.  This is equivalant to
         * scaling n-1 by pressure.
         */
        *n  = layer_lo->nbase - 1.0;
        *n *= P / layer_lo->Pbase;
        *n += 1.0;
    } else if (lnum_hi < lnum_lo && lnum_lo < model->nlayers) {
        /*
         * z is bracketed by two model definition levels.
         * Compute n by interpolating in P/T.
         */
        double T;
        if (model->PTmode & PTMODE_TBASE) {
            T = log_x_interp(
                    layer_lo->Pbase,
                    layer_lo->Tbase,
                    layer_hi->Pbase,
                    layer_hi->Tbase,
                    P);
        } else {
            T = layer_lo->T;
        }
        *n = lin_interp(
                layer_lo->Pbase / layer_lo->Tbase,
                layer_lo->nbase,
                layer_hi->Pbase / layer_hi->Tbase,
                layer_hi->nbase,
                P / T);
    } else {
        /*
         * z is below the lowest model definition level.  Assume
         * constant molecular refractivity and scale n-1 by P/T.
         */
        double gamma = layer_base_lapse_rate(model, lnum_hi);
        double T = layer_hi->Tbase + gamma * (layer_hi->zbase - z);
        *n  = layer_hi->nbase - 1.0;
        *n *= (P / layer_hi->Pbase) * (layer_hi->Tbase / T);
        *n += 1.0;
    }
    return 0;
}   /* find_n_from_z() */


/***********************************************************
* static int find_P_from_z(
*         double *P, const model_t *model, const double z)
*
* Purpose:
*   Given a height z, finds the corresponding pressure level
*   in a hydrostatic model, by interpolating log(P) linearly
*   in z.  This function assumes that layer base heights are
*   up to date on model definition layers, and that model
*   pressure limits have been set.
*
*   If z is not bracketed by model definition levels, P is
*   extrapolated isothermally at the top of the atmosphere,
*   and extrapolated according to the extrapolation mode
*   setting at the base of the atmosphere.
*
*   If the pressure corresponding to z is lower than the
*   nominal top-of-atmosphere pressure model->Ptoa, then the
*   height z is regarded as lying outside the atmosphere,
*   and this function returns P = 0.
*
*   If the pressure corresponding to z is higher than the
*   base pressure extrapolation limit model->Pmax, P is
*   clamped at model->Pmax.
*
* Arguments:
*   double *P      - pointer to receive computed pressure
*   model_t *model - pointer to model data structure
*   const double z - height
*
* Return:
*   0 on success, 1 on failure or if P is clamped at Pmax
************************************************************/

static int find_P_from_z(
        double *P, const model_t *model, const double z)
{
    int lnum_lo, lnum_hi;
    layer_t *layer_lo = NULL;
    layer_t *layer_hi = NULL;

    if (model->Pmax == 0.0) {
        /*
         * There are no model definition layers at non-zero
         * pressure.
         */
        *P = 0.0;
        return 0;
    }
    /*
     * Scan for two model definition layers with nonzero base
     * pressures that have zbase values that bracket z.  Here,
     * lnum_hi and lnum_low refer to indices of model definition
     * levels above and below the interpolated level,
     * respectively.
     */
    lnum_hi = -1;
    for (lnum_lo = 0; lnum_lo < model->nlayers; ++lnum_lo) {
        layer_lo = model->layer[lnum_lo];
        if (layer_lo->type != LAYER_TYPE_DEF) {
            /* skip if not a model definition layer */
            continue;
        } else if (layer_lo->Pbase == 0.0) {
            /* skip zero-pressure model definition layer */
            continue;
        } else if (layer_lo->dP == 0.0 && layer_lo->h >= 0.0) {
            /* skip plume layer */
            continue;
        } else if (layer_lo->zbase > z) {
            /* update upper z boundary */
            lnum_hi  = lnum_lo;
            layer_hi = model->layer[lnum_hi];
            continue;
        } else { 
            /* found lower z boundary */
            break;
        }
    }
    if (lnum_hi < 0 && lnum_lo == model->nlayers) {
        /*
         * No model definition layers above or below z with
         * nonzero base pressure.  Not a normal case.
         */
        *P = 0.0;
    } else if (lnum_hi < 0) {
        /*
         * z is in the highest model definition layer with
         * nonzero base pressure.  Extrapolate isothermally from
         * the base of this layer upwards.  If the pressure is
         * below the top-of-atmosphere extrapolation limit, set
         * it to zero.
         */
        double R, T, x;
        T = model->PTmode & PTMODE_TBASE ? layer_lo->Tbase : layer_lo->T;
        if (T == 0.0) {
            errlog(188, 0);
            *P = 0.0;
            return 1;
        }
        R  = PCONST_KB / (layer_lo->Mair * PCONST_AMU);
        R *= unit_tab[UNIT_M].factor * unit_tab[UNIT_M].factor;
        x  = (layer_lo->zbase - z) * layer_lo->gbase / (R * T);
        *P = layer_lo->Pbase * exp(x);
        if (*P < model->Ptoa)
            *P = 0.0;
    } else if (lnum_hi < lnum_lo && lnum_lo < model->nlayers) {
        /*
         * z is bracketed by two model definition levels.
         * Compute P by linearly interpolating log(P) in z.
         */
        *P = log_y_interp(
                layer_lo->zbase,
                layer_lo->Pbase,
                layer_hi->zbase,
                layer_hi->Pbase,
                z);
    } else {
        /*
         * z is below the lowest model definition level.  Get the
         * lapse rate gamma corresponding to the base extension
         * mode setting and extrapolate.
         */
        double gamma, dz, Tbar, R, x;
        gamma = layer_base_lapse_rate(model, lnum_hi);
        dz   = layer_hi->zbase - z;
        Tbar = layer_hi->Tbase + gamma * 0.5 * dz;
        if (Tbar <= 0.0) {
            errlog(188, 0);
            *P = 0.0;
            return 1;
        }
        R  = PCONST_KB / (layer_hi->Mair * PCONST_AMU);
        R *= unit_tab[UNIT_M].factor * unit_tab[UNIT_M].factor;
        x  = (layer_hi->zbase - z) * layer_hi->gbase / (R * Tbar);
        *P = layer_hi->Pbase * exp(x);
        if (*P > model->Pmax) {
            *P = model->Pmax;
            return 1;
        }
    }
    return 0;
}   /* find_P_from_z() */


/***********************************************************
* static int find_z_from_P(
*        double *z, const model_t *model, const double P)
*
* Purpose:
*   Given a pressure level P, finds the corresponding height
*   z by interpolating linearly in log(P) between model
*   definition levels.  This function assumes that base
*   level heights have been computed on all model definition
*   layers, and model atmosphere pressure limits have been
*   set.
*
*   If P is not bracketed by model definition levels, z is
*   found by isothermal extrapolation.
*
* Return:
*   0 on success, 1 on failure
************************************************************/

static int find_z_from_P(
        double *z, const model_t *model, const double P)
{
    int lnum_lo, lnum_hi;
    layer_t *layer_lo = NULL;
    layer_t *layer_hi = NULL;

    if (P <= 0.0) {
        /*
         * Fail for non-positive P
         */
        errlog(163, 0);
        return 1;
    } else if (P < model->Ptoa) {
        /*
         * If P is more than two scale heights above the
         * highest model definition level, log a warning.
         */
        errlog(164, 0);
    } else if (P > model->Pmax) {
        /*
         * Fail if P would extrapolate more than one scale
         * height below the lowest model definition level.
         * This also catches the case when there are no model
         * definition levels at non-zero pressure.
         */
        errlog(171, 0);
        *z = 0.0;
        return 1;
    }
    /*
     * Scan for two model definition layers having Pbase
     * values that bracket P.
     */
    lnum_hi = -1;
    for (lnum_lo = 0; lnum_lo < model->nlayers; ++lnum_lo) {
        layer_lo = model->layer[lnum_lo];
        if (layer_lo->type != LAYER_TYPE_DEF) {
            /* skip if not a model definition layer */
            continue;
        } else if (layer_lo->Pbase <  P) {
            /* update upper P boundary */
            lnum_hi  = lnum_lo;
            layer_hi = model->layer[lnum_hi];
        } else if (layer_lo->Pbase >= P) {
            /* found lower P boundary */
            break;
        }
    }
    if (lnum_hi < 0 && lnum_lo == model->nlayers) {
        /*
         * No model definition layers above or below P.  This
         * case should already have been trapped by the tests
         * above.
         */
        errlog(205, __LINE__);
        return 1;
    } else if (lnum_hi < 0) {
        /*
         * P is in the highest model definition layer in the
         * atmosphere.  Extrapolate isothermally from the base of
         * this layer upwards to find z.
         */
        double R, T;
        T  = model->PTmode & PTMODE_TBASE ? layer_lo->Tbase : layer_lo->T;
        R  = PCONST_KB / (layer_lo->Mair * PCONST_AMU);
        R *= unit_tab[UNIT_M].factor * unit_tab[UNIT_M].factor;
        *z = layer_lo->zbase -
            (R * T / layer_lo->gbase) * log(P / layer_lo->Pbase);
    } else if (lnum_hi < lnum_lo && lnum_lo < model->nlayers) {
        /*
         * P is bracketed by two model definition levels.
         * Compute z by linearly interpolating z in log(P).
         */
        *z = log_x_interp(
                layer_lo->Pbase,
                layer_lo->zbase,
                layer_hi->Pbase,
                layer_hi->zbase,
                P);
    } else {
        /*
         * P is below the lowest model definition level.  Get the
         * lapse rate gamma corresponding to the base extension
         * mode setting and extrapolate.
         */
        double gamma, R, x, dz;
        gamma = layer_base_lapse_rate(model, lnum_hi);
        R   = PCONST_KB / (layer_hi->Mair * PCONST_AMU);
        R  *= unit_tab[UNIT_M].factor * unit_tab[UNIT_M].factor;
        x   = log(P / layer_hi->Pbase);
        dz  = x * R * layer_hi->Tbase;
        dz /= layer_hi->gbase - 0.5 * R * gamma * x;
        *z  = layer_hi->zbase - dz;
    }
    return 0;
}   /* find_z_from_P() */


/***********************************************************
* static void set_hydrostatic_pressure_limits(model_t *model)
*
* Purpose:
*   Set pressure limits for hydrostatic models.  These are
*   the top of atmosphere pressure Ptoa, defined as the
*   product of TOA_EPSILON and the base pressure of the
*   highest model definition layer; and the pressure
*   extrapolation limit Pmax, defined as the product of
*   P_OVERSHOOT_LIMIT and the base pressure of the lowest
*   model layer.  See the definitions of these constants
*   near the top of this file for details.
*
* Arguments:
*   model_t *model - pointer to model data structure
************************************************************/

static void set_hydrostatic_pressure_limits(model_t *model)
{
    int lnum;
    /*
     * If this is not a hydrostatic model, do nothing.
     */
    if (!(model->PTmode & PTMODE_HYDROSTATIC))
        return;
    /*
     * Scan for the highest model definition layer having a
     * non-zero base pressure.
     */
    model->Ptoa = 0.0;
    for (lnum = 0; lnum < model->nlayers; ++lnum) {
        layer_t *layer = model->layer[lnum];
        if (layer->type == LAYER_TYPE_DEF && layer->Pbase > 0.0) {
            model->Ptoa = layer->Pbase;
            break;
        }
    }
    model->Ptoa *= TOA_EPSILON;
    /*
     * Scan for the lowest model definition layer having a
     * non-zero base pressure.
     */
    model->Pmax = 0.0;
    for (lnum = model->nlayers - 1; lnum >= 0; --lnum) {
        layer_t *layer = model->layer[lnum];
        if (layer->type == LAYER_TYPE_DEF && layer->Pbase > 0.0) {
            model->Pmax = layer->Pbase;
            break;
        }
    }
    model->Pmax *= P_OVERSHOOT_LIMIT;
    return;
} /* set_hydrostatic_pressure_limits() */


/***********************************************************
* static int set_user_defined_levels(model_t *model, model_t *lmodel)
*
* Purpose:
*   Inserts any user-defined interpolated source, observing,
*   and tangent levels into the model, inserting the
*   additional layers needed to split the associated parent
*   model definition layers.  If any of these layers were
*   carried over from a prior computation, they will be
*   adjusted as needed.  This adjustment will include
*   deletion and re-insertion if the level has moved to a
*   different parent layer, changing the layer ordering.
*
* Arguments:
*   model_t *model  - pointer to model data structure
*   model_t *lmodel - pointer to a model data structure
*                     containing scalar data from a prior
*                     computation
* Return:
*   0 on success, 1 on error.
************************************************************/

static int set_user_defined_levels(model_t *model, model_t *lmodel)
{
    set_hydrostatic_pressure_limits(model);
    if (model->geometry & GEOMETRY_POBS_USER_DEFINED) {
        double zobs, nobs;
        if (    find_z_from_P(&zobs, model, model->Pobs) ||
                find_n_from_z(&nobs, model, zobs) ||
                insert_interpolated_level(
                    model,
                    lmodel,
                    LAYER_TYPE_OBS,
                    zobs,
                    model->Pobs,
                    nobs)) {
            errlog(190, 0);
            return 1;
        }
    }
    if (model->geometry & GEOMETRY_ZOBS_USER_DEFINED) {
        double Pobs, nobs;
        if (find_P_from_z(&Pobs, model, model->zobs)) {
            errlog(171, 0);
            errlog(190, 0);
            return 1;
        }
        if (    find_n_from_z(&nobs, model, model->zobs) ||
                insert_interpolated_level(
                    model,
                    lmodel,
                    LAYER_TYPE_OBS,
                    model->zobs,
                    Pobs,
                    nobs)) {
            errlog(190, 0);
            return 1;
        }
    }
    if (model->geometry & GEOMETRY_PSOURCE_USER_DEFINED) {
        double zsource, nsource;
        if (    find_z_from_P(&zsource, model, model->Psource) ||
                find_n_from_z(&nsource, model, zsource) ||
                insert_interpolated_level(
                    model,
                    lmodel,
                    LAYER_TYPE_SOURCE,
                    zsource,
                    model->Psource,
                    nsource)) {
            errlog(191, 0);
            return 1;
        }
    }
    if (model->geometry & GEOMETRY_ZSOURCE_USER_DEFINED) {
        double Psource, nsource;
        if (find_P_from_z(&Psource, model, model->zsource)) {
            errlog(171, 0);
            errlog(191, 0);
            return 1;
        }
        if (    find_n_from_z(&nsource, model, model->zsource) ||
                insert_interpolated_level(
                    model,
                    lmodel,
                    LAYER_TYPE_SOURCE,
                    model->zsource,
                    Psource,
                    nsource)) {
            errlog(191, 0);
            return 1;
        }
    }
    if (model->geometry & GEOMETRY_PTAN_USER_DEFINED) {
        double ztan, ntan;
        if (    find_z_from_P(&ztan, model, model->Ptan) ||
                find_n_from_z(&ntan, model, ztan) ||
                insert_interpolated_level(
                    model,
                    lmodel,
                    LAYER_TYPE_TAN,
                    ztan,
                    model->Ptan,
                    ntan)) {
            errlog(192, 0);
            return 1;
        }
    }
    if (model->geometry & GEOMETRY_ZTAN_USER_DEFINED) {
        double Ptan, ntan;
        if (find_P_from_z(&Ptan, model, model->ztan)) {
            errlog(171, 0);
            errlog(192, 0);
            return 1;
        }
        if (    find_n_from_z(&ntan, model, model->ztan) ||
                insert_interpolated_level(
                    model,
                    lmodel,
                    LAYER_TYPE_TAN,
                    model->ztan,
                    Ptan,
                    ntan)) {
            errlog(192, 0);
            return 1;
        }
    }
    return 0;
}   /* set_user_defined_levels() */


/***********************************************************
* static int set_propagation_path_indices(model_t *model)
*
* Purpose:
*   Sets up the layer index ranges for the propagation path.
*
*   The propagation path is divided across the tangent point
*   into upward and downward segments.  Both of these
*   segments occur for limb paths and paths below the
*   astronomical horizon through the tangent point.
*
*   Paths above the horizon (which includes all plane-
*   parallel paths), and paths below the horizon from a
*   source point on the near side of the tangent point, have
*   a single upward or downward segment.  In these cases,
*   the unused segment is assigned reversed index ranges
*   relative to the propagation direction so it will always
*   be skipped.
*
*   An error can occur if the propagation path does not
*   intersect a user-defined source level.
*
* Arguments:
*   model_t *model  - pointer to model data structure
*
* Return:
*   0 on success, 1 on error.
************************************************************/

static int set_propagation_path_indices(model_t *model)
{
    int obs_lnum    = get_obs_lnum(model);
    int source_lnum = get_source_lnum(model);
    int tan_lnum    = get_tan_lnum(model);
    int base_lnum   = get_lnum_by_type(model, LAYER_TYPE_DEF);

    /*
     * Moderate extrapolation of observing, source, or tangent
     * levels below the base level of the lowest defined model is
     * allowed.  Nevertheless, log a warning to alert the user,
     * unless an extrapolation mode has been explicitly specified
     * in an "extend" statement.
     */
    if (!(model->PTmode & PTMODE_EXTEND_USER_DEFINED)) {
        if (obs_lnum > base_lnum)
            errlog(165, 0);
        if (source_lnum > base_lnum)
            errlog(177, 0);
        if (tan_lnum > base_lnum)
            errlog(166, 0);
    }
    if (model->geometry & GEOMETRY_LIMB) {
        /*
         * Limb path.  The tangent point is either user-defined,
         * or defaults to the bottom model layer.
         */
        model->path_begin = 0;
        model->path_mid   = tan_lnum;
        model->path_end   = 0;
    } else if (fabs(model->za) <= PI_ON_TWO) {
        /*
         * Path on or above the astronomical horizon.  Fail if
         * the source level is below the astronomical horizon.
         */
        if (source_lnum > obs_lnum) {
            errlog(178, 0);
            return 1;
        } else {
            model->path_begin = source_lnum + 1;
            model->path_mid   = obs_lnum;
            model->path_end   = obs_lnum + 1;
        }
    } else if (!(model->geometry & GEOMETRY_SOURCE_NEAR)) {
        /*
         * Path through tangent point below the astronomical
         * horizon to a source on the far side of the tangent
         * point.  Fail if the source level is below the tangent
         * level.
         */
        if (source_lnum > tan_lnum) {
            errlog(179, 0);
            return 1;
        } else {
            model->path_begin = source_lnum + 1;
            model->path_mid   = tan_lnum;
            model->path_end   = obs_lnum + 1;
        }
    } else {
        /*
         * Path below the astronomical horizon to a source on the
         * near side of the tangent point.  Fail if the source
         * level isn't between the tangent and observing levels.
         *
         * A tangent level in this case could be below the
         * planetary surface, and will not have been inserted in
         * the model, so check heights instead of layer indices.
         */
        double zobs    = model->layer[obs_lnum]->zbase;
        double zsource = model->layer[source_lnum]->zbase;
        if (zsource > zobs || zsource < model->ztan) {
            errlog(180, 0);
            return 1;
        } else {
            model->path_begin = source_lnum + 1;
            model->path_mid   = source_lnum;
            model->path_end   = obs_lnum + 1;
        }
    }
    if (model->geometry & GEOMETRY_REVERSE) {
        /*
         * For reverse geometry, swap path_begin and path_end.
         * This has the effect of reversing both the order and
         * the direction of the path.
         */
        int temp          = model->path_end;
        model->path_end   = model->path_begin;
        model->path_begin = temp;
    }
    /*
     * Set the minimum layer index for the path.
     */
    if (model->path_begin < model->path_end)
        model->path_min = model->path_begin;
    else
        model->path_min = model->path_end;
    return 0;
}   /* set_propagation_path_indices() */


/***********************************************************
* int setup_atmospheric_model(model_t *model, model_t *lmodel)
*
* Purpose:
*   Sets up the atmospheric model and path geometry.
*
* Arguments:
*   model_t *model  - pointer to model data structure
*   model_t *lmodel - pointer to a model data structure
*                     containing scalar data from a prior
*                     computation
*
* Return:
*   0 on success, 1 on error.
************************************************************/

int setup_atmospheric_model(model_t *model, model_t *lmodel)
{
    /*
     * If there are no model layers, set up path ranges that
     * cover no layers and return.
     */
    if (model->nlayers <= 0) {
        model->path_begin =  0;
        model->path_mid   = -1;
        model->path_end   =  0;
        return 0;
    }
    /*
     * Set the observer zenith angle za and sec(za).
     */
    if (model->geometry & GEOMETRY_SEC_ZA_USER_DEFINED) {
        model->za = acos(1.0 / model->sec_za);
    } else {
        if (model->geometry & GEOMETRY_PLANE_PARALLEL) {
            /*
             * Plane-parallel models are restricted to zenith
             * angles less than 90 degrees.
             */
            if ((PI_ON_TWO - fabs(model->za)) < FLT_EPSILON) {
                errlog(161, 0);
                return 1;
            }
        }
        model->sec_za = 1.0 / cos(model->za);
    }
    /*
     * The atmospheric model is anchored to the layers as defined
     * in the model configuration file (the "model definition
     * layers").  In hydrostatic models, additional interpolated
     * or extrapolated levels may subsequently be added for
     * user-defined source, observing, or tangent levels.  These
     * require the addition of new model layers created by
     * splitting a model definition layer (a "parent layer") into
     * parts, or by extrapolating a parent layer below the layer
     * stack.
     *
     * The first pass through the layer stack is carried out on
     * model definition layers only.  Any added layers carried
     * over from a prior computation are ignored, and parent
     * layers are restored to their unmodified state.  This pass
     * computes scalar parameters including mixing ratios,
     * geometric heights, and base level refractive indicies,
     * which are inherited by or interpolated on user-defined
     * levels and their associated layers.
     */
    restore_def_layer_explicit_col_densities(model);
    set_def_layer_pressures(model);
    set_def_layer_temperatures(model);
    set_def_layer_vmrs_and_apply_Nscale_factors(model);
    set_def_layer_heights_and_col_densities(model);
    set_def_layer_refractive_indices(model);
    /*
     * Add user-defined interpolated source, observing, or
     * tangent levels and associated layers.  Existing
     * interpolated layers carried over from a prior computation
     * are deleted and reinserted if the layer ordering will
     * change, otherwise layers are adjusted in place.
     * 
     * Because these operations affect the number of layers or
     * their ordering, they are carried out on both *model and
     * *lmodel, but level positions and other variables are
     * updated on *model only.  On initial layer insertion, or
     * in the case of a layer deletion and re-insertion, the
     * default initialization of parameters on *lmodel will
     * trigger (re)computation of all spectra on that layer.
     */
    if (set_user_defined_levels(model, lmodel))
        return 1;
    /*
     * Now that we know the observer level height and zenith
     * angle (or the tangent level height in limb mode) we have
     * what we need to compute the refractive invariant.
     */
    if (compute_refractive_invariant(model))
        return 1;
    /*
     * If we're in spherical geometry mode, and looking below the
     * astronomical horizon, iteratively solve for and insert the
     * tangent level.
     */
    if ((model->geometry & GEOMETRY_SPHERICAL) &&
            (fabs(model->za) > PI_ON_TWO)) {
        if (solve_for_tan_level(model, lmodel))
            return 1;
    }
    /*
     * With all interpolated levels now in place, set and adjust
     * remaining variables on interpolated and parent layers.
     */
    set_interpolated_layer_variables(model);
    compute_adjusted_selfbroad_vmrs(model);
    /*
     * Compute the propagation path geometry.
     */
    if (set_propagation_path_indices(model))
        return 1;
    if (compute_airmass_and_refraction(model))
        return 1;
    return 0;
}   /* setup_atmospheric_model() */


/***********************************************************
* static int solve_for_tan_level(
*       model_t *model, model_t *lmodel)
*
* Purpose:
*   In spherical geometry mode, for paths below the
*   astronomical horizon, this function solves for the
*   tangent height and inserts a tangent level.  Because the
*   refractive index n varies with altitude, this has to be
*   done iteratively.
*
* Arguments:
*   model_t *model  - pointer to model data structure
*
* Return:
*   0 on success, 1 on error or failure to converge.
************************************************************/

static int solve_for_tan_level(model_t *model, model_t *lmodel)
{
    double ntan, ztan;
    int    iter, lnum;
    /*
     * If a tangent layer already exists from a prior
     * computation, start the iteration at the base of that
     * layer.  Otherwise, start at the base of the lowest model
     * layer.
     */
    if ((lnum = get_lnum_by_type(model, LAYER_TYPE_TAN)) < 0)
        lnum = model->nlayers - 1;
    ztan = model->layer[lnum]->zbase;
    for (iter = 0; iter < ZTAN_ITER_MAX; ++iter) {
        double n, ztan_last;
        if (find_n_from_z(&n, model, ztan)) {
            errlog(168, 0);
            return 1;
        }
        ztan_last = ztan;
        ztan      = (model->p / n) - model->R0;
        if (fabs((ztan - ztan_last) / ztan) < ZTAN_TOL)
            break;
    }
    /*
     * The iteration may fail to converge to a solution due to
     * ducting, extrapolation beyond defined layers, or for some
     * other reason.  If so, log an error and give up.
     */
    if (iter == ZTAN_ITER_MAX) {
        errlog(168, 0);
        return 1;
    }
    /*
     * If the tangent level pressure is higher than the
     * extrapolation limit, and the source is on the far side of
     * the tangent point from the observer, log an error and
     * fail.
     */
    if (find_P_from_z(&model->Ptan, model, ztan) &&
            !(model->geometry & GEOMETRY_SOURCE_NEAR)) {
        errlog(196, 0);
        return 1;
    }
    /*
     * If the tangent level height is below the planetary surface
     * and the source is on the far side of the tangent point
     * from the observer, log a warning, but continue
     * nonetheless.
     */
    if (ztan < 0.0 && !(model->geometry & GEOMETRY_SOURCE_NEAR))
        errlog(167, 0);
    model->ztan = ztan;
    if (find_n_from_z(&ntan, model, ztan))
        return 1;
    /*
     * Insert a tangent level if the source level is on the far
     * side of the tangent point.
     */
    if (!(model->geometry & GEOMETRY_SOURCE_NEAR) &&
            insert_interpolated_level(
                model,
                lmodel,
                LAYER_TYPE_TAN,
                model->ztan,
                model->Ptan,
                ntan)) {
            errlog(182, 0);
            return 1;
    }
    return 0;
}   /* solve_for_tan_level() */


/***********************************************************
* double source_to_obs_geometric_distance(model_t *model)
*
* Purpose:
*   Computes the geometric (i.e. unrefracted) distance
*   between source and observer.  If either path endpoint
*   lies at infinity, the return value is DBL_MAX.
*
* Arguments:
*   model_t *model - pointer to model data structure
*
* Return:
*   geometric distance
************************************************************/

double source_to_obs_geometric_distance(model_t *model)
{
    if (model->path_begin == 0 || model->path_end == 0)
        return DBL_MAX;
    if (model->geometry & GEOMETRY_PLANE_PARALLEL) {
        double x = 0.0;
        double z = 0.0;
        int lnum;
        for (lnum = model->path_begin; lnum <= model->path_mid; ++lnum) {
            layer_t *layer = model->layer[lnum];
            double dz;
            if (layer->h > 0.0)
                dz = layer->h;
            else
                dz = model->layer[lnum - 1]->zbase - layer->zbase;
            x += dz * tan(layer->za_base);
            z += dz;
        }
        for (lnum = model->path_mid; lnum >= model->path_end; --lnum) {
            layer_t *layer = model->layer[lnum];
            double dz;
            if (layer->h > 0.0)
                dz = layer->h;
            else
                dz = model->layer[lnum - 1]->zbase - layer->zbase;
            x += dz * tan(layer->za_base);
            z += dz;
        }
        return sqrt(x * x + z * z);
    } else {
        double phi = 0.0;
        double Robs, Rsource;
        int lnum;
        for (lnum = model->path_begin; lnum <= model->path_mid; ++lnum)
            phi += model->layer[lnum]->dphi;
        for (lnum = model->path_mid  ; lnum >= model->path_end; --lnum)
            phi += model->layer[lnum]->dphi;
        lnum    = get_obs_lnum(model);
        Robs    = model->R0 + model->layer[lnum]->zbase;
        lnum    = get_source_lnum(model);
        Rsource = model->R0 + model->layer[lnum]->zbase;
        return sqrt(Robs * Robs + Rsource * Rsource - 
                2.0 * Robs * Rsource * cos(phi));
    }
}   /* source_to_obs_geometric_distance() */


/***********************************************************
* double source_to_obs_path_distance(model_t *model)
*
* Purpose:
*   Computes the total refracted path distance from source
*   to observer.  If either path endpoint lies at infinity,
*   the return value is DBL_MAX.
*
* Arguments:
*   model_t *model - pointer to model data structure
*
* Return:
*   path distance
************************************************************/

double source_to_obs_path_distance(model_t *model)
{
    double d = 0.0;
    int lnum;

    if (model->path_begin == 0 || model->path_end == 0)
        return DBL_MAX;
    for (lnum = model->path_begin; lnum <= model->path_mid; ++lnum) {
        layer_t *layer = model->layer[lnum];
        if (layer->h > 0.0)
            d += layer->h;
        else
            d += layer->m * (model->layer[lnum - 1]->zbase - layer->zbase);
    }
    for (lnum = model->path_mid; lnum >= model->path_end; --lnum) {
        layer_t *layer = model->layer[lnum];
        if (layer->h > 0.0)
            d += layer->h;
        else
            d += layer->m * (model->layer[lnum - 1]->zbase - layer->zbase);
    }
    return d;
}   /* source_to_obs_path_distance() */
        

/***********************************************************
* double source_to_obs_projected_distance(model_t *model)
*
* Purpose:
*   In plane-parallel models, computes the length of the
*   horizontal plane projection of the path from source to
*   observer.
*
*   In spherical models, computes the projected great circle
*   map distance along the planetary surface between source
*   and observer nadir points.  
*
* Arguments:
*   model_t *model - pointer to model data structure
*
* Return:
*   Projected map distance
************************************************************/

double source_to_obs_projected_distance(model_t *model)
{
    if (model->geometry & GEOMETRY_PLANE_PARALLEL) {
        double x = 0.0;
        int lnum;
        /*
         * If either path endpoint is at infinity in a
         * plane-parallel model, the projected path length
         * is infinite unless the zenith angle is zero.
         */
        if (model->path_begin == 0 || model->path_end == 0)
            return fabs(model->za) < DBL_MIN ? 0.0 : DBL_MAX; 
        for (lnum = model->path_begin; lnum <= model->path_mid; ++lnum) {
            layer_t *layer = model->layer[lnum];
            double dz;
            if (layer->h > 0.0)
                dz = layer->h;
            else
                dz = model->layer[lnum - 1]->zbase - layer->zbase;
            x += dz * tan(layer->za_base);
        }
        for (lnum = model->path_mid; lnum >= model->path_end; --lnum) {
            layer_t *layer = model->layer[lnum];
            double dz;
            if (layer->h > 0.0)
                dz = layer->h;
            else
                dz = model->layer[lnum - 1]->zbase - layer->zbase;
            x += dz * tan(layer->za_base);
        }
        return x;
    } else {
        double phi = 0.0;
        int lnum;
        for (lnum = model->path_begin; lnum <= model->path_mid; ++lnum)
            phi += model->layer[lnum]->dphi;
        for (lnum = model->path_mid  ; lnum >= model->path_end; --lnum)
            phi += model->layer[lnum]->dphi;
        return phi * model->R0;
    }
}   /* source_to_obs_projected_distance() */


/***********************************************************
* double total_airmass(model_t *model)
*
* Purpose:
*   Computes the total relative airmass by column density,
*   irrespective of species, for the propagation path.
*
*   Airmass is normalized to the sum of the zenith column
*   densities of all layers in the path from source to
*   observer.  If the path traverses a tangent point, layers
*   that are traversed twice are only counted once in the
*   zenith column density sum.
*
* Arguments:
*   model_t *model - pointer to model data structure
*
* Return:
*   Total nominal airmass, or 0.0 if undefined.
************************************************************/

double total_airmass(model_t *model)
{
    double N     = 0.0;
    double N_los = 0.0;
    int lnum;
    int path_begin = model->path_begin;
    int path_mid   = model->path_mid;
    int path_end   = model->path_end;
    int path_min   = path_begin < path_end ? path_begin : path_end;

    /*
     * Zenith total column density.
     */
    for (lnum = path_mid; lnum >= path_min; --lnum) {
        layer_t *layer = model->layer[lnum];
        int cnum;
        for (cnum = 0; cnum < layer->ncols; ++cnum) {
            column_t *column = layer->column[cnum];
            int ctype = column->col_typenum;
            if (!(col_type[ctype].flags & COL_PARAMETRIC))
                N += column->N_scaled;
        }
    }
    /*
     * Line-of-sight total column density.
     */
    for (lnum = model->path_begin; lnum <= model->path_mid; ++lnum) {
        layer_t *layer = model->layer[lnum];
        int cnum;
        for (cnum = 0; cnum < layer->ncols; ++cnum) {
            column_t *column = layer->column[cnum];
            int ctype = column->col_typenum;
            if (!(col_type[ctype].flags & COL_PARAMETRIC))
                N_los += column->N_scaled * layer->m;
        }
    }
    for (lnum = model->path_mid; lnum >= model->path_end; --lnum) {
        layer_t *layer = model->layer[lnum];
        int cnum;
        for (cnum = 0; cnum < layer->ncols; ++cnum) {
            column_t *column = layer->column[cnum];
            int ctype = column->col_typenum;
            if (!(col_type[ctype].flags & COL_PARAMETRIC))
                N_los += column->N_scaled * layer->m;
        }
    }
    /*
     * The airmass is the ratio of line-of-sight to zenith column
     * density.
     */
    return N > 0.0 ? N_los / N : 0.0;
}   /* total_airmass() */


/***********************************************************
* double total_refraction(model_t *model)
*
* Purpose:
*   Computes the total refraction through the model
*   propagation path.  For spherical paths, the sign of the
*   refractive bending is considered positive when the line
*   of sight is concave towards the planetary center,
*   irrespective of propagation direction.  Similarly, for
*   plane-parallel geometry, the sign is positive for
*   paths for which the zenith angle is smaller at lower
*   altitude, irrespective of propagation direction.
*
*
* Arguments:
*   model_t *model  - pointer to model data structure
*
* Return:
*   Total refraction in radians.
************************************************************/

double total_refraction(model_t *model)
{
    double R = 0.0;
    int lnum;

    if (model->geometry & GEOMETRY_PLANE_PARALLEL) {
        /*
         * In the plane parallel case, the total refraction is
         * just the difference in zenith angles between the two
         * ends of a path segment.
         */
        if (model->path_begin <= model->path_mid) {
            R += model->layer[model->path_begin]->za_base;
            R -= model->layer[model->path_mid]->za_base;
        }
        if (model->path_mid >= model->path_end) {
            R += model->layer[model->path_end]->za_base;
            R -= model->layer[model->path_mid]->za_base;
        }
    } else {
        /*
         * In the spherical case, the refractive bending is zero
         * for a layer with an upper boundary pressure of zero,
         * so all such layers are skipped.  (One way to think of
         * this is that, formally, such a layer has a pressure
         * gradient asymptotically equal to zero, and hence no
         * refractive index gradient.)
         */
        for (lnum = model->path_begin; lnum <= model->path_mid; ++lnum) {
            if (lnum != 0 && model->layer[lnum - 1]->Pbase > 0.0)
                R += model->layer[lnum - 1]->za_base
                    - model->layer[lnum]->za_base
                    + model->layer[lnum]->dphi;
        }
        for (lnum = model->path_mid; lnum >= model->path_end; --lnum) {
            if (lnum != 0 && model->layer[lnum - 1]->Pbase > 0.0)
                R += model->layer[lnum - 1]->za_base
                    - model->layer[lnum]->za_base
                    + model->layer[lnum]->dphi;
        }
    }
    return R;
}   /* total_refraction() */
