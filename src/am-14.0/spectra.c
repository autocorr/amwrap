/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* spectra.c                      S. Paine rev. 2024 April 17
*
* Functions for computing various output spectra including
* those derived from the optical depth and radiance spectra.
************************************************************/

#include <math.h>
#include <float.h>
#include <stdlib.h>

#include "abscoeff.h"
#include "am_types.h"
#include "column.h"
#include "errlog.h"
#include "math_const.h"
#include "model.h"
#include "output.h"
#include "phys_const.h"
#include "planck.h"
#include "specfunc.h"
#include "transform.h"
#include "units.h"


/***********************************************************
* int compare_spectral_subgrid_ranges(model_t *model, model_t *lmodel)
*
* Purpose:
*   Checks whether the spectral subgrid ranges being
*   computed have changed from the last model computation
*   to the current one.  This might happen when computing
*   IF spectra if the LO frequency has changed.
*
* Arguments:
*   model_t *model  - pointer to model data structure
*   model_t *lmodel - pointer to a model data structure
*       containing scalar data from a prior computation.
*
* Return:
*   1 if the subgrid range has changed or if lmodel == NULL
*   0 if there have been no changes
************************************************************/

int compare_spectral_subgrid_ranges(model_t *model, model_t *lmodel)
{
    return
        (lmodel == NULL) ||
        (model->isub[0].min != lmodel->isub[0].min) ||
        (model->isub[0].max != lmodel->isub[0].max) ||
        (model->isub[1].min != lmodel->isub[1].min) ||
        (model->isub[1].max != lmodel->isub[1].max);
}   /* compare_spectral_subgrid_ranges() */


/***********************************************************
* void compute_k_out(model_t *model)
*
* Purpose:
*   Copies the first spectral absorption coefficient from
*   the model to the array model.k_out.  A warning is logged
*   if the model contains zero or more than one absorption
*   coefficient.  Appropriate units are set according to
*   whether the absorption coefficient is a molecular one or
*   a binary one.
*
* Arguments:
*   model_t *model - pointer to model structure
************************************************************/

void compute_k_out(model_t *model)
{
    int lnum, cnum, knum;
    int subgrid;
    int abscoeff_count = 0;

    for (lnum = 0; lnum < model->nlayers; ++lnum) {
        layer_t *layer = model->layer[lnum];
        for (cnum = 0; cnum < layer->ncols; ++cnum) {
            column_t *column = layer->column[cnum];
            for (knum = 0; knum < column->n_abscoeffs; ++knum) {
                abscoeff_t *abscoeff = column->abscoeff[knum];
                if (abscoeff->k == NULL)
                    continue;
                if (abscoeff_count == 0) {
                    output[OUTPUT_K].default_unitnum =
                        k_type[abscoeff->k_typenum].unitnum;
                    /*
                     * If no output unit was set by the user, assign
                     * it here.  Otherwise, if a user unit was
                     * assigned, correct with a warning if
                     * needed.
                     */
                    if (output[OUTPUT_K].unitnum == UNIT_AUTO) {
                        output[OUTPUT_K].unitnum =
                            k_type[abscoeff->k_typenum].unitnum;
                    } else if (
                            unit_tab[output[OUTPUT_K].unitnum].group !=
                            unit_tab[output[OUTPUT_K].default_unitnum].group) {
                        errlog(200,0);
                        output[OUTPUT_K].unitnum =
                            output[OUTPUT_K].default_unitnum;
                    }
                    output[OUTPUT_K].k_typenum = abscoeff->k_typenum;
                    for (subgrid = 0; subgrid <= 1; ++subgrid) {
                        gridsize_t i;
                        gridsize_t imin = model->isub[subgrid].min;
                        gridsize_t imax = model->isub[subgrid].max;
                        for (i = imin; i <= imax; ++i)
                            model->k_out[i] = abscoeff->k[i];
                    }
                }
                ++abscoeff_count;
            }
        }
    }
    /*
     * If the model contains no spectral absorption coefficients,
     * set k_out to zero and log a warning.
     */
    if (abscoeff_count == 0) {
        for (subgrid = 0; subgrid <= 1; ++subgrid) {
            gridsize_t i;
            gridsize_t imin = model->isub[subgrid].min;
            gridsize_t imax = model->isub[subgrid].max;
            for (i = imin; i <= imax; ++i)
                model->k_out[i] = 0.0;
        }
        errlog(198, 0);
    }
    /*
     * If the model contains more than one spectral absorption
     * coefficent, log a warning
     */
    if (abscoeff_count > 1)
        errlog(199, abscoeff_count);
    return;
}   /* compute_k_out() */


/***********************************************************
* int compute_delay_spectrum(model_t *model)
*
* Purpose:
*   Computes the excess delay spectrum, which is divided
*   into two contributions.  Optical transitions contribute
*   an essentially frequency-independent refractivity, which
*   produces a delay proportional to column density only.
*   Lower-lying (mainly rotational) transitions contribute a
*   dispersive delay which depends on all the factors which
*   influence these spectra, including P, T, mixing ratios
*   and column densities.
*
*   In this function, a first pass through the model
*   computes the total optical delay.  In am this is
*   associated with particular h2o and dry_air column types.
*   Then, the dispersive delay associated with all types of
*   absorption which have been included in the model is
*   computed by Hilbert transformation of the opacity
*   spectrum.  Both contributions are then added together
*   to produce the total spectral delay L.
*
* Arguments:
*   model_t *model - pointer to model structure
*
* Return:
*   0 if successful, 1 otherwise
************************************************************/

int compute_delay_spectrum(model_t *model)
{
    gridsize_t i, if0, imax, jr, ji, np2;
    double Lopt;
    double *z;
    int lnum, cnum;

    if (model->L == NULL) {
        errlog(71, 0);
        return 1;
    }
    /*
     * Compute the optical delay.  This is the constant part of the
     * delay spectrum arising from optical transitions.
     */
    Lopt = 0.0;
    for (lnum = model->path_begin; lnum <= model->path_mid; ++lnum) {
        for (cnum = 0; cnum < model->layer[lnum]->ncols; ++cnum) {
            Lopt += column_optical_delay(model, lnum, cnum);
        }
    }
    for (lnum = model->path_mid; lnum >= model->path_end; --lnum) {
        for (cnum = 0; cnum < model->layer[lnum]->ncols; ++cnum) {
            Lopt += column_optical_delay(model, lnum, cnum);
        }
    }
    if0 = (gridsize_t)(0.5 + model->f[0] / model->df);
    np2 = 2 * model->nLpad;
    /*
     * Allocate temporary space for a complex array z[0..2*nLpad-1]
     */
    if ((z = (double*)malloc(np2 * sizeof(double))) == NULL) {
        errlog(68, 0);
        return 1;
    }
    /*
     * The imaginary part ni of the refractive index n = nr + i * ni,
     * integrated over the propagation path, is
     *
     *   integral(ni * ds) = c * tau / (4 * pi * f).
     *
     * Put this into the real part of z, antisymmetrically about f = 0.
     * Zero fill the rest of the array.
     */
    for (i = 0; i < if0; ++i) {
        jr = i << 1;
        ji = jr + 1;
        z[jr] = 0.0;
        z[ji] = 0.0;
    }
    jr = if0 << 1;
    ji = jr + 1;
    z[jr] = (model->f[0] < F_EPSILON) ?
        0.0 : (model->tau[0] * C_ON_FOURPI / model->f[0]);
    z[ji] = 0.0;
    for (i = 1; i < model->ngrid; ++i) {
        jr = (i + if0) << 1;
        ji = jr + 1;
        z[jr] = model->tau[i] * C_ON_FOURPI / model->f[i];
        z[ji] = 0.0;
    }
    imax = model->nLpad - if0 - model->ngrid;
    for (i = if0 + model->ngrid; i <= imax; ++i) {
        jr = i << 1;
        ji = jr + 1;
        z[jr] = 0.0;
        z[ji] = 0.0;
    }
    for (i = model->ngrid - 1; i > 0; --i) {
        jr = (model->nLpad - if0 - i) << 1;
        ji = jr + 1;
        z[jr] = -model->tau[i] * C_ON_FOURPI / model->f[i];
        z[ji] = 0.0;
    }
    for (i = model->nLpad - if0; i < model->nLpad; ++i) {
        jr = i << 1;
        ji = jr + 1;
        z[jr] = 0.0;
        z[ji] = 0.0;
    }
    /*
     * Compute the delay spectrum L = integral[(nr-1) * dz], which
     * is a Hilbert transform of the integrated imaginary part.
     */
    hilbert(z, (unsigned long)model->nLpad);
    /*
     * Combine the delay spectrum and the optical delay.
     */
    for (i = 0; i < model->ngrid; ++i) {
        jr = (i + if0) << 1;
        ji = jr + 1;
        model->L[i] = z[ji] + Lopt;
    }
    free(z);
    return 0;
}   /* compute_delay_spectrum() */


/***********************************************************
* void compute_free_space_loss(model_t *model)
*
* Purpose:
*   Computes the natural log of the free space path loss.
*
* Arguments:
*   model_t *model - pointer to model structure
************************************************************/

void compute_free_space_loss(model_t *model)
{
    if (model->path_begin == 0 || model->path_end == 0) {
        /*
         * If one path endpoint is at infinity, the free space
         * path loss is not meaningful.  Log an error and write
         * zeros into tau_fsl.
         */
        int subgrid;
        errlog(194, 0);
        for (subgrid = 0; subgrid <= 1; ++subgrid) {
            gridsize_t i;
            gridsize_t imin = model->isub[subgrid].min;
            gridsize_t imax = model->isub[subgrid].max;
            for (i = imin; i <= imax; ++i)
                model->tau_fsl[i] = 0.0;
        }
    } else {
        int subgrid;
        double d = source_to_obs_path_distance(model);
        for (subgrid = 0; subgrid <= 1; ++subgrid) {
            gridsize_t i;
            gridsize_t imin = model->isub[subgrid].min;
            gridsize_t imax = model->isub[subgrid].max;
            for (i = imin; i <= imax; ++i) {
                double f   = model->f[i];
                double fsl = FOURPI_ON_C * f * d;
                fsl *= fsl;
                if (fsl < 1.0) {
                    /*
                     * (4 pi d) < lambda.  Clamp fsl at 1
                     * (i.e. tau_fsl = 0), and log a warning
                     * unless f == 0.0.
                     */
                    model->tau_fsl[i] = 0.0;
                    if (f != 0.0)
                        errlog(195, 0);
                } else {
                    model->tau_fsl[i] = log(fsl);
                }
            }
        }
    }
    return;
}   /* compute_free_space_loss() */


/***********************************************************
* int compute_if_spectra(model_t *model)
*
* Purpose:
*   Computes IF spectra from spectra computed on the model
*   frequency grid.
*
* Arguments:
*   model_t *model - pointer to model structure
*
* Return:
*   0 if successful, 1 otherwise
************************************************************/

int compute_if_spectra(model_t *model)
{
    gridsize_t i;

    if (!model->ifmode)
        return 0;
    for (i = 0; i < OUTPUT_END_OF_TABLE; ++i) {
        if ((output[i].flags & OUTPUT_ACTIVE) && (i != OUTPUT_FREQUENCY)) {
            double *s_mod;  /* model spectrum, on model frequency grid */
            double *s_if;   /* if spectrum, on if frequency grid       */
            gridsize_t j;
            s_mod = output[i].spectrum;
            s_if = model->fif;  /* use fif array for scratch space */
            if (model->ifmode & IFMODE_DSB) {
                double su, sl, wtu, wtl;
                gridsize_t ju, jl;
                wtu = model->dsb_utol_ratio / (1.0 + model->dsb_utol_ratio);
                wtl = 1.0 - wtu;
                for (j = 0,
                    ju = model->iusb.min, jl = model->ilsb.max;
                    j < model->nif;
                    ++j, ++ju, --jl) {
                    su = s_mod[ju]
                        + model->fif_delta * (s_mod[ju+1] - s_mod[ju]);
                    sl = s_mod[jl]
                        + model->fif_delta * (s_mod[jl+1] - s_mod[jl]);
                    s_if[j] = wtu * su + wtl * sl;
                }
            } else if (model->ifmode & IFMODE_USB) {
                gridsize_t ju;
                for (j = 0, ju = model->iusb.min; j < model->nif; ++j, ++ju) {
                    s_if[j] = s_mod[ju]
                        + model->fif_delta * (s_mod[ju+1] - s_mod[ju]);
                }
            } else if (model->ifmode & IFMODE_LSB) {
                gridsize_t jl;
                for (j = 0, jl = model->ilsb.max; j < model->nif; ++j, --jl) {
                    s_if[j] = s_mod[jl]
                        + model->fif_delta * (s_mod[jl+1] - s_mod[jl]);
                }
            } else {
                errlog(86, __LINE__);
                return 1;
            }
            /* copy IF spectrum array back to model spectrum array */
            for (j = 0; j < model->nif; ++j)
                s_mod[j] = s_if[j];
        }
    }
    /* fill in the fif array */
    for (i = 0; i < model->nif; ++i)
        model->fif[i] = model->fif_0 + i * model->df;
    return 0;
}   /* compute_if_spectra() */


/***********************************************************
* void compute_radiance_difference_spectrum(model_t *model, model_t *lmodel)
*
* Purpose:
*   Computes the radiance difference spectrum, defined as
*  
*    I_diff = I - I_ref
*  
*   where
*  
*    I_ref = B(T_ref)
*  
*   I_diff models the output of spectrometers which
*   effectively measure the difference in radiance between
*   the scene and a reference blackbody load at T = T_ref.
*   Examples are instruments that chop between reference
*   load and scene, and rapid-scan Fourier spectrometers
*   that have one input port of the interferometer
*   terminated by a reference load.
*
* Arguments:
*   model_t *model - pointer to model structure
*   model_t *lmodel - pointer to a model data structure
*       containing scalar data from a prior computation.
************************************************************/

void compute_radiance_difference_spectrum(model_t *model, model_t *lmodel)
{
    int subgrid;

    /*
     * Compute I_ref if needed.
     */
    if (lmodel == NULL ||
        fabs(model->Tref - lmodel->Tref) > (DBL_EPSILON * model->Tref) ||
        compare_spectral_subgrid_ranges(model, lmodel)) {
        for (subgrid = 0; subgrid <= 1; ++subgrid) {
            gridsize_t imin = model->isub[subgrid].min;
            gridsize_t imax = model->isub[subgrid].max;
            B_Planck(
                &model->I_ref[imin],
                &model->f[imin],
                model->df,
                imax - imin + 1,
                model->Tref);
        }
    }
    for (subgrid = 0; subgrid <= 1; ++subgrid) {
        gridsize_t i;
        gridsize_t imin = model->isub[subgrid].min;
        gridsize_t imax = model->isub[subgrid].max;
        for (i = imin; i <= imax; ++i)
            model->I_diff[i] = model->I[i] - model->I_ref[i];
    }
    return;
}   /* compute_radiance_difference_spectrum() */


/***********************************************************
* void compute_spectral_Tsys(model_t *model)
*
* Purpose:
*   Computes a system temperature spectrum.  This is a
*   shifted and scaled version of the Rayleigh-Jeans
*   temperature spectrum, defined as
*
*     Tsys[i] = (Trj[i] + Trx) * rx_gain_factor
*
*   where rx_gain_factor is intended to be used for modeling
*   fluctuations in receiver gain when fitting series of
*   spectra.
*
* Arguments:
*   model_t *model - pointer to model structure
************************************************************/

void compute_spectral_Tsys(model_t *model)
{
    int subgrid;
    for (subgrid = 0; subgrid <= 1; ++subgrid) {
        gridsize_t i;
        gridsize_t imin = model->isub[subgrid].min;
        gridsize_t imax = model->isub[subgrid].max;
        for (i = imin; i <= imax; ++i)
            model->Tsys[i] =
                (model->Trj[i] + model->Trx) * model->rx_gain_factor;
    }
    return;
}   /* compute_spectral_Tsys() */


/***********************************************************
* void compute_spectral_Y_factor(model_t *model)
*
* Purpose:
*   Compute the spectral Y-factor, defined as
*
*     Y[i] = {(Trj[i] + Trx) / [Trj(Tref) + Trx]} * rx_gain_factor
*
*   where Trj(Tref) is a Rayleigh-Jeans brightness
*   temperature computed for a physical reference
*   temperature Tref.  Note that if Y is convolved with an
*   ILS, we are making the implicit assumption here that the
*   frequency dependence of Trj(Tref) is negligible across
*   the width of the ILS.  The receiver noise temperature
*   Trx is assumed to be already expressed in terms of an
*   equivalent Rayleigh-Jeans source brightness temperature.
*
*   As in the case of Tsys, rx_gain_factor is intended to be
*   used for modeling fluctuations in receiver gain when
*   fitting series of spectra.
*
* Arguments:
*   model_t *model - pointer to model structure
************************************************************/

void compute_spectral_Y_factor(model_t *model)
{
    /*
     * If the denominator in the Y-factor is too small, log an error
     * and set the spectrum to zero.
     */
    if ((model->Tref + model->Trx) < DBL_EPSILON) {
        int subgrid;
        errlog(110, 0);
        for (subgrid = 0; subgrid <= 1; ++subgrid) {
            gridsize_t i;
            gridsize_t imin = model->isub[subgrid].min;
            gridsize_t imax = model->isub[subgrid].max;
            for (i = imin; i <= imax; ++i)
                model->Y[i] = 0.0;
        }
    } else {
        int subgrid;
        for (subgrid = 0; subgrid <= 1; ++subgrid) {
            gridsize_t i;
            gridsize_t imin = model->isub[subgrid].min;
            gridsize_t imax = model->isub[subgrid].max;
            for (i = imin; i <= imax; ++i)
                model->Y[i] = model->rx_gain_factor *
                    (model->Trj[i] + model->Trx) /
                    (T_Rayleigh_Jeans(model->f[i], model->Tref) + model->Trx);
        }
    }
    return;
}   /* compute_spectral_Y_factor() */


/***********************************************************
* void compute_Tb(model_t *model)
*
* Purpose:
*   Compute Planck brightness temperature spectrum from the
*   spectral radiance.
*
* Arguments:
*   model_t *model - pointer to model structure
************************************************************/

void compute_Tb(model_t *model)
{
    int subgrid;

    for (subgrid = 0; subgrid <= 1; ++subgrid) {
        gridsize_t i;
        gridsize_t imin = model->isub[subgrid].min;
        gridsize_t imax = model->isub[subgrid].max;
        /*
         * The first point is handled specially, in case the frequency
         * grid starts at zero.  lim Tb, f->0 = T0, so set Tb(f=0) = T0.
         * This prevents ringing if Tb gets convolved with an ILS.
         */
        if (imin == 0) {
            model->Tb[0] = model->f[0] < F_EPSILON ?
                model->T0 : T_Planck(model->f[0], model->I[0]);
            ++imin;
        }
        #ifdef _OPENMP 
        #pragma omp parallel for schedule(guided, 128)\
            if (model->isub[0].max - model->isub[0].min >= 16384)
        #endif
        for (i = imin; i <= imax; ++i)
            model->Tb[i] = T_Planck(model->f[i], model->I[i]);
    }
    return;
}   /* compute_Tb() */


/***********************************************************
* void compute_transmittance(model_t *model)
*
* Purpose:
*   Compute spectral transmittance from opacity.
*
* Arguments:
*   model_t *model - pointer to model structure
************************************************************/

void compute_transmittance(model_t *model)
{
    int subgrid;

    for (subgrid = 0; subgrid <= 1; ++subgrid) {
        gridsize_t i;
        gridsize_t imin = model->isub[subgrid].min;
        gridsize_t imax = model->isub[subgrid].max;
        #ifdef _OPENMP
        #pragma omp parallel for schedule(guided, 128)\
            if (imax - imin >= 16384)
        #endif
        for (i = imin; i <=imax; ++i)
            model->tx[i] = am_exp(-model->tau[i]);
    }
    return;
}   /* compute_transmittance() */


/***********************************************************
* void compute_Trj(model_t *model)
*
* Purpose:
*   Compute Rayleigh-Jeans brightness temperature spectrum
*   from the spectral radiance.
*
* Arguments:
*   model_t *model - pointer to model structure
************************************************************/

void compute_Trj(model_t *model)
{
    int subgrid;

    for (subgrid = 0; subgrid <= 1; ++subgrid) {
        gridsize_t i;
        gridsize_t imin = model->isub[subgrid].min;
        gridsize_t imax = model->isub[subgrid].max;
        /*
         * The first point is handled specially, in case the frequency
         * grid starts at zero.  lim Trj, f->0 = T0, so set Tb(f=0) = T0.
         * This prevents ringing if Trj gets convolved with an ILS.
         */
        if (imin == 0) {
            model->Trj[0] = model->f[0] < F_EPSILON ?
                model->T0 : model->I[0] / (TWOKB_ON_CSQUARED * model->f2[0]);
            ++imin;
        }
        /*
         * This loop was found not to gain significantly from being
         * parallelized.
         */
        for (i = imin; i <= imax; ++i)
            model->Trj[i] = model->I[i] / (TWOKB_ON_CSQUARED * model->f2[i]);
    }
    return;
}   /* compute_Trj() */


/***********************************************************
* static int set_IF_spectrum_subgrid_ranges(model_t *model)
*
* Purpose:
*
*   If IF spectra are being computed, this function sets the
*   ranges of grid points used to compute LSB and/or USB
*   spectra.
*
*   If no IF range has been specified, the IF grid range is
*   set to the maximum range supported by the model grid.
*   Otherwise, it is set to the intersection of the
*   specified IF grid and the range supported by the model
*   grid.
*
*   The model frequency grid is always aligned to f = 0,
*   whereas the IF frequency grid will be aligned to
*   fif = 0.  That is, the origin of the IF frequency grid
*   is flo, with the same grid interval df as the model
*   grid:
*
*                                           f = flo
*                ->|       |<- df             |
*   f   <- |       |       |       |       |       |       | ->
*   fif <-    |       |       |       |       |       |      ->
*                ->|  |<- df * fif_delta      |
*                                           fif = 0
*
*   Here, positive (USB) and negative (LSB) IF frequencies
*   are shown.  However, the IF frequency grid is always
*   output as positive.
*
*   Besides the grid range, this function sets three
*   additional parameters in the model structure:
*
*     model->fif_delta - normalized offset between a point
*       on the IF grid, and the model grid point at or below
*       the IF grid point.
*     model->fif_0 -  the lowest IF frequency
*     model->nif - the number of IF frequency grid points.
*
* Arguments:
*   model_t *model - pointer to model structure
*
* Return:
*   0 if successful, 1 otherwise
************************************************************/

int set_IF_spectrum_subgrid_ranges(model_t *model)
{
    double lsb_fmin = 0.0, lsb_fmax = 0.0;
    double usb_fmin = 0.0, usb_fmax = 0.0;
    int fif_defined = (model->fif_max >= 0.0);

    /*
     * At least two grid points are needed for IF spectra
     */
    if (model->ngrid < 2) {
        errlog(88, 0);
        return 1;
    }
    /*
     * LSB frequency range
     */
    if (model->ifmode & (IFMODE_LSB | IFMODE_DSB)) {
        lsb_fmax = model->flo - model->f[0];
        if (lsb_fmax < 0.0) {
            errlog(85, 0);
            return 1;
        }
        if (fif_defined && (lsb_fmax > model->fif_max * (1.0 + DBL_EPSILON)))
            lsb_fmax = model->fif_max * (1.0 + DBL_EPSILON);
        lsb_fmin = model->flo - model->f[model->ngrid - 1];
        if (lsb_fmin < 0.0) 
            lsb_fmin = 0.0;
        if (fif_defined && (lsb_fmin < model->fif_min)) {
            if (model->fif_min > lsb_fmax) {
                errlog(85, 0);
                return 1;
            }
            lsb_fmin = model->fif_min;
        }
    }
    /*
     * USB frequency range.
     */
    if (model->ifmode & (IFMODE_USB | IFMODE_DSB)) {
        usb_fmax = model->f[model->ngrid - 1] - model->flo;
        if (usb_fmax < 0.0) {
            errlog(85, 0);
            return 1;
        }
        if (fif_defined && (usb_fmax > model->fif_max * (1.0 + DBL_EPSILON)))
            usb_fmax = model->fif_max * (1.0 + DBL_EPSILON);
        usb_fmin = model->f[0] - model->flo;
        if (usb_fmin < 0.0)
            usb_fmin = 0.0;
        if (fif_defined && (usb_fmin < model->fif_min)) {
            if (model->fif_min > usb_fmax) {
                errlog(85, 0);
                return 1;
            }
            usb_fmin = model->fif_min;
        }
    }
    /*
     * For DSB, find the overlap between the LSB and USB ranges.
     */
    if (model->ifmode & (IFMODE_DSB)) {
        if (lsb_fmax < usb_fmax)
            usb_fmax = lsb_fmax;
        else
            lsb_fmax = usb_fmax;
        if (lsb_fmin > usb_fmin)
            usb_fmin = lsb_fmin;
        else
            lsb_fmin = usb_fmin;
    }
    /*
     * Model grid indices covering the LSB and USB spectra.  The
     * ranges ilsb and iusb cover the model grid points just below
     * each corresponding IF grid point.  Note that, in principle,
     * there could be an index rounding inconsistency, in the DSB
     * case, between LSB and USB, but this would occur near
     * fif_delta = 0.0, so the resulting interpolation error would be
     * negligible.
     */
    if (model->ifmode & (IFMODE_LSB | IFMODE_DSB)) {
        double f;
        model->fif_0 =
            model->df * ceil((lsb_fmin * (1.0 - DBL_EPSILON)) / model->df);
        model->nif = (gridsize_t)ceil(((lsb_fmax - model->fif_0) *
            (1.0 - DBL_EPSILON)) / model->df);
        f = model->flo - model->fif_0;
        model->ilsb.max = (gridsize_t)floor(((f - model->f[0]) *
            (1.0 + DBL_EPSILON)) / model->df);
        model->fif_delta = (f - model->f[model->ilsb.max]) / model->df;
        model->ilsb.min = model->ilsb.max - (model->nif - 1);
    }
    if (model->ifmode & (IFMODE_USB | IFMODE_DSB)) {
        double f;
        model->fif_0 = model->df *
            ceil((usb_fmin * (1.0 - DBL_EPSILON)) / model->df);
        model->nif = (gridsize_t)ceil(((usb_fmax - model->fif_0) *
            (1.0 - DBL_EPSILON)) / model->df);
        f = model->flo + model->fif_0;
        model->iusb.min = (gridsize_t)floor(((f - model->f[0]) *
            (1.0 + DBL_EPSILON)) / model->df);
        model->fif_delta = (f - model->f[model->iusb.min]) / model->df;
        model->iusb.max = model->iusb.min + (model->nif - 1);
    }
    return 0;
}   /* set_IF_spectrum_subgrid_ranges() */


/***********************************************************
* int set_spectral_subgrid_ranges(model_t *model, model_t *lmodel)
*
* Purpose:
*
*   This function sets ranges for two subgrids of the model
*   frequency grid.  For normal, LSB, and USB spectra, one
*   of these subgrid ranges is empty.  For DSB spectra, if
*   the USB and LSB ranges do not overlap, the two subgrid
*   ranges correspond to USB and LSB.  If they do overlap,
*   they are consolidated into a single subgrid.
*
*   Note that the IF sideband ranges computed in
*   set_IF_spectrum_subgrid_ranges() are the range of model
*   grid indices for the model grid points just at or below
*   each corresponding IF frequency point, so the ranges set
*   here include an extra point at the top for
*   interpolation.
*
*   These abbreviated sub ranges are used to accelerate
*   computations at the layer and column levels, which
*   typically occur repetitively during fits or Jacobian
*   computations.  (This optimization is most important for
*   DSB spectra.)
*
*   If an ILS is being applied, or if a delay spectrum is
*   being computed, the model will always be computed over
*   the entire frequency grid; in these cases any difference
*   between the model and IF grid ranges is assumed to be
*   spectral padding desired by the user.
*
* Arguments:
*   model_t *model  - pointer to model data structure
*   model_t *lmodel - pointer to a model data structure
*       containing scalar data from a prior computation.
*
* Return:
*   0 if successful, 1 otherwise
************************************************************/

int set_spectral_subgrid_ranges(model_t *model, model_t *lmodel)
{
    /*
     * For IF spectra, initialize or adjust subgrid ranges as needed.
     */
    if (model->ifmode) {
        if ((lmodel == NULL) ||
            (fabs(model->flo - lmodel->flo)
                > DBL_EPSILON * model->flo) ||
            (fabs(model->fif_min - lmodel->fif_min)
                > DBL_EPSILON * model->fif_min) ||
            (fabs(model->fif_max - lmodel->fif_max)
                > DBL_EPSILON * model->fif_max)) {
            if (set_IF_spectrum_subgrid_ranges(model))
                return 1;
        }
    }
    if ((model->ifmode == 0) || (model->ils != NULL) || (model->L != NULL)) {
        model->isub[0].min = 0;
        model->isub[0].max = model->ngrid - 1;
        model->isub[1].min = 0;
        model->isub[1].max = -1;
    } else if (model->ifmode & IFMODE_LSB) {
        model->isub[0].min = model->ilsb.min;
        model->isub[0].max = model->ilsb.max + 1;
        model->isub[1].min = 0;
        model->isub[1].max = -1;
    } else if (model->ifmode & IFMODE_USB) {
        model->isub[0].min = model->iusb.min;
        model->isub[0].max = model->iusb.max + 1;
        model->isub[1].min = 0;
        model->isub[1].max = -1;
    } else if (model->ifmode & IFMODE_DSB) {
        if (model->ilsb.max + 1 >= model->iusb.min) {
            model->isub[0].min = model->ilsb.min;
            model->isub[0].max = model->iusb.max + 1;
            model->isub[1].min = 0;
            model->isub[1].max = -1;
        } else {
            model->isub[0].min = model->ilsb.min;
            model->isub[0].max = model->ilsb.max + 1;
            model->isub[1].min = model->iusb.min;
            model->isub[1].max = model->iusb.max + 1;
        }
    } else {
        errlog(86, __LINE__);
        return 1;
    }
    return 0;
}   /* set_spectral_subgrid_ranges() */
