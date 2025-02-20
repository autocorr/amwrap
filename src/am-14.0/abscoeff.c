/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* abscoeff.c                     S. Paine rev. 2024 April 20
*
* Absorption coefficient computation and cache lookups.
************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>

#include "abscoeff.h"
#include "am_alloc.h"
#include "am_types.h"
#include "cia.h"
#include "dcache.h"
#include "errlog.h"
#include "h2o_continuum.h"
#include "h2o_ice.h"
#include "h2o_liquid.h"
#include "kcache.h"
#include "linesum.h"
#include "output.h"
#include "phys_const.h"
#include "units.h"

#include "ch4.h"
#include "ch3cn.h"
#include "ch3oh.h"
#include "clo.h"
#include "co.h"
#include "co2.h"
#include "hbr.h"
#include "hcl.h"
#include "hcn.h"
#include "hf.h"
#include "h2co.h"
#include "hno3.h"
#include "h2o.h"
#include "h2o2.h"
#include "ho2.h"
#include "hocl.h"
#include "h2s.h"
#include "nh3.h"
#include "n2o.h"
#include "no.h"
#include "no2.h"
#include "o.h"
#include "o2.h"
#include "o3.h"
#include "ocs.h"
#include "oh.h"
#include "so2.h"
#include "oneline.h"

/*
 * Any  additions to  this table  between program  version number
 * changes  must  be  made  at  the  end.  This  is because  the
 * enum  value  in  abscoeff.h associated  with each  absorption
 * coefficient type  is part  of the header data  for absorption
 * coefficient arrays that are cached to disk.
 */
struct k_type_tabentry k_type[] = {
    {"none",
        UNIT_NONE, 0, LINESHAPE_NONE, 0, 0},
    {"ch4",
        UNIT_CM2,  0, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"12ch4",
        UNIT_CM2,  1, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"13ch4",    
        UNIT_CM2,  2, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"12ch3d",  
        UNIT_CM2,  3, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"ch3cn",   
        UNIT_CM2,  0, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"12ch3_12c14n",
        UNIT_CM2,  1, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"ch3oh",   
        UNIT_CM2,  0, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"12ch3_16oh",
        UNIT_CM2,  1, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"co",      
        UNIT_CM2,  0, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"12c_16o", 
        UNIT_CM2,  1, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"13c_16o", 
        UNIT_CM2,  2, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"12c_18o", 
        UNIT_CM2,  3, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"12c_17o", 
        UNIT_CM2,  4, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"13c_18o", 
        UNIT_CM2,  5, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"13c_17o", 
        UNIT_CM2,  6, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"co2",     
        UNIT_CM2,  0, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"12c_16o2",
        UNIT_CM2,  1, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"13c_16o2",
        UNIT_CM2,  2, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"16o_12c_18o",
        UNIT_CM2,  3, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"16o_12c_17o",
        UNIT_CM2,  4, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"16o_13c_18o",
        UNIT_CM2,  5, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"16o_13c_17o",
        UNIT_CM2,  6, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"12c_18o2",
        UNIT_CM2,  7, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"clo",     
        UNIT_CM2,  0, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"35cl_16o",
        UNIT_CM2,  1, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"37cl_16o",
        UNIT_CM2,  2, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"hbr",     
        UNIT_CM2,  0, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"h_79br",  
        UNIT_CM2,  1, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"h_81br",  
        UNIT_CM2,  2, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"hcn",     
        UNIT_CM2,  0, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"h_12c_14n",
        UNIT_CM2,  1, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"h_13c_14n",
        UNIT_CM2,  2, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"h_12c_15n",
        UNIT_CM2,  3, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"h2co",    
        UNIT_CM2,  0, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"h2_12c_16o",
        UNIT_CM2,  1, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"h2_13c_16o",
        UNIT_CM2,  2, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"h2_12c_18o",
        UNIT_CM2,  3, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"hcl",     
        UNIT_CM2,  0, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"h_35cl",  
        UNIT_CM2,  1, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"h_37cl",  
        UNIT_CM2,  2, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"hf",      
        UNIT_CM2,  0, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"h_19f",   
        UNIT_CM2,  1, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"hno3",    
        UNIT_CM2,  0, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"h_14n_16o3",
        UNIT_CM2,  1, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"h2o_lines",
        UNIT_CM2,  0, LINESHAPE_VVH_750, CACHEABLE, LINESUM_DEP_FLAGS},
    {"h2_16o",  
        UNIT_CM2,  1, LINESHAPE_VVH_750, CACHEABLE, LINESUM_DEP_FLAGS},
    {"h2_18o",  
        UNIT_CM2,  2, LINESHAPE_VVH_750, CACHEABLE, LINESUM_DEP_FLAGS},
    {"h2_17o",  
        UNIT_CM2,  3, LINESHAPE_VVH_750, CACHEABLE, LINESUM_DEP_FLAGS},
    {"hd_16o",  
        UNIT_CM2,  4, LINESHAPE_VVH_750, CACHEABLE, LINESUM_DEP_FLAGS},
    {"hd_18o",  
        UNIT_CM2,  5, LINESHAPE_VVH_750, CACHEABLE, LINESUM_DEP_FLAGS},
    {"hd_17o",  
        UNIT_CM2,  6, LINESHAPE_VVH_750, CACHEABLE, LINESUM_DEP_FLAGS},
    {"h2o_air_continuum",
        UNIT_CM5,  0, LINESHAPE_NONE, 0, DEP_ON_T},
    {"h2o_self_continuum",
        UNIT_CM5,  0, LINESHAPE_NONE, 0, DEP_ON_T},
    {"h2o2",    
        UNIT_CM2,  0, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"h2_16o2", 
        UNIT_CM2,  1, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"ho2",     
        UNIT_CM2,  0, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"h_16o2",  
        UNIT_CM2,  1, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"hocl",    
        UNIT_CM2,  0, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"h_16o_35cl",
        UNIT_CM2,  1, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"h_16o_37cl",
        UNIT_CM2,  2, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"h2s",     
        UNIT_CM2,  0, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"h2_32s",  
        UNIT_CM2,  1, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"h2_34s",  
        UNIT_CM2,  2, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"h2_33s",  
        UNIT_CM2,  3, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"iwp_abs_Rayleigh",
        UNIT_CM2,  0, LINESHAPE_NONE, 0, DEP_ON_T},
    {"lwp_abs_Rayleigh",
        UNIT_CM2,  0, LINESHAPE_NONE, 0, DEP_ON_T},
    {"nh3",     
        UNIT_CM2,  0, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"14nh3",   
        UNIT_CM2,  1, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"15nh3",   
        UNIT_CM2,  2, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"n2n2",    
        UNIT_CM5,  0, LINESHAPE_NONE, CACHEABLE, DEP_ON_T},
    {"n2o",     
        UNIT_CM2,  0, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"14n2_16o",
        UNIT_CM2,  1, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"14n_15n_16o",
        UNIT_CM2,  2, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"15n_14n_16o",
        UNIT_CM2,  3, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"14n2_18o",
        UNIT_CM2,  4, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"14n2_17o",
        UNIT_CM2,  5, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"no",      
        UNIT_CM2,  0, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"14n_16o", 
        UNIT_CM2,  1, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"no2",     
        UNIT_CM2,  0, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"14n_16o2",
        UNIT_CM2,  1, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"o",       
        UNIT_CM2,  0, LINESHAPE_VOIGT_KIELKOPF, CACHEABLE, LINESUM_DEP_FLAGS},
    {"16o",     
        UNIT_CM2,  1, LINESHAPE_VOIGT_KIELKOPF, CACHEABLE, LINESUM_DEP_FLAGS},
    {"o2_coupled",
        UNIT_CM2,  0, LINESHAPE_VVW_COUPLED, CACHEABLE, LINESUM_DEP_FLAGS},
    {"16o2_coupled",
        UNIT_CM2,  1, LINESHAPE_VVW_COUPLED, CACHEABLE, LINESUM_DEP_FLAGS},
    {"o2_uncoupled",
        UNIT_CM2,  0, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"16o2_uncoupled",
        UNIT_CM2,  1, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"16o_18o_uncoupled",
        UNIT_CM2,  2, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"16o_17o_uncoupled",
        UNIT_CM2,  3, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"o2o2",    
        UNIT_CM5,  0, LINESHAPE_NONE, CACHEABLE, DEP_ON_T},
    {"o3",      
        UNIT_CM2,  0, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"16o3",    
        UNIT_CM2,  1, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"16o2_18o",
        UNIT_CM2,  2, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"16o_18o_16o",
        UNIT_CM2,  3, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"16o2_17o",
        UNIT_CM2,  4, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"16o_17o_16o",
        UNIT_CM2,  5, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"ocs",     
        UNIT_CM2,  0, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"16o_12c_32s",
        UNIT_CM2,  1, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"16o_12c_34s",
        UNIT_CM2,  2, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"16o_13c_32s",
        UNIT_CM2,  3, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"16o_12c_33s",
        UNIT_CM2,  4, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"18o_12c_32s",
        UNIT_CM2,  5, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"oh",      
        UNIT_CM2,  0, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"16o_h",   
        UNIT_CM2,  1, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"18o_h",   
        UNIT_CM2,  2, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"16o_d",   
        UNIT_CM2,  3, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"so2",     
        UNIT_CM2,  0, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"32s_16o2",
        UNIT_CM2,  1, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"34s_16o2",
        UNIT_CM2,  2, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"oneline", 
        UNIT_CM2,  0, LINESHAPE_GROSS, CACHEABLE, LINESUM_DEP_FLAGS},
    {"END",     
        UNIT_NONE, 0, LINESHAPE_NONE, 0, 0},
};

static int compute_absorption_coefficient(
        model_t*,
        int,
        int,
        int,
        double*,
        double);


/***********************************************************
* int get_absorption_coefficient(
*         model_t *model,
*         model_t *lmodel,
*         int lnum,
*         int cnum,
*         int knum)
*
* Purpose:
*   Fills a request for absorption coefficient knum, of
*   column cnum, of layer lnum of the model, by either
*
*     (1) finding it in the disk cache (dcache),
*
*     (2) interpolating between previously computed values
*         at two bracketing temperatures in the memory cache
*         (kcache), or
*
*     (3) calling compute_absorption_coefficient().
*
*     In case (3), the disk cache, memory cache, or both are
*     updated with the computed result.
*
* Arguments:
*   model_t *model  - pointer to model structure
*   model_t *lmodel - pointer to a model structure
*                     containing scalar data from a prior
*                     computation, if any.  NULL otherwise.
*   int lnum        - layer number
*   int cnum        - column number
*   int knum        - absorption coefficient number
*
* Return:
*   0 if OK, 1 if an error occurred
************************************************************/

int get_absorption_coefficient(
        model_t *model,
        model_t *lmodel,
        int lnum,
        int cnum,
        int knum)
{
    int bypass_dcache;
    int no_kcache;
    int outside_kcache_Trange;
    layer_t    *layer    = model->layer[lnum];
    column_t   *column   = layer->column[cnum];
    abscoeff_t *abscoeff = column->abscoeff[knum];

    /*
     * If k[] for this column is not being cached in memory
     * (kcache), or if the temperature is outside the cachable
     * range, then look it up in the disk cache (dcache).  If it
     * isn't in the dcache, compute it and save it in the dcache.
     * Note that the disk cache lookup and save functions do
     * nothing if no disk cache directory has been specified in
     * the environment.
     *
     * During a fit, if this is an update computation of a
     * previously-computed model (i.e. lmodel != NULL), then
     * absorption coefficients are only cached to disk if they
     * are kcache entries, which are computed on a discrete T
     * grid.  Otherwise the dcache is bypassed.  The rationale
     * for this is that caching absorption coefficients which
     * depend on continuously-varying fit variables would produce
     * a large number of cache misses, evicting useful cache data
     * and replacing it with data unlikely to be needed again.  A
     * fit in progress is recognized by a 0 in the OUTPUT_ACTIVE
     * bit of output[ALL_OUTPUTS].flags.
     */
    if (abscoeff->kcache == NULL) {
        no_kcache = 1;
    } else {
        no_kcache = 0;
        /*
         * The absorption coefficients cached in the kcache
         * depend on pressure and mixing ratio.  If these have
         * changed, any affected data needs to be flushed.
         * In most cases, a kcache will not have been created
         * when these variables are subject to change; exceptions
         * are when an interpolated level position is a fit
         * variable, and when changes are made during during
         * streaming fits using the "set" facility.)
         */
        if (lmodel != NULL) {
            layer_t    *llayer    = lmodel->layer[lnum];
            column_t   *lcolumn   = llayer->column[cnum];
            abscoeff_t *labscoeff = lcolumn->abscoeff[knum];
            int        k_typenum  = abscoeff->k_typenum;
            /*
             * If the layer pressure has changed, flush all
             * kcache data on this layer.  Otherwise, if the
             * mixing ratio associated with this particular
             * absorption coefficient has changed, flush the
             * kcache for this absorption coefficient only.
             */
            if ((k_type[k_typenum].dep_flags & DEP_ON_P) &&
                    (fabs(layer->P -
                          llayer->P) >
                     (DBL_EPSILON * layer->P))) {
                clear_layer_kcache_entries(model, lnum);
            } else if ((k_type[k_typenum].dep_flags & DEP_ON_VMR_SELFBROAD) &&
                    (fabs(abscoeff->vmr_selfbroad -
                          labscoeff->vmr_selfbroad) >
                     (DBL_EPSILON * abscoeff->vmr_selfbroad))) {
                int j;
                for (j = 0; j < model->nkcache; ++j) {
                    kcache_free(abscoeff->kcache[j]);
                    abscoeff->kcache[j] = NULL;
                }
            }
        }
    }
    outside_kcache_Trange =
        (layer->T < model->kcache_Tmin) ||
        (layer->T
         > (model->kcache_Tmin + (model->nkcache - 1) * model->kcache_dT));
    bypass_dcache =
        (lmodel != NULL) &&
        (abscoeff->kcache == NULL) &&
        !(output[ALL_OUTPUTS].flags & OUTPUT_ACTIVE);
    if (bypass_dcache) {
        if (compute_absorption_coefficient(
                    model,
                    lnum,
                    cnum,
                    knum,
                    abscoeff->k,
                    layer->T)) {
            return 1;
        }
    } else if (no_kcache || outside_kcache_Trange) {
        if (dcache_lookup_absorption_coefficient(
                    model,
                    lnum,
                    cnum,
                    knum,
                    abscoeff->k,
                    layer->T)) {
            if (compute_absorption_coefficient(
                        model,
                        lnum,
                        cnum,
                        knum,
                        abscoeff->k,
                        layer->T)) {
                return 1;
            }
            dcache_save_absorption_coefficient(
                    model,
                    lnum,
                    cnum,
                    knum,
                    abscoeff->k,
                    layer->T);
        }
    } else {
        /*
         * Interpolate the absorption coefficient between a pair
         * of adjacent entries in the kcache at the two grid
         * temperatures bracketing the requested one.  If either
         * or both of these has not yet been computed, allocate a
         * kcache entry and either get the coefficient from the
         * disk cache or compute it.
         *
         * Note that during a kcache_alloc(), the other
         * coefficient of the pair is locked in memory by a call
         * to kcache_lock_entry(). This is in case kcache memory
         * has run out, forcing kcache_alloc() to throw out old
         * entries to make space for new ones.
         */
        double u;
        double Tlo;
        double *k, *klo, *khi;
        gridsize_t i;
        unsigned int m;
        u = 1.0 / model->kcache_dT;
        m = (unsigned int)(u * (layer->T - model->kcache_Tmin));
        Tlo = model->kcache_Tmin + (double)m * model->kcache_dT;
        if (abscoeff->kcache[m] == NULL) {
            kcache_log(KCACHE_MISS);
            kcache_lock_entry(abscoeff->kcache[m+1]);
            if ((abscoeff->kcache[m] = kcache_alloc(model)) == NULL) {
                errlog(54, 0);
                return 1;
            }
            kcache_unlock_entry(abscoeff->kcache[m+1]);
            if (dcache_lookup_absorption_coefficient(
                        model,
                        lnum,
                        cnum,
                        knum,
                        abscoeff->kcache[m],
                        Tlo)) {
                if (compute_absorption_coefficient(
                            model,
                            lnum,
                            cnum,
                            knum,
                            abscoeff->kcache[m],
                            Tlo)) {
                    return 1;
                }
                dcache_save_absorption_coefficient(
                        model,
                        lnum,
                        cnum,
                        knum,
                        abscoeff->kcache[m],
                        Tlo);
            }
        } else {
            kcache_log(KCACHE_HIT);
        }
        if (abscoeff->kcache[m+1] == NULL) {
            kcache_log(KCACHE_MISS);
            kcache_lock_entry(abscoeff->kcache[m]);
            if ((abscoeff->kcache[m+1] = kcache_alloc(model)) == NULL) {
                errlog(54, 0);
                return 1;
            }
            kcache_unlock_entry(abscoeff->kcache[m]);
            if (dcache_lookup_absorption_coefficient(
                        model,
                        lnum,
                        cnum,
                        knum,
                        abscoeff->kcache[m+1],
                        Tlo + model->kcache_dT)) {
                if (compute_absorption_coefficient(
                            model,
                            lnum,
                            cnum,
                            knum,
                            abscoeff->kcache[m+1],
                            Tlo + model->kcache_dT)) {
                    return 1;
                }
                dcache_save_absorption_coefficient(
                        model,
                        lnum,
                        cnum,
                        knum,
                        abscoeff->kcache[m+1],
                        Tlo + model->kcache_dT);
            }
        } else {
            kcache_log(KCACHE_HIT);
        }
        /*
         * Now do the interpolation.
         */
        k   = abscoeff->k;
        klo = abscoeff->kcache[m];
        khi = abscoeff->kcache[m+1];
        u  *= layer->T - Tlo;
        for (i = 0; i < model->ngrid; ++i)
            k[i] = klo[i] + u * (khi[i] - klo[i]);
    }
    if (abscoeff->unres_lines != 0)
        errlog(5, 0);
    return 0;
}   /* get_absorption_coefficient() */


/***********************************************************
* static int compute_absorption_coefficient(
*         model_t *model,
*         int lnum,
*         int cnum,
*         int knum,
*         double *k,
*         double T)
*
* Purpose:
*   Computes absorption coefficient knum, of column cnum, of
*   layer lnum of the model.  The address k of the
*   absorption coefficient array, and the temperature T, are
*   passed as separate parameters to accommodate computation
*   of cached absorption coefficients on the fixed T grid
*   used for the kcache.
*
* Arguments:
*   model_t *model - pointer to model structure
*   int lnum       - layer number
*   int cnum       - column number
*   int knum       - absorption coefficient number
*   double *k      - pointer to absorption coefficient array
*                    to be computed
*   double T       - temperature [K]
*
* Return:
*   0 if OK, 1 if an error occurred
************************************************************/

static int compute_absorption_coefficient(
        model_t *model,
        int lnum,
        int cnum,
        int knum,
        double *k,
        double T)
{
    layer_t *layer;
    abscoeff_t *abscoeff;
    static double Qrat[MAX_ISO];

    layer    = model->layer[lnum];
    abscoeff = layer->column[cnum]->abscoeff[knum];

    switch (abscoeff->k_typenum) {
        case K_TYPE_CH4:
        case K_TYPE_12CH4:
        case K_TYPE_13CH4:
        case K_TYPE_12CH3D:
            if (Qratio( T,
                        ch4_Tref,
                        ch4_Qtab,
                        ch4_Qtab_rows,
                        ch4_Qtab_cols,
                        Qrat)) {
                errlog(7, 0);
                return 1;
            }
            if (linesum(k,
                        model->f,
                        model->f2,
                        model->df,
                        model->ngrid,
                        ch4_cat,
                        NULL,
                        ch4_num_lines,
                        layer->lineshape[abscoeff->k_typenum],
                        k_type[abscoeff->k_typenum].iso,
                        layer->P,
                        abscoeff->vmr_selfbroad,
                        T,
                        ch4_Tref,
                        Qrat,
                        ch4_abundance_tab,
                        ch4_mass_tab,
                        model->tol,
                        &(abscoeff->unres_lines))) {
                return 1;
            }
            break;
        case K_TYPE_CH3CN:
        case K_TYPE_12CH3_12C14N:
            if (Qratio( T,
                        ch3cn_Tref,
                        ch3cn_Qtab,
                        ch3cn_Qtab_rows,
                        ch3cn_Qtab_cols,
                        Qrat)) {
                errlog(203, 0);
                return 1;
            }
            if (linesum(k,
                        model->f,
                        model->f2,
                        model->df,
                        model->ngrid,
                        ch3cn_cat,
                        NULL,
                        ch3cn_num_lines,
                        layer->lineshape[abscoeff->k_typenum],
                        k_type[abscoeff->k_typenum].iso,
                        layer->P,
                        abscoeff->vmr_selfbroad,
                        T,
                        ch3cn_Tref,
                        Qrat,
                        ch3cn_abundance_tab,
                        ch3cn_mass_tab,
                        model->tol,
                        &(abscoeff->unres_lines))) {
                return 1;
            }
            break;
        case K_TYPE_CH3OH:
        case K_TYPE_12CH3_16OH:
            if (Qratio( T,
                        ch3oh_Tref,
                        ch3oh_Qtab,
                        ch3oh_Qtab_rows,
                        ch3oh_Qtab_cols,
                        Qrat)) {
                errlog(154, 0);
                return 1;
            }
            if (linesum(k,
                        model->f,
                        model->f2,
                        model->df,
                        model->ngrid,
                        ch3oh_cat,
                        NULL,
                        ch3oh_num_lines,
                        layer->lineshape[abscoeff->k_typenum],
                        k_type[abscoeff->k_typenum].iso,
                        layer->P,
                        abscoeff->vmr_selfbroad,
                        T,
                        ch3oh_Tref,
                        Qrat,
                        ch3oh_abundance_tab,
                        ch3oh_mass_tab,
                        model->tol,
                        &(abscoeff->unres_lines))) {
                return 1;
            }
            break;
        case K_TYPE_CO:
        case K_TYPE_12C_16O:
        case K_TYPE_13C_16O:
        case K_TYPE_12C_18O:
        case K_TYPE_12C_17O:
        case K_TYPE_13C_18O:
        case K_TYPE_13C_17O:
            if (Qratio( T,
                        co_Tref,
                        co_Qtab,
                        co_Qtab_rows,
                        co_Qtab_cols,
                        Qrat)) {
                errlog(9, 0);
                return 1;
            }
            if (linesum(k,
                        model->f,
                        model->f2,
                        model->df,
                        model->ngrid,
                        co_cat,
                        NULL,
                        co_num_lines,
                        layer->lineshape[abscoeff->k_typenum],
                        k_type[abscoeff->k_typenum].iso,
                        layer->P,
                        abscoeff->vmr_selfbroad,
                        T,
                        co_Tref,
                        Qrat,
                        co_abundance_tab,
                        co_mass_tab,
                        model->tol,
                        &(abscoeff->unres_lines))) {
                return 1;
            }
            break;
        case K_TYPE_CO2:
        case K_TYPE_12C_16O2:
        case K_TYPE_13C_16O2:
        case K_TYPE_16O_12C_18O:
        case K_TYPE_16O_12C_17O:
        case K_TYPE_16O_13C_18O:
        case K_TYPE_16O_13C_17O:
        case K_TYPE_12C_18O2:
            if (Qratio( T,
                        co2_Tref,
                        co2_Qtab,
                        co2_Qtab_rows,
                        co2_Qtab_cols,
                        Qrat)) {
                errlog(132, 0);
                return 1;
            }
            if (linesum(k,
                        model->f,
                        model->f2,
                        model->df,
                        model->ngrid,
                        co2_cat,
                        NULL,
                        co2_num_lines,
                        layer->lineshape[abscoeff->k_typenum],
                        k_type[abscoeff->k_typenum].iso,
                        layer->P,
                        abscoeff->vmr_selfbroad,
                        T,
                        co2_Tref,
                        Qrat,
                        co2_abundance_tab,
                        co2_mass_tab,
                        model->tol,
                        &(abscoeff->unres_lines))) {
                return 1;
            }
            break;
        case K_TYPE_ClO:
        case K_TYPE_35Cl_16O:
        case K_TYPE_37Cl_16O:
            if (Qratio( T,
                        clo_Tref,
                        clo_Qtab,
                        clo_Qtab_rows,
                        clo_Qtab_cols,
                        Qrat)) {
                errlog(150, 0);
                return 1;
            }
            if (linesum(k,
                        model->f,
                        model->f2,
                        model->df,
                        model->ngrid,
                        clo_cat,
                        NULL,
                        clo_num_lines,
                        layer->lineshape[abscoeff->k_typenum],
                        k_type[abscoeff->k_typenum].iso,
                        layer->P,
                        abscoeff->vmr_selfbroad,
                        T,
                        clo_Tref,
                        Qrat,
                        clo_abundance_tab,
                        clo_mass_tab,
                        model->tol,
                        &(abscoeff->unres_lines))) {
                return 1;
            }
            break;
        case K_TYPE_HBr:
        case K_TYPE_H_79Br:
        case K_TYPE_H_81Br:
            if (Qratio( T,
                        hbr_Tref,
                        hbr_Qtab,
                        hbr_Qtab_rows,
                        hbr_Qtab_cols,
                        Qrat)) {
                errlog(147, 0);
                return 1;
            }
            if (linesum(k,
                        model->f,
                        model->f2,
                        model->df,
                        model->ngrid,
                        hbr_cat,
                        NULL,
                        hbr_num_lines,
                        layer->lineshape[abscoeff->k_typenum],
                        k_type[abscoeff->k_typenum].iso,
                        layer->P,
                        abscoeff->vmr_selfbroad,
                        T,
                        hbr_Tref,
                        Qrat,
                        hbr_abundance_tab,
                        hbr_mass_tab,
                        model->tol,
                        &(abscoeff->unres_lines))) {
                return 1;
            }
            break;
        case K_TYPE_HCN:
        case K_TYPE_H_12C_14N:
        case K_TYPE_H_13C_14N:
        case K_TYPE_H_12C_15N:
            if (Qratio( T,
                        hcn_Tref,
                        hcn_Qtab,
                        hcn_Qtab_rows,
                        hcn_Qtab_cols,
                        Qrat)) {
                errlog(151, 0);
                return 1;
            }
            if (linesum(k,
                        model->f,
                        model->f2,
                        model->df,
                        model->ngrid,
                        hcn_cat,
                        NULL,
                        hcn_num_lines,
                        layer->lineshape[abscoeff->k_typenum],
                        k_type[abscoeff->k_typenum].iso,
                        layer->P,
                        abscoeff->vmr_selfbroad,
                        T,
                        hcn_Tref,
                        Qrat,
                        hcn_abundance_tab,
                        hcn_mass_tab,
                        model->tol,
                        &(abscoeff->unres_lines))) {
                return 1;
            }
            break;
        case K_TYPE_H2CO:
        case K_TYPE_H2_12C_16O:
        case K_TYPE_H2_13C_16O:
        case K_TYPE_H2_12C_18O:
            if (Qratio( T,
                        h2co_Tref,
                        h2co_Qtab,
                        h2co_Qtab_rows,
                        h2co_Qtab_cols,
                        Qrat)) {
                errlog(155, 0);
                return 1;
            }
            if (linesum(k,
                        model->f,
                        model->f2,
                        model->df,
                        model->ngrid,
                        h2co_cat,
                        NULL,
                        h2co_num_lines,
                        layer->lineshape[abscoeff->k_typenum],
                        k_type[abscoeff->k_typenum].iso,
                        layer->P,
                        abscoeff->vmr_selfbroad,
                        T,
                        h2co_Tref,
                        Qrat,
                        h2co_abundance_tab,
                        h2co_mass_tab,
                        model->tol,
                        &(abscoeff->unres_lines))) {
                return 1;
            }
            break;
        case K_TYPE_HCl:
        case K_TYPE_H_35Cl:
        case K_TYPE_H_37Cl:
            if (Qratio( T,
                        hcl_Tref,
                        hcl_Qtab,
                        hcl_Qtab_rows,
                        hcl_Qtab_cols,
                        Qrat)) {
                errlog(148, 0);
                return 1;
            }
            if (linesum(k,
                        model->f,
                        model->f2,
                        model->df,
                        model->ngrid,
                        hcl_cat,
                        NULL,
                        hcl_num_lines,
                        layer->lineshape[abscoeff->k_typenum],
                        k_type[abscoeff->k_typenum].iso,
                        layer->P,
                        abscoeff->vmr_selfbroad,
                        T,
                        hcl_Tref,
                        Qrat,
                        hcl_abundance_tab,
                        hcl_mass_tab,
                        model->tol,
                        &(abscoeff->unres_lines))) {
                return 1;
            }
            break;
        case K_TYPE_HF:
        case K_TYPE_H_19F:
            if (Qratio( T,
                        hf_Tref,
                        hf_Qtab,
                        hf_Qtab_rows,
                        hf_Qtab_cols,
                        Qrat)) {
                errlog(149, 0);
                return 1;
            }
            if (linesum(k,
                        model->f,
                        model->f2,
                        model->df,
                        model->ngrid,
                        hf_cat,
                        NULL,
                        hf_num_lines,
                        layer->lineshape[abscoeff->k_typenum],
                        k_type[abscoeff->k_typenum].iso,
                        layer->P,
                        abscoeff->vmr_selfbroad,
                        T,
                        hf_Tref,
                        Qrat,
                        hf_abundance_tab,
                        hf_mass_tab,
                        model->tol,
                        &(abscoeff->unres_lines))) {
                return 1;
            }
            break;
        case K_TYPE_HNO3:
        case K_TYPE_H_14N_16O3:
            if (Qratio( T,
                        hno3_Tref,
                        hno3_Qtab,
                        hno3_Qtab_rows,
                        hno3_Qtab_cols,
                        Qrat)) {
                errlog(144, 0);
                return 1;
            }
            if (linesum(k,
                        model->f,
                        model->f2,
                        model->df,
                        model->ngrid,
                        hno3_cat,
                        NULL,
                        hno3_num_lines,
                        layer->lineshape[abscoeff->k_typenum],
                        k_type[abscoeff->k_typenum].iso,
                        layer->P,
                        abscoeff->vmr_selfbroad,
                        T,
                        hno3_Tref,
                        Qrat,
                        hno3_abundance_tab,
                        hno3_mass_tab,
                        model->tol,
                        &(abscoeff->unres_lines))) {
                return 1;
            }
            break;
        case K_TYPE_H2O_LINES:
        case K_TYPE_H2_16O:
        case K_TYPE_H2_18O:
        case K_TYPE_H2_17O:
        case K_TYPE_HD_16O:
        case K_TYPE_HD_18O:
        case K_TYPE_HD_17O:
            if (Qratio( T,
                        h2o_Tref,
                        h2o_Qtab,
                        h2o_Qtab_rows,
                        h2o_Qtab_cols,
                        Qrat)) {
                errlog(4, 0);
                return 1;
            }
            if (linesum(k,
                        model->f,
                        model->f2,
                        model->df,
                        model->ngrid,
                        h2o_cat,
                        NULL,
                        h2o_num_lines,
                        layer->lineshape[abscoeff->k_typenum],
                        k_type[abscoeff->k_typenum].iso,
                        layer->P,
                        abscoeff->vmr_selfbroad,
                        T,
                        h2o_Tref,
                        Qrat,
                        h2o_abundance_tab,
                        h2o_mass_tab,
                        model->tol,
                        &(abscoeff->unres_lines))) {
                return 1;
            }
            break;
        case K_TYPE_H2O_AIR_CONTINUUM:
            if (H2O_air_continuum(k, model->f, model->ngrid, T))
                return 1;
            break;
        case K_TYPE_H2O_SELF_CONTINUUM:
            if (H2O_self_continuum(k, model->f, model->ngrid, T))
                return 1;
            break;
        case K_TYPE_H2O2:
        case K_TYPE_H2_16O2:
            if (Qratio( T,
                        h2o2_Tref,
                        h2o2_Qtab,
                        h2o2_Qtab_rows,
                        h2o2_Qtab_cols,
                        Qrat)) {
                errlog(134, 0);
                return 1;
            }
            if (linesum(k,
                        model->f,
                        model->f2,
                        model->df,
                        model->ngrid,
                        h2o2_cat,
                        NULL,
                        h2o2_num_lines,
                        layer->lineshape[abscoeff->k_typenum],
                        k_type[abscoeff->k_typenum].iso,
                        layer->P,
                        abscoeff->vmr_selfbroad,
                        T,
                        h2o2_Tref,
                        Qrat,
                        h2o2_abundance_tab,
                        h2o2_mass_tab,
                        model->tol,
                        &(abscoeff->unres_lines))) {
                return 1;
            }
            break;
        case K_TYPE_HO2:
        case K_TYPE_H_16O2:
            if (Qratio( T,
                        ho2_Tref,
                        ho2_Qtab,
                        ho2_Qtab_rows,
                        ho2_Qtab_cols,
                        Qrat)) {
                errlog(135, 0);
                return 1;
            }
            if (linesum(k,
                        model->f,
                        model->f2,
                        model->df,
                        model->ngrid,
                        ho2_cat,
                        NULL,
                        ho2_num_lines,
                        layer->lineshape[abscoeff->k_typenum],
                        k_type[abscoeff->k_typenum].iso,
                        layer->P,
                        abscoeff->vmr_selfbroad,
                        T,
                        ho2_Tref,
                        Qrat,
                        ho2_abundance_tab,
                        ho2_mass_tab,
                        model->tol,
                        &(abscoeff->unres_lines))) {
                return 1;
            }
            break;
        case K_TYPE_HOCl:
        case K_TYPE_H_16O_35Cl:
        case K_TYPE_H_16O_37Cl:
            if (Qratio( T,
                        hocl_Tref,
                        hocl_Qtab,
                        hocl_Qtab_rows,
                        hocl_Qtab_cols,
                        Qrat)) {
                errlog(152, 0);
                return 1;
            }
            if (linesum(k,
                        model->f,
                        model->f2,
                        model->df,
                        model->ngrid,
                        hocl_cat,
                        NULL,
                        hocl_num_lines,
                        layer->lineshape[abscoeff->k_typenum],
                        k_type[abscoeff->k_typenum].iso,
                        layer->P,
                        abscoeff->vmr_selfbroad,
                        T,
                        hocl_Tref,
                        Qrat,
                        hocl_abundance_tab,
                        hocl_mass_tab,
                        model->tol,
                        &(abscoeff->unres_lines))) {
                return 1;
            }
            break;
        case K_TYPE_H2S:
        case K_TYPE_H2_32S:
        case K_TYPE_H2_34S:
        case K_TYPE_H2_33S:
            if (Qratio( T,
                        h2s_Tref,
                        h2s_Qtab,
                        h2s_Qtab_rows,
                        h2s_Qtab_cols,
                        Qrat)) {
                errlog(153, 0);
                return 1;
            }
            if (linesum(k,
                        model->f,
                        model->f2,
                        model->df,
                        model->ngrid,
                        h2s_cat,
                        NULL,
                        h2s_num_lines,
                        layer->lineshape[abscoeff->k_typenum],
                        k_type[abscoeff->k_typenum].iso,
                        layer->P,
                        abscoeff->vmr_selfbroad,
                        T,
                        h2s_Tref,
                        Qrat,
                        h2s_abundance_tab,
                        h2s_mass_tab,
                        model->tol,
                        &(abscoeff->unres_lines))) {
                return 1;
            }
            break;
        case K_TYPE_IWP_ABS_RAYLEIGH:
            if (iwp_abs_Rayleigh(k, model->f, model->ngrid, T))
                return 1;
            break;
        case K_TYPE_LWP_ABS_RAYLEIGH:
            if (lwp_abs_Rayleigh(k, model->f, model->ngrid, T))
                return 1;
            break;
        case K_TYPE_NH3:
        case K_TYPE_14NH3:
        case K_TYPE_15NH3:
            if (Qratio( T,
                        nh3_Tref,
                        nh3_Qtab,
                        nh3_Qtab_rows,
                        nh3_Qtab_cols,
                        Qrat)) {
                errlog(156, 0);
                return 1;
            }
            if (linesum(k,
                        model->f,
                        model->f2,
                        model->df,
                        model->ngrid,
                        nh3_cat,
                        NULL,
                        nh3_num_lines,
                        layer->lineshape[abscoeff->k_typenum],
                        k_type[abscoeff->k_typenum].iso,
                        layer->P,
                        abscoeff->vmr_selfbroad,
                        T,
                        nh3_Tref,
                        Qrat,
                        nh3_abundance_tab,
                        nh3_mass_tab,
                        model->tol,
                        &(abscoeff->unres_lines))) {
                return 1;
            }
            break;
        case K_TYPE_N2N2:
            if (N2N2_cia(k, model->f, model->ngrid, T))
                return 1;
            break;
        case K_TYPE_N2O:
        case K_TYPE_14N2_16O:
        case K_TYPE_14N_15N_16O:
        case K_TYPE_15N_14N_16O:
        case K_TYPE_14N2_18O:
        case K_TYPE_14N2_17O:
            if (Qratio( T,
                        n2o_Tref,
                        n2o_Qtab,
                        n2o_Qtab_rows,
                        n2o_Qtab_cols,
                        Qrat)) {
                errlog(10, 0);
                return 1;
            }
            if (linesum(k,
                        model->f,
                        model->f2,
                        model->df,
                        model->ngrid,
                        n2o_cat,
                        NULL,
                        n2o_num_lines,
                        layer->lineshape[abscoeff->k_typenum],
                        k_type[abscoeff->k_typenum].iso,
                        layer->P,
                        abscoeff->vmr_selfbroad,
                        T,
                        n2o_Tref,
                        Qrat,
                        n2o_abundance_tab,
                        n2o_mass_tab,
                        model->tol,
                        &(abscoeff->unres_lines))) {
                return 1;
            }
            break;
        case K_TYPE_NO:
        case K_TYPE_14N_16O:
            if (Qratio( T,
                        no_Tref,
                        no_Qtab,
                        no_Qtab_rows,
                        no_Qtab_cols,
                        Qrat)) {
                errlog(145, 0);
                return 1;
            }
            if (linesum(k,
                        model->f,
                        model->f2,
                        model->df,
                        model->ngrid,
                        no_cat,
                        NULL,
                        no_num_lines,
                        layer->lineshape[abscoeff->k_typenum],
                        k_type[abscoeff->k_typenum].iso,
                        layer->P,
                        abscoeff->vmr_selfbroad,
                        T,
                        no_Tref,
                        Qrat,
                        no_abundance_tab,
                        no_mass_tab,
                        model->tol,
                        &(abscoeff->unres_lines))) {
                return 1;
            }
            break;
        case K_TYPE_NO2:
        case K_TYPE_14N_16O2:
            if (Qratio( T,
                        no2_Tref,
                        no2_Qtab,
                        no2_Qtab_rows,
                        no2_Qtab_cols,
                        Qrat)) {
                errlog(146, 0);
                return 1;
            }
            if (linesum(k,
                        model->f,
                        model->f2,
                        model->df,
                        model->ngrid,
                        no2_cat,
                        NULL,
                        no2_num_lines,
                        layer->lineshape[abscoeff->k_typenum],
                        k_type[abscoeff->k_typenum].iso,
                        layer->P,
                        abscoeff->vmr_selfbroad,
                        T,
                        no2_Tref,
                        Qrat,
                        no2_abundance_tab,
                        no2_mass_tab,
                        model->tol,
                        &(abscoeff->unres_lines))) {
                return 1;
            }
            break;
        case K_TYPE_O:
        case K_TYPE_16O:
            if (Qratio( T,
                        o_Tref,
                        o_Qtab,
                        o_Qtab_rows,
                        o_Qtab_cols,
                        Qrat)) {
                errlog(137, 0);
                return 1;
            }
            if (linesum(k,
                        model->f,
                        model->f2,
                        model->df,
                        model->ngrid,
                        o_cat,
                        NULL,
                        o_num_lines,
                        layer->lineshape[abscoeff->k_typenum],
                        k_type[abscoeff->k_typenum].iso,
                        layer->P,
                        abscoeff->vmr_selfbroad,
                        T,
                        o_Tref,
                        Qrat,
                        o_abundance_tab,
                        o_mass_tab,
                        model->tol,
                        &(abscoeff->unres_lines))) {
                return 1;
            }
            break;
        case K_TYPE_O2_COUPLED:
        case K_TYPE_16O2_COUPLED:
            if (Qratio( T,
                        o2_Tref,
                        o2_Qtab,
                        o2_Qtab_rows,
                        o2_Qtab_cols,
                        Qrat)) {
                errlog(6, 0);
                return 1;
            }
            if (linesum(k,
                        model->f,
                        model->f2,
                        model->df,
                        model->ngrid,
                        o2_coupled_cat,
                        o2_line_coupling_coeffs,
                        o2_num_coupled_lines,
                        layer->lineshape[abscoeff->k_typenum],
                        k_type[abscoeff->k_typenum].iso,
                        layer->P,
                        abscoeff->vmr_selfbroad,
                        T,
                        o2_Tref,
                        Qrat,
                        o2_abundance_tab,
                        o2_mass_tab,
                        model->tol,
                        &(abscoeff->unres_lines))) {
                return 1;
            }
            break;
        case K_TYPE_O2_UNCOUPLED:
        case K_TYPE_16O2_UNCOUPLED:
        case K_TYPE_16O_18O_UNCOUPLED:
        case K_TYPE_16O_17O_UNCOUPLED:
            if (Qratio( T,
                        o2_Tref,
                        o2_Qtab,
                        o2_Qtab_rows,
                        o2_Qtab_cols,
                        Qrat)) {
                errlog(6, 0);
                return 1;
            }
            if (linesum(k,
                        model->f,
                        model->f2,
                        model->df,
                        model->ngrid,
                        o2_uncoupled_cat,
                        NULL,
                        o2_num_uncoupled_lines,
                        layer->lineshape[abscoeff->k_typenum],
                        k_type[abscoeff->k_typenum].iso,
                        layer->P,
                        abscoeff->vmr_selfbroad,
                        T,
                        o2_Tref,
                        Qrat,
                        o2_abundance_tab,
                        o2_mass_tab,
                        model->tol,
                        &(abscoeff->unres_lines))) {
                return 1;
            }
            break;
        case K_TYPE_O2O2:
            if (O2O2_cia(k, model->f, model->ngrid, T))
                return 1;
            break;
        case K_TYPE_O3:
        case K_TYPE_16O3:
        case K_TYPE_16O2_18O:
        case K_TYPE_16O_18O_16O:
        case K_TYPE_16O2_17O:
        case K_TYPE_16O_17O_16O:
            if (Qratio( T,
                        o3_Tref,
                        o3_Qtab,
                        o3_Qtab_rows,
                        o3_Qtab_cols,
                        Qrat)) {
                errlog(8, 0);
                return 1;
            }
            if (linesum(k,
                        model->f,
                        model->f2,
                        model->df,
                        model->ngrid,
                        o3_cat,
                        NULL,
                        o3_num_lines,
                        layer->lineshape[abscoeff->k_typenum],
                        k_type[abscoeff->k_typenum].iso,
                        layer->P,
                        abscoeff->vmr_selfbroad,
                        T,
                        o3_Tref,
                        Qrat,
                        o3_abundance_tab,
                        o3_mass_tab,
                        model->tol,
                        &(abscoeff->unres_lines))) {
                return 1;
            }
            break;
        case K_TYPE_OCS:
        case K_TYPE_16O_12C_32S:
        case K_TYPE_16O_12C_34S:
        case K_TYPE_16O_13C_32S:
        case K_TYPE_16O_12C_33S:
        case K_TYPE_18O_12C_32S:
            if (Qratio( T,
                        ocs_Tref,
                        ocs_Qtab,
                        ocs_Qtab_rows,
                        ocs_Qtab_cols,
                        Qrat)) {
                errlog(104, 0);
                return 1;
            }
            if (linesum(k,
                        model->f,
                        model->f2,
                        model->df,
                        model->ngrid,
                        ocs_cat,
                        NULL,
                        ocs_num_lines,
                        layer->lineshape[abscoeff->k_typenum],
                        k_type[abscoeff->k_typenum].iso,
                        layer->P,
                        abscoeff->vmr_selfbroad,
                        T,
                        ocs_Tref,
                        Qrat,
                        ocs_abundance_tab,
                        ocs_mass_tab,
                        model->tol,
                        &(abscoeff->unres_lines))) {
                return 1;
            }
            break;
        case K_TYPE_OH:
        case K_TYPE_16O_H:
        case K_TYPE_18O_H:
        case K_TYPE_16O_D:
            if (Qratio( T,
                        oh_Tref,
                        oh_Qtab,
                        oh_Qtab_rows,
                        oh_Qtab_cols,
                        Qrat)) {
                errlog(136, 0);
                return 1;
            }
            if (linesum(k,
                        model->f,
                        model->f2,
                        model->df,
                        model->ngrid,
                        oh_cat,
                        NULL,
                        oh_num_lines,
                        layer->lineshape[abscoeff->k_typenum],
                        k_type[abscoeff->k_typenum].iso,
                        layer->P,
                        abscoeff->vmr_selfbroad,
                        T,
                        oh_Tref,
                        Qrat,
                        oh_abundance_tab,
                        oh_mass_tab,
                        model->tol,
                        &(abscoeff->unres_lines))) {
                return 1;
            }
            break;
        case K_TYPE_SO2:
        case K_TYPE_32S_16O2:
        case K_TYPE_34S_16O2:
            if (Qratio( T,
                        so2_Tref,
                        so2_Qtab,
                        so2_Qtab_rows,
                        so2_Qtab_cols,
                        Qrat)) {
                errlog(138, 0);
                return 1;
            }
            if (linesum(k,
                        model->f,
                        model->f2,
                        model->df,
                        model->ngrid,
                        so2_cat,
                        NULL,
                        so2_num_lines,
                        layer->lineshape[abscoeff->k_typenum],
                        k_type[abscoeff->k_typenum].iso,
                        layer->P,
                        abscoeff->vmr_selfbroad,
                        T,
                        so2_Tref,
                        Qrat,
                        so2_abundance_tab,
                        so2_mass_tab,
                        model->tol,
                        &(abscoeff->unres_lines))) {
                return 1;
            }
            break;
        case K_TYPE_ONELINE:
            if (Qratio( T,
                        oneline_Tref,
                        oneline_Qtab,
                        oneline_Qtab_rows,
                        oneline_Qtab_cols,
                        Qrat)) {
                errlog(11, 0);
                return 1;
            }
            if (linesum(k,
                        model->f,
                        model->f2,
                        model->df,
                        model->ngrid,
                        oneline_cat,
                        NULL,
                        oneline_num_lines,
                        layer->lineshape[abscoeff->k_typenum],
                        k_type[abscoeff->k_typenum].iso,
                        layer->P,
                        abscoeff->vmr_selfbroad,
                        T,
                        oneline_Tref,
                        Qrat,
                        oneline_abundance_tab,
                        oneline_mass_tab,
                        model->tol,
                        &(abscoeff->unres_lines))) {
                return 1;
            }
            break;
        default:
            errlog(24, abscoeff->k_typenum);
            return 1;
    }
    return 0;
}   /* compute_absorption_coefficient() */
