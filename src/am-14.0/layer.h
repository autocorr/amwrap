/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* layer.h                         S. Paine rev. 2024 July 16
*
* Constants and declarations for layer.c
************************************************************/

#ifndef AM_LAYER_H
#define AM_LAYER_H

#include "am_types.h"

/*
 * Layer types
 */
enum {
    LAYER_TYPE_NONE,
    LAYER_TYPE_DEF,    /* model definition layer                   */
    LAYER_TYPE_OBS,    /* interpolated layer above observing level */
    LAYER_TYPE_SOURCE, /* interpolated layer above source level    */
    LAYER_TYPE_TAN,    /* interpolated layer above tangent level   */
    LAYER_TYPE_END
};

/*
 * The following constants correspond to bits in layer_t.vmr_stat and
 * column_t.vmr_stat.  Some bits track error conditions associated
 * with volume mixing ratio computations.  If one of these error bits
 * is set during computation of a mixing ratio, and that mixing ratio
 * turns out to be needed in a subsequent computation, then the master
 * VMR_ERROR bit will also be set.
 */
enum {
    VMR_ERROR                        = 0x1,
    VMR_ZERO_DEFAULT                 = 0x2,
    VMR_SUM_EXCEEDS_UNITY            = 0x4,
    VMR_DEFAULT_ON_SCALED_COLUMN     = 0x8,
    VMR_USER_DEFINED                 = 0x10,
    VMR_BY_RH_LIQUID                 = 0x20,
    VMR_BY_RH_ICE                    = 0x40,
    VMR_WEIGHTED_MEAN_MASS_UNDEFINED = 0x80
};

enum {
    VMR_BY_RH = (VMR_BY_RH_LIQUID | VMR_BY_RH_ICE)
};

enum {
    VMR_ALL_ERRS = (
        VMR_ERROR |
        VMR_ZERO_DEFAULT |
        VMR_SUM_EXCEEDS_UNITY |
        VMR_DEFAULT_ON_SCALED_COLUMN |
        VMR_WEIGHTED_MEAN_MASS_UNDEFINED
        )
};

/*
 * VMR_SUM_WARNING_TOL is the threshold for generating a warning
 * that the sum of mixing ratios on a layer exceeds unity.
 */
#define VMR_SUM_WARNING_TOL 3.0e-3

void   compute_adjusted_selfbroad_vmrs(model_t*);
int    get_lnum_by_type(const model_t*, const int);
int    get_obs_lnum(model_t*);
int    get_source_lnum(model_t*);
int    get_tan_lnum(model_t*);
int    insert_interpolated_level(
        model_t*,
        model_t*,
        const int,
        const double,
        const double,
        const double);
double layer_base_lapse_rate(const model_t*, const int);
void   reset_layer_performance_timers(model_t*);
void   restore_def_layer_explicit_col_densities(model_t*);
void   set_def_layer_heights_and_col_densities(model_t*);
void   set_def_layer_temperatures(model_t*);
void   set_def_layer_pressures(model_t*);
void   set_def_layer_refractive_indices(model_t*);
void   set_def_layer_vmrs_and_apply_Nscale_factors(model_t*);
int    set_interpolated_layer_variables(model_t*);

#endif /* AM_LAYER_H */
