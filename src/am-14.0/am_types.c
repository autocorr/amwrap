/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* am_types.c                      S. Paine rev. 2024 July 22
*
* Initializers for certain defined types.
************************************************************/

#include <stdio.h>

#include "abscoeff.h"
#include "am_types.h"
#include "column.h"
#include "fit.h"
#include "ils.h"
#include "layer.h"
#include "model.h"
#include "molecules.h"
#include "output.h"
#include "phys_const.h"
#include "units.h"


/*
 * Initializers for dimensioned quantities are in am internal units.
 * These may differ from the initializers for the corresponding user
 * units, which are used for output.
 */
const model_t MODEL_INIT = {
    0.0,                            /* fmin                 */
    0.0,                            /* fmax                 */
    0.0,                            /* df                   */
    -1.0,                           /* fout_min             */
    -1.0,                           /* fout_max             */
    0.0,                            /* flo                  */
    0.0,                            /* fif_0                */
    0.0,                            /* fif_delta            */
    -1.0,                           /* fif_min              */
    -1.0,                           /* fif_max              */
    0.0,                            /* ils_fwhm             */
    0.0,                            /* ils_fif              */
    1.0,                            /* dsb_utol_ratio       */
    1.0,                            /* rx_gain_factor       */
    1.0,                            /* RH_offset_exp        */
    1.0,                            /* RH_scale             */
    G_STD,                          /* g                    */
    0.0,                            /* dg_dz                */
    R_EARTH,                        /* R0                   */
    0.0,                            /* z0                   */
    0.0,                            /* p                    */
    0.0,                            /* Ptoa                 */
    0.0,                            /* Pmax                 */
    0.0,                            /* Pobs                 */
    0.0,                            /* zobs                 */
    0.0,                            /* Psource              */
    0.0,                            /* zsource              */
    0.0,                            /* Ptan                 */
    0.0,                            /* ztan                 */
    0.0,                            /* tol                  */
    3.0e-3,                         /* selfbroad_vmr_tol    */
    0.0,                            /* T0                   */
    0.0,                            /* Trx                  */
    0.0,                            /* Tref                 */
    190.0,                          /* kcache_Tmin          */
    320.0,                          /* kcache_Tmax          */
    1.0,                            /* kcache_dT            */
    1.0,                            /* sec_za               */
    0.0,                            /* za                   */
    0.0,                            /* am_runtime           */
    0.0,                            /* od_runtime           */
    0.0,                            /* rt_runtime           */
    0.0,                            /* spec_runtime         */
    0.0,                            /* runtime              */
    NULL,                           /* f                    */
    NULL,                           /* f2                   */
    NULL,                           /* fif                  */
    NULL,                           /* ils                  */
    NULL,                           /* ilsworkspace         */
    NULL,                           /* tau                  */
    NULL,                           /* tx                   */
    NULL,                           /* I0                   */
    NULL,                           /* I                    */
    NULL,                           /* I_ref                */
    NULL,                           /* I_diff               */
    NULL,                           /* Tb                   */
    NULL,                           /* Trj                  */
    NULL,                           /* Tsys                 */
    NULL,                           /* Y                    */
    NULL,                           /* L                    */
    NULL,                           /* tau_fsl              */
    NULL,                           /* k_out                */
    NULL,                           /* layer                */
    0,                              /* ngrid                */
    0,                              /* nif                  */
    0,                              /* npad                 */
    0,                              /* nLpad                */
    {0, 0},                         /* ilsb                 */
    {0, 0},                         /* iusb                 */
    {{0, 0}, {0, 0}},               /* isub[]               */
    0,                              /* nlayers              */
    0,                              /* path_begin           */
    -1,                             /* path_mid             */
    0,                              /* path_end             */
    0,                              /* path_min             */
    131,                            /* nkcache              */
    AM_UNIT_FREQUENCY,              /* fmin_unitnum         */
    AM_UNIT_FREQUENCY,              /* fmax_unitnum         */
    AM_UNIT_FREQUENCY,              /* df_unitnum           */
    AM_UNIT_FREQUENCY,              /* fout_min_unitnum     */
    AM_UNIT_FREQUENCY,              /* fout_max_unitnum     */
    AM_UNIT_FREQUENCY,              /* flo_unitnum          */
    AM_UNIT_FREQUENCY,              /* fif_min_unitnum      */
    AM_UNIT_FREQUENCY,              /* fif_max_unitnum      */
    AM_UNIT_ACCEL,                  /* g_unitnum            */
    AM_UNIT_ACCEL_GRADIENT,         /* dg_dz_unitnum        */
    UNIT_KM,                        /* R0_unitnum           */
    UNIT_M,                         /* z0_unitnum           */
    AM_UNIT_PRESSURE,               /* Pobs_unitnum         */
    UNIT_M,                         /* zobs_unitnum         */
    AM_UNIT_PRESSURE,               /* Psource_unitnum      */
    UNIT_M,                         /* zsource_unitnum      */
    AM_UNIT_PRESSURE,               /* Ptan_unitnum         */
    UNIT_M,                         /* ztan_unitnum         */
    GEOMETRY_DEFAULTS,              /* geometry             */
    0,                              /* ifmode               */
    AM_UNIT_FREQUENCY,              /* ils_fwhm_unitnum     */
    AM_UNIT_FREQUENCY,              /* ils_fif_unitnum      */
    ILSMODE_NORMAL,                 /* ilsmode              */
    ILS_NONE,                       /* ils_typenum          */
    0,                              /* RH_adj_flag          */
    AM_UNIT_TEMPERATURE,            /* T0_unitnum           */
    AM_UNIT_TEMPERATURE,            /* Tref_unitnum         */
    AM_UNIT_TEMPERATURE,            /* Trx_unitnum          */
    PTMODE_DEFAULTS,                /* PTmode               */
    UNIT_DEG,                       /* za_unitnum           */
    0                               /* log_runtimes         */
};

const layer_t LAYER_INIT = {
    0.0,                    /* P                        */
    0.0,                    /* T                        */
    0.0,                    /* dP                       */
    0.0,                    /* dP_def                   */
    0.0,                    /* Pbase                    */
    0.0,                    /* dphi                     */
    0.0,                    /* Tbase                    */
    -1.0,                   /* h (-1.0 means undefined) */
    0.0,                    /* m (airmass)              */
    MASS_DRY_AIR,           /* Mair                     */
    MASS_DRY_AIR,           /* M0                       */
    G_STD,                  /* gbase                    */
    1.0,                    /* nbase                    */
    0.0,                    /* zbase                    */
    0.0,                    /* za_base                  */
    NULL,                   /* B                        */
    NULL,                   /* tau                      */
    NULL,                   /* tx                       */
    NULL,                   /* column                   */
    NULL,                   /* lineshape                */
    NULL,                   /* strict_selfbroad         */
    NULL,                   /* Mair_flag                */
    0.0,                    /* runtime                  */
    0,                      /* tagnum                   */
    LAYER_TYPE_NONE,        /* type                     */
    0,                      /* ncols                    */
    UNIT_NONE,              /* P_unitnum                */
    UNIT_NONE,              /* T_unitnum                */
    UNIT_NONE,              /* h_unitnum                */
    0,                      /* updateflag               */
    0                       /* vmr_stat                 */
};

const column_t COLUMN_INIT = {
    0.0,                  /* N           */
    0.0,                  /* N_def       */
    0.0,                  /* N_scaled    */
    0.0,                  /* xvmr        */
    0.0,                  /* vmr_scaled  */
    0.0,                  /* RH          */
    0.0,                  /* RH_adj      */
    NULL,                 /* ztau        */
    NULL,                 /* abscoeff    */
    0.0,                  /* runtime     */
    0,                    /* col_typenum */
    0,                    /* n_abscoeffs */
    0,                    /* N_mode      */
    AM_UNIT_COL_DENSITY,  /* N_unitnum   */
    0,                    /* vmr_stat    */
    0                     /* unres_lines */
};

const abscoeff_t ABSCOEFF_INIT = {
    0.0,          /* vmr_selfbroad */
    NULL,         /* k             */
    NULL,         /* kcache        */
    K_TYPE_NONE,  /* k_type        */
    0             /* unres_lines   */
};

const simplex_t SIMPLEX_INIT = {
    0.0,    /* E1                                */
    0.0,    /* E2                                */
    1e-4,   /* tol                               */
    NULL,   /* name                              */
    NULL,   /* init                              */
    NULL,   /* scale                             */
    NULL,   /* unitnum                           */
    NULL,   /* mapping                           */
    NULL,   /* varptr                            */
    NULL,   /* vertex                            */
    NULL,   /* E                                 */
    NULL,   /* pbar                              */
    NULL,   /* p1                                */
    NULL,   /* p2                                */
    NULL,   /* lowest vertex at last convergence */
    0,      /* ihigh                             */
    0,      /* ilow                              */
    0,      /* n                                 */
    0,      /* iter                              */
    0,      /* reinit                            */
    0,      /* restart_count                     */
    1       /* logarithmic                       */
};

const fit_data_t FIT_DATA_INIT = {
    NULL,                               /* filename              */
    NULL,                               /* fp                    */
    NULL,                               /* f (frequency)         */
    NULL,                               /* s (spectral data)     */
    NULL,                               /* s_mod                 */
    NULL,                               /* res                   */
    NULL,                               /* res_est               */
    NULL,                               /* b (channel bandwidth) */
    NULL,                               /* w (weight factor)     */
    0.0,                                /* mean_var              */
    0.0,                                /* mean_var_tracked      */
    -1.0,                               /* res_track_gain        */
    0.0,                                /* runtime               */
    0,                                  /* npts                  */
    0,                                  /* nalloc                */
    OUTPUT_NONE,                        /* data_type             */
    FIT_ESTIMATOR_ABSOLUTE_RESIDUALS,   /* estimator_type        */
    AM_UNIT_FREQUENCY,                  /* f_unitnum             */
    UNIT_NONE,                          /* s_unitnum             */
    1,                                  /* f_col                 */
    2,                                  /* s_col                 */
    0,                                  /* b_col                 */
    0,                                  /* w_col                 */
    1000,                               /* max_iter              */
    5,                                  /* max_restarts          */
    0,                                  /* nfiles                */
    0,                                  /* open_filenum          */
    0,                                  /* blocknum              */
    FIT_OUTPUT_CONFIG |
        FIT_OUTPUT_SPECTRUM |
        FIT_OUTPUT_RESIDUAL             /* output_mode           */
};

