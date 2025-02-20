/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* am_alloc.c                        S. Paine rev. 2023 May 5
*
* Memory allocation functions.
************************************************************/

#include <float.h>
#include <math.h>
#include <stdlib.h>

#include "abscoeff.h"
#include "am_alloc.h"
#include "am_types.h"
#include "column.h"
#include "errlog.h"
#include "kcache.h"
#include "layer.h"
#include "model.h"
#include "molecules.h"
#include "nscale.h"
#include "output.h"
#include "simplex.h"
#include "units.h"

/*
 * Constants controlling grow_fit_data_arrays().
 */
enum {
    FIT_DATA_NALLOC_INIT = 512,
    FIT_DATA_NALLOC_GROW = 2
};

static void free_abscoeffs(column_t*, int);
static void free_columns(layer_t*, int);
static void free_layer(model_t*, int);
static void free_model_layers(model_t*);


/***********************************************************
* int add_column(layer_t *layer, int col_typenum)
*
* Purpose:
*   Adds a new column structure to a layer, and initializes
*   it.
*
* Arguments:
*   layer_t *layer  - pointer to model structure
*   int col_typenum - type number of column
*
* Return:
*   0 if OK
*   1 on error
************************************************************/

int add_column(layer_t *layer, int col_typenum)
{
    int i, n, n_abs;
    column_t **tptr, *column;

    n = layer->ncols;
    /*
     * Enlarge the table of pointers to column structures
     */
    if ((tptr = (column_t**)realloc(
                    layer->column, (n + 1) * sizeof(column_t*))) == NULL) {
        errlog(45, 0);
        return 1;
    }
    /*
     * realloc() succeeded, OK to update pointer.
     */
    layer->column = tptr;
    ++layer->ncols;
    /*
     * Allocate memory for the new column structure and
     * initialize.
     */
    if ((layer->column[n] = (column_t*)malloc(sizeof(column_t))) == NULL) {
        errlog(15, 0);
        return 1;
    }
    column  = layer->column[n];
    *column = COLUMN_INIT;
    column->col_typenum = col_typenum;
    n_abs = col_type[col_typenum].n_abscoeffs;
    column->n_abscoeffs = n_abs;
    /*
     * Allocate memory for the table of absorption coefficient
     * structures for this column.
     */
    if ((column->abscoeff =
                (abscoeff_t**)malloc(n_abs * sizeof(abscoeff_t*))) == NULL) {
        errlog(99, 0);
        return 1;
    }
    /*
     * Allocate and initialize the absorption coefficient
     * structures.
     */
    for (i = 0; i < n_abs; ++i) {
        if ((column->abscoeff[i] =
                    (abscoeff_t*)malloc(sizeof(abscoeff_t))) == NULL) {
            errlog(100, 0);
            return 1;
        }
        *(column->abscoeff[i]) = ABSCOEFF_INIT;
        column->abscoeff[i]->k_typenum = col_type[col_typenum].k_type[i];
    }
    return 0;
}   /* add_column() */


/***********************************************************
* int add_layer(model_t *model, int type)
*
* Purpose:
*   Adds a new layer to an existing model structure.
*
* Arguments:
*   model_t *model - pointer to model structure
*   type           - layer type, defined in am_types.h.
*
* Return:
*   0 if OK
*   1 on error
************************************************************/

int add_layer(model_t *model, int type)
{
    int i, n;
    layer_t **tptr;

    if (type <= LAYER_TYPE_NONE || type >= LAYER_TYPE_END) {
        errlog(162, type);
        return 1;
    }
    n = model->nlayers;
    /*
     * Enlarge the table of pointers to layer structures.
     */
    if ((tptr = (layer_t**)realloc(model->layer, (n + 1) * sizeof(layer_t*)))
            == NULL) {
        errlog(44, 0);
        return 1;
    }
    /*
     * realloc() succeeded, OK to update pointer.
     */
    model->layer = tptr;
    ++model->nlayers;
    /*
     * Allocate space for the new layer structure and initialize.
     */
    if ((model->layer[n] = (layer_t*)malloc(sizeof(layer_t))) == NULL) {
        errlog(13, 0);
        return 1;
    }
    *(model->layer[n])    = LAYER_INIT;
    model->layer[n]->type = type;
    /*
     * Allocate space for lineshape type array, and set defaults.
     */
    if ((model->layer[n]->lineshape =
                (int*)malloc(K_TYPE_END_OF_TABLE * sizeof(int))) == NULL) {
        errlog(129, 0);
        return 1;
    }
    for (i = 0; i < K_TYPE_END_OF_TABLE; ++i)
        model->layer[n]->lineshape[i] = k_type[i].default_lineshape;
    /*
     * Allocate space for strict_selfbroad flag array, and set
     * defaults.
     */
    if ((model->layer[n]->strict_selfbroad =
                (int*)malloc(K_TYPE_END_OF_TABLE * sizeof(int))) == NULL) {
        errlog(130, 0);
        return 1;
    }
    for (i = 0; i < K_TYPE_END_OF_TABLE; ++i)
        model->layer[n]->strict_selfbroad[i] = 0;
    /*
     * Allocate space for Mair_flag array, and set defaults.
     */
    if ((model->layer[n]->Mair_flag =
                (int*)malloc(COL_TYPE_END_OF_TABLE * sizeof(int))) == NULL) {
        errlog(131,0);
        return 1;
    }
    for (i = 0; i < COL_TYPE_END_OF_TABLE; ++i)
        model->layer[n]->Mair_flag[i] = 0;
    model->layer[n]->Mair_flag[COL_TYPE_H2O] = 1;
    model->layer[n]->Mair_flag[COL_TYPE_O3] = 1;
    return 0;
}   /* add_layer() */


/***********************************************************
* int alloc_layer_arrays(model_t *model, int lnum)
*
* Purpose:
*   Allocates space for spectral arrays in layer lnum of
*   *model, and in all of the column and absorption
*   coefficient structures contained in the layer.
*
*   Array allocation is skipped when certain arrays will
*   never need to be computed.  Besides saving memory, the
*   NULL pointer is used elsewhere in the program as a flag
*   to skip the corresponding computation.  In columns with
*   a fixed column density equal to zero, the k[] and ztau[]
*   arrays are not allocated, unless this is a fit or a
*   Jacobian computation, or if a spectral absorption
*   coefficient is a requested output.  In layers with no
*   non-empty columns, tau[] and tx[] are not allocated.
*   B[] is allocated on an empty layer when layer
*   temperatures are specified by Tbase.
*
* Arguments:
*   model_t *model - pointer to model structure
*   int lnum       - layer number
*
* Return:
*   0 if successful, 1 otherwise
************************************************************/

int alloc_layer_arrays(model_t *model, int lnum)
{
    int cnum;
    int layer_not_empty = 0;
    layer_t *layer = model->layer[lnum];

    for (cnum = 0; cnum < layer->ncols; ++cnum) {
        column_t *column = layer->column[cnum];
        int col_typenum = column->col_typenum;
        int knum;
        /*
         * If this is not a fit or jacobian computation, and the
         * column density for this column will always be zero,
         * skip array allocation unless k has been requested as a
         * program output.
         */
        if (
                !(output[ALL_OUTPUTS].flags & OUTPUT_FITTED) &&
                !(output[ALL_OUTPUTS].flags & OUTPUT_JACOBIAN) &&
                !(output[OUTPUT_K].flags & OUTPUT_USER)
                ) {
            /*
             * Skip array allocation if Nscale is zero.
             */
            if (lookup_Nscale(col_typenum, layer->tagnum) < DBL_MIN)
                continue;
            /*
             * Skip array allocation if the column density is set
             * by mixing ratio, and the mixing ratio is fixed at
             * zero either directly or as zero relative humidity.
             */
            if (column->N_mode & N_BY_VMR) {
                if (column->vmr_stat & VMR_BY_RH) {
                    if (column->RH < DBL_MIN  &&
                            fabs(model->RH_offset_exp - 1.0) < DBL_MIN) {
                        continue;
                    }
                } else {
                    /*
                     * xvmr->0 as vmr->0, so just test xvmr.
                     */
                    if (column->xvmr < DBL_MIN) {
                        continue;
                    }
                }
            }
            /*
             * Skip array allocation if the column has been
             * defined with an explicit column density equal to
             * zero.
             */
            if ((column->N_mode & N_EXPLICIT) && (column->N_def < DBL_MIN))
                continue;
        }
        /*
         * Skip array allocation if arrays are not needed for
         * this column type, as for delay-only column types.
         */
        if (col_type[col_typenum].flags & COL_NO_ARRAYS)
            continue;
        /*
         * OK, allocate arrays
         */
        layer_not_empty = 1;
        for (knum = 0; knum < column->n_abscoeffs; ++knum) {
            if ((column->abscoeff[knum]->k =
                        (double*)malloc(model->ngrid * sizeof(double)))
                    == NULL) {
                errlog(16, 0);
                return 1;
            }
        }
        if ((column->ztau = (double*)malloc(model->ngrid * sizeof(double)))
                == NULL) {
            errlog(17, 0);
            return 1;
        }
    }
    if (layer_not_empty || (model->PTmode & PTMODE_TBASE)) {
        if ((layer->B = (double*)malloc(model->ngrid * sizeof(double)))
                == NULL) {
            errlog(93, lnum);
            return 1;
        }
    }
    if (layer_not_empty) {
        if ((layer->tau = (double*)malloc(model->ngrid * sizeof(double)))
                == NULL) {
            errlog(14, lnum);
            return 1;
        }
        if ((layer->tx = (double*)malloc(model->ngrid * sizeof(double)))
                == NULL) {
            errlog(94, lnum);
            return 1;
        }
    }
    return 0;
}   /* alloc_layer_arrays() */


/***********************************************************
* int alloc_model_arrays(model_t *model)
*
* Purpose:
*   Allocates space for the frequency grid and spectral
*   arrays in the top level model structure, and in all
*   layers, columns, and absorption coefficients.  The
*   addresses of those arrays which can be model outputs are
*   copied into the global output table.
*
* Arguments:
*   model_t *model - pointer to model structure
*
* Return:
*   0 if successful, 1 otherwise
************************************************************/

int alloc_model_arrays(model_t *model)
{
    int lnum;

    if (model->ngrid == 0) {
        errlog(67, 0);
        return 1;
    }
    if ((model->f =
                (double*)malloc(model->ngrid * sizeof(double))) == NULL) {
        errlog(21, 0);
        return 1;
    }
    if ((model->f2 =
                (double*)malloc(model->ngrid * sizeof(double))) == NULL) {
        errlog(21, 0);
        return 1;
    }
    if (model->ifmode) {
        if ((model->fif =
                    (double*)malloc(model->ngrid * sizeof(double))) == NULL) {
            errlog(58, 0);
            return 1;
        }
        output[OUTPUT_FREQUENCY].spectrum = model->fif;
    } else {
        output[OUTPUT_FREQUENCY].spectrum = model->f;
    }
    if ((model->tau =
                (double*)malloc(model->ngrid * sizeof(double))) == NULL) {
        errlog(63, 0);
        return 1;
    }
    output[OUTPUT_OPACITY].spectrum = model->tau;
    if ((model->tx =
                (double*)malloc(model->ngrid * sizeof(double))) == NULL) {
        errlog(28, 0);
        return 1;
    }
    output[OUTPUT_TRANSMITTANCE].spectrum = model->tx;
    if ((model->I0 =
                (double*)malloc(model->ngrid * sizeof(double))) == NULL) {
        errlog(92, 0);
        return 1;
    }
    if ((model->I =
                (double*)malloc(model->ngrid * sizeof(double))) == NULL) {
        errlog(29, 0);
        return 1;
    }
    output[OUTPUT_RADIANCE].spectrum = model->I;
    if ((model->I_ref =
                (double*)malloc(model->ngrid * sizeof(double))) == NULL) {
        errlog(119, 0);
        return 1;
    }
    if ((model->I_diff =
                (double*)malloc(model->ngrid * sizeof(double))) == NULL) {
        errlog(120, 0);
        return 1;
    }
    output[OUTPUT_RADIANCE_DIFF].spectrum = model->I_diff;
    if ((model->Tb =
                (double*)malloc(model->ngrid * sizeof(double))) == NULL) {
        errlog(64, 0);
        return 1;
    }
    output[OUTPUT_TB_PLANCK].spectrum = model->Tb;
    if ((model->Trj =
                (double*)malloc(model->ngrid * sizeof(double))) == NULL) {
        errlog(53, 0);
        return 1;
    }
    output[OUTPUT_TB_RAYLEIGH_JEANS].spectrum = model->Trj;
    if ((model->Tsys =
                (double*)malloc(model->ngrid * sizeof(double))) == NULL) {
        errlog(108, 0);
        return 1;
    }
    output[OUTPUT_TSYS].spectrum = model->Tsys;
    if ((model->Y =
                (double*)malloc(model->ngrid * sizeof(double))) == NULL) {
        errlog(109, 0);
        return 1;
    }
    output[OUTPUT_Y].spectrum = model->Y;
    /*
     * Allocation of L is a special case because the padded
     * length just might exceed the maximum capacity of
     * gridsize_t.
     */
    if (output[OUTPUT_DELAY].flags & (OUTPUT_USER | OUTPUT_FITTED)) {
        if (model->nLpad == 0) {
            errlog(40, 0);
            return 1;
        }
        if ((model->L =
                    (double*)malloc(model->ngrid * sizeof(double))) == NULL) {
            errlog(66, 0);
            return 1;
        }
        output[OUTPUT_DELAY].spectrum = model->L;
    }
    if ((model->tau_fsl =
                (double*)malloc(model->ngrid * sizeof(double))) == NULL) {
        errlog(193, 0);
        return 1;
    }
    output[OUTPUT_FREE_SPACE_LOSS].spectrum = model->tau_fsl;
    if ((model->k_out =
                (double*)malloc(model->ngrid * sizeof(double))) == NULL) {
        errlog(197, 0);
        return 1;
    }
    output[OUTPUT_K].spectrum = model->k_out;
    for (lnum = 0; lnum < model->nlayers; ++lnum) {
        if (alloc_layer_arrays(model, lnum))
                return 1;
    }
    return 0;
}   /* alloc_model_arrays() */


/***********************************************************
* int alloc_jacobians(model_t *model, simplex_t *simplex)
*
* Purpose:
*   For each output array for which the OUTPUT_JACOBIAN flag
*   bit is set, this function allocates an array for the
*   Jacobian with respect to each model variable named in
*   simplex, and an array for the numerical derivative error
*   estimates.
*
* Arguments:
*   model_t *model - pointer to model structure
*   simplex_t *simplex - pointer to simplex structure
* Return:
*   0 on success, 1 otherwise
************************************************************/

int alloc_jacobians(model_t *model, simplex_t *simplex)
{
    unsigned int k;
    
    for (k = 1; k < OUTPUT_END_OF_TABLE; ++k) {
        if (output[k].flags & OUTPUT_JACOBIAN) {
            unsigned int i;
            if (((output[k].jacobian = (double**)malloc(
                        simplex->n * sizeof(double*))) == NULL) ||
                ((output[k].jacobian_round_err = (double**)malloc(
                        simplex->n * sizeof(double*))) == NULL) ||
                ((output[k].jacobian_trunc_err = (double**)malloc(
                        simplex->n * sizeof(double*))) == NULL) ||
                ((output[k].jacobian_sorted_err = (double**)malloc(
                        simplex->n * sizeof(double*))) == NULL)) {
                errlog(97, k);
                return 1;
            }
            for (i = 0; i < simplex->n; ++i) {
                if (((output[k].jacobian[i] = (double*)malloc(
                            model->ngrid * sizeof(double))) == NULL) ||
                    ((output[k].jacobian_round_err[i] = (double*)malloc(
                            model->ngrid * sizeof(double))) == NULL) ||
                    ((output[k].jacobian_trunc_err[i] = (double*)malloc(
                            model->ngrid * sizeof(double))) == NULL) ||
                    ((output[k].jacobian_sorted_err[i] = (double*)malloc(
                            model->ngrid * sizeof(double))) == NULL)) {
                    errlog(98, i);
                    return 1;
                }
            }
        }
    }
    return 0;
}   /* alloc_jacobians() */


/***********************************************************
* int clear_layer_kcache_entries(model_t *model, int lnum)
*
* Purpose:
*   This function clears all the kcache entries for all
*   absorption coefficients on layer lnum of the model.
*   This is done when the kcache data have become invalid,
*   for example if the layer pressure has been changed with
*   the 'set' keyword in a fit data stream.
*
* Arguments:
*   model_t *model - pointer to model structure
*   int lnum - layer number
*
* Return:
*   0 if OK, 1 on error
************************************************************/

int clear_layer_kcache_entries(model_t *model, int lnum)
{
    layer_t *layer;
    int cnum;
    if (lnum < 0 || lnum >= model->nlayers) {
        errlog(105, lnum);
        return 1;
    }
    layer = model->layer[lnum];
    for (cnum = 0; cnum < layer->ncols; ++cnum) {
        column_t *column = layer->column[cnum];
        int knum;
        for (knum = 0; knum < column->n_abscoeffs; ++knum) {
            abscoeff_t *abscoeff = column->abscoeff[knum];
            int j;
            if (abscoeff->kcache == NULL)
                continue; /* no kcache for this absorption coefficient */
            for (j = 0; j < model->nkcache; ++j) {
                kcache_free(abscoeff->kcache[j]);
                abscoeff->kcache[j] = NULL;
            }
        }
    }
    return 0;
}   /* clear_layer_kcache_entries() */


/***********************************************************
* int copy_layer_allocations(model_t *model, layer_t *player, layer_t *layer)
*
* Purpose:
*   Allocates spectral arrays on *layer, duplicating the
*   array allocation pattern of a parent layer *player
*   having the same dimensions.
*
* Arguments:
*   model_t *model  - pointer to model structure
*   layer_t *player - pointer to parent layer
*   layer_t *layer  - pointer to a newly-created layer
*
* Return:
*   0 if successful, 1 otherwise
************************************************************/

int copy_layer_allocations(model_t *model, layer_t *player, layer_t *layer)
{
    int cnum;

    if (layer->ncols != player->ncols) {
        errlog(173, 0);
        return 1;
    }
    for (cnum = 0; cnum < player->ncols; ++cnum) {
        int knum;
        column_t  *column =  layer->column[cnum];
        column_t *pcolumn = player->column[cnum];
        if (column->n_abscoeffs != pcolumn->n_abscoeffs) {
            errlog(173, 0);
            return 1;
        }
        for (knum = 0; knum < column->n_abscoeffs; ++knum) {
            if ((pcolumn->abscoeff[knum]->k != NULL) &&
                    ((column->abscoeff[knum]->k =
                      (double*)malloc(model->ngrid * sizeof(double)))
                     == NULL)) {
                errlog(181, 0);
                return 1;
            }
        }
        if ((pcolumn->ztau != NULL) && 
                ((column->ztau =
                  (double*)malloc(model->ngrid * sizeof(double)))
                 == NULL)) {
            errlog(181, 0);
            return 1;
        }
    }
    if ((player->B != NULL) &&
            ((layer->B = (double*)malloc(model->ngrid * sizeof(double)))
             == NULL)) {
        errlog(181, 0);
        return 1;
    } else {
    }
    if ((player->tau != NULL) &&
            ((layer->tau = (double*)malloc(model->ngrid * sizeof(double)))
             == NULL)) {
        errlog(181, 0);
        return 1;
    }
    if ((player->tx != NULL) &&
            ((layer->tx = (double*)malloc(model->ngrid * sizeof(double)))
             == NULL)) {
        errlog(181, 0);
        return 1;
    }
    return 0;
}   /* copy_layer_allocations() */


/***********************************************************
* int copy_layer_dimensions(layer_t *player, layer_t *layer)
*
* Purpose:
*   Duplicates the column structure of parent layer on a
*   newly-created layer.
*
* Arguments:
*   layer_t *player - pointer to parent layer
*   layer_t *layer  - pointer to a newly-created layer
*
* Return:
*   0 if successful, 1 otherwise
************************************************************/

int copy_layer_dimensions(layer_t *player, layer_t *layer)
{
    int cnum;

    for (cnum = 0; cnum < player->ncols; ++cnum) {
        int col_typenum = player->column[cnum]->col_typenum;
        if (add_column(layer, col_typenum))
            return 1;
    }
    return 0;
}   /* copy_layer_dimensions() */


/***********************************************************
* int copy_model_dimensions(model_t *model, model_t *lmodel)
*
* Purpose:
*   Duplicates the layer and column structure of model in
*   the initially empty model structure lmodel.  lmodel is
*   assumed to have been initialized to MODEL_INIT, which
*   sets the array sizes of lmodel to zero.
*
* Arguments:
*   model_t *model, *lmodel - pointers to model structures
*
* Return:
*   0 if successful, 1 otherwise
************************************************************/

int copy_model_dimensions(model_t *model, model_t *lmodel)
{
    int lnum;

    for (lnum = 0; lnum < model->nlayers; ++lnum) {
        if (add_layer(lmodel, model->layer[lnum]->type))
            return 1;
        if (copy_layer_dimensions(model->layer[lnum], lmodel->layer[lnum]))
            return 1;
    }
    return 0;
}   /* copy_model_dimensions() */


/***********************************************************
* int delete_layer(model_t *model, int lnum)
*
* Purpose:
*   Frees memory associated with model layer lnum and
*   removes the corresponding entry from the layer table.
*   The layer table entries are shifted to fill in the gap
*   left by the removed layer, and the number of layers is
*   decremented.  The layer table size is left unchanged,
*   and the extra slot at the end is set to NULL.
*
*   This function only deletes the layer, it does not handle
*   any needed adjustments to adjacent layer parameters,
*   which is the responsibility of the calling function.
*   The calling function is also responsible for keeping
*   model structures expected to have identical dimensions
*   (i.e. model, lmodel) in sync through consistent calls
*   to this function.
*
* Arguments:
*   model_t *model - pointer to model structure
*   int lnum       - index of layer to be deleted
*
* Return:
*   0 if successful, 1 otherwise
************************************************************/

int delete_layer(model_t *model, int lnum)
{
    if (lnum < 0 || lnum > model->nlayers - 1)
        return 1;
    free_layer(model, lnum);
    --model->nlayers;
    for (;lnum < model->nlayers; ++lnum)
        model->layer[lnum] = model->layer[lnum+1];
    model->layer[model->nlayers] = NULL;
    return 0;
}   /* delete_layer() */


/***********************************************************
* static void free_abscoeffs(column_t *column, int nkcache)
*
* Purpose:
*   Frees the memory allocated to all the absorption
*   coefficients associated with a column, and sets
*   the abscoeff array pointer to NULL.
*
* Arguments:
*   column_t *column - pointer to column structure
*   int nkcache - number of kcache entries
************************************************************/

static void free_abscoeffs(column_t *column, int nkcache)
{
    int i, j;

    /*
     * Free the allocated memory associated with each absorption
     * coefficient structure.
     */
    for (i = 0; i < column->n_abscoeffs; ++i) {
        free(column->abscoeff[i]->k);
        if (column->abscoeff[i]->kcache != NULL) {
            for (j = 0; j < nkcache; ++j)
                kcache_free(column->abscoeff[i]->kcache[j]);
            free(column->abscoeff[i]->kcache);
        }
        free(column->abscoeff[i]);
    }
    /*
     * Free the abscoeff pointer array.
     */
    free(column->abscoeff);
    column->abscoeff = NULL;
    return;
}   /* free_abscoeffs() */


/***********************************************************
* static void free_columns(layer_t *layer, int nkcache)
*
* Purpose:
*   Frees the memory allocated to all the columns within a
*   layer, and sets the column array pointer to NULL.
*
* Arguments:
*   layer_t *layer - pointer to layer structure
*   int lnum - layer number
*   int nkcache - number of kcache entries
************************************************************/

static void free_columns(layer_t *layer, int nkcache)
{
    int cnum;

    /*
     * Free the allocated memory associated with each column.
     */
    for (cnum = 0; cnum < layer->ncols; ++cnum) {
        free_abscoeffs(layer->column[cnum], nkcache);
        free(layer->column[cnum]->ztau);
        free(layer->column[cnum]);
    }
    /*
     * Free the column pointer array.
     */
    free(layer->column);
    layer->column = NULL;
    return;
}   /* free_columns() */


/***********************************************************
* void free_fit_data_entities(fit_data_t *fit_data)
*
* Purpose:
*   Frees the memory allocated to all the entities within
*   a fit_data structure and NULLs the pointers.
*
* Arguments:
*   fit_data_t *fit_data pointer to fit_data structure.
************************************************************/

void free_fit_data_entities(fit_data_t *fit_data)
{
    int i;
    for (i = 0; i < fit_data->nfiles; ++i)
        free(fit_data->filename[i]);
    free(fit_data->filename);
    fit_data->filename = NULL;
    free(fit_data->f);
    fit_data->f = NULL;
    free(fit_data->s);
    fit_data->s = NULL;
    free(fit_data->s_mod);
    fit_data->s_mod = NULL;
    free(fit_data->res);
    fit_data->res = NULL;
    free(fit_data->res_est);
    fit_data->res_est = NULL;
    free(fit_data->b);
    fit_data->b = NULL;
    free(fit_data->w);
    fit_data->w = NULL;
    return;
}   /* free_fit_data_entities() */


/***********************************************************
* void free_jacobians(simplex_t *simplex)
*
* Purpose:
*   Frees memory allocated for Jacobian arrays and array
*   tables in the global output table.
************************************************************/
void free_jacobians(simplex_t *simplex)
{
    unsigned int k;

    for (k = 1; k < OUTPUT_END_OF_TABLE; ++k) {
        if (output[k].jacobian != NULL) {
            unsigned int i;
            for (i = 0; i < simplex->n; ++i) {
                free(output[k].jacobian[i]);
                free(output[k].jacobian_round_err[i]);
                free(output[k].jacobian_trunc_err[i]);
                free(output[k].jacobian_sorted_err[i]);
            }
            free(output[k].jacobian);
            free(output[k].jacobian_round_err);
            free(output[k].jacobian_trunc_err);
            free(output[k].jacobian_sorted_err);
        }
    }
    return;
} /* free_jacobians() */


/***********************************************************
* static void free_layer(model_t *model, int lnum)
*
* Purpose:
*   Frees the memory allocated to a layer.
*
* Arguments:
*   model_t *model - pointer to model structure
*   int lnum - layer number to be freed
************************************************************/

static void free_layer(model_t *model, int lnum)
{
    layer_t *layer = model->layer[lnum];

    free_columns(layer, model->nkcache);
    free(layer->B);
    free(layer->tau);
    free(layer->tx);
    free(layer->lineshape);
    free(layer->strict_selfbroad);
    free(layer->Mair_flag);
    free(layer);
    return;
}   /* free_layer() */


/***********************************************************
* static void free_model_layers(model_t *model)
*
* Purpose:
*   Frees the memory allocated to each of the layers in a
*   model, and sets the layer array pointer to NULL.
*
* Arguments:
*   model_t *model - pointer to model structure
************************************************************/

static void free_model_layers(model_t *model)
{
    int lnum;

    for (lnum = 0; lnum < model->nlayers; ++lnum)
        free_layer(model, lnum);
    free(model->layer);
    model->layer = NULL;
    return;
}   /* free_model_layers() */


/***********************************************************
* void free_model_entities(model_t *model)
*
* Purpose:
*   Frees the memory allocated to all the entities within
*   a model data structure and NULLs the pointers.
*
* Arguments:
*   model_t *model - pointer to model structure.
************************************************************/

void free_model_entities(model_t *model)
{
    free_model_layers(model);
    free(model->f);
    model->f = NULL;
    free(model->f2);
    model->f2 = NULL;
    free(model->fif);
    model->fif = NULL;
    output[OUTPUT_FREQUENCY].spectrum = NULL;
    free(model->ils);
    model->ils = NULL;
    free(model->ilsworkspace);
    model->ilsworkspace = NULL;
    free(model->tau);
    model->tau = NULL;
    output[OUTPUT_OPACITY].spectrum = NULL;
    free(model->tx);
    model->tx = NULL;
    output[OUTPUT_TRANSMITTANCE].spectrum = NULL;
    free(model->I0);
    model->I0 = NULL;
    free(model->I);
    model->I = NULL;
    output[OUTPUT_RADIANCE].spectrum = NULL;
    free(model->I_ref);
    model->I_ref = NULL;
    free(model->I_diff);
    model->I_diff = NULL;
    output[OUTPUT_RADIANCE_DIFF].spectrum = NULL;
    free(model->Tb);
    model->Tb = NULL;
    output[OUTPUT_TB_PLANCK].spectrum = NULL;
    free(model->Trj);
    model->Trj = NULL;
    output[OUTPUT_TB_RAYLEIGH_JEANS].spectrum = NULL;
    free(model->Tsys);
    model->Tsys = NULL;
    output[OUTPUT_TSYS].spectrum = NULL;
    free(model->Y);
    model->Y = NULL;
    output[OUTPUT_Y].spectrum = NULL;
    free(model->L);
    model->L = NULL;
    output[OUTPUT_DELAY].spectrum = NULL;
    free(model->tau_fsl);
    model->tau_fsl = NULL;
    output[OUTPUT_FREE_SPACE_LOSS].spectrum = NULL;
    return;
}   /* free_model_entities() */


/***********************************************************
* int grow_fit_data_arrays(fit_data_t *fit_data)
*
* Purpose:
*   Increases the memory allocated to all the arrays within
*   a fit_data structure, or makes an initial allocation.
*
* Arguments:
*   fit_data_t *fit_data pointer to fit_data structure.
*
* Return:
*   0 if OK, 1 on realloc() error.
************************************************************/

int grow_fit_data_arrays(fit_data_t *fit_data)
{
    void *tptr;

    if (fit_data->nalloc == 0) {
        fit_data->nalloc = FIT_DATA_NALLOC_INIT;
    } else {
        fit_data->nalloc *= FIT_DATA_NALLOC_GROW;
    }
    if ((tptr = realloc(fit_data->f, fit_data->nalloc * sizeof(double)))
            == NULL)
        return 1;
    fit_data->f = (double *)tptr;
    if ((tptr = realloc(fit_data->s, fit_data->nalloc * sizeof(double)))
            == NULL)
        return 1;
    fit_data->s = (double *)tptr;
    if ((tptr = realloc(fit_data->s_mod, fit_data->nalloc * sizeof(double)))
            == NULL)
        return 1;
    fit_data->s_mod = (double *)tptr;
    if ((tptr = realloc(fit_data->res, fit_data->nalloc * sizeof(double)))
            == NULL)
        return 1;
    fit_data->res = (double *)tptr;
    if ((tptr = realloc(fit_data->res_est, fit_data->nalloc * sizeof(double)))
            == NULL)
        return 1;
    fit_data->res_est = (double *)tptr;
    if ((tptr = realloc(fit_data->b, fit_data->nalloc * sizeof(double)))
            == NULL)
        return 1;
    fit_data->b = (double *)tptr;
    if ((tptr = realloc(fit_data->w, fit_data->nalloc * sizeof(double)))
            == NULL)
        return 1;
    fit_data->w = (double *)tptr;
    return 0;
}   /* grow_fit_data_arrays() */


/***********************************************************
* int init_kcache(abscoeff_t *abscoeff, int nkcache)
*
* Purpose:
*   This function checks whether the absorption coefficient
*   type of abscoeff is cacheable.  If it is, it allocates
*   memory for nkcache entries in the kcache table, and
*   initializes all the table entries to NULL.
*
* Arguments:
*   abscoeff_t *abscoeff - pointer to abscoeff structure
*   int nkcache          - number of kcache entries
*
* Return:
*   1 on error, 0 otherwise.
************************************************************/

int init_kcache(abscoeff_t *abscoeff, int nkcache)
{
    int i;

    if (!(k_type[abscoeff->k_typenum].comp_flags & CACHEABLE))
        return 0;
    if (abscoeff->kcache != NULL) {
        errlog(50, 0);
        return 1;
    }
    if ((abscoeff->kcache = (double**)malloc(nkcache * sizeof(double*)))
            == NULL) {
        errlog(51, 0);
        return 1;
    }
    for (i = 0; i < nkcache; ++i)
        abscoeff->kcache[i] = NULL;
    return 0;
}   /* init_kcache() */
