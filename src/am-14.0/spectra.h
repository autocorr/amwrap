/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* spectra.h                         S. Paine rev. 2023 May 5
*
* Declarations for spectra.c
************************************************************/

#ifndef AM_SPECTRA_H
#define AM_SPECTRA_H

#include "am_types.h"

int  compare_spectral_subgrid_ranges(model_t*, model_t*);
void compute_k_out(model_t*);
int  compute_delay_spectrum(model_t*);
void compute_free_space_loss(model_t*);
int  compute_if_spectra(model_t*);
void compute_radiance_difference_spectrum(model_t*, model_t*);
void compute_spectral_Tsys(model_t*);
void compute_spectral_Y_factor(model_t*);
void compute_Tb(model_t*);
void compute_transmittance(model_t*);
void compute_Trj(model_t*);
int  set_IF_spectrum_subgrid_ranges(model_t*);
int  set_spectral_subgrid_ranges(model_t*, model_t*);

#endif /* AM_SPECTRA_H */
