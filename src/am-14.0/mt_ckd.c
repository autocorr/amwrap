/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* mt_ckd.c                      S. Paine rev. 2024 August 06
*
* Water vapor continuum data from the MT_CKD continuum
* model from AER, Inc.  See comments below for copyright,
* version, and usage information.  Note that data in this
* file are only a subset of MT_CKD covering am's frequency
* range.
************************************************************/

#include "mt_ckd.h"

/*
 * Spectral absorption coefficients for water vapor self- and
 * foreign-induced continuum absorption are from the MT_CKD
 * (Mlawer, Tobin, Clough, Kneizys, Davies) continuum model.
 * MT_CKD is maintained by the Radiative Transfer Working Group
 * at AER, Inc., and is distributed with the following copyright,
 * reference, and version information:
 *
 *   The MT_CKD Water Vapor Continuum - 4.2
 *
 *   Copyright ©, Atmospheric and Environmental Research, Inc.,
 *   2022  All rights reserved. This source code was developed as
 *   part of the LBLRTM software and is designed for scientific
 *   and research purposes. Atmospheric and Environmental
 *   Research Inc. (AER) grants USER the right to download,
 *   install, use and copy this software and data for scientific
 *   and research purposes only. This software and data may be
 *   redistributed as long as this copyright notice is reproduced
 *   on any copy made and appropriate acknowledgment is given to
 *   AER. This software and data or any modified version of this
 *   software and data may not be incorporated into proprietary
 *   software or data or commercial software or data offered for
 *   sale without the express written consent of AER. This
 *   software and data are provided as is without any expressed
 *   or implied warranties.
 *
 *   Address questions to: aer_contnm@aer.com
 *
 *   General reference: Mlawer et al. (2012),
 *   doi:10.1098/rsta.2011.0295
 *
 *   MT_CKD_4.2: Changes to self, foreign continuum, and self
 *   temperature dependance in IR window (590-1400 cm-1);
 *   Mlawer/Mascio/Turner. For previous version information see
 *   GitHub repository Whats New
 *   Wiki(http://github.com/AER-RC/MT_CKD).
 *
 * See also:
 *
 *   E.J. Mlawer, K.E.Cady-Periera, J. Mascio, and I.E. Gordon
 *   2023, "The inclusion of the MT_CKD water vapor continuum
 *   model in the HITRAN molecular spectroscopic database."
 *   JQSRT 306:108645.
 *   (https://doi.org/10.1016/j.jqsrt.2023.108645)
 *
 * and references therein.
 */

const double mt_ckd_ref_press = 1013.0;
const double mt_ckd_ref_temp  = 296.0;

const double mt_ckd_for_absco_ref[] = {
 8.562900e-23, 8.827840e-23, 9.182100e-23, 9.796020e-23, 1.139710e-22,
 1.376484e-22, 1.544904e-22, 1.507922e-22, 1.386438e-22, 1.215105e-22,
 1.062604e-22, 9.151380e-23, 8.008000e-23, 6.970302e-23, 6.054750e-23,
 5.243129e-23, 4.559360e-23, 3.762800e-23, 3.054901e-23, 2.562640e-23,
 2.197152e-23, 1.809267e-23, 1.483440e-23, 1.167269e-23, 8.809307e-24,
 7.130614e-24, 5.932630e-24, 5.171887e-24, 4.220168e-24, 3.427428e-24,
 2.685861e-24, 2.095871e-24, 1.858127e-24, 1.475934e-24, 1.100543e-24,
 8.588344e-25, 6.498760e-25, 4.874671e-25, 3.860516e-25, 3.000811e-25,
 2.343264e-25, 2.087365e-25, 1.771056e-25, 1.588460e-25, 1.408050e-25,
 1.180898e-25, 1.024386e-25, 9.223488e-26, 7.915559e-26, 6.692405e-26,
 5.720098e-26, 5.124893e-26, 4.477977e-26,
};

const double mt_ckd_self_absco_ref[] = {
 2.877000e-21, 2.855000e-21, 2.731000e-21, 2.490000e-21, 2.178000e-21,
 1.863000e-21, 1.595000e-21, 1.381000e-21, 1.221000e-21, 1.095000e-21,
 9.882000e-22, 8.926000e-22, 8.041000e-22, 7.217000e-22, 6.450000e-22,
 5.740000e-22, 5.088000e-22, 4.492000e-22, 3.952000e-22, 3.466000e-22,
 3.031000e-22, 2.644000e-22, 2.302000e-22, 2.001000e-22, 1.738000e-22,
 1.509000e-22, 1.311000e-22, 1.139000e-22, 9.913000e-23, 8.639000e-23,
 7.544000e-23, 6.603000e-23, 5.794000e-23, 5.096000e-23, 4.495000e-23,
 3.974000e-23, 3.523000e-23, 3.130000e-23, 2.786000e-23, 2.485000e-23,
 2.221000e-23, 1.987000e-23, 1.780000e-23, 1.596000e-23, 1.433000e-23,
 1.287000e-23, 1.156000e-23, 1.040000e-23, 9.349000e-24, 8.410000e-24,
 7.568000e-24, 6.810000e-24, 6.130000e-24,
};

const double mt_ckd_self_texp[] = {
 6.413000e+00, 6.414000e+00, 6.275000e+00, 6.049000e+00, 5.789000e+00,
 5.557000e+00, 5.330000e+00, 5.107000e+00, 4.893000e+00, 4.685000e+00,
 4.481000e+00, 4.282000e+00, 4.090000e+00, 3.905000e+00, 3.720000e+00,
 3.548000e+00, 3.379000e+00, 3.215000e+00, 3.056000e+00, 2.903000e+00,
 2.756000e+00, 2.614000e+00, 2.475000e+00, 2.344000e+00, 2.219000e+00,
 2.100000e+00, 1.983000e+00, 1.873000e+00, 1.768000e+00, 1.670000e+00,
 1.575000e+00, 1.486000e+00, 1.450000e+00, 1.432000e+00, 1.429000e+00,
 1.431000e+00, 1.428000e+00, 1.458000e+00, 1.517000e+00, 1.570000e+00,
 1.628000e+00, 1.689000e+00, 1.753000e+00, 1.819000e+00, 1.886000e+00,
 1.976000e+00, 2.032000e+00, 2.106000e+00, 2.190000e+00, 2.259000e+00,
 2.330000e+00, 2.397000e+00, 2.466000e+00,
};

const double mt_ckd_wavenumbers[] = {
 0.000000e+00, 1.000000e+01, 2.000000e+01, 3.000000e+01, 4.000000e+01,
 5.000000e+01, 6.000000e+01, 7.000000e+01, 8.000000e+01, 9.000000e+01,
 1.000000e+02, 1.100000e+02, 1.200000e+02, 1.300000e+02, 1.400000e+02,
 1.500000e+02, 1.600000e+02, 1.700000e+02, 1.800000e+02, 1.900000e+02,
 2.000000e+02, 2.100000e+02, 2.200000e+02, 2.300000e+02, 2.400000e+02,
 2.500000e+02, 2.600000e+02, 2.700000e+02, 2.800000e+02, 2.900000e+02,
 3.000000e+02, 3.100000e+02, 3.200000e+02, 3.300000e+02, 3.400000e+02,
 3.500000e+02, 3.600000e+02, 3.700000e+02, 3.800000e+02, 3.900000e+02,
 4.000000e+02, 4.100000e+02, 4.200000e+02, 4.300000e+02, 4.400000e+02,
 4.500000e+02, 4.600000e+02, 4.700000e+02, 4.800000e+02, 4.900000e+02,
 5.000000e+02, 5.100000e+02, 5.200000e+02,
};

const int mt_ckd_ntab =
    (int)(sizeof(mt_ckd_wavenumbers) / sizeof(mt_ckd_wavenumbers[0]));

