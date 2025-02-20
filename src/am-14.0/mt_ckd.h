/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* mt_ckd.h                        S. Paine rev. 2023 June 27
*
* External declarations for water vapor continuum data from
* the MT_CKD continuum model, defined in mt_ckd.c.  MT_CKD
* is maintained by the Radiative Transfer Working Group at
* AER, Inc.  See mt_ckd.c for copyright, reference, and
* version information.
************************************************************/

#ifndef AM_MT_CKD_H
#define AM_MT_CKD_H

extern const double mt_ckd_ref_press;
extern const double mt_ckd_ref_temp;
extern const double mt_ckd_for_absco_ref[];
extern const double mt_ckd_self_absco_ref[];
extern const double mt_ckd_self_texp[];
extern const double mt_ckd_wavenumbers[];
extern const int    mt_ckd_ntab;

#endif /* AM_MT_CKD_H */
