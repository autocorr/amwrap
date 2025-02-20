/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* molecules.h                   S. Paine rev. 2024 August 06
*
* Constants for molecular masses, isotopic abundances, and
* default atmospheric mixing ratios
************************************************************/


#ifndef AM_MOLECULES_H
#define AM_MOLECULES_H

/*
 * Dry air
 *
 * Average molar mass from
 *
 *  R.M. Goody and Y.L. Yung 1989, "Atmospheric Radiation" 2nd ed.
 *  Oxford University Press.
 *
 * This is a standard value.  It needn't agree exactly with the
 * weighted average molar mass using the default mixing ratios defined
 * below, but it does so within their realistic variability in time
 * and space.
 */
#define MASS_DRY_AIR          2.896400e+01
#define FRAC_DRY_AIR          1.000000e+00


/* 
 * N2
 *
 * Dry air mixing ratio from Table 1.1 in
 *
 *  R.M. Goody and Y.L. Yung 1989, "Atmospheric Radiation" 2nd ed.
 *  Oxford University Press.
 * 
 * Average molar mass from NIST Chemistry WebBook, http://webbook.nist.gov/ 
 */
#define MASS_N2               2.801340e+01
#define FRAC_N2               7.808400e-01


/*
 * O2
 *
 * Dry air mixing ratio from Table 1.1 in
 *
 *  R.M. Goody and Y.L. Yung 1989, "Atmospheric Radiation" 2nd ed.
 *  Oxford University Press.
 * 
 * See also:
 * 
 *  L. Machta and E. Hughes 1970, "Atmospheric Oxygen in 1967 to
 *  1970." Science 168:1582.
 *
 * Relative to this reference value, the O2 concentration is subject
 * to an annual cycle and secular drift; presently, O2/N2 varies with
 * an amplitude of about 100 per meg annually and drifts downward by
 * about 20 per meg annually.  See:
 *
 *  L.N.T. Nguyen et al. 2022, "Two decades of flask observations of
 *  atmospheric delta(O2/N2), CO2, and APO at stations Lutjewad (the
 *  Netherlands) and Mace Head (Ireland), and three years from Halley
 *  station (Antarctica)."  Earth Syst. Sci. Data 14:991.
 *
 * Molar mass and abundance data are from the HITRAN molparam.txt file
 */
#define MASS_O2               3.199858e+01
#define MASS_16O2             3.198983e+01
#define MASS_16O_18O          3.399408e+01
#define MASS_16O_17O          3.299404e+01

#define ABUN_16O2             9.952620e-01
#define ABUN_16O_18O          3.991410e-03
#define ABUN_16O_17O          7.422350e-04

#define FRAC_O2               2.094600e-01
#define FRAC_16O2             (FRAC_O2 * ABUN_16O2)
#define FRAC_16O_18O          (FRAC_O2 * ABUN_16O_18O)
#define FRAC_16O_17O          (FRAC_O2 * ABUN_16O_17O)


/* 
 * Ar
 *
 * Dry air mixing ratio from Table 1.1 in
 *
 *  R.M. Goody and Y.L. Yung 1989, "Atmospheric Radiation" 2nd ed.
 *  Oxford University Press.
 * 
 * Average molar mass from NIST Chemistry WebBook, http://webbook.nist.gov/ 
 */
#define MASS_AR               3.994800e+01
#define FRAC_AR               9.340000e-03


/*
 * CO2
 *
 * The default dry air mixing ratio for co2 is based on the 2019
 * globally-averaged surface value from Table 2.2 in
 *
 *  Gulev, S.K., P.W. Thorne, J. Ahn, F.J. Dentener, C.M. Domingues,
 *  S. Gerland, D. Gong, D.S. Kaufman, H.C. Nnamchi, J. Quaas, J.A.
 *  Rivera, S. Sathyendranath, S.L. Smith, B. Trewin, K. von
 *  Schuckmann, and R.S. Vose, 2021, "Changing State of the Climate
 *  System." In "Climate Change 2021: The Physical Science Basis.
 *  Contribution of Working Group I to the Sixth Assessment Report of
 *  the Intergovernmental Panel on Climate Change" [Masson-Delmotte, V.,
 *  P. Zhai, A. Pirani, S.L. Connors, C. Péan, S. Berger, N. Caud, Y.
 *  Chen, L. Goldfarb, M.I. Gomis, M. Huang, K. Leitzell, E. Lonnoy,
 *  J.B.R. Matthews, T.K. Maycock, T. Waterfield, O. Yelekçi, R. Yu,
 *  and B. Zhou (eds.)]. Cambridge University Press, Cambridge, United
 *  Kingdom and New York, NY, USA, pp. 287–422,
 *  doi:10.1017/9781009157896.004.
 *
 * which is 409.9 ppm, a secular increase of about 20 ppm since 2011,
 * on top of a latitude-dependent seasonal cycle with amplitude of
 * order 10 ppm.  The value here is rounded accordingly.  CO2 is
 * well-mixed in the terrestrial atmosphere up to 100 km.
 *
 * Molar mass and abundance data are from the HITRAN molparam.txt file
 */
#define MASS_CO2              4.400974e+01
#define MASS_12C_16O2         4.398983e+01
#define MASS_13C_16O2         4.499318e+01
#define MASS_16O_12C_18O      4.599408e+01
#define MASS_16O_12C_17O      4.499404e+01
#define MASS_16O_13C_18O      4.699743e+01
#define MASS_16O_13C_17O      4.599740e+01
#define MASS_12C_18O2         4.799832e+01

#define ABUN_12C_16O2         9.842040e-01
#define ABUN_13C_16O2         1.105740e-02
#define ABUN_16O_12C_18O      3.947070e-03
#define ABUN_16O_12C_17O      7.339890e-04
#define ABUN_16O_13C_18O      4.434460e-05
#define ABUN_16O_13C_17O      8.246230e-06
#define ABUN_12C_18O2         3.957340e-06

#define FRAC_CO2              4.100000e-04
#define FRAC_12C_16O2         (FRAC_CO2 * ABUN_12C_16O2)
#define FRAC_13C_16O2         (FRAC_CO2 * ABUN_13C_16O2)
#define FRAC_16O_12C_18O      (FRAC_CO2 * ABUN_16O_12C_18O)
#define FRAC_16O_12C_17O      (FRAC_CO2 * ABUN_16O_12C_17O)
#define FRAC_16O_13C_18O      (FRAC_CO2 * ABUN_16O_13C_18O)
#define FRAC_16O_13C_17O      (FRAC_CO2 * ABUN_16O_13C_17O)
#define FRAC_12C_18O2         (FRAC_CO2 * ABUN_12C_18O2)


/*
 * The next three trace species are included in the dry_air column
 * type, but not in the dry_air_std column type.
 */

/*
 * CH4
 *
 * The default dry air mixing ratio for ch4 is based on the 2019
 * globally-averaged surface value from Table 2.2 in
 *
 *  Gulev, S.K., P.W. Thorne, J. Ahn, F.J. Dentener, C.M. Domingues,
 *  S. Gerland, D. Gong, D.S. Kaufman, H.C. Nnamchi, J. Quaas, J.A.
 *  Rivera, S. Sathyendranath, S.L. Smith, B. Trewin, K. von
 *  Schuckmann, and R.S. Vose, 2021, "Changing State of the Climate
 *  System." In "Climate Change 2021: The Physical Science Basis.
 *  Contribution of Working Group I to the Sixth Assessment Report of
 *  the Intergovernmental Panel on Climate Change" [Masson-Delmotte, V.,
 *  P. Zhai, A. Pirani, S.L. Connors, C. Péan, S. Berger, N. Caud, Y.
 *  Chen, L. Goldfarb, M.I. Gomis, M. Huang, K. Leitzell, E. Lonnoy,
 *  J.B.R. Matthews, T.K. Maycock, T. Waterfield, O. Yelekçi, R. Yu,
 *  and B. Zhou (eds.)]. Cambridge University Press, Cambridge, United
 *  Kingdom and New York, NY, USA, pp. 287–422,
 *  doi:10.1017/9781009157896.004.
 *
 * which is 1866.3 ppb, a secular increase of about 60 ppb since 2011.
 * The value here is rounded accordingly.  CH4 is well-mixed in the
 * troposphere, then the mixing ratio falls steadily above the
 * tropopause.
 * 
 * Molar mass and abundance data are from the HITRAN molparam.txt file
 */
#define MASS_CH4              1.604306e+01
#define MASS_12CH4            1.603130e+01
#define MASS_13CH4            1.703466e+01
#define MASS_12CH3D           1.703748e+01

#define ABUN_12CH4            9.882740e-01
#define ABUN_13CH4            1.110310e-02
#define ABUN_12CH3D           6.157510e-04

#define FRAC_CH4              1.900000e-06
#define FRAC_12CH4            (FRAC_CH4 * ABUN_12CH4)
#define FRAC_13CH4            (FRAC_CH4 * ABUN_13CH4)
#define FRAC_12CH3D           (FRAC_CH4 * ABUN_12CH3D)


/*
 * N2O
 *
 * The default dry air mixing ratio for n2o is based on the 2019
 * globally-averaged surface value from Table 2.2 in
 *
 *  Gulev, S.K., P.W. Thorne, J. Ahn, F.J. Dentener, C.M. Domingues,
 *  S. Gerland, D. Gong, D.S. Kaufman, H.C. Nnamchi, J. Quaas, J.A.
 *  Rivera, S. Sathyendranath, S.L. Smith, B. Trewin, K. von
 *  Schuckmann, and R.S. Vose, 2021, "Changing State of the Climate
 *  System." In "Climate Change 2021: The Physical Science Basis.
 *  Contribution of Working Group I to the Sixth Assessment Report of
 *  the Intergovernmental Panel on Climate Change" [Masson-Delmotte, V.,
 *  P. Zhai, A. Pirani, S.L. Connors, C. Péan, S. Berger, N. Caud, Y.
 *  Chen, L. Goldfarb, M.I. Gomis, M. Huang, K. Leitzell, E. Lonnoy,
 *  J.B.R. Matthews, T.K. Maycock, T. Waterfield, O. Yelekçi, R. Yu,
 *  and B. Zhou (eds.)]. Cambridge University Press, Cambridge, United
 *  Kingdom and New York, NY, USA, pp. 287–422,
 *  doi:10.1017/9781009157896.004.
 *
 * which is 332.1 ppb, a secular increase of about 8 ppb since 2011.
 * The value here is rounded accordingly.  N2O is well-mixed in the
 * troposphere, then the mixing ratio falls steadily above the
 * tropopause.
 *
 * Molar mass and abundance data are from the HITRAN molparam.txt file
 */
#define MASS_N2O              4.401267e+01
#define MASS_14N2_16O         4.400106e+01
#define MASS_14N_15N_16O      4.499810e+01
#define MASS_15N_14N_16O      4.499810e+01
#define MASS_14N2_18O         4.600531e+01
#define MASS_14N2_17O         4.500528e+01

#define ABUN_14N2_16O         9.903330e-01
#define ABUN_14N_15N_16O      3.640930e-03
#define ABUN_15N_14N_16O      3.640930e-03
#define ABUN_14N2_18O         1.985820e-03
#define ABUN_14N2_17O         3.692800e-04

#define FRAC_N2O              3.300000e-07
#define FRAC_14N2_16O         (FRAC_N2O * ABUN_14N2_16O)
#define FRAC_14N_15N_16O      (FRAC_N2O * ABUN_14N_15N_16O)
#define FRAC_15N_14N_16O      (FRAC_N2O * ABUN_15N_14N_16O)
#define FRAC_14N2_18O         (FRAC_N2O * ABUN_14N2_18O)
#define FRAC_14N2_17O         (FRAC_N2O * ABUN_14N2_17O)


/*
 * CO
 *
 * Dry air mixing ratio from Table 1.1 in
 *
 *  R.M. Goody and Y.L. Yung 1989, "Atmospheric Radiation" 2nd ed.
 *  Oxford University Press.
 *
 * which is an average value across the height of the troposphere.
 * The mixing ratio declines from about 150 ppb at the surface to a
 * minimum of about 10 ppb in the stratosphere, then rises again.
 *
 * Molar mass and abundance data are from the HITRAN molparam.txt file
 */
#define MASS_CO               2.801045e+01
#define MASS_12C_16O          2.799491e+01
#define MASS_13C_16O          2.899827e+01
#define MASS_12C_18O          2.999916e+01
#define MASS_12C_17O          2.899913e+01
#define MASS_13C_18O          3.100252e+01
#define MASS_13C_17O          3.000249e+01

#define ABUN_12C_16O          9.865440e-01
#define ABUN_13C_16O          1.108360e-02
#define ABUN_12C_18O          1.978220e-03
#define ABUN_12C_17O          3.678670e-04
#define ABUN_13C_18O          2.222500e-05
#define ABUN_13C_17O          4.132920e-06

#define FRAC_CO               7.000000e-08
#define FRAC_12C_16O          (FRAC_CO * ABUN_12C_16O)
#define FRAC_13C_16O          (FRAC_CO * ABUN_13C_16O)
#define FRAC_12C_18O          (FRAC_CO * ABUN_12C_18O)
#define FRAC_12C_17O          (FRAC_CO * ABUN_12C_17O)
#define FRAC_13C_18O          (FRAC_CO * ABUN_13C_18O)
#define FRAC_13C_17O          (FRAC_CO * ABUN_13C_17O)


/*
 * The trace species beyond this point are not included in the
 * dry_air column type.
 */

/*
 * CH3CN
 *
 * Molar mass and abundance data are from the HITRAN molparam.txt file
 */
#define MASS_CH3CN            4.102655e+01
#define MASS_12CH3_12C14N     4.102655e+01

#define ABUN_12CH3_12C14N     9.738660e-01


/*
 * CH3OH
 *
 * Molar mass and abundance data are from the HITRAN molparam.txt file
 */
#define MASS_CH3OH            3.202622e+01
#define MASS_12CH3_16OH       3.202622e+01

#define ABUN_12CH3_16OH       9.859300e-01


/*
 * ClO
 *
 * Molar mass and abundance data are from the HITRAN molparam.txt file
 */
#define MASS_ClO              5.144764e+01
#define MASS_35Cl_16O         5.096377e+01
#define MASS_37Cl_16O         5.296082e+01

#define ABUN_35Cl_16O         7.559080e-01
#define ABUN_37Cl_16O         2.417200e-01


/*
 * HBr
 *
 * Molar mass and abundance data are from the HITRAN molparam.txt file
 */
#define MASS_HBr              8.091143e+01
#define MASS_H_79Br           7.992616e+01
#define MASS_H_81Br           8.192412e+01

#define ABUN_H_79Br           5.067810e-01
#define ABUN_H_81Br           4.930630e-01


/*
 * HCN
 *
 * Molar mass and abundance data are from the HITRAN molparam.txt file
 */
#define MASS_HCN              2.702562e+01
#define MASS_H_12C_14N        2.701090e+01
#define MASS_H_13C_14N        2.801425e+01
#define MASS_H_12C_15N        2.800793e+01

#define ABUN_H_12C_14N        9.851140e-01
#define ABUN_H_13C_14N        1.106760e-02
#define ABUN_H_12C_15N        3.621740e-03


/*
 * H2CO
 *
 * Molar mass and abundance data are from the HITRAN molparam.txt file
 */
#define MASS_H2CO             3.002566e+01
#define MASS_H2_12C_16O       3.001056e+01
#define MASS_H2_13C_16O       3.101392e+01
#define MASS_H2_12C_18O       3.201481e+01

#define ABUN_H2_12C_16O       9.862370e-01
#define ABUN_H2_13C_16O       1.108020e-02
#define ABUN_H2_12C_18O       1.977610e-03


/*
 * HCl
 *
 * Molar mass and abundance data are from the HITRAN molparam.txt file
 */
#define MASS_HCl              3.646055e+01
#define MASS_H_35Cl           3.597668e+01
#define MASS_H_37Cl           3.797373e+01

#define ABUN_H_35Cl           7.575870e-01
#define ABUN_H_37Cl           2.422570e-01


/*
 * HF
 *
 * Molar mass and abundance data are from the HITRAN molparam.txt file
 */
#define MASS_HF               2.000623e+01
#define MASS_H_19F            2.000623e+01

#define ABUN_H_19F            9.998440e-01


/*
 * HNO3
 *
 * Molar mass and abundance data are from the HITRAN molparam.txt file
 */
#define MASS_HNO3             6.299564e+01
#define MASS_H_14N_16O3       6.299564e+01

#define ABUN_H_14N_16O3       9.891100e-01


/*
 * H2O
 *
 * Molar mass and abundance data are from the HITRAN molparam.txt file
 */
#define MASS_H2O              1.801526e+01
#define MASS_H2_16O           1.801056e+01
#define MASS_H2_18O           2.001481e+01
#define MASS_H2_17O           1.901478e+01
#define MASS_HD_16O           1.901674e+01
#define MASS_HD_18O           2.102098e+01
#define MASS_HD_17O           2.002096e+01

#define ABUN_H2_16O           9.973170e-01
#define ABUN_H2_18O           1.999830e-03
#define ABUN_H2_17O           3.718840e-04
#define ABUN_HD_16O           3.106930e-04
#define ABUN_HD_18O           6.230030e-07
#define ABUN_HD_17O           1.158530e-07


/*
 * H2O2
 *
 * Molar mass and abundance data are from the HITRAN molparam.txt file
 */
#define MASS_H2O2             3.400548e+01
#define MASS_H2_16O2          3.400548e+01

#define ABUN_H2_16O2          9.949520e-01


/*
 * HO2
 *
 * Molar mass and abundance data are from the HITRAN molparam.txt file
 */
#define MASS_HO2              3.299766e+01
#define MASS_H_16O2           3.299766e+01

#define ABUN_H_16O2           9.951070e-01


/*
 * HOCl
 *
 * Molar mass and abundance data are from the HITRAN molparam.txt file
 */
#define MASS_HOCl             5.245547e+01
#define MASS_H_16O_35Cl       5.197159e+01
#define MASS_H_16O_37Cl       5.396864e+01

#define ABUN_H_16O_35Cl       7.557900e-01
#define ABUN_H_16O_37Cl       2.416830e-01


/*
 * H2S
 *
 * Molar mass and abundance data are from the HITRAN molparam.txt file
 */
#define MASS_H2S              3.407935e+01
#define MASS_H2_32S           3.398772e+01
#define MASS_H2_34S           3.598351e+01
#define MASS_H2_33S           3.498710e+01

#define ABUN_H2_32S           9.498840e-01
#define ABUN_H2_34S           4.213690e-02
#define ABUN_H2_33S           7.497660e-03


/*
 * NH3
 *
 * Molar mass and abundance data are from the HITRAN molparam.txt file
 */
#define MASS_NH3              1.703020e+01
#define MASS_14NH3            1.702655e+01
#define MASS_15NH3            1.802358e+01

#define ABUN_14NH3            9.958720e-01
#define ABUN_15NH3            3.661290e-03


/*
 * NO
 *
 * Molar mass and abundance data are from the HITRAN molparam.txt file
 */
#define MASS_NO               2.999799e+01
#define MASS_14N_16O          2.999799e+01

#define ABUN_14N_16O          9.939740e-01


/*
 * NO2
 *
 * Molar mass and abundance data are from the HITRAN molparam.txt file
 */
#define MASS_NO2              4.599290e+01
#define MASS_14N_16O2         4.599290e+01

#define ABUN_14N_16O2         9.916160e-01


/*
 * O
 *
 * Molar mass and abundance data are from the HITRAN molparam.txt file
 */
#define MASS_O                1.599492e+01
#define MASS_16O              1.599492e+01

#define ABUN_16O              9.976280e-01


/*
 * O3
 *
 * Molar mass and abundance data are from the HITRAN molparam.txt file
 */
#define MASS_O3               4.799783e+01
#define MASS_16O3             4.798474e+01
#define MASS_16O2_18O         4.998899e+01
#define MASS_16O_18O_16O      4.998899e+01
#define MASS_16O2_17O         4.898896e+01
#define MASS_16O_17O_16O      4.898896e+01

#define ABUN_16O3             9.929010e-01
#define ABUN_16O2_18O         3.981940e-03
#define ABUN_16O_18O_16O      1.990970e-03
#define ABUN_16O2_17O         7.404750e-04
#define ABUN_16O_17O_16O      3.702370e-04


/*
 * OCS
 *
 * Molar mass and abundance data are from the HITRAN molparam.txt file
 */
#define MASS_OCS              6.007183e+01
#define MASS_16O_12C_32S      5.996699e+01
#define MASS_16O_12C_34S      6.196278e+01
#define MASS_16O_13C_32S      6.097034e+01
#define MASS_16O_12C_33S      6.096637e+01
#define MASS_18O_12C_32S      6.197123e+01

#define ABUN_16O_12C_32S      9.373950e-01
#define ABUN_16O_12C_34S      4.158280e-02
#define ABUN_16O_13C_32S      1.053150e-02
#define ABUN_16O_12C_33S      7.399080e-03
#define ABUN_18O_12C_32S      1.879670e-03


/*
 * OH
 *
 * Molar mass and abundance data are from the HITRAN molparam.txt file
 */
#define MASS_OH               1.700691e+01
#define MASS_16O_H            1.700274e+01
#define MASS_18O_H            1.900699e+01
#define MASS_16O_D            1.800891e+01

#define ABUN_16O_H            9.974730e-01
#define ABUN_18O_H            2.000140e-03
#define ABUN_16O_D            1.553710e-04


/*
 * SO2
 *
 * Molar mass and abundance data are from the HITRAN molparam.txt file
 */
#define MASS_SO2              6.404667e+01
#define MASS_32S_16O2         6.396190e+01
#define MASS_34S_16O2         6.595770e+01

#define ABUN_32S_16O2         9.456780e-01
#define ABUN_34S_16O2         4.195030e-02

/*
 * oneline is a test molecule having a single line at 500 GHz, and
 * a very small mass such that the Doppler- and pressure-broadened
 * widths at 500 mbar and 296 K are equal.
 */
#define MASS_ONELINE          3.796124e-5

#define ABUN_ONELINE          1.0


#endif /* AM_MOLECULES_H */
