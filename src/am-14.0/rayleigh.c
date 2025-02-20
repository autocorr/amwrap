/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* rayleigh.c                      S. Paine rev. 2023 June 26 
*
* Absorption by particles in the Rayleigh limit.
************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "am_types.h"
#include "errlog.h"
#include "math_const.h"
#include "phys_const.h"
#include "rayleigh.h"

/*
 * Notes on absorption and scattering in the Rayleigh limit:
 *
 * For spherical particles of radius r, define the size parameter
 *
 *      x = 2 pi r / lambda = 2 pi r nu / c
 *
 * This is the optical phase delay corresponding to the particle
 * radius.  The Rayleigh (small particle) approximation holds
 * when the phase shift across the particle is negligible both
 * inside and outside the particle.  That is,
 *
 *      x << 1, |nx| << 1
 *
 * where n is the complex refractive index.  In this limit, the
 * scattering and absorption efficiencies are [*]
 *
 *      Qsca = (8/3) x^4 |(eps - 1) / (eps + 2)|^2 
 *
 * and
 *
 *      Qabs = 4x Im[(eps - 1) / (eps + 2)]
 *           = 12x epsi / |eps + 2|^2
 *           = 12x epsi / [(epsr+2)^2 + epsi^2]
 *
 * where eps = epsr + i epsi is the complex permittivity of the 
 * particle.
 *
 * (For the important case of water droplets in a non-precipitating
 * cloud, a typical particle radius is ~10 microns.  For water,
 * |n| ~ 10 at microwave frequencies, but falls to |n| ~ 2 by 1 THz.
 * The Rayleigh approximation will hold to accuracy better than 1% for
 * frequencies up to about 300 GHz, with Qsca << Qabs over this range.
 * See J. C. Ku and J. D. Felske 1984, JQSRT 31:569.)
 *
 * The absorption cross-section per particle is
 *
 *      Cabs = pi r^2 Qabs
 *           = 18 pi vp (nu / c) {epsi / [(epsr+2)^2 + epsi^2]}
 * 
 * where vp = (4/3) pi r^3 is the volume per particle.  For a column
 * density np particles per unit area, the optical depth is
 *
 *      tau = np Cabs
 *
 * In am, for consistency with molecular absorption coefficients,
 * particle column density is represented internally in the same
 * units [molecule cm^-2] that are used for gases.  The molecular
 * column density is
 *
 *      N = rho np vp
 *
 * where rho is the molecular density per unit volume.  Then
 *
 *      tau = (18 pi N / rho) (nu / c) {epsi / [(epsr+2)^2 + epsi^2]}
 * 
 * and the molecular absorption coefficient is
 *
 *      k = tau / N
 *        = (18 pi / rho) (nu / c) {epsi / [(epsr+2)^2 + epsi^2]}
 * 
 * [*] See, for example
 *
 *  M.I. Mischenko, L. D. Travis, and A. A. Lacis 2002,
 *  "Scattering, Absorption, and Emission of Light by Small Particles."
 *  Cambridge U. Press and NASA Goddard Institute for Space Studies, NY.
 *
 * The original (Cambridge U. Press) edition is out of print, but an
 * updated electronic version is maintained at
 *
 *  http://www.giss.nasa.gov/staff/mmishchenko/books.html
 */

/***********************************************************
* int Rayleigh_mol_abscoeff(
*   double *k,
*   const double *f,
*   const gridsize_t ngrid,
*   const double *epsr,
*   const double *epsi,
*   const double rho)
*
* Purpose:
*   Computes the molecular absorption coefficient [cm^2]
*   for particles in the Rayleigh limit, given the complex
*   permittivity and the molecular density of the particle
*   material.
*
* Arguments:
*   double *k - absorption coefficient [cm^2]
*   const double *f - frequency grid [GHz]
*   const gridsize_t ngrid - number of grid points
*   const double *epsr - real part of permittivity
*   const double *epsi - imaginary part of permittivity
*   const double rho - particle molecular density [cm^-3]
*
* Return:
*   0 if OK, 1 if an error occurred
************************************************************/

int Rayleigh_mol_abscoeff(
    double *k,
    const double *f,
    const gridsize_t ngrid,
    const double *epsr,
    const double *epsi,
    const double rho)
{
    gridsize_t i;
    double r1;

    if (rho <= 0.0) {
        errlog(144, 0);
        return 1;
    }
    r1 = 18.0 * PI / (rho * C_CM_GHZ);
    for (i = 0; i < ngrid; ++i) {
        double r2 = epsr[i] + 2.0;
        k[i] = r1 * f[i] * epsi[i] / (r2 * r2 + epsi[i] * epsi[i]);
    }
    return 0;
}   /* Rayleigh_mol_abscoeff() */
