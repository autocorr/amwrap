/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* errlog.c                   S. Paine rev. 2024 September 26
*
* The global error log.
************************************************************/

#include <stdio.h>

#include "errlog.h"

static int ERRFLAG = 0; /* flag that gets set if any error occurs */

static struct err {
    int num;            /* index number */
    int count;          /* number of times this error occurred */
    int data;           /* optional descriptive data for last occurance */
    const char *text;   /* descriptive text string (may include %d for data) */
} errtab[] = {
    {   0, 0, 0,
        "! Error: Column definition encountered on layer %d with P <= 0."
    },
    {   1, 0, 0,
        "! Error: Column definition encountered on layer %d with T <= 0."
    },
    {   2, 0, 0,
        "! Error: N (scaled) < 0 on layer %d."
    },
    {   3, 0, 0,
        "! Error: Mixing ratio < 0 on layer %d."
    },
    {   4, 0, 0,
        "! Error: T outside tabulated range for h2o partition sum."
    },
    {   5, 0, 0,
        "! Warning: Encountered in-band lines narrower than the frequency\n"
        "!          grid spacing.  The output configuration data includes\n"
        "!          the unresolved line count after each column definition\n"
        "!          for which this occurred."
    },
    {   6, 0, 0,
        "! Error: T outside tabulated range for o2 partition sum."
    },
    {   7, 0, 0,
        "! Error: T outside tabulated range for ch4 partition sum."
    },
    {   8, 0, 0,
        "! Error: T outside tabulated range for o3 partition sum."
    },
    {   9, 0, 0,
        "! Error: T outside tabulated range for co partition sum."
    },
    {  10, 0, 0,
        "! Error: T outside tabulated range for n2o partition sum."
    },
    {  11, 0, 0,
        "! Error: T outside tabulated range for oneline partition sum."
    },
    {  12, 0, 0,
        "! Error: Unrecognized unit in model configuration file."
    },
    {  13, 0, 0,
        "! Error: Allocation of memory for new layer structure failed."
    },
    {  14, 0, 0,
        "! Error: Allocation of memory for layer %d opacity failed."
    },
    {  15, 0, 0,
        "! Error: Allocation of memory for new column structure failed."
    },
    {  16, 0, 0,
        "! Error: Allocation of memory for column absorption coefficient\n"
        "!        failed."
    },
    {  17, 0, 0,
        "! Error: Allocation of memory for column opacity failed."
    },
    {  18, 0, 0,
        "! Error: Allocation of memory for simplex failed at %d dimensions."
    },
    {  19, 0, 0,
        "! Error: Could not open model configuration file."
    },
    {  20, 0, 0,
        "! Error: parse error."
    },
    {  21, 0, 0,
        "! Error: Allocation of memory for frequency grid failed."
    },
    {  22, 0, 0,
        "! Warning: An error occurred while writing a temporary cache file."
    },
    {  23, 0, 0,
        "! Warning: Unable to create a file in the cache directory named\n"
        "!          in the environment."
    },
    {  24, 0, 0,
        "! Error: Unrecognized absorption coefficient type (%d) in\n"
        "!        compute_absorption_coefficient()."
    },
    {  25, 0, 0,
        "! Error: Unrecognized mapping (%d) in map_differential()"}, 
    {  26, 0, 0,
        "! Error: Simplex variable index out of range in\n"
        "!        simplex_variable_range()."
    },
    {  27, 0, 0,
        "! Error: Unrecognized ILS type number (%d) in ilsfunc()."
    },
    {  28, 0, 0,
        "! Error: Allocation of memory for tx (transmittance) array failed."
    },
    {  29, 0, 0,
        "! Error: Allocation of memory for I (radiance) array failed."
    },
    {  30, 0, 0,
        "! Error: Allocation of memory failed in initialize_ils()."
    },
    {  31, 0, 0,
        "! Error: The environment string for the fit data directory path\n"
        "!        is too long."
    },
    {  32, 0, 0,
        "! Warning: T outside fitted range for EBC lineshape in cia.c"
    },
    {  33, 0, 0,
        "! Error: Allocation of memory for fit data file table failed."
    },
    {  34, 0, 0,
        "! Warning: No fit variables were specified."
    },
    {  35, 0, 0,
        "! Error: Model structure dimensions do not match in\n"
        "!        copy_model_scalars()."
    },
    {  36, 0, 0,
        "! Error: Could not open fit data file for reading."
    },
    {  37, 0, 0,
        "! Error: Memory reallocation failure in read_fit_data_file()."
    },
    {  38, 0, 0,
        "! Error: The environment string for the fit output directory path\n"
        "!        is too long."
    },
    {  39, 0, 0,
        "! Error: Undefined mode in set_active_outputs()"
    },
    {  40, 0, 0,
        "! Error: The frequency grid, extended to f = 0 for delay spectrum\n"
        "!        computation, is too large."
    },
    {  41, 0, 0,
        "! Error: Too few fit data points."
    },
    {  42, 0, 0,
        "! Warning: Negative channel bandwidth in fit data; absolute value\n"
        "!          assumed."
    },
    {  43, 0, 0,
        "! Warning: The fit output directory named in the environment may\n"
        "!          be invalid."
    },
    {  44, 0, 0,
        "! Error: Reallocation of memory for layer pointer table failed."
    },
    {  45, 0, 0,
        "! Error: Reallocation of memory for column pointer table failed."
    },
    {  46, 0, 0,
        "! Error: Fit variables or Jacobian differentiation variables must\n"
        "!        be positive and nonzero when log mapping of simplex\n"
        "!        coordinates is in effect (this is the default).  Use the\n"
        "!        configuration file statement\n"
        "!\n"
        "!          simplex_log 0\n"
        "!\n"
        "!        to turn off log mapping in simplex space."
    },
    {  47, 0, 0,
        "! Error: ilsworkspace == NULL in ils_convolve()"
    },
    {  48, 0, 0,
        "! Error: Invalid estimator type in update_simplex_vertex()."
    },
    {  49, 0, 0,
        "! Error: Unrecognized column type number (%d) in column_opacity()."
    },
    {  50, 0, 0,
        "! Error: Attempted to reallocate a kcache table which was already\n"
        "!        allocated."
    },
    {  51, 0, 0,
        "! Error: Allocation of memory failed in init_column_kcache()."
    },
    {  52, 0, 0,
        "! Error: init_kcache() failed in init_fit_kcache() on layer %d."
    },
    {  53, 0, 0,
        "! Error: Allocation of memory for Trj (Rayleigh-Jeans brightness\n"
        "!        temperature) array failed."
    },
    {  54, 0, 0,
        "! Error: Allocation of a kcache entry failed."
    },
    {  55, 0, 0,
        "! Warning: Unable to rename file in promote_to_mru()."
    },
    {  56, 0, 0,
        "! Warning: Unable to rename file in insert_as_mru()."
    },
    {  57, 0, 0,
        "! Warning: The h2o continuum is being computed at an effective\n"
        "!          water vapor mixing ratio of 0.0.  The affected column\n"
        "!          is indicated in the output configuration data."
    },
    {  58, 0, 0,
        "! Error: Allocation of memory for fif array failed."
    },
    {  59, 0, 0,
        "! Error: Allocation of memory for line data failed in linesum()."
    },
    {  60, 0, 0,
        "! Error: Unrecognized lineshape type (%d) in linesum(),"
    },
    {  61, 0, 0,
        "! Error: The initial simplex scale must be non-zero."
    },
    {  62, 0, 0,
        "! Error: The frequency grid needed for convolution with the ILS\n"
        "!        is too large."
    },
    {  63, 0, 0,
        "! Error: Allocation of memory for tau (opacity) array failed."
    },
    {  64, 0, 0,
        "! Error: Allocation of memory for Tb (Planck brightness\n"
        "!        temperature) array failed."
    },
    {  65, 0, 0,
        "! Error: Unrecognized mapping (%d) in unmap_differential()"}, 
    {  66, 0, 0,
        "! Error: Allocation of memory for L (delay) array failed."
    },
    {  67, 0, 0,
        "! Error: No frequency grid has been defined."
    },
    {  68, 0, 0,
        "! Error: Allocation of memory failed in compute_delay_spectrum()."
    },
    {  69, 0, 0,
        "! Warning: Found a corrupt or binary-incompatable disk cache file."
    },
    {  70, 0, 0,
        "! Error: The environment string for the disk cache directory path\n"
        "!        is too long."
    },
    {  71, 0, 0,
        "! Error: model->L == NULL in compute_delay_spectrum()."
    },
    {  72, 0, 0,
        "! Error: compute_model() failed in compute_jacobians()."
    },
    {  73, 0, 0,
        "! Error: The hash modulus specified in the environment was not\n"
        "!        valid."
    },
    {  74, 0, 0,
        "! Error: The hash modulus specified in the environment cannot\n"
        "!        exceed %d."
    },
    {  75, 0, 0,
        "! Warning: remove_with_retry() failed in dcache.c, line %d"
    },
    {  76, 0, 0,
        "! Error: model == NULL in compute_model()"
    },
    {  77, 0, 0,
        "! Warning: The ILS is too broad relative to the range of the\n"
        "!          frequency grid."
    },
    {  78, 0, 0,
        "! Error: Unrecognized mapping (%d) in compute_jacobians()"}, 
    {  79, 0, 0,
        "! Error: lambda = %d in clebsch_squared()."
    },
    {  80, 0, 0,
        "! Error: negative j value in clebsch_squared()."
    },
    {  81, 0, 0,
        "! Error: zero Doppler linewidth due to T = 0."
    },
    {  82, 0, 0,
        "! Error: zero collisional linewidth due to P = 0."
    },
    {  83, 0, 0,
        "! Warning: The ILS is not resolved by the frequency grid.  It has\n"
        "!          been replaced by an impulse function."
    },
    {  84, 0, 0,
        "! Warning: Encountered a column for which strict self-broadening\n"
        "!          was in effect, and for which there was a problem with\n"
        "!          the volume mixing ratio computation.  The affected\n"
        "!          column is indicated in the output configuration data."
    },
    {  85, 0, 0,
        "! Error: IF spectra could not be computed for the specified LO\n"
        "!        frequency and frequency grid range."
    },
    {  86, 0, 0,
        "! Error: bad IF spectrum mode in model.c, line %d"
    },
    {  87, 0, 0,
        "! Warning: The model frequency grid and output frequency ranges\n"
        "!          do not overlap."
    },
    {  88, 0, 0,
        "! Error: IF spectra cannot be computed for model frequency grids\n"
        "!        with fewer than two points."
    },
    {  89, 0, 0,
        "! Warning: Negative weight factor found in fit data; absolute\n"
        "!          value assumed."
    },
    {  90, 0, 0,
        "! Error: The fit data format string must match between two and\n"
        "!        four data fields."
    },
    {  91, 0, 0,
        "! Error: Invalid fit data format specification.  Valid\n"
        "!        specifications for fields to be read are %%le, %%lf,\n"
        "!        or %%lg, which are all equivalent."
    },
    {  92, 0, 0,
        "! Error: Allocation of memory for I0 (background radiance) array\n"
        "!        failed."
    },
    {  93, 0, 0,
        "! Error: Allocation of memory for layer %d Planck function array\n"
        "!        failed."
    },
    {  94, 0, 0,
        "! Error: Allocation of memory for layer %d transmittance failed."
    },
    {  95, 0, 0,
        "! Error: The simplex vertex table was already allocated in\n"
        "!        create_null_simplex()."
    },
    {  96, 0, 0,
        "! Error: Unrecognized zenith angle dependence (%d) in\n"
        "!        compute_model()."
    },
    {  97, 0, 0,
        "! Error: Allocation of Jacobian array table for output %d failed."
    },
    {  98, 0, 0,
        "! Error: Allocation of Jacobian array for variable %d failed."
    },
    {  99, 0, 0,
        "! Error: Allocation of memory for absorption coefficient pointer\n"
        "!        table failed."
    },
    { 100, 0, 0,
        "! Error: Allocation of memory for new absorption coefficient\n"
        "!        structure failed."
    },
    { 101, 0, 0,
        "! Warning: Encountered a layer for which the sum of the volume\n"
        "!          mixing ratios for all columns exceeds unity.  The\n"
        "!          affected layer is indicated in the output configuration\n"
        "!          data."
    },
    { 102, 0, 0,
        "! Error: The frequency grid lies outside the tabulated frequency\n"
        "!        range in interpolate_continuum_table()."
    },
    { 103, 0, 0,
        "! Error: Not a delay column type in column_optical_delay()."
    },
    { 104, 0, 0,
        "! Error: T outside tabulated range for ocs partition sum."
    },
    { 105, 0, 0,
        "! Error: layer number (%d) out of bounds in\n"
        "!        clear_layer_kcache_entries()."
    },
    { 106, 0, 0,
        "! Error: Encountered a redundant differentiation or fit variable.\n"
        "!        This might have been caused by a duplicate line in the\n"
        "!        configuration file."
    },
    { 107, 0, 0,
        "! Warning: The user-specified characteristic scale for computing\n"
        "!          a Jacobian with respect to mixing ratio exceeded 1.0,\n"
        "!          and has been clamped at 1.0."
    },
    { 108, 0, 0,
        "! Error: Allocation of memory for Tsys array failed."
    },
    { 109, 0, 0,
        "! Error: Allocation of memory for Y array failed."
    },
    { 110, 0, 0,
        "! Error: Denominator (Tref + Trx) too small for Y computation."
    },
    { 111, 0, 0,
        "! Warning: Ignoring model configuration command found after first\n"
        "!          data point in block."
    },
    { 112, 0, 0,
        "! Warning: An automatically-computed column density has been\n"
        "!          multiplied by an Nscale factor."
    },
    { 113, 0, 0,
        "! Warning: The kcache memory limit specified in the environment\n"
        "!          was not valid."
    },
    { 114, 0, 0,
        "! Error: The kcache memory limit specified in the environment is\n"
        "!        too small.  Increase the limit, or use\n"
        "!\n"
        "!          kcache off\n"
        "!\n"
        "!        in the model configuration file to disable the kcache."
    },
    { 115, 0, 0,
        "! Warning: Ignoring fit restart command encountered within a\n"
        "!          data block."
    },
    { 116, 0, 0,
        "! Warning: RH has been used on a layer for which T is outside\n"
        "!          the applicable range (123 K < T < 332 K) of the\n"
        "!          Murphy-Koop formula for saturated vapor pressure over\n"
        "!          liquid water."
    },
    { 117, 0, 0,
        "! Warning: Encountered unexpected column type number %d in\n"
        "!          H2O_liquid_Psat()."
    },
    { 118, 0, 0,
        "! Warning: RH_offset would have produced a negative relative\n"
        "!          humidity.  Instead, the relative humidity was\n"
        "!          clamped at 0%."
    },
    { 119, 0, 0,
        "! Error: Allocation of memory for I_ref array failed."
    },
    { 120, 0, 0,
        "! Error: Allocation of memory for I_diff array failed."
    },
    { 121, 0, 0,
        "! Error: Allocation of memory for variable name string failed."
    },
    { 122, 0, 0,
        "! Warning: A mixing ratio computed from relative humidity was\n"
        "!          clamped at 1.0."
    },
    { 123, 0, 0,
        "! Error: Unrecognized mapping (%d) in map_variable()"}, 
    { 124, 0, 0,
        "! Error: Unrecognized mapping (%d) in unmap_variable()"}, 
    { 125, 0, 0,
        "! Error: collisional linewidth error due to T = 0."
    },
    { 126, 0, 0,
        "! Warning: RHi has been used on a layer for which T is outside\n"
        "!          the applicable range (111 K < T < 273.16 K) of the\n"
        "!          Murphy-Koop forumula for saturated vapor pressure over\n"
        "!          ice."
    },
    { 127, 0, 0,
        "! Warning: Encountered unexpected column type number %d in\n"
        "!          H2O_ice_Psat()."
    },
    { 128, 0, 0,
        "! Error: Invalid diff_type in compute_jacobians()."
    },
    { 129, 0, 0,
        "! Error: allocation of memory for layer lineshape table failed"
    },
    { 130, 0, 0,
        "! Error: allocation of memory for layer strict_selfbroad table\n"
        "!        failed"
    },
    { 131, 0, 0,
        "! Error: allocation of memory for layer Mair_flag table failed"
    },
    { 132, 0, 0,
        "! Error: T outside tabulated range for co2 partition sum."
    },
    { 133, 0, 0,
        "! Warning: The user-specified scale for the differentiation\n"
        "!          variable flo would have produced a differentiation\n"
        "!          step larger than the frequency grid interval df,\n"
        "!          resulting in a non-constant IF grid size.  The\n"
        "!          outermost differentiation step has been clamped at df."
    },
    { 134, 0, 0,
        "! Error: T outside tabulated range for h2o2 partition sum."
    },
    { 135, 0, 0,
        "! Error: T outside tabulated range for ho2 partition sum."
    },
    { 136, 0, 0,
        "! Error: T outside tabulated range for oh partition sum."
    },
    { 137, 0, 0,
        "! Error: T outside tabulated range for o partition sum."
    },
    { 138, 0, 0,
        "! Error: T outside tabulated range for so2 partition sum."
    },
    { 139, 0, 0,
        "! Error: T outside valid range for lwp_abs_Rayleigh."
    },
    { 140, 0, 0,
        "! Warning: lwp_abs_Rayleigh has been found on a layer for which\n"
        "!          T is below the typical lower limit (-35 C) for liquid\n"
        "!          water clouds."
    },
    { 141, 0, 0,
        "! Error: memory allocation failure in lwp_abs_Rayleigh()."
    },
    { 142, 0, 0,
        "! Error: frequency outside tabulated range in\n"
        "!        H2O_ice_permittivity_Warren_Brandt()."
    },
    { 143, 0, 0,
        "! Error: NaN or bad number on line %d of fit data."
    },
    { 144, 0, 0,
        "! Error: T outside tabulated range for hno3 partition sum."
    },
    { 145, 0, 0,
        "! Error: T outside tabulated range for no partition sum."
    },
    { 146, 0, 0,
        "! Error: T outside tabulated range for no2 partition sum."
    },
    { 147, 0, 0,
        "! Error: T outside tabulated range for hbr partition sum."
    },
    { 148, 0, 0,
        "! Error: T outside tabulated range for hcl partition sum."
    },
    { 149, 0, 0,
        "! Error: T outside tabulated range for hf partition sum."
    },
    { 150, 0, 0,
        "! Error: T outside tabulated range for clo partition sum."
    },
    { 151, 0, 0,
        "! Error: T outside tabulated range for hcn partition sum."
    },
    { 152, 0, 0,
        "! Error: T outside tabulated range for hocl partition sum."
    },
    { 153, 0, 0,
        "! Error: T outside tabulated range for h2s partition sum."
    },
    { 154, 0, 0,
        "! Error: T outside tabulated range for ch3oh partition sum."
    },
    { 155, 0, 0,
        "! Error: T outside tabulated range for h2co partition sum."
    },
    { 156, 0, 0,
        "! Error: T outside tabulated range for nh3 partition sum."
    },
    { 157, 0, 0,
        "! Warning: Encounted a layer on which the weighted mean molecular\n"
        "!          mass was undefined.  The default value was used instead.\n"
        "!          The affected layer is indicated in the output\n"
        "!          configuration data."
    },
    { 158, 0, 0,
        "! Error: Water ice was encountered on a layer for which the\n"
        "!        dielectric properties would have been extrapolated to\n"
        "!        an unphysically high temperature (T > 300 K)."
    },
    { 159, 0, 0,
        "! Warning: Water ice was encountered on a layer for which the\n"
        "!          dielectric properties were extrapolated beyond the\n"
        "!          triple point (T > 273.16 K)"
    },
    { 160, 0, 0,
        "! Error: Allocation of memory for scratch space dk failed in\n"
        "!        linesum()."
    },
    { 161, 0, 0,
        "! Error: zenith angle must be less than 90 degrees for plane\n"
        "!        parallel geometry."
    },
    { 162, 0, 0,
        "! Error: bad layer type (%d) in add_layer()."
    },
    { 163, 0, 0,
        "! Error: Attempted to find a level height at P <= 0 in\n"
        "         find_z_from_P()."
    },
    { 164, 0, 0,
        "! Warning: A level height is being extrapolated more than two\n"
        "           pressure scale heights above the highest defined level\n"
        "           in the atmospheric model."
    },
    { 165, 0, 0,
        "! Warning: The observing level has been extrapolated below the\n"
        "!          lowest defined model level.  By default, this is done\n"
        "!          by extending the base layer isothermally.  Setting\n"
        "!          the default mode explicity with the statement\n"
        "!            PTmode extend isothermal\n"
        "!          will suppress this warning, as will setting any other\n"
        "!          extrapolation mode.\n"
    },
    { 166, 0, 0,
        "! Warning: The tangent level has been extrapolated isothermally\n"
        "!          below the lowest defined model level.  Use an explicit\n"
        "!          \"extend\" statement to set the extrapolation mode and\n"
        "!          suppress this warning."
    },
    { 167, 0, 0,
        "! Warning: The propagation path passes through a tangent point\n"
        "!          extrapolated below the planetary radius R0."
    },
    { 168, 0, 0,
        "! Error: Iterative solution for the tangent level failed to\n"
        "!        converge."
    },
    { 169, 0, 0,
        "! Warning: The entire propagation path lies within a model layer\n"
        "!          with an upper boundary at infinity (P = 0).  Note that\n"
        "!          for such a layer, the air mass is asymptotically equal\n"
        "!          to unity."
    },
    { 170, 0, 0,
        "! Error: The -a option requires a model file (or stdin)."
    },
    { 171, 0, 0,
        "! Error: A level height is being extrapolated more than one\n"
        "         pressure scale height below the lowest defined level in\n"
        "         the atmospheric model."
    },
    { 172, 0, 0,
        "! Error: The line of sight is trapped in a duct at layer %d."
    },
    { 173, 0, 0,
        "! Error: unmatched layer dimensions in copy_layer_allocations().\n"
    },
    { 174, 0, 0,
        "! Warning: Refraction is turned off by default.  Use the statements\n"
        "!            refract radio\n"
        "!          or\n"
        "!            refract optical\n"
        "!          to select radio or optical refractivity models.  To\n"
        "!          suppress this warning in spherical models, use\n"
        "!            refract none\n"
        "!          to turn off refraction explicitly."
    },
    { 175, 0, 0,
        "! Warning: R0 has no effect in plane-parallel models."
    },
    { 176, 0, 0,
        "! Warning: The vertical gravity gradient dg_dz was assigned a\n"
        "!          positive value.  The normal sign convention is\n"
        "!          dg/dz < 0."
    },
    { 177, 0, 0,
        "! Warning: The source level has been extrapolated isothermally\n"
        "!          below the lowest defined model level.  Use an explicit\n"
        "!          \"extend\" statement to set the extrapolation mode and\n"
        "!          suppress this warning."
    },
    { 178, 0, 0,
        "! Error: The source level is below the observer level, but the view\n"
        "!        direction is at or above the astronomical horizon."
    },
    { 179, 0, 0,
        "! Error: The source level is below the tangent level."
    },
    { 180, 0, 0,
        "! Error: A sub-horizon source level defined to be on the near\n"
        "!        side of the tangent point must lie between the tangent\n"
        "!        level and the observing level."
    },
    { 181, 0, 0,
        "! Error: memory allocation failed in copy_layer_allocations().\n"
    },
    { 182, 0, 0,
        "! Error: Failed to insert tangent level in solve_for_tan_level().\n"
    },
    { 183, 0, 0,
        "! Warning: Double-sideband opacity spectra are only meaningful in\n"
        "!          the optically-thin limit, or if tau(usb) and tau(lsb)\n"
        "!          are approximately equal."
    },
    { 184, 0, 0,
        "! Warning: None of the specified output spectra were eligible\n"
        "!          for convolution with an ILS."
    },
    { 185, 0, 0,
        "! Warning: Not all of the specified output spectra were eligible\n"
        "!          for convolution with an ILS.  See the output config\n"
        "!          data for details."
    },
    { 186, 0, 0,
        "! Error: bad layer type (%d) in insert_interpolated_level()."
    },
    { 187, 0, 0,
        "! Error: delete_layer() failed in insert_interpolated_level()\n"
        "!        with lnum = %d"
    },
    { 188, 0, 0,
        "! Error: T <= 0 in find_P_from_z().  A possible cause is aa\n"
        "!        model layer having T = 0 with P > 0."
    },
    { 189, 0, 0,
        "! Error: Attempted to extrapolate too far below the lowest model\n"
        "!        layer in find_P_from_z()."
    },
    { 190, 0, 0,
        "! Error: Failed to set observing level."
    },
    { 191, 0, 0,
        "! Error: Failed to set source level."
    },
    { 192, 0, 0,
        "! Error: Failed to set tangent level."
    },
    { 193, 0, 0,
        "! Error: Allocation of memory for tau_fsl array failed."
    },
    { 194, 0, 0,
        "! Warning: Free space path loss is not meaningful for paths\n"
        "!          having one endpoint at infinity."
    },
    { 195, 0, 0,
        "! Warning: tau_fsl was clamped at zero for frequencies f\n"
        "!          for which (4 pi d) < (c / f)"
    },
    { 196, 0, 0, 
        "! Error: The source is on the far side of the tangent point from\n"
        "!        the observer, and the tangent level is more than one\n"
        "!        scale height below the lowest model definition level.\n"
        "!\n"
        "!  Hint: If a source on the near side of the tangent point is\n"
        "!        intended, as for a nadir-viewing satellite, use the\n"
        "!        \"near\" keyword in a Psource or zsource statement."
    },
    { 197, 0, 0,
        "! Error: Allocation of memory for k_out array failed."
    },
    { 198, 0, 0,
        "! Warning: The output spectrum \"k\" was requested, but no spectral\n"
        "!          absorption coefficient was computed for this model."},    
    { 199, 0, 0,
        "! Warning: The output spectrum \"k\" is normally associated with a\n"
        "!          minimal model configuration invoking a single spectral\n"
        "!          absorption coefficient computation.  Here, the output\n"
        "!          is the first of %d spectral absorption coefficients\n"
        "!          that were computed for this model."
    },
    { 200, 0, 0,
        "! Warning: The user-specified units for k in the output statement\n"
        "!          did not match the dimensions of the computed spectral\n"
        "!          absorption coefficient.  They have been replaced with\n"
        "!          correct default units."
    },
    { 201, 0, 0,
        "! Error: Memory allocation failure in interpolate_continuum_table()."
    },
    { 202, 0, 0,
        "! Error: A \"set\" command embedded in the fit data stream\n"
        "         contained an error, and the fit was terminated."
    },
    { 203, 0, 0,
        "! Error: T outside tabulated range for ch3cn partition sum."
    },
    { 204, 0, 0,
        "! Error: Unsupported output format (%d).\n"
    },
    { 205, 0, 0,
        "! Error: Unexpected failure to locate model definition level(s)\n"
        "!        in function find_z_from_P() in model.c line %d."    
    },
    { 206, 0, 0,
        "! Error: Bad PTmode (%d) in layer_base_lapse_rate().\n"
    },
    { 207, 0, 0,
        "! Warning: The midpoint temperature of layer %d, which appears\n"
        "!          to be a plume or instrument layer, is affected by the\n"
        "!          base temperature of the layer immediately above.\n"
        "!          If needed, the plume temperature can be decoupled from\n"
        "!          the atmosphere by adding adjacent empty layers having\n"
        "!          zero pressure thickness and suitably-chosen base\n"
        "!          temperatures."
    },
    { 208, 0, 0,
        "! Warning: The base temperature of layer %d, which appears to be\n"
        "!          a plume or instrument layer, is affecting the midpoint\n"
        "!          temperature of the layer immediately below.\n"
        "!          If needed, the plume temperature can be decoupled from\n"
        "!          the atmosphere by adding adjacent empty layers having\n"
        "!          zero pressure thickness and suitably-chosen base\n"
        "!          temperatures."
    },
    { 209, 0, 0,
        "! Warning: For models with layers defined by midpoint temperatures,\n"
        "!          base extrapolation is always isothermal.  A PTmode\n"
        "!          statement other than\n"
        "!            PTmode extend isothermal\n"
        "!          will have no effect."
    },
    { 210, 0, 0,
        "! Warning: Failed to set output stream to binary mode in\n"
        "!          write_model_spectra_as_npy()."
    },
};


/***********************************************************
* int errlog(const int errnum, const int data)
*
* Purpose:
*   Increments the error count for error type errnum, and
*   sets errtab[errnum].data equal to data.  The global
*   error flag is set.
*
* Arguments:
*   int errnum - error number in the error table errtab[]
*   int data   - integer data for optional data field
*
* Return:
*   new value of errtab[errnum].count
************************************************************/

int errlog(const int errnum, const int data)
{
    ERRFLAG = 1;
    ++errtab[errnum].count;
    errtab[errnum].data = data;
    return errtab[errnum].count;
}   /* errlog() */


/***********************************************************
* int errstat(void)
*
* Purpose:
*   Reports error log status
*
* Return:
*   1 if any errors or warnings have been logged
*   0 otherwise
************************************************************/

int errstat(void)
{
    return ERRFLAG;
}   /* errstat() */


/***********************************************************
* int errtest(int errnum)
*
* Purpose:
*   Tests whether error number errnum has occurred
*
* Return:
    Current value of errtab[errnum].count
************************************************************/

int errtest(const int errnum)
{
    return errtab[errnum].count;
}   /* errtest() */


/***********************************************************
* int print_errlog(void)
*
* Purpose:
*   Prints the text string and count for every error in the
*   error log table having a non-zero count.
*
* Return:
*   0 if there were no errors, 1 if there were.
************************************************************/

int print_errlog(void)
{
    unsigned int i;

    if (!ERRFLAG)
        return 0;
    fprintf(stderr, "\n");
    for (i = 0; i < sizeof(errtab) / sizeof(struct err); ++i) {
        if (errtab[i].count > 0) {
            fprintf(stderr, errtab[i].text, errtab[i].data);
            if (errtab[i].count > 1)
                fprintf(stderr, "  Count: %d", errtab[i].count);
            fprintf(stderr, "\n");
        }
    }
    fprintf(stderr, "\n");
    fflush(stderr);
    return 1;
}   /* print_errlog() */
