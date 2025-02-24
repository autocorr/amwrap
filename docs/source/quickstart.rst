Getting started
===============
Atmospheric models can be configured directly using the :class:`amwrap.Model`
class or from standard climatologies using the
:class:`amwrap.Model.from_climatology` constructor, which takes the name of the
climatology to be used or an instance of the :class:`amwrap.Climatology` class.
The standard climatologies for the continental United States, i.e., the "US
Standard", "US Midlatitude Winter", etc. are included from Anderson et al.
(1986). The results are returned in a :class:`pandas.DataFrame`:

.. code-block:: python

   import amwrap
   m = amwrap.Model.from_climatology("midlatitude_winter")
   df = m.run()
   df.head()

The atmospheric data in a climatology can be accessed directly from a
:class:`amwrap.Climatology` class instance:

.. code-block:: python

   # List available climatologies:
   amwrap.Climatology.names
   cl = amwrap.Climatology("midlatitude_winter")
   dir(cl)      # list attributes
   cl.pressure  # access specific quantities

Models can also be constructed directly from :class:`astropy.units.Quantity`
arrays as well. See the documentation for the :class:`amwrap.Model` class for
details on initialization.

.. code-block:: python

   m = amwrap.Model(cl.pressure, cl.temperature, {"h2o": cl.mixing_ratio["h2o"]})
   df = m.run()  # default is single-threaded with `parallel=False`
   df = m.run(parallel=True)  # use number of threads equal to number of CPUs.

Other arguments that can be set include the zenith viewing angle, frequency
start, stop, and width, the output quantities, and a scaling factor on the
water vapor present in the troposphere.

Output column types for AM include frequency, opacity, transmittance, radiance,
radiance difference (against the radiance of the background, by default the
CMB), brightness temperature in the full Planck definition, Rayleigh-Jeans
brightness temperature, delay, and absorption coeffience. The default outputs
are ``('frequency', 'brightness temperature', 'opacity', 'delay')``. A full
list of the output descriptor names and their units can be found in the
:attr:`amwrap.Model.valid_output_descriptors`.

Note that input arrays *must* have approrpiate units associated with them from
the AstroPy units package.  The following shows an example of using
:class:`astropy.units.Quantity`, for a more complete description of the units
system, see the `units module` documentation. The units can be any particular
dimensionally equivalent value. Volumetric mixing ratios are dimensionless
and expressed using the :class:`astropy.units.dimensionless` quantity.

.. code-block:: python

   import numpy as np
   from astropy import units as u
   cl = amwrap.Climatology("midlatitude_winter")
   temperature = 5 * u.K * np.ones_like(cl.temperature)
   m = amwrap.Model(cl.pressure, temperature, {"h2o": cl.mixing_ratio["h2o"]})
   df = m.run()

.. _units module: https://docs.astropy.org/en/stable/units/index.html

