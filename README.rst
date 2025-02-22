Python wrapper for AM
=====================
A simple Python wrapper for Scott Paine's AM atmospheric radiative transfer
code.  For details on AM see the official Zenodo links to the `source`_ and
`documentation`_. Further information can be found on Scott Paine's `website`_.

Installing ``amwrap`` will automatically compile AM from source files
distributed with the package.  If copies of AM are found in the user's
``PATH``, then these will be used instead.  Documentation on how to use this
package may be found at https://amwrap.readthedocs.io .

Installation
------------
Currently only Unix-like operating systems (i.e., Linux and macOS) are
supported. Building AM depends on GNU Make and a C compiler, such as GCC. The
parallel version of AM requires a C compiler with OpenMP support.

To install ``amwrap``, run the following from the command line:

.. code-block:: bash

   pip install git+https://github.com/autocorr/amwrap.git

or alternatively:

.. code-block:: bash

   git clone https://github.com/autocorr/amwrap.git
   cd amwrap
   pip install .

References
----------
The climatologies in the ``amwrap/climatology`` directory are taken from
``pyrtlib.climatologies`` module which themselves are taken from Anderson et
al. (1986) "AFGL Atmospheric Constituent Profiles (0-120km)", AFGL-TR-86-0110.

License
-------
This wrapper is authored by Brian Svoboda copyright 2025 and released under the
GNU General Public License Agreement Version 3 (GPLv3). The full text of the
license is supplied in the ``LICENSE`` file included with the software. Portions
of this wrapper are adapted or copied from the ``pyrtlib`` Python library
written by Salvatore Larosa that are themselves licensed under the GPLv3.

AM is authored by Scott Paine of the Smithsonian Astrophysical Observatory.
The AM software is a work of the United States and may be used freely, with
attribution and credit to the Smithsonian Astrophysical Observatory. The
program is intended for educational, scholarly or research purposes.

.. _source: https://zenodo.org/records/13748403
.. _documentation: https://zenodo.org/records/13748391
.. _website: https://lweb.cfa.harvard.edu/~spaine/am/index.html
