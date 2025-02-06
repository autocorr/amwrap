Python wrapper for AM
=====================
A simple Python wrapper for Scott Paine's ``am`` atmospheric radiative transfer
code.  For details on ``am`` see the official Zenodo links to the `source`_ and
`documentation`_. Further information can be found on Scott Paine's `website`_.
This wrapper requires the ``am`` executable to be available on the host
operating system and in the user's path.  Documentation on how to use this
wrapper may be found at https://amwrap.readthedocs.io .

References
----------
The climatologies in the ``amwrap.data`` directory are taken from
``pyrtlib.climatologies`` module which themselves are taken from Anderson et
al. (1986) "AFGL Atmospheric Constituent Profiles (0-120km)", AFGL-TR-86-0110.

License
-------
This wrapper is authored by Brian Svoboda copyright 2025 and released under the
GNU General Public License Agreement Version 3 (GPLv3). The full text of the
license is supplied in the ``LICENSE`` file included with the software. Portions
of this wrapper are adapted or copied from the ``pyrtlib`` Python library
written by Salvatore Larosa that are themselves licensed under the GPLv3.

The ``am`` code is a work of the United States and may be used freely, with
attribution and credit to the Smithsonian Astrophysical Observatory. The
program is intended for educational, scholarly or research purposes.

.. _source: https://zenodo.org/records/13748403
.. _documentation: https://zenodo.org/records/13748391
.. _website: https://lweb.cfa.harvard.edu/~spaine/am/index.html
