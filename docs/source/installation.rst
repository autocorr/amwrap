Installation and setup
======================
AM is built and compiled automatically when installing ``amwrap`` and does not
need to be installed beforehand.  Currently only Unix-like operating systems
(i.e., Linux and macOS) are supported.  Building AM depends on GNU Make and a C
compiler, such as GCC. The parallel version of AM requires a C compiler with
OpenMP support.

To install ``amwrap``, run the following from the command line:

.. code-block:: bash

   pip install git+https://github.com/autocorr/amwrap.git

or alternatively:

.. code-block:: bash

   git clone https://github.com/autocorr/amwrap.git
   cd amwrap
   pip install .

If binary executables for AM are available in the user's path, this version
will be used instead of the automatically installed version. To see which
version of AM is being used run:

.. code-block:: python

   import amwrap
   print(amwrap.AM_PARALLEL.exec_name)
   print(amwrap.AM_SERIAL.exec_name)

If the above show ``"am"`` or ``"am-serial"``, the user compiled versions
are being used, otherwise the absolute path to the executable distributed
with ``amwrap`` will be shown.

Please report any installation errors to the GitHub `issues`_ page.

.. _issues: https://github.com/autocorr/amwrap/issues
