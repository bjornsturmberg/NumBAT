.. _chap-install-label:

Installation
================

The source code for NumBAT is hosted `here on Github <https://github.com/bjornsturmberg/NumBAT>`_. Please download the latest release from here.

NumBAT has been developed on Ubuntu 16.04 with the following package versions: Python 3.5.3, Numpy 1.11.0, Suitesparse 4.4.6, and Gmsh 2.10.1.
It has also been successfully installed by users on Debian, RedHat and on Windows 10 (installing Ubuntu after enabling the Windows Subystem for Linux - steps 3 here https://msdn.microsoft.com/en-au/commandline/wsl/install_guide) and with different versions of packages, but these installations have not been as thoroughly documented so may require user testing.

In general, you can simply run the setup script ::

    $ ./setup.sh

or, depending on your system configuration as ::

    $ sudo ./setup.sh

from the ``NumBAT/`` directory.

Before doing so you may wish to update your system ::

    $ sudo apt-get update
    $ sudo apt-get upgrade

Or, if you prefer to do things manually, this is equivalent to ::

    $ sudo apt-get install -y <dependencies>
    $ cd backend/fortran/
    $ make
    $ cd ../../tests/
    $ nosetests3

where the <dependencies> packages are listed dependencies.txt. Note that it is safer to pip install matplotlib than apt-get'ing as will install matplotlib 2.0 without conflicting older versions.

For optimal results ::

    $ cp NumBAT/backend/NumBATstyle.mplstyle ~/.config/matplotlib/stylelib/

or replace plt.style.use('NumBATstyle') in NumBAT/backend/plotting.py with your own prefered matplotlib style file.

**This is all there is, there isn't any more.**

Well there's more if you want to change it up.

The Fortran components (NumBAT source code and libraries) have been successfully compiled with intel's ``ifortran`` as well as open-source ``gfortran``. In this documentation we use ``gfortran``, but this can be easily adjusted in ``NumBAT/backend/fortran/Makefile``

On non-ubuntu OSes you may also need to compile a local version of Suitesparse, which is described in the next section.

Manual installation of SuiteSparse
----------------------------------

The FEM routine used in NumBAT makes use of the highly optimised `UMFPACK <https://www.cise.ufl.edu/research/sparse/umfpack/>`_ (Unsymmetric MultiFrontal Package) direct solver for sparse matrices developed by Prof. Timothy A. Davis. This is distributed as part of the  SuiteSparse libraries under a GPL license. It can be downloaded from `https://www.cise.ufl.edu/research/sparse/SuiteSparse/ <https://www.cise.ufl.edu/research/sparse/SuiteSparse/>`_

This is the process we have used in the past, however this was some years ago and may need to be modified.

Unpack SuiteSparse into ``NumBAT/backend/fortran/``, it should create a directory there; ``SuiteSparse/``
Make a directory where you want SuiteSparse installed, in my case SS_installed ::

    $ mkdir SS_installed/

edit SuiteSparse/SuiteSparse\_config/SuiteSparse\_config.mk for consistency across the whole build; i.e. if using intel fortran compiler ::

    line 75 F77 = gfortran --> ifort

set path to install folder::

    line 85 INSTALL_LIB = /$Path_to_EMustack/NumBAT/backend/fortran/SS_installed/lib
    line 86 INSTALL_INCLUDE = /$Path_to_EMustack/NumBAT/backend/fortran/SS_installed/include

line 290ish commenting out all other references to these::

    F77 = ifort
    CC = icc
    BLAS   = -L/apps/intel-ct/12.1.9.293/mkl/lib/intel64 -lmkl_rt
    LAPACK = -L/apps/intel-ct/12.1.9.293/mkl/lib/intel64 -lmkl_rt

Now make new directories for the paths you gave 2 steps back::

    $ mkdir SS_installed/lib SS_installed/include

Download `metis-4.0 <http://glaros.dtc.umn.edu/gkhome/fsroot/sw/metis/OLD>`_ and unpack metis into SuiteSparse/ Now move to the metis directory::

    $ cd SuiteSparse/metis-4.0

Optionally edit ``metis-4.0/Makefile.in`` as per ``SuiteSparse/README.txt`` plus with ``-fPIC``::

    CC = gcc
    or
    CC = icc
    OPTFLAGS = -O3 -fPIC

Now make ``metis`` (still in SuiteSparse/metis-4.0/)::

    $ make

Now move back to ``NumBAT/backend/fortran/`` ::

    $ cp SuiteSparse/metis-4.0/libmetis.a SS_installed/lib/

and then move to ``SuiteSparse/`` and execute the following::

    $ make library
    $ make install
    $ cd SuiteSparse/UMFPACK/Demo
    $ make fortran64
    $ cp SuiteSparse/UMFPACK/Demo/umf4_f77zwrapper64.o into SS_installed/lib/

Copy the libraries into ``NumBAT/backend/fortran/Lib/`` so that ``NumBAT/`` is a complete package that can be moved across machine without alteration. This will override the pre-compiled libraries from the release (you may wish to save these somewhere).::

    $ cp SS_installed/lib/*.a NumBAT/backend/fortran/Lib/
    $ cp SS_installed/lib/umf4_f77zwrapper64.o NumBAT/backend/fortran/Lib/


NumBAT Makefile

Edit ``NumBAT/backend/fortran/Makefile`` to reflect what compiler you are using and how you installed the libraries. The Makefile has further details.

Then finally run the setup.sh script!

.. _sec-contribute-label:

Contributing to NumBAT
----------------------------------

NumBAT is open source software licensed under the GPL with all source and documentation available
at `github.com <https://github.com/bjornsturmberg/NumBAT.git>`_. We welcome additions to NumBAT code, documentation and the materials library. Interested users should fork the standard release from github and make a pull request when ready.  For major changes, we strongly suggest contacting the NumBAT team before starting work at ``michael.steel@mq.edu.au``.
