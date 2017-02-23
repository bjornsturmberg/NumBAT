Installation
================

The source code for NumBAT is hosted `here on Github <https://github.com/bjornsturmberg/NumBAT>`_. Please download the latest release from here.

NumBAT has been developed on Ubuntu and is easiest to install on this platform. Simply run the setup script ::

    $ sudo /setup.sh

Or, if you prefer to do things manually, this is equivalent to ::

    $ sudo apt-get update
    $ sudo apt-get upgrade
    $ sudo apt-get install -y <dependencies>
    $ cd backend/fortran/
    $ make
    $ cd ../../tests/
    $ nosetests

where the <dependencies> packages are listed dependencies.txt.

**This is all there is, there isn't any more.**

Well there's more if you want to change it up.

The Fortran components (NumBAT source code and libraries) have been successfully compiled with intel's ifortran as well as open-source gfortran. In this documentation we use gfortran, but this can be easily adjusted in NumBAT/backend/fortran/Makefile

On non-ubuntu OS you may also need to compile a local version of Suitesparse, which is described in the next section.

Manual installation of SuiteSparse
----------------------------------

The FEM routine used in NumBAT makes use of the highly optimised `UMFPACK <https://www.cise.ufl.edu/research/sparse/umfpack/>`_ (Unsymmetric MultiFrontal Package) direct solver for sparse matrices developed by Prof. Timothy A. Davis. This is distributed as part of the  SuiteSparse libraries under a GPL license. It can be downloaded from `https://www.cise.ufl.edu/research/sparse/SuiteSparse/ <https://www.cise.ufl.edu/research/sparse/SuiteSparse/>`_

This is the process I followed in my installations, however this was some years ago and may need to be modified.

Unpack SuiteSparse into NumBAT/backend/fortran/, it should create a directory there; SuiteSparse/
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

Optionally edit metis-4.0/Makefile.in as per SuiteSparse/README.txt plus with -fPIC::

    CC = gcc
    or
    CC = icc
    OPTFLAGS = -O3 -fPIC

Now make metis (still in SuiteSparse/metis-4.0/)::

    $ make

Now move back to NumBAT/backend/fortran/ ::

    $ cp SuiteSparse/metis-4.0/libmetis.a SS_installed/lib/

and then move to SuiteSparse/ and execute the following::

    $ make library
    $ make install
    $ cd SuiteSparse/UMFPACK/Demo
    $ make fortran64
    $ cp SuiteSparse/UMFPACK/Demo/umf4_f77zwrapper64.o into SS_installed/lib/

Copy the libraries into NumBAT/backend/fortran/Lib/ so that NumBAT/ is a complete package that can be moved across machine without alteration. This will override the pre-compiled libraries from the release (you may wish to save these somewhere).::

    $ cp SS_installed/lib/*.a NumBAT/backend/fortran/Lib/
    $ cp SS_installed/lib/umf4_f77zwrapper64.o NumBAT/backend/fortran/Lib/


NumBAT Makefile

Edit NumBAT/backend/fortran/Makefile to reflect what compiler you are using and how you installed the libraries. The Makefile has further details.

Then finally run the setup.sh script!
