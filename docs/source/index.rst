.. NumBAT documentation master file, created by
   sphinx-quickstart on Wed Oct 19 15:23:46 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to NumBAT's documentation!
==================================

Contents:

.. toctree::
   :maxdepth: 2


==================
Introduction
==================

.. toctree::
    :maxdepth: 4

    intro


==================
Installation
==================

.. toctree::
    :maxdepth: 4

    install

====================================
Basic Usage and Tutorials
====================================

.. toctree::
    :maxdepth: 4

    tutorial


==================
Python Backend
==================

.. toctree::
    :maxdepth: 4

    objects
    materials
    mode_calcs
    integration
    plotting

==================
Fortran Backends
==================

The intention of NumBAT is that the Fortran FEM routines are essentially black boxes. They are called from mode_calcs.py and return the modal fields. However, there are a few important things to know about the workings of these routines.

.. toctree::
    :maxdepth: 4

    fem





Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

