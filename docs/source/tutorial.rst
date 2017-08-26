Simulation Structure
------------------------------------------------

Simulations with NumBAT are generally carried out using a python script file.
This file is kept in its own directory which is placed in the NumBAT directory.
All results of the simulation are automatically created within this directory. This directory then serves as a complete record of the calculation. Often, we will also save the simulation objects (scattering matrices, propagation constants etc.) within this folder for future inspection, manipulation, plotting, etc.

Traditionally the name of the python script file begins with simo\-. This is convenient for setting terminal alias' for running the script.
Throughout the tutorial the script file will be called simo.py.

To start a simulation open a terminal and change into the directory containing the simo.py file.
To run this script::

    $ python3 simo.py

To have direct access to the simulation objects upon the completion of the script use,::

    $ python3 -i simo.py

This will return you into an interactive python session in which all simulation objects are accessible. In this session you can access the docstrings of objects, classes and methods. For example::

    >>> from pydoc import help
    >>> help(objects.Struct)

where we have accessed the docstring of the Struct class from objects.py

To review how material types and physical dimensions are represented in the mesh geometries go to NumBAT/docs/msh_type_lib and view the relevant .png file.

The parameters lc_bkg, lc2, lc3 set the fineness of the FEM mesh. lc_bkg sets the reference background mesh size, larger lc_bkg = larger (more coarse) mesh. In NumBAT the x-dimension of the unit cell is traditionally normalised to unity, in which case there will be lc_bkg mesh elements along the horizontal outside edge; in other words the outside edge is divided into lc_bkg elements. At the interface between materials the mesh is refined to be lc_bkg/lc2, therefore larger lc2 = finer mesh at these interfaces. The meshing program automatically adjusts the mesh size to smoothly transition from a point that has one mesh parameter to points that have other meshing parameters. The mesh it typically also refined at the centers of important regions, eg in the center of a waveguide, which is done with lc3, which just like lc2 refines the mesh size at these points as lc_bkg/lc3.

Choosing appropriate values of lc_bkg, lc2, lc3 is crucial NumBAT to give accurate results. The values depend strongly on the type of structure being studied, and so it is recommended to carry out a convergence test before delving into new structures (see Tutorial 5) starting from similar parameters as used in the tutoarial simulations. You can also visually check the resolution of your mesh by setting check_msh=True when you define your objects.Struct (see Tutorial 1), or by running the following command ::
    
    NumBAT/backend/fortran/msh$ gmsh <msh_name>.msh


In the remainder of this chapter we go through a number of example simo.py files. But before we do, another quick tip about running simulations within screen sessions, which allow you to disconnect from servers leaving them to continue your processes.

.. raw:: latex

    \clearpage

Screen Sessions
------------------------------------------------
::

    screen

is an extremely useful little linux command. In the context of long-ish calculations it has two important applications; ensuring your calculation is unaffected if your connection to a remote machine breaks, and terminating calculations that have hung without closing the terminal.
For more information see the manual::

    $ man screen

or see online discussions `here <http://www.howtoforge.com/linux_screen>`_, `and here <http://www.rackaid.com/blog/linux-screen-tutorial-and-how-to/>`_.


The screen session or also called screen instance looks just like your regular terminal/putty, but you can disconnect from it (close putty, turn off your computer etc.) and later reconnect to the screen session and everything inside of this will have kept running. You can also reconnect to the session from a different computer via ssh.

Basic Usage
,,,,,,,,,,,,,,,,,,,,,

To install screen::

    $ sudo apt-get install screen

To open a new screen session::

    $ screen

We can start a new calculation here::

    $ cd NumBAT/tutorials/
    $ python simo-tut_01-first_calc.py

We can then detach from the session (leaving everything in the screen running) by typing::

    Ctrl +a
    Ctrl +d

We can now monitor the processes in that session::

    $ top

Where we note the numerous running python processes that NumBAT has started. Watching the number of processes is useful for checking if a long simulation is near completion (which is indicated by the number of processes dropping to less than the specified num_cores).

We could now start another screen and run some more calculations in this terminal (or do anything else).
If we want to access the first session we 'reattach' by typing::

    Ctrl +a +r

Or entering the following into the terminal::

    $ screen -r

If there are multiple sessions use::

    $ screen -ls

to get a listing of the sessions and their ID numbers. To reattach to a particular screen, with ID 1221::

    $ screen -r 1221

To terminate a screen from within type::

    Ctrl+d

Or, taking the session ID from the previous example::

    screen -X -S 1221 kill



Terminating NumBAT simos
,,,,,,,,,,,,,,,,,,,,,,,,

If a simulation hangs, we can kill all python instances upon the machine::

    $ pkill python3

If a calculation hangs from within a screen session one must first detach from that session then kill python, or if it affects multiple instances, you can kill screen. A more targeted way to kill processes is using their PID::

    $ kill PID

Or if this does not suffice be a little more forceful::

    $ kill -9 PID

The PID is found from one of two ways::

    $ top
    $ ps -fe | grep username


.. raw:: latex

    \clearpage


Tutorial
--------

In this section we go through a number of simple simulations that demonstrate the basic use of NumBAT.


Basic SBS Gain Calculation
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
This example, contained in ``tutorials/simo-tut_01-first_calc.py`` calculates the backward SBS gain for a rectangular silicon waveguide surrounded by air.

The sequence of operations is

  #. Import NumBAT modules
  #. Define the structure shape and dimensions
  #. Specify the electromagnetic and acoustic modes to be solved for
  #. Construct the waveguide with ``objects.Struct``
  #. Solve the electromagnetic problem. ``mode_calcs.calc_EM_modes`` returns an object containing modes and their propagation constants as ``Eig_values`` in 1/m.
  #. Convert the EM eigenvalue to an effective index
  #. Identify the desired acoustic wavenumber and solve the acoustic problem. ``mode_calcs.calc_AC_modes`` returns an object containing the modes for propagation constant ``k_AC`` and acoustic frequencies as ``Eig_values`` in Hz.
  #. Calculate the total SBS gain, contributions from photoelasticity and moving boundary effects, and the acoustic loss


.. literalinclude:: ../../tutorials/simo-tut_01-first_calc.py
    :lines: 0-


SBS Gain Spectra
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
This example, contained in ``tutorials/simo-tut_02-gain_spectra-npsave.py`` considers the same structure
but adds plotting of fields, gain spectra and techniques for saving and reusing data from earlier
calculations. 

Elements to note:
  #. Both electric and magnetic fields can be selected using ``'EM_E'`` or ``'EM_H'`` as the value of ``EM_AC` in 
       ``plotting.mode_fields``.
  #. ``np.savez`` and ``np.load`` allow storage of arbitrary data between simulations.

.. literalinclude:: ../../tutorials/simo-tut_02-gain_spectra-npsave.py
    :lines: 0-


.. figure:: ../../tutorials/tut_02-gain_spectra-MB_PE_comps.png
   :scale: 70 %
   
   Gain spectra showing gain due to electrostriction, gain due to radiation pressure, and total gain.

.. figure:: ../../tutorials/tut_02-gain_spectra-MB_PE_comps_zoom.png
   :scale: 70 %
   
   Zoomed in gain spectra.


Investigating Dispersion
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

.. literalinclude:: ../../tutorials/simo-tut_03_1-dispersion-npload.py
    :lines: 0-

Investigating Dispersion
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

.. literalinclude:: ../../tutorials/simo-tut_03_2-dispersion-multicore.py
    :lines: 0-


Parameter Scan of Widths
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

.. literalinclude:: ../../tutorials/simo-tut_04-scan_widths.py
    :lines: 0-


Convergence Study
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

.. literalinclude:: ../../tutorials/simo-tut_05-convergence_study.py
    :lines: 0-


Silica Nanowire 
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

.. literalinclude:: ../../tutorials/simo-tut_06-silica_nanowire.py
    :lines: 0-


Slot Waveguide
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

.. literalinclude:: ../../tutorials/simo-tut_07-slot.py
    :lines: 0-


Covered Slot Waveguide
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

.. literalinclude:: ../../tutorials/simo-tut_08-slot_coated-scan.py
    :lines: 0-


Chalcogenide Rib Waveguide
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

.. literalinclude:: ../../tutorials/simo-tut_09-chalc_rib.py
    :lines: 0-

Anisotropic Elastic Materials 
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

.. literalinclude:: ../../tutorials/simo-tut_10-anisotropy.py


.. raw:: latex

    \clearpage



Literature Examples
---------------------

Having gotten familiar with NumBAT, we now set out to replicate a number of examples from the literature.


Rectangular Waveguide - Silica
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

.. literalinclude:: ../../lit_examples/simo-lit_01-Laude-AIPAdv_2013-silica.py
    :lines: 0-


Rectangular Waveguide - Silicon
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

.. literalinclude:: ../../lit_examples/simo-lit_02-Laude-AIPAdv_2013-silicon.py
    :lines: 0-


Tapered Fibre
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

.. literalinclude:: ../../lit_examples/simo-lit_03-Florez-NatComm_2016.py
    :lines: 0-


Tapered Fibre - Scanning Widths
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

.. literalinclude:: ../../lit_examples/simo-lit_04-Beugnot-NatComm_2014.py
    :lines: 0-


Waveguide on a Pedestal - Fields
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

.. literalinclude:: ../../lit_examples/simo-lit_05-Van_Laer-NatPhot_2015.py
    :lines: 0-


FSBS - Waveguide on a Pedestal NJP
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

.. literalinclude:: ../../lit_examples/simo-lit_06-Van_Laer-NJP_2015.py
    :lines: 0-


FSBS - Rib Waveguide
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

.. literalinclude:: ../../lit_examples/simo-lit_07-Rakich_rib.py
    :lines: 0-


FSBS - Coated Rib Waveguide
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

.. literalinclude:: ../../lit_examples/simo-lit_08-Rakich_rib-coated.py
    :lines: 0-

