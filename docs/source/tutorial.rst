.. _chap-tutorial-label:

Simulation Procedure
------------------------------------------------

Simulations with NumBAT are generally carried out using a python script file.
This file is kept in its own directory which is placed in the NumBAT directory.
All results of the simulation are automatically created within this directory. This directory then serves as a complete record of the calculation. Often, we will also save the simulation objects within this folder for future inspection, manipulation, plotting, etc.

Traditionally the name of the python script file begins with simo\-. This is convenient for setting terminal alias' for running the script.
Throughout the tutorial the script file will be called simo.py.

To start a simulation open a terminal and change into the directory containing the ``simo.py`` file.
To run this script::

    $ python3 simo.py

To have direct access to the simulation objects upon the completion of the script use::

    $ python3 -i simo.py

This will return you into an interactive python session in which all simulation objects are accessible. In this session you can access the docstrings of objects, classes and methods. For example::

    >>> from pydoc import help
    >>> help(objects.Struct)

where we have accessed the docstring of the Struct class from ``objects.py``.


Script Structure
----------------------------

As will be seen in the tutorials below, most NumBAT scripts proceed with a standard
structure: 
  - defining materials
  - defining waveguide geometries and associating them with material properties
  - solving electromagnetic and acoustic modes 
  - calculating gain and other derived quantities

The following section provides some information about the pre-defined range of waveguide
structures and the key parameters controlling finite-element meshing.
Information on how to add new structures to NumBAT is provided in :ref:`sec-newmesh-label`.

Waveguide Geometries
----------------------

The following figures give some examples of how material types and physical 
dimensions are represented in the mesh geometries. These can also be found in the directory::

    >>>  NumBAT/docs/msh_type_lib 

as a series of ``.png`` file.

.. figure:: ../msh_type_lib/1.png
   :scale: 30 %

   Rectangular waveguide.

.. figure:: ../msh_type_lib/1_circular.png
   :scale: 15 %

   Elliptical waveguide.

.. figure:: ../msh_type_lib/2.png
   :scale: 30 %

   Coupled rectangular waveguides.

.. figure:: ../msh_type_lib/rib.png
   :scale: 30 %

   A conventional rib waveguide.

.. figure:: ../msh_type_lib/rib_coated.png
   :scale: 30 %

   A coated rib waveguide.

.. figure:: ../msh_type_lib/rib_double_coated.png
   :scale: 30 %

   A rib waveguide on two substrates.

.. figure:: ../msh_type_lib/slot.png
   :scale: 30 %

   A slot waveguide (``material_a`` is low index).

.. figure:: ../msh_type_lib/slot_coated.png
   :scale: 30 %

   A coated slot waveguide (``material_a`` is low index).

.. figure:: ../msh_type_lib/onion.png
   :scale: 30 %

   A concentric layered structure.

.. raw:: latex

    \clearpage



The parameters ``lc_bkg``, ``lc2``, ``lc3``  to be encountered below set the fineness of the FEM mesh. ``lc_bkg`` sets the reference background mesh size, larger ``lc_bkg`` = larger (more coarse) mesh. In NumBAT the x-dimension of the unit cell is traditionally normalised to unity, in which case there will be ``lc_bkg`` mesh elements along the horizontal outside edge; in other words the outside edge is divided into ``lc_bkg`` elements. At the interface between materials the mesh is refined to be ``lc_bkg/lc2``, therefore larger ``lc2`` = finer mesh at these interfaces. The meshing program automatically adjusts the mesh size to smoothly transition from a point that has one mesh parameter to points that have other meshing parameters. The mesh is typically also refined at the centers of important regions, such as in the center of a waveguide, which is done with ``lc3``, which analogously to ``lc2``, refines the mesh size at these points as ``lc_bkg/lc3``. For definition of lc3+ parameters see the particular .geo file.

Choosing appropriate values of ``lc_bkg``, ``lc2``, ``lc3`` is crucial NumBAT to give accurate results. The values depend strongly on the type of structure being studied, and so it is recommended to carry out a convergence test before delving into new structures (see Tutorial 5) starting from similar parameters as used in the tutoarial simulations. You can also visually check the resolution of your mesh by setting ``plt_mesh=True`` or ``check_mesh=True`` when you define your ``objects.Struct`` - the first saves a png of the mesh the second opens mesh in gmsh - (see Tutorial 1), or by running the following command ::
    
    NumBAT/backend/fortran/msh$ gmsh <msh_name>.msh


In the remainder of this chapter we go through a number of example ``simo.py`` files. But before we do, another quick tip about running simulations within screen sessions, which allow you to disconnect from servers leaving them to continue your processes.

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



Terminating NumBAT simulations
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

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

In this section we walk through a number of simple simulations that demonstrate the basic use of NumBAT.
:ref:`sec-literature-label` looks at a number of literature examples taken from many of
the well-known groups in this field.
The full Python interface is documented in :ref:`chap-pythonbackend-label`.



Basic SBS Gain Calculation
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
This example, contained in ``tutorials/simo-tut_01-first_calc.py`` calculates the backward SBS gain for a rectangular silicon waveguide surrounded by air.

The sequence of operations (annotated in the source code below as Step 1, Step 2, etc) is

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

.. raw:: latex

    \clearpage


SBS Gain Spectra
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
This example, contained in ``tutorials/simo-tut_02-gain_spectra-npsave.py`` considers the same structure
but adds plotting of fields, gain spectra and techniques for saving and reusing data from earlier
calculations. 

Elements to note:
  #. Both electric and magnetic fields can be selected using ``EM_E`` or ``EM_H`` as the value of ``EM_AC`` in 
       ``plotting.mode_fields``.
  #. ``np.savez`` and ``np.load`` allow storage of arbitrary data between simulations.

.. literalinclude:: ../../tutorials/simo-tut_02-gain_spectra-npsave.py
    :lines: 0-


The following figures show a selection of electromagnetic and acoustic mode profiles produced
in this example.

.. figure:: ../../tutorials/tut_02-fields/EM_E_field_0.png
   :scale: 40 %
   
   Fundamental optical mode fields.


.. figure:: ../../tutorials/tut_02-fields/AC_field_2.png
   :scale: 40 %
   
   Acoustic mode with high gain due to moving boundary effect.


.. figure:: ../../tutorials/tut_02-fields/AC_field_4.png
   :scale: 40 %
   
   Acoustic mode with high gain due to moving boundary effect.

.. raw:: latex

    \clearpage

This example also generates gain spectra.

.. _fig-gainspec1-label:

.. figure:: ../../tutorials/tut_02-gain_spectra-MB_PE_comps.png
   :scale: 35 %
   
   Gain spectra showing gain due to the photoelastic effect, gain due to moving boundary effect, and the total gain.


.. figure:: ../../tutorials/tut_02-gain_spectra-MB_PE_comps_zoom.png
   :scale: 35 %
   
   Zoomed-in gain spectra from :ref:`fig-gainspec1-label`.

.. raw:: latex

    \clearpage


Investigating Dispersion and np.save/np.load
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

.. literalinclude:: ../../tutorials/simo-tut_03_1-dispersion-npload.py
    :lines: 0-


.. figure:: ../../tutorials/tut_03_1-dispersion_npload_symmetrised.png
   :scale: 70 %
   
   Acoustic dispersion diagram with modes categorised by symmetry.

.. raw:: latex

    \clearpage



Investigating Dispersion and multiprocessing
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

.. literalinclude:: ../../tutorials/simo-tut_03_2-dispersion-multicore.py
    :lines: 0-


.. figure:: ../../tutorials/tut_03_2-dispersion_multicore.png
   :scale: 70 %
   
   Acoustic dispersion diagram ploted as lines.

.. raw:: latex

    \clearpage


Parameter Scan of Widths
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

.. literalinclude:: ../../tutorials/simo-tut_04-scan_widths.py
    :lines: 0-


.. figure:: ../../tutorials/tut_04-gain_spectra-waterfall.png
   :scale: 70 %
   
   Gain spectra as function of waveguide width.

.. raw:: latex

    \clearpage


Convergence Study
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

.. literalinclude:: ../../tutorials/simo-tut_05-convergence_study.py
    :lines: 0-


.. figure:: ../../tutorials/tut_05-convergence-freq_EM.png
   :scale: 50 %
   
   Convergence of optical mode frequencies.


.. figure:: ../../tutorials/tut_05-convergence-freq_AC.png
   :scale: 50 %
   
   Convergence of acoustic mode frequencies.


.. figure:: ../../tutorials/tut_05-convergence-Gain_PE.png
   :scale: 50 %
   
   Convergence of photoelastic gain.


.. figure:: ../../tutorials/tut_05-convergence-Gain_MB.png
   :scale: 50 %
   
   Convergence of moving boundary gain.


.. figure:: ../../tutorials/tut_05-convergence-Gain.png
   :scale: 50 %
   
   Convergence of total gain.

.. raw:: latex

    \clearpage


Silica Nanowire 
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

.. literalinclude:: ../../tutorials/simo-tut_06-silica_nanowire.py
    :lines: 0-


.. figure:: ../../tutorials/tut_06-gain_spectra-MB_PE_comps_SiO2_NW.png
   :scale: 50 %
   
   Gain spectra showing gain due to photoelastic effect, gain due to moving boundary effect, and total gain.

.. raw:: latex

    \clearpage


Slot Waveguide
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

.. literalinclude:: ../../tutorials/simo-tut_07-slot.py
    :lines: 0-


.. figure:: ../../tutorials/tut_07-gain_spectra-MB_PE_comps_slot.png
   :scale: 50 %
   
   Gain spectra showing gain due to photoelastic effect, gain due to moving boundary effect, and total gain.

.. raw:: latex

    \clearpage


Slot Waveguide Scan Covering
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

.. literalinclude:: ../../tutorials/simo-tut_08-slot_coated-scan.py
    :lines: 0-


.. figure:: ../../tutorials/tut_08-freq_changes.png
   :scale: 50 %
   
   Acoustic frequencies as function of covering layer thickness.

.. raw:: latex

    \clearpage


Anisotropic Elastic Materials 
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

.. literalinclude:: ../../tutorials/simo-tut_09-anisotropy.py


.. raw:: latex

    \clearpage


Multilayered 'Onion'
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

.. literalinclude:: ../../tutorials/simo-tut_10-onion.py


.. raw:: latex

    \clearpage




.. _sec-literature-label:

Literature Examples
---------------------

Having become somewhat familiar with NumBAT, we now set out to replicate a number of examples 
from the recent literature.
The examples are presented in chronological order. 
We note the particular importance of examples 5-8 which include experimental and numerical results that are in good agreement.


2013 - Laude - AIP Adv - BSBS - Rectangular Waveguide - Silica
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

This example is based on the calculation of backward SBS
in a small rectangular silica waveguide described in V. Laude and J.-C. Beugnot, 
`Generation of phonons from electrostriction in small-core optical waveguides 
<http://dx.doi.org/10.1063/1.4801936>`_, *AIP Advances* **3**, 042109 (2013).

Observe the use of a material named ``materials.SiO2_2013_Laude`` 
specifically modelled on the parameters in this paper.
This technique allows users to easily compare exactly to other authors
without changing their preferred material values for their own samples and experiments.

.. literalinclude:: ../../lit_examples/simo-lit_01-Laude-AIPAdv_2013-silica.py
    :lines: 0-


.. figure:: ../../lit_examples/lit_01-fields/EM_E_field_0.png
   :scale: 50 %
   
   Fundamental optical mode fields.


.. figure:: ../../lit_examples/lit_01-fields/AC_field_4.png
   :scale: 50 %
   
   High gain acoustic mode, marked as C in paper.


.. figure:: ../../lit_examples/lit_01-fields/AC_field_55.png
   :scale: 50 %
   
   High gain acoustic mode, marked as D in paper.


.. figure:: ../../lit_examples/lit_01-gain_spectra-MB_PE_comps-logy.png
   :scale: 50 %
   
   Gain spectra on semilogy axis.
   

.. figure:: ../../lit_examples/lit_01-gain_spectra-MB_PE_comps_zoom.png
   :scale: 50 %
   
   Gain spectra zoomed in on mode D.

.. raw:: latex

    \clearpage

2013 - Laude - AIP Adv - BSBS - Rectangular Waveguide - Silicon
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

This example again follows the paper of
V. Laude and J.-C. Beugnot, 
`Generation of phonons from electrostriction in small-core optical waveguides 
<http://dx.doi.org/10.1063/1.4801936>`_, *AIP Advances* **3**, 042109 (2013),
but this time looks at the *silicon* waveguide case.

.. literalinclude:: ../../lit_examples/simo-lit_02-Laude-AIPAdv_2013-silicon.py
    :lines: 0-


.. figure:: ../../lit_examples/lit_02-fields/AC_field_4.png
   :scale: 50 %
   
   High gain acoustic mode, marked as G in paper.


.. figure:: ../../lit_examples/lit_02-gain_spectra-MB_PE_comps-logy.png
   :scale: 50 %
   
   Gain spectra on semilogy axis.

.. raw:: latex

    \clearpage

2014 - Beugnot - Nat Comm - BSBS - Tapered Fibre - Scanning Widths
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

This example is based on the calculation of backward SBS
in a micron scale optical fibre described in J.-C. Beugnot *et al.*, 
`Brillouin light scattering from surface acoustic waves in a subwavelength-diameter optical fibre
<http://dx.doi.org/10.1038/ncomms6242>`_, *Nature Communications* **5**, 5242 (2014).

.. literalinclude:: ../../lit_examples/simo-lit_03-Beugnot-NatComm_2014.py
    :lines: 0-


.. figure:: ../../lit_examples/lit_03-gain-width_scan.png
   :scale: 50 %
   
   Full acoustic wave spectrum for silica microwire, as per Fig. 4a in paper.

.. raw:: latex

    \clearpage

2015 - Van Laer - Nat Phot - FSBF - Waveguide on a Pedestal
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
This example is based on the calculation of forward SBS
in a pedestal silicon waveguide described in R. Van Laer *et al.*, 
`Interaction between light and highly confined hypersound in a silicon photonic nanowire 
<http://dx.doi.org/10.1038/ncomms6242>`_, *Nature Photonics* **9**, 199 (2015).

Note that the absence of an absorptive boundary in the acoustic model 
causes a problem where the slab layer significantly distorting acoustic modes.
Adding this feature is a priority for the next release of NumBAT.
The following example shows an approximate way to avoid the problem for now.

.. literalinclude:: ../../lit_examples/simo-lit_04-pillar-Van_Laer-NatPhot_2015.py
    :lines: 0-


.. figure:: ../../lit_examples/lit_04-pillar-fields/EM_E_field_0.png
   :scale: 50 %
   
   Fundamental optical mode fields.


.. figure:: ../../lit_examples/lit_04-pillar-fields/AC_field_38.png
   :scale: 50 %
   
   Dominant high gain acoustic mode.
   Note how the absence of an absorptive boundary on the SiO2 slab causes this layer to significantly distorted the acoustic modes.


We may also choose to study the simplified situation where the pedestal is removed.


.. literalinclude:: ../../lit_examples/simo-lit_04-no_pillar-Van_Laer-NatPhot_2015.py
    :lines: 0-

Which gives good agreement for the gain spectrum.

.. figure:: ../../lit_examples/lit_04-pillar-fields/lit_04-no_pillar-gain_spectra-MB_PE_comps.png
   :scale: 50 %
   
   Gain spectrum for the simplified case of a waveguide surrounded by vacuum.


.. raw:: latex

    \clearpage

2015 - Van Laer - New J Phys - FSBF - Waveguide without Pedestal
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

This example continues  the study of forward SBS
in a pedestal silicon waveguide described in R. Van Laer *et al.*, 
`Interaction between light and highly confined hypersound in a silicon photonic nanowire 
<http://dx.doi.org/10.1038/ncomms6242>`_, *Nature Photonics* **9**, 199 (2015).

In this case, we simply remove the pedestal and model the main rectangular waveguide.
This makes the acoustic loss calculation incorrect but avoids the problem of acoustic
energy being excessively concentrated in the substrate.

.. literalinclude:: ../../lit_examples/simo-lit_05-Van_Laer-NJP_2015.py
    :lines: 0-


.. figure:: ../../lit_examples/lit_05-fields/EM_E_field_0.png
   :scale: 50 %
   
   Fundamental optical mode fields.


.. figure:: ../../lit_examples/lit_05-fields/AC_field_6.png
   :scale: 50 %
   
   Dominant high gain acoustic mode.


.. raw:: latex

    \clearpage

2016 - Florez - Nat Comm - BSBS - Tapered Fibre - Self Cancel - d = 550 nm
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

This example looks at the phenomenon of Brillouin "self-cancellation" due to 
the electrostrictive and radiation pressure effects acting with opposite sign. 
This was described in O. Florez *et al.*, `Brillouin self-cancellation 
<http://dx.doi.org/10.1038/ncomms11759>`_, *Nature Communications* **7**, 11759 (2016).

.. literalinclude:: ../../lit_examples/simo-lit_06_1-Florez-NatComm_2016-d550nm.py
    :lines: 0-


.. figure:: ../../lit_examples/lit_06_1-fields/AC_field_4.png
   :scale: 50 %
   
   :math:`TR_{21}` acoustic mode fields of a nanowire with diameter 550 nm.


.. figure:: ../../lit_examples/lit_06_1-fields/AC_field_5.png
   :scale: 50 %
   
   :math:`R_{01}` acoustic mode fields of a nanowire with diameter 550 nm.


.. figure:: ../../lit_examples/lit_06_1-gain_spectra-MB_PE_comps.png
   :scale: 50 %
   
   Gain spectra of a nanowire with diameter 550 nm, matching blue curve of Fig. 3b in paper.


.. figure:: ../../lit_examples/lit_06_1-gain_spectra-MB_PE_comps-5.png
   :scale: 50 %

.. figure:: ../../lit_examples/lit_06_1-gain_spectra-MB_PE_comps-6.png
   :scale: 50 %

.. figure:: ../../lit_examples/lit_06_1-gain_spectra-MB_PE_comps-8.png
   :scale: 50 %

.. figure:: ../../lit_examples/lit_06_1-gain_spectra-MB_PE_comps-11.png
   :scale: 50 %
   
   Zoomed in gain spectra around gaint peaks of 550 nm diameter NW.

.. raw:: latex

    \clearpage

2016 - Florez - Nat Comm - BSBS - Tapered Fibre - Self Cancel - d = 1160 nm
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

This example again looks at the paper 
 O. Florez *et al.*, `Brillouin self-cancellation <http://dx.doi.org/10.1038/ncomms11759>`_, *Nature Communications* **7**, 11759 (2016),
but now for a wider core.

.. literalinclude:: ../../lit_examples/simo-lit_06_2-Florez-NatComm_2016-1160nm.py
    :lines: 0-


.. figure:: ../../lit_examples/lit_06_2-gain_spectra-MB_PE_comps.png
   :scale: 50 %
   
   Gain spectra of a nanowire with diameter 1160 nm, as in Fig. 4 of Florez, showing near perfect cancellation at 5.4 GHz.


.. figure:: ../../lit_examples/lit_06_2-gain_spectra-MB_PE_comps-logy.png
   :scale: 50 %
   
   Gain spectra of a nanowire with diameter 1160 nm, as in Fig. 4 of paper, showing near perfect cancellation at 5.4 GHz.


.. raw:: latex

    \clearpage

2016 - Kittlaus - Nat Phot - FSBF - Rib Waveguide
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

This examples explores a first geometry showing large forward SBS in silicon
as described in E. Kittlaus *et al.*, `Large Brillouin amplification in silicon 
<http://dx.doi.org/10.1038/nphoton.2016.112>`_, *Nature Photonics* **10**, 463 (2016).



.. literalinclude:: ../../lit_examples/simo-lit_07-Kittlaus-NatPhot_2016.py
    :lines: 0-


.. figure:: ../../lit_examples/lit_07-fields/EM_E_field_0.png
   :scale: 50 %
   
   Fundamental optical mode fields.


.. figure:: ../../lit_examples/lit_07-fields/AC_field_19.png
   :scale: 50 %
   
   Dominant high gain acoustic mode.


.. figure:: ../../lit_examples/lit_07-gain_spectra-MB_PE_comps.png
   :scale: 50 %
   
   Gain spectra showing gain due to photoelastic effect, gain due to moving boundary effect, and total gain.


.. raw:: latex

    \clearpage

2017 - Kittlaus - Nat Comm - FSBF - Intermode
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

.. literalinclude:: ../../lit_examples/simo-lit_08-Kittlaus-NatComm_2017.py
    :lines: 0-


.. figure:: ../../lit_examples/lit_08-fields/EM_E_field_0.png
   :scale: 50 %
   
   Fundamental (symmetric TE-like) optical mode fields.


.. figure:: ../../lit_examples/lit_08-fields/EM_E_field_1.png
   :scale: 50 %
   
   2nd lowest order (anti-symmetric TE-like) optical mode fields.


.. figure:: ../../lit_examples/lit_08-fields/AC_field_23.png
   :scale: 50 %
   
   Dominant high gain acoustic mode.


.. figure:: ../../lit_examples/lit_08-gain_spectra-MB_PE_comps.png
   :scale: 50 %
   
   Gain spectra showing gain due to photoelastic effect, gain due to moving boundary effect, and total gain.

.. raw:: latex

    \clearpage


2017 - Morrison - Optica - BSBS - Chalcogenide Rib Waveguide
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

.. literalinclude:: ../../lit_examples/simo-lit_09-Morrison-Optica_2017.py
    :lines: 0-


.. figure:: ../../lit_examples/lit_09-gain_spectra-MB_PE_comps.png
   :scale: 50 %
   
   Gain spectra showing gain due to photoelastic effect, gain due to moving boundary effect, and total gain.
