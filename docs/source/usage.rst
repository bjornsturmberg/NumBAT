.. _chap-usage-label:

Simulation Procedure
------------------------------------------------

Simulations with NumBAT are generally carried out using a python script file.
This file is kept in its own directory which is placed in the NumBAT directory.
All results of the simulation are automatically created within this directory. This directory then serves as a complete record of the calculation. Often, we will also save the simulation objects within this folder for future inspection, manipulation, plotting, etc.

Throughout the tutorial the script file will be called simo.py.

These files can be edited using your choice of text editor (for instance running the following in the terminal ``$ nano simo.py``) or an IDE (for instance pycharm) which allow you to run and debug code within the IDE.

To start a simulation open a terminal and change into the directory containing the ``simo.py`` file.

To start we run an example simulation from the tutorials directory. To move to this directory in the terminal enter::

    $ cd <path to installation>/NumBAT/tutorials

To run this script execute::

    $ python3 simo.py

To save the results from the simulation that are displayed upon execution (the print statements in simo.py) use::

    $ python3 ./simo.py | tee log-simo.log

This may require you to update the permissions for the simo.py file to make it executable. This is done in the terminal as::

    $ chmod +x simo.py

To have direct access to the simulation objects upon the completion of the script use::

    $ python3 -i simo.py

This will execute the simo.py script and then return you into an interactive python session within the terminal. This terminal session provides the user experience of an ipython type shell where the python environment and all the simulation objects are as in the simo.py script. In this session you can access the docstrings of objects, classes and methods. For example::

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

The following section provides some information about specifying material properties and waveguide
structures, as well as the key parameters for controlling the finite-element meshing.
Information on how to add new structures to NumBAT is provided in :ref:`sec-newmesh-label`.


Materials
----------------------

In order to calculate the modes of a structure we must specify the acoustic and optical properties of all constituent materials.

In NumBAT, this data is read in from json files, which are stored in /NumBAT/backend/material_data

These files not only provide the numerical values for optical and acoustic variables, but record how these variables have been arrived at. Often they are taken from the literature.

The intention of this arrangement is to create a library of materials that can we hope can form a standard amongst the research community. 
They also allow users to check the sensitivity of their results on particular parameters for a given material.

At present, the material library contains:
  - Vacuum
  - As2S3_2016_Smith
  - As2S3_2017_Morrison
  - GaAs_2016_Smith
  - Si_2013_Laude
  - Si_2015_Van_Laer
  - Si_2016_Smith
  - SiO2_2013_Laude
  - SiO2_2015_Van_Laer
  - SiO2_2016_Smith
  - Si_test_anisotropic

All available materials are loaded into NumBAT into the materials.materials_dict dictionary, 
whose keys are the json file names. 
Materials can easily be added to this by copying any of these files as a template and 
modifying the properties to suit. The Si_test_anisotropic file contains all the variables
that NumBAT is setup to read. We ask that stable parameters (particularly those used
for published results) be added to the NumBAT git repository using the same naming convention.


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



The parameters ``lc_bkg``, ``lc_refine_1``, ``lc_refine_2``  to be encountered below set the fineness of the FEM mesh. ``lc_bkg`` sets the reference background mesh size, larger ``lc_bkg`` = larger (more coarse) mesh. In NumBAT it is also possible to refine the mesh near interfaces and near select points in the domain, as highlighted in the figures above. This is done using the ``lc_refine_`` commands, which we now discuss. At the interface between materials the mesh is refined to be ``lc_bkg/lc_refine_1``, therefore larger ``lc_refine_1`` = finer mesh at these interfaces. The meshing program automatically adjusts the mesh size to smoothly transition from a point that has one mesh parameter to points that have other meshing parameters. The mesh is typically also refined at the centers of important regions, such as in the center of a waveguide, which is done with ``lc_refine_2``, which analogously to ``lc_refine_1``, refines the mesh size at these points as ``lc_bkg/lc_refine_2``. For definition of ``lc_refine_3+`` parameters see the particular .geo file.

Choosing appropriate values of ``lc_bkg``, ``lc_refine_1``, ``lc_refine_2`` is crucial NumBAT to give accurate results. The values depend strongly on the type of structure being studied, and so it is recommended to carry out a convergence test before delving into new structures (see Tutorial 5) starting from similar parameters as used in the tutorial simulations. In NumBAT the x-dimension of the unit cell is traditionally normalised to unity, in which case there will be ``lc_bkg`` mesh elements along the horizontal outside edge; in other words the outside edge is divided into ``lc_bkg`` elements. 

You can also visually check the resolution of your mesh by setting ``plt_mesh=True`` or ``check_mesh=True`` when you define your ``objects.Struct`` - the first saves a png of the mesh (in NumBAT/backend/fortran/msh/) the second opens mesh in gmsh - (see Tutorial 1). The NumBAT generated .msh file is stored in NumBAT/backend/fortran/msh/ which can be viewed by running the following command ::
    
    NumBAT/backend/fortran/msh$ gmsh <msh_name>.msh

Users on WSL will need to first run an X listener (such as XMING) in Windows in order for the "plt_mesh=True" feature to work.
Once the X listener is running, execute the following in the terminal::

    $ sudo apt-get install x11-apps
    $ export DISPLAY=:0
    $ xclock

where the last command is simply to check the setup. Once this is confirmed to be operating smoothly, the "plt_mesh=True" command will then run as anticipated and generate two png files (one for the geometry and one for the mesh) in NumBAT/backend/fortran/msh/. Note the X windows that open must be manually closed for the calculation to continue, and after unexpected restarts the X window may no longer display output but the png files will contain the necessary features.

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


