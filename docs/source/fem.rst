FEM Mode Solvers
====================

Making New Mesh
------------------------------------------------

At some point you may well wish to study a structure that is not described by an existing NumBAT mesh template. In this section we provide an example of how to create a new mesh. In this case we will create a rib waveguide that is has a coating surrounding the guiding region.

Creating a mesh is typically a two step process: first we define the points that define the outline of the structures, then we define the lines connecting the points and the surfaces formed out of the lines. The first step is best done in a text editor in direct code, while the second can be done using the open source program `gmsh <http://geuz.org/gmsh/>`_ GUI. 

To start we are going to make a copy of NumBAT/backend/fortran/msh/empty_msh_template.geo ::

    $ cd NumBAT/backend/fortran/msh/
    $ cp empty_msh_template.geo rib_coated_msh_template.geo

Step 1
~~~~~~~~~~~~~~~~~~~~

Opening the new file in a text editor you see it contains points defining the unit cell. The points are defined as ::

    Point(1) = {x, y, z, meshing_value}

We start by adding the two points that define the top of the substrate (the bottom will be the bottom edge of the unit cell at {0, -h} and {d, -h}). We use a placeholder slab thickness of 100 nm, which is normalised by the width of the unit cell. ::

    slab1 = 100;
    s1 = slab1/d_in_nm;
    Point(5) = {0, -h+s1, 0, lc};
    Point(6) = {d, -h+s1, 0, lc};

We then add a further layer on top of the bottom slab, this time using a placeholder thickness of 50 nm. Note that each point must be labeled by a unique number.::

	slab2 = 50;
	s2 = slab2/d_in_nm;
	Point(7) = {0, -h+s1+s2, 0, lc};
	Point(8) = {d, -h+s1+s2, 0, lc};
 
We next define the peak of the rib, which involves a width and a height, ::

	ribx = 200;
	riby = 30;
	rx = ribx/d_in_nm;
	ry = riby/d_in_nm;
	Point(9) = {d/2-rx/2, -h+s1+s2, 0, lc2};
	Point(10) = {d/2+rx/2, -h+s1+s2, 0, lc2};
	Point(11) = {d/2-rx/2, -h+s1+s2+ry, 0, lc2};
	Point(12) = {d/2+rx/2, -h+s1+s2+ry, 0, lc2};

Lastly we coat the whole structure with a conformal layer. ::

	coatx = 20;
	coaty = 20;
	cx = coatx/d_in_nm;
	cy = coaty/d_in_nm;
	Point(13) = {0, -h+s1+s2+cy, 0, lc};
	Point(14) = {d, -h+s1+s2+cy, 0, lc};
	Point(15) = {d/2-rx/2-cx, -h+s1+s2+cy, 0, lc};
	Point(16) = {d/2+rx/2+cx, -h+s1+s2+cy, 0, lc};
	Point(17) = {d/2-rx/2-cx, -h+s1+s2+2*cy, 0, lc};
	Point(18) = {d/2+rx/2+cx, -h+s1+s2+2*cy, 0, lc};


Step 2
~~~~~~~~~~~~~~~~~~~~

To create the lines that connect the points, and the mesh surfaces it is easiest to use gmsh (although it can also be written directly in code). Open your geometry file in gmsh::
	
	NumBAT/backend/fortran/msh$ gmsh rib_coated_msh_template.geo

Navigate through the side menu to Modules/Geometry/Elementary Entities/Add and click "Straight Line". Now click consecutively on the point you wish to connect.

Next navigate through the side menu to Modules/Geometry/Elementary Entities/Add and click "Plane Surface".


FEM Errors
-----------

There are 2 errors that can be easily triggered within the Fortran FEM routines. These both cause them to simulation to abort and the terminal to be unresponsive (until you kill python or the screen session).

The first of these is ::

	Error with _naupd, info_32 =           -3
	Check the documentation in _naupd.
	Aborting...

Long story short, this indicates that the FEM mesh is too coarse for solutions for higher order Bloch modes (Eigenvaules) to converge. To see this run the simulation with FEM_debug = 1 (in mode_calcs.py) and it will print the number of converged Eigenvalues nconv != nval.
This error is easily fixed by increasing the mesh resolution. Decrease 'lc_bkg' and/or increase 'lc2' etc.


The second error is :: 

	Error with _naupd, info_32 =           -8
	Check the documentation in _naupd.
	Aborting...

This is the opposite problem, when the mesh is so fine that the simulation is overloading the memory of the machine. More accurately the memory depends on the number of Eigenvalues being calculated as well as the number of FEM mesh points.
The best solution to this is to increase 'lc_bkg' and/or decrease 'lc2' etc.