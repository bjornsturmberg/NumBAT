 

FEM Mode Solvers
====================

.. _sec-newmesh-label:

Making New Mesh
------------------------------------------------

At some point you may well wish to study a structure that is not described by an existing NumBAT mesh template. In this section we provide an example of how to create a new mesh. In this case we will create a rib waveguide that is has a coating surrounding the guiding region.

Creating a mesh is typically a three step process: first we define the points that define the outline of the structures, then we define the lines connecting the points and the surfaces formed out of the lines. The first step is best done in a text editor in direct code, while the second can be done using the open source program `gmsh <http://geuz.org/gmsh/>`_ GUI. The third step involves adding some lines to the NumBAT backend.

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
    Point(17) = {d/2-rx/2-cx, -h+s1+s2+2*cy+ry, 0, lc};
    Point(18) = {d/2+rx/2+cx, -h+s1+s2+2*cy+ry, 0, lc};


Step 2
~~~~~~~~~~~~~~~~~~~~

To create the lines that connect the points, and the mesh surfaces it is easiest to use gmsh (although it can also be written directly in code). Open your geometry file in gmsh::
    
    NumBAT/backend/fortran/msh$ gmsh rib_coated_msh_template.geo

Navigate through the side menu to Modules/Geometry/Elementary entities/Add and click "Straight line". Now click consecutively on the point you wish to connect.

Navigate through the side menu to Modules/Geometry/Elementary entities/Add and click "Plane surface". Now click on the boundary of each enclosed area. Remember to separate out your inclusion from the background by highlighting it when asked for “hole boundaries". If the inclusion is complicated it is best to carve up the background area into smaller simpler areas that don't have any inclusions ("holes"), for example see slot coated.

Navigate through the side menu to Modules/Geometry/Physical groups/Add and click "Line". Now click on the lines that make up each side of the unit cell boundary, pressing the "e" key to end your selection once the each side is fully highlighted. 

Navigate through the side menu to Modules/Geometry/Physical groups/Add and click "Surface". Now click on all the surfaces of a given material type (in this example there is only one surface per material). It is crucial to remember the order you defined the physical surfaces in. Now open the .geo file in your favorite text editor, scroll to the bottom, and change the numbering of the physical surfaces to start at 1, and to increase by one per surface type. Eg. by tradition 1 is the background material, 2 is the waveguide, 3 is the bottom substrate, and 4 is the cladding. ::

    Physical Surface(1) = {24};
    Physical Surface(2) = {28};
    Physical Surface(3) = {30};
    Physical Surface(4) = {26};

The important thing is to make a note of the chosen labeling! This is best done by taking a screen-shot of the geometry in gmsh, labeling this with material types and physical dimensions, and then adding this file to the NumBAT/docs/msh_type_lib folder.


Step 3
~~~~~~~~~~~~~~~~~~~~

The last step is to add your geometry to the make_mesh function in NumBAT/backend/objects.py.

This involves adding a new elif statement for the inc_shape, in this case 'rib_coated', and then adding lines that define how the final mesh will be created based on the template. This involves giving the mesh a name, specifying the number of element types, and modifying the template geometric parameters. See objects.py for details.

One last thing, if the geometry contains only rectangular shapes, and all elements are therefore linear (rather than curvi-linear), you should also add the inc_shape name to the self.linear_element_shapes list in objects.py. This will ensure that the most efficient semi-analytic integration routines are used. If NumBAT is not told that the mesh is linear it will default to using numerical quadrature.









FEM Errors
-----------

There are 2 main errors that can be easily triggered within the Fortran FEM routines. These cause them to simulation to abort and the terminal to be unresponsive (until you kill python or the screen session).

The first of these is ::

    VALPR_64: info_32 != 0 : 
    VALPR_64: iparam_32(5) = 
    VALPR_64: number of converged values =    
    py_calc_modes.f: convergence problem with valpr_64
    py_calc_modes.f: You should probably increase resolution of mesh!
    py_calc_modes.f: n_conv != nval :

Long story short, this indicates that the FEM mesh is too coarse for solutions for higher order Bloch modes (Eigenvaules) to converge. 
This error is easily fixed by increasing the mesh resolution. Decrease 'lc_bkg' and/or increase 'lc2' etc.


The second error is :: 

    Error with _naupd, info_32 =           -8
    Check the documentation in _naupd.
    Aborting...

This is the opposite problem, when the mesh is so fine that the simulation is overloading the memory of the machine. More accurately the memory depends on the number of Eigenvalues being calculated as well as the number of FEM mesh points.
The best solution to this is to increase 'lc_bkg' and/or decrease 'lc2' etc.
