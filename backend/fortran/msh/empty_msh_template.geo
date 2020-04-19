// Template mesh geometry file from which to build future template mesh.

// Meshing parameters
lc = 0; // base level meshing parameter
lc_refine_1 = lc/1; // mesh refinement parameters
lc_refine_2 = lc/1; // mesh refinement parameters

// Unit cell
d_in_nm = 1000;
dy_in_nm = 600;
d = 1; // normalised unitcell limits x
h = dy_in_nm/d_in_nm; // unitcell limit y

Point(1) = {0, 0, 0, lc};
Point(2) = {0, -h, 0, lc};
Point(3) = {d, -h, 0, lc};
Point(4) = {d, 0, 0,lc};