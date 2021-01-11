// Template mesh geometry file for a rib waveguide.

// Force Gmsh to use legacy msh file format v2
Mesh.MshFileVersion = 2.2;

d = 1; // grating period
d_in_nm = 100;
dy_in_nm = 50;
dy = dy_in_nm/d_in_nm;
a1 = 20;
a1y = 10;
radius1 = (a1/(2*d_in_nm))*d;
radius1y = (a1y/(2*d_in_nm))*d;

slabx = 80;
slaby = 10;
slab_w = slabx/d_in_nm;
slab_h = slaby/d_in_nm;

lc = 0; // background and unitcell edge
lc_refine_1 = lc/1; // rib
lc_refine_2 = lc/1; // slab

hy = dy/2 + (slab_h/2) + radius1y; // 
hx = 0.;


Point(1) = {0, 0, 0, lc};
Point(2) = {-hx, -dy, 0, lc};
Point(3) = {-hx+d, -dy, 0, lc};
Point(4) = {d, 0, 0,lc};

// Slab
Point(5) = {d/2-slab_w/2, -hy+slab_h, 0, lc_refine_2};
Point(6) = {d/2+slab_w/2, -hy+slab_h, 0, lc_refine_2};
Point(13) = {d/2-slab_w/2, -hy, 0, lc_refine_2};
Point(14) = {d/2+slab_w/2, -hy, 0, lc_refine_2};

// Rib
Point(7) = {-hx+d/2-radius1, -hy+slab_h, 0, lc_refine_1};
Point(8) = {-hx+d/2+radius1, -hy+slab_h, 0, lc_refine_1};
Point(9) = {-hx+d/2-radius1, -hy+2*radius1y+slab_h, 0, lc_refine_1};
Point(10) = {-hx+d/2+radius1, -hy+2*radius1y+slab_h, 0, lc_refine_1};

Point(11) = {0, -hy+slab_h, 0, lc};
Point(12) = {d, -hy+slab_h, 0, lc};
Point(15) = {-hx+d/2+radius1, 0, 0, lc};
Point(16) = {-hx+d/2-radius1, 0, 0, lc};
Point(17) = {-hx+d/2+radius1, -dy, 0, lc};
Point(18) = {-hx+d/2-radius1, -dy, 0, lc};

Line(5) = {5, 7};
Line(6) = {7, 9};
Line(7) = {9, 10};
Line(8) = {10, 8};
Line(9) = {8, 7};
Line(10) = {8, 6};
Line(11) = {6, 14};
Line(12) = {14, 13};
Line(13) = {13, 5};
Line(14) = {1, 11};
Line(15) = {11, 2};
Line(16) = {4, 12};
Line(17) = {12, 3};
Line Loop(18) = {14, 15, 2, -17, -16, 4};
Line(19) = {11, 5};
Line(20) = {6, 12};
Line Loop(25) = {6, 7, 8, 9};
Plane Surface(26) = {25};
Line Loop(27) = {9, -5, -13, -12, -11, -10};
Plane Surface(28) = {27};
Line(32) = {1, 16};
Line(33) = {16, 9};
Line(34) = {10, 15};
Line(35) = {15, 16};
Line(36) = {15, 4};
Line Loop(37) = {5, 6, -33, -32, 14, 19};
Plane Surface(38) = {37};
Line Loop(39) = {7, 34, 35, 33};
Plane Surface(40) = {39};
Line Loop(41) = {34, 36, 16, -20, -10, -8};
Plane Surface(42) = {41};
Line(44) = {2, 18};
Line(45) = {18, 17};
Line(46) = {17, 3};
Line Loop(48) = {19, -13, -12, -11, 20, 17, -46, -45, -44, -15};
Plane Surface(49) = {48};


Physical Line(29) = {14, 15};
Physical Line(31) = {17, 16};
Physical Line(43) = {32, 35, 36};
Physical Line(47) = {44, 45, 46};

Physical Surface(1) = {38, 40, 42, 49};
Physical Surface(2) = {26};
Physical Surface(3) = {28};