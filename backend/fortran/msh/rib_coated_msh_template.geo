// Coated rib waveguide mesh template.

// Meshing parameters
lc = 0; // base level meshing parameter
lc2 = lc/1; // mesh refinement parameters
lc3 = lc/1; // mesh refinement parameters

// Unit cell
d_in_nm = 1000;
dy_in_nm = 600;
d = 1; // normalised unitcell limits x
h = dy_in_nm/d_in_nm; // unitcell limit y

Point(1) = {0, 0, 0, lc};
Point(2) = {0, -h, 0, lc};
Point(3) = {d, -h, 0, lc};
Point(4) = {d, 0, 0,lc};

// Bottom slab
slab1 = 100;
s1 = slab1/d_in_nm;
Point(5) = {0, -h+s1, 0, lc};
Point(6) = {d, -h+s1, 0, lc};

// Waveguide slab
slab2 = 50;
s2 = slab2/d_in_nm;
Point(7) = {0, -h+s1+s2, 0, lc};
Point(8) = {d, -h+s1+s2, 0, lc};

// Rib
ribx = 200;
riby = 30;
rx = ribx/d_in_nm;
ry = riby/d_in_nm;
Point(9) = {d/2-rx/2, -h+s1+s2, 0, lc2};
Point(10) = {d/2+rx/2, -h+s1+s2, 0, lc2};
Point(11) = {d/2-rx/2, -h+s1+s2+ry, 0, lc2};
Point(12) = {d/2+rx/2, -h+s1+s2+ry, 0, lc2};

// Coating
coatx = 20; // around rib
coaty = 20;
cx = coatx/d_in_nm;
cy = coaty/d_in_nm;
Point(13) = {0, -h+s1+s2+cy, 0, lc};
Point(14) = {d, -h+s1+s2+cy, 0, lc};
Point(15) = {d/2-rx/2-cx, -h+s1+s2+cy, 0, lc};
Point(16) = {d/2+rx/2+cx, -h+s1+s2+cy, 0, lc};
Point(17) = {d/2-rx/2-cx, -h+s1+s2+2*cy, 0, lc};
Point(18) = {d/2+rx/2+cx, -h+s1+s2+2*cy, 0, lc};Line(1) = {1, 4};


Line(2) = {4, 14};
Line(3) = {14, 8};
Line(4) = {8, 6};
Line(5) = {6, 3};
Line(6) = {3, 2};
Line(7) = {2, 5};
Line(8) = {5, 7};
Line(9) = {7, 13};
Line(10) = {13, 1};
Line(11) = {5, 6};
Line(12) = {8, 10};
Line(13) = {10, 9};
Line(14) = {9, 7};
Line(15) = {13, 15};
Line(16) = {15, 17};
Line(17) = {17, 18};
Line(18) = {18, 16};
Line(19) = {14, 16};
Line(20) = {10, 12};
Line(21) = {12, 11};
Line(22) = {11, 9};
