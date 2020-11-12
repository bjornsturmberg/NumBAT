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

coatx = 2;
coaty = 2;
coat_w = coatx/d_in_nm;
coat_h = coaty/d_in_nm;

lc = 0; // background and unitcell edge
lc_refine_1 = lc/1; // rib
lc_refine_2 = lc/1; // slab
lc_refine_3 = lc/1; // coat

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

// Coat
Point(15) = {d/2-slab_w/2, -hy+slab_h+coat_h, 0, lc_refine_3};
Point(16) = {d/2+slab_w/2, -hy+slab_h+coat_h, 0, lc_refine_3};
Point(17) = {d/2-(radius1+coat_w), -hy+slab_h+coat_h, 0, lc_refine_3};
Point(18) = {d/2+(radius1+coat_w), -hy+slab_h+coat_h, 0, lc_refine_3};
Point(19) = {d/2-(radius1+coat_w), -hy+slab_h+coat_h+2*radius1y, 0, lc_refine_3};
Point(20) = {d/2+(radius1+coat_w), -hy+slab_h+coat_h+2*radius1y, 0, lc_refine_3};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 7};
Line(6) = {7, 9};
Line(7) = {9, 10};
Line(8) = {10, 8};
Line(9) = {8, 7};
Line(10) = {8, 6};
Line(11) = {6, 14};
Line(12) = {14, 13};
Line(13) = {13, 5};

Physical Line(21) = {1};
Physical Line(22) = {2};
Physical Line(23) = {3};
Physical Line(24) = {4};

Line(25) = {5, 15};
Line(26) = {15, 17};
Line(27) = {17, 19};
Line(28) = {19, 20};
Line(29) = {20, 18};
Line(30) = {18, 16};
Line(31) = {16, 6};
Line Loop(32) = {4, 1, 2, 3};
Line Loop(33) = {28, 29, 30, 31, 11, 12, 13, 25, 26, 27};
Plane Surface(34) = {32, 33};
Line Loop(35) = {6, 7, 8, 9};
Plane Surface(36) = {35};
Line Loop(37) = {9, -5, -13, -12, -11, -10};
Plane Surface(38) = {37};
Line Loop(39) = {26, 27, 28, 29, 30, 31, -10, -8, -7, -6, -5, 25};
Plane Surface(40) = {39};
Physical Surface(1) = {34};
Physical Surface(2) = {36};
Physical Surface(3) = {38};
Physical Surface(4) = {40};
