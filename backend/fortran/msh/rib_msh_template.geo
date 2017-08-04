// Template mesh geometry file for a rib waveguide.

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
lc2 = lc/1; // rib
lc3 = lc/1; // slab

hy = dy/2 + (slab_h/2) + radius1y; // 
hx = 0.;


Point(1) = {0, 0, 0, lc};
Point(2) = {-hx, -dy, 0, lc};
Point(3) = {-hx+d, -dy, 0, lc};
Point(4) = {d, 0, 0,lc};

// Slab
Point(5) = {d/2-slab_w/2, -hy+slab_h, 0, lc3};
Point(6) = {d/2+slab_w/2, -hy+slab_h, 0, lc3};
Point(13) = {d/2-slab_w/2, -hy, 0, lc3};
Point(14) = {d/2+slab_w/2, -hy, 0, lc3};

// Rib
Point(7) = {-hx+d/2-radius1, -hy+slab_h, 0, lc2};
Point(8) = {-hx+d/2+radius1, -hy+slab_h, 0, lc2};
Point(9) = {-hx+d/2-radius1, -hy+2*radius1y+slab_h, 0, lc2};
Point(10) = {-hx+d/2+radius1, -hy+2*radius1y+slab_h, 0, lc2};

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
Line Loop(14) = {1, 2, 3, 4};
Line Loop(15) = {5, 6, 7, 8, 10, 11, 12, 13};
Plane Surface(16) = {14, 15};
Line Loop(17) = {7, 8, 9, 6};
Plane Surface(18) = {17};
Line Loop(19) = {9, -5, -13, -12, -11, -10};
Plane Surface(20) = {19};

Physical Line(21) = {1};
Physical Line(22) = {2};
Physical Line(23) = {3};
Physical Line(24) = {4};

Physical Surface(1) = {16};
Physical Surface(2) = {18};
Physical Surface(3) = {20};
