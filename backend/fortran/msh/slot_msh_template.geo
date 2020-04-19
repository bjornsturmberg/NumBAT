// Template mesh geometry file for a slot waveguide.

d = 1; // grating period
d_in_nm = 100;
dy_in_nm = 50;
dy = dy_in_nm/d_in_nm;
a1 = 20;
a1y = 10;
radius1 = (a1/(2*d_in_nm))*d;
radius1y = (a1y/(2*d_in_nm))*d;
a2 = 20;
radius2 = (a2/(2*d_in_nm))*d;
radius2y = radius1y;

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

// Slot centre
Point(7) = {-hx+d/2-radius1, -hy+slab_h, 0, lc_refine_1};
Point(8) = {-hx+d/2+radius1, -hy+slab_h, 0, lc_refine_1};
Point(9) = {-hx+d/2-radius1, -hy+2*radius1y+slab_h, 0, lc_refine_1};
Point(10) = {-hx+d/2+radius1, -hy+2*radius1y+slab_h, 0, lc_refine_1};

// Slot left wall
Point(20) = {-hx+d/2-radius1-radius2, -hy+slab_h, 0, lc_refine_1};
Point(21) = {-hx+d/2+radius1+radius2, -hy+slab_h, 0, lc_refine_1};
Point(22) = {-hx+d/2-radius1-radius2, -hy+2*radius1y+slab_h, 0, lc_refine_1};
Point(23) = {-hx+d/2+radius1+radius2, -hy+2*radius1y+slab_h, 0, lc_refine_1};

Point(11) = {0, -hy+slab_h, 0, lc};
Point(12) = {d, -hy+slab_h, 0, lc};
Point(15) = {-hx+d/2+radius1, 0, 0, lc};
Point(16) = {-hx+d/2-radius1, 0, 0, lc};
Point(17) = {-hx+d/2+radius1, -dy, 0, lc};
Point(18) = {-hx+d/2-radius1, -dy, 0, lc};
Point(24) = {-hx+d/2+radius1+radius2, 0, 0, lc};
Point(25) = {-hx+d/2-radius1-radius2, 0, 0, lc};
Point(26) = {-hx+d/2+radius1+radius2, -dy, 0, lc};
Point(27) = {-hx+d/2-radius1-radius2, -dy, 0, lc};


Line(6) = {7, 9};
Line(7) = {9, 10};
Line(8) = {10, 8};
Line(9) = {8, 7};
Line(11) = {6, 14};
Line(12) = {14, 13};
Line(13) = {13, 5};
Line(14) = {1, 11};
Line(15) = {11, 2};
Line(16) = {4, 12};
Line(17) = {12, 3};
Line(19) = {11, 5};
Line(20) = {6, 12};
Line(33) = {16, 9};
Line(34) = {10, 15};
Line(35) = {15, 16};
Line(45) = {18, 17};
Line(47) = {1, 25};
Line(48) = {25, 16};
Line(49) = {25, 22};
Line(50) = {22, 9};
Line(51) = {22, 20};
Line(52) = {20, 7};
Line(53) = {20, 5};
Line(54) = {8, 21};
Line(55) = {21, 6};
Line(56) = {21, 23};
Line(57) = {23, 24};
Line(58) = {15, 24};
Line(59) = {24, 4};
Line(60) = {10, 23};
Line(85) = {2, 27};
Line(86) = {27, 18};
Line(87) = {17, 26};
Line(88) = {26, 3};
Line Loop(61) = {47, 49, 51, 53, -19, -14};
Plane Surface(62) = {61};
Line Loop(89) = {19, -13, -12, -11, 20, 17, -88, -87, -45, -86, -85, -15};
Plane Surface(90) = {89};
Line Loop(65) = {53, -13, -12, -11, -55, -54, 9, -52};
Plane Surface(66) = {65};
Line Loop(67) = {52, 6, -50, 51};
Plane Surface(68) = {67};
Line Loop(69) = {50, -33, -48, 49};
Plane Surface(70) = {69};
Line Loop(71) = {33, 7, 34, 35};
Plane Surface(72) = {71};
Line Loop(73) = {7, 8, 9, 6};
Plane Surface(74) = {73};
Line Loop(75) = {8, 54, 56, -60};
Plane Surface(76) = {75};
Line Loop(77) = {55, 20, -16, -59, -57, -56};
Plane Surface(78) = {77};
Line Loop(79) = {57, -58, -34, 60};
Plane Surface(80) = {79};
Physical Line(81) = {47, 48, 35, 58, 59};
Physical Line(82) = {16, 17};
Physical Line(84) = {14, 15};
Physical Line(91) = {85, 86, 45, 87, 88};
Physical Surface(1) = {62, 70, 72, 80, 78, 64, 90};
Physical Surface(2) = {74};
Physical Surface(3) = {66};
Physical Surface(4) = {68, 76};
