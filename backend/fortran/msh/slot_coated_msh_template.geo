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

c1y = 10;
cov1y = (c1y/(2*d_in_nm))*d;

lc = 0; // background and unitcell edge
lc_refine_1 = lc/1; // rib
lc_refine_2 = lc/1; // slab
lc_refine_3 = lc/1; // coating

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

// Coating
Point(30) = {-hx+d/2-radius1, -hy+2*radius1y+slab_h+cov1y, 0, lc_refine_3};
Point(31) = {-hx+d/2+radius1, -hy+2*radius1y+slab_h+cov1y, 0, lc_refine_3};
Point(32) = {-hx+d/2-radius1-radius2, -hy+2*radius1y+slab_h+cov1y, 0, lc_refine_3};
Point(33) = {-hx+d/2+radius1+radius2, -hy+2*radius1y+slab_h+cov1y, 0, lc_refine_3};

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
Line(35) = {15, 16};
Line(45) = {18, 17};
Line(47) = {1, 25};
Line(48) = {25, 16};
Line(50) = {22, 9};
Line(51) = {22, 20};
Line(52) = {20, 7};
Line(53) = {20, 5};
Line(54) = {8, 21};
Line(55) = {21, 6};
Line(56) = {21, 23};
Line(58) = {15, 24};
Line(59) = {24, 4};
Line(60) = {10, 23};
Line(85) = {2, 27};
Line(86) = {27, 18};
Line(87) = {17, 26};
Line(88) = {26, 3};
Physical Line(81) = {47, 48, 35, 58, 59};
Physical Line(82) = {16, 17};
Physical Line(84) = {14, 15};
Physical Line(91) = {85, 86, 45, 87, 88};
Line(92) = {32, 22};
Line(93) = {32, 30};
Line(94) = {30, 9};
Line(95) = {30, 31};
Line(96) = {31, 33};
Line(97) = {33, 23};
Line(98) = {10, 31};
Line(99) = {31, 15};
Line(100) = {24, 33};
Line(101) = {30, 16};
Line(102) = {25, 32};
Line Loop(103) = {102, 93, 101, -48};
Plane Surface(104) = {103};
Line Loop(105) = {101, -35, -99, -95};
Plane Surface(106) = {105};
Line Loop(107) = {58, 100, -96, 99};
Plane Surface(108) = {107};
Line Loop(109) = {96, 97, -60, 98};
Plane Surface(110) = {109};
Line Loop(111) = {98, -95, 94, 7};
Plane Surface(112) = {111};
Line Loop(113) = {94, -50, -92, 93};
Plane Surface(114) = {113};
Line Loop(115) = {102, 92, 51, 53, -19, -14, 47};
Plane Surface(116) = {115};
Line Loop(117) = {50, -6, -52, -51};
Plane Surface(118) = {117};
Line Loop(119) = {53, -13, -12, -11, -55, -54, 9, -52};
Plane Surface(120) = {119};
Line Loop(121) = {8, 9, 6, 7};
Plane Surface(122) = {121};
Line Loop(123) = {8, 54, 56, -60};
Plane Surface(124) = {123};
Line Loop(125) = {56, -97, -100, 59, 16, -20, -55};
Plane Surface(126) = {125};
Line Loop(127) = {11, 12, 13, -19, 15, 85, 86, 45, 87, 88, -17, -20};
Plane Surface(128) = {127};
Physical Surface(1) = {116, 128, 126, 108, 106, 104};
Physical Surface(2) = {122};
Physical Surface(3) = {120};
Physical Surface(4) = {124, 118};
Physical Surface(5) = {114, 112, 110};
