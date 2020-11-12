// Template mesh geometry file for a single inclusion on a slab.
// Inclusion can be circular/elliptical (default), or square/rectangular.

// Force Gmsh to use legacy msh file format v2
Mesh.MshFileVersion = 2.2;

d = 1; // grating period
ff = 0;
d_in_nm = 100;
dy_in_nm = 50;
dy = dy_in_nm/d_in_nm;
a1 = 20;
a1y = 10;
radius1 = (a1/(2*d_in_nm))*d;
radius1y = (a1y/(2*d_in_nm))*d;

a2 = 10;
a2y = 20;
radius2 = (a2/(2*d_in_nm))*d;
radius2y = (a2y/(2*d_in_nm))*d;
sep = 10;
b = sep/(2*d_in_nm);

rect = 1;

slab_width = d_in_nm;
slab_height = 10;
slab_w = slab_width/d_in_nm;
slab_h = slab_height/d_in_nm;
slab_w_full = 0;
If(slab_w == 1)
    slab_w_full = 1;
EndIf

slab2_width = d_in_nm;
slab2_height = 5;
slab2_w = slab2_width/d_in_nm;
slab2_h = slab2_height/d_in_nm;
slab2_w_full = 0;
If(slab2_w == 1)
    slab2_w_full = 1;
EndIf

lc = 0; // 0.501 0.201 0.0701;
lc_refine_1 = lc/1; // on cylinder surfaces
lc_refine_2 = lc/1; // cylinder1 centres
lc_refine_3 = lc/1; // centres of top and bottom
lc_refine_4 = lc/1; // slab

hy = dy; // Thickness: square profile => hy=d
hx = 0.;


Point(1) = {0, 0, 0, lc};
Point(2) = {-hx, -hy, 0, lc};
Point(3) = {-hx+d, -hy, 0, lc};
Point(4) = {d, 0, 0,lc};

// Slab
Point(252) = {0, -hy+slab_h+slab2_h, 0, lc_refine_4};
Point(253) = {d, -hy+slab_h+slab2_h, 0, lc_refine_4};
Point(350) = {0, -hy+slab2_h, 0, lc_refine_4};
Point(351) = {d, -hy+slab2_h, 0, lc_refine_4};

// Inclusion
Point(5) = {-radius1-b-hx+d/2, -hy+slab_h+slab2_h+radius1y, 0, lc_refine_2};
Point(6) = {-radius1-b-hx+d/2, -hy+slab_h+slab2_h+2*radius1y, 0, lc_refine_1};
Point(7) = {-radius1-b-hx+d/2-radius1, -hy+slab_h+slab2_h+radius1y, 0, lc_refine_1};
Point(8) = {-radius1-b-hx+d/2, -hy+slab_h+slab2_h, 0, lc_refine_1};
Point(9) = {-radius1-b-hx+d/2+radius1, -hy+slab_h+slab2_h+radius1y, 0, lc_refine_1};

Point(10) = {-2*radius1-b-hx+d/2, 0, 0, lc_refine_3};
Point(11) = {-b-hx+d/2, 0, 0, lc_refine_3};
Point(12) = {b-hx+d/2, 0, 0, lc_refine_3};
Point(13) = {2*radius2+b-hx+d/2, 0, 0, lc_refine_3};

Point(18) = {radius2+b-hx+d/2, -hy+slab_h+slab2_h+radius2y, 0, lc_refine_2};
Point(19) = {radius2+b-hx+d/2, -hy+slab_h+slab2_h+2*radius2y, 0, lc_refine_1};
Point(20) = {radius2+b-hx+d/2-radius2, -hy+slab_h+slab2_h+radius2y, 0, lc_refine_1};
Point(21) = {radius2+b-hx+d/2, -hy+slab_h+slab2_h, 0, lc_refine_1};
Point(22) = {radius2+b-hx+d/2+radius2, -hy+slab_h+slab2_h+radius2y, 0, lc_refine_1};

Line(1) = {1, 10};
Line(2) = {10, 11};
Line(3) = {11, 12};
Line(4) = {12, 13};
Line(5) = {13, 4};
Line(6) = {1, 252};
Line(7) = {252, 350};
Line(107) = {2, 350};
Line(9) = {351, 253};
Line(109) = {351, 3};
Line(10) = {253, 4};
Physical Line(110) = {6, 7, 107};
Physical Line(111) = {10, 9, 109};
Physical Line(14) = {5, 4, 3, 2, 1};
Line(15) = {6, 5};
Line(16) = {5, 8};
Line(17) = {5, 7};
Line(18) = {5, 9};
Line(19) = {19, 18};
Line(20) = {18, 21};
Line(21) = {20, 18};
Line(22) = {18, 22};

Line(112) = {350, 351};

Line(23) = {2, 3};
Physical Line(24) = {23};


If(rect == 0)
    Ellipsis(27) = {9,5,6,6};
    Ellipsis(28) = {6,5,7,7};
    Ellipsis(29) = {7,5,8,8};
    Ellipsis(30) = {8,5,9,9};

    Ellipsis(31) = {22,18,19,19};
    Ellipsis(32) = {19,18,20,20};
    Ellipsis(33) = {20,18,21,21};
    Ellipsis(34) = {21,18,22,22};

    Line(35) = {10, 7};
    Line(36) = {11, 9};
    Line(37) = {12, 20};
    Line(38) = {13, 22};

    Line(39) = {252, 8};
    Line(40) = {8, 21};
    Line(41) = {21, 253};
    Line Loop(113) = {6, 39, -29, -35, -1};
    Plane Surface(114) = {113};
    Line Loop(115) = {39, 40, 41, -9, -112, -7};
    Plane Surface(116) = {115};
    Line Loop(117) = {38, -34, 41, 10, -5};
    Plane Surface(118) = {117};
    Line Loop(119) = {112, 109, -23, 107};
    Plane Surface(120) = {119};

    Line Loop(121) = {35, -28, -27, -36, -2};
    Plane Surface(122) = {121};
    Line Loop(123) = {36, -30, 40, -33, -37, -3};
    Plane Surface(124) = {123};
    Line Loop(125) = {37, -32, -31, -38, -4};
    Plane Surface(126) = {125};
    Line Loop(127) = {31, 19, 22};
    Plane Surface(128) = {127};
    Line Loop(129) = {22, -34, -20};
    Plane Surface(130) = {129};
    Line Loop(131) = {20, -33, 21};
    Plane Surface(132) = {131};
    Line Loop(133) = {21, -19, 32};
    Plane Surface(134) = {133};
    Line Loop(135) = {27, 15, 18};
    Plane Surface(136) = {135};
    Line Loop(137) = {15, 17, -28};
    Plane Surface(138) = {137};
    Line Loop(139) = {17, 29, -16};
    Plane Surface(140) = {139};
    Line Loop(141) = {16, 30, -18};
    Plane Surface(142) = {141};

    Physical Surface(1) = {114, 122, 124, 126, 118};
    Physical Surface(2) = {134, 128, 130, 132};
    Physical Surface(5) = {136, 142, 140, 138};
    Physical Surface(3) = {116};
    Physical Surface(4) = {120};

EndIf


If(rect == 1)
    Point(150) = {-radius1-b-hx+d/2+radius1, -hy+slab_h+slab2_h+2*radius1y, 0,lc_refine_2};
    Point(151) = {-radius1-b-hx+d/2-radius1, -hy+slab_h+slab2_h+2*radius1y, 0,lc_refine_2};
    Point(152) = {-radius1-b-hx+d/2+radius1, -hy+slab_h+slab2_h, 0,lc_refine_2};
    Point(153) = {-radius1-b-hx+d/2-radius1, -hy+slab_h+slab2_h, 0,lc_refine_2};

    Point(154) = {radius2+b-hx+d/2+radius2, -hy+slab_h+slab2_h+2*radius2y, 0,lc_refine_2};
    Point(155) = {radius2+b-hx+d/2-radius2, -hy+slab_h+slab2_h+2*radius2y, 0,lc_refine_2};
    Point(156) = {radius2+b-hx+d/2+radius2, -hy+slab_h+slab2_h, 0,lc_refine_2};
    Point(157) = {radius2+b-hx+d/2-radius2, -hy+slab_h+slab2_h, 0,lc_refine_2};

    Line(48) = {10, 151};
    Line(26) = {151, 6};
    Line(27) = {6, 150};
    Line(28) = {150, 11};
    Line(29) = {150, 9};
    Line(30) = {152, 9};
    Line(31) = {152, 8};
    Line(32) = {8, 153};
    Line(33) = {153, 7};
    Line(34) = {7, 151};
    Line(36) = {12, 155};
    Line(37) = {155, 19};
    Line(38) = {19, 154};
    Line(39) = {154, 13};
    Line(40) = {154, 22};
    Line(41) = {22, 156};
    Line(43) = {156, 21};
    Line(44) = {21, 157};
    Line(45) = {157, 20};
    Line(46) = {20, 155};
    Line(47) = {157, 152};

    Line(35) = {252, 153};
    Line(42) = {156, 253};

    Line Loop(113) = {6, 35, 33, 34, -48, -1};
    Plane Surface(114) = {113};
    Line Loop(115) = {35, -32, -31, -47, -44, -43, 42, -9, -112, -7};
    Plane Surface(116) = {115};
    Line Loop(117) = {107, 112, 109, -23};
    Plane Surface(118) = {117};
    Line Loop(119) = {48, 26, 27, 28, -2};
    Plane Surface(120) = {119};
    Line Loop(121) = {26, 15, 17, 34};
    Plane Surface(122) = {121};
    Line Loop(123) = {15, 18, -29, -27};
    Plane Surface(124) = {123};
    Line Loop(125) = {30, -18, 16, -31};
    Plane Surface(126) = {125};
    Line Loop(127) = {17, -33, -32, -16};
    Plane Surface(128) = {127};
    Line Loop(129) = {28, 3, 36, -46, -45, 47, 30, -29};
    Plane Surface(130) = {129};
    Line Loop(131) = {36, 37, 38, 39, -4};
    Plane Surface(132) = {131};
    Line Loop(133) = {37, 19, -21, 46};
    Plane Surface(134) = {133};
    Line Loop(135) = {19, 22, -40, -38};
    Plane Surface(136) = {135};
    Line Loop(137) = {22, 41, 43, -20};
    Plane Surface(138) = {137};
    Line Loop(139) = {20, 44, 45, 21};
    Plane Surface(140) = {139};
    Line Loop(141) = {41, 42, 10, -5, -39, 40};
    Plane Surface(142) = {141};
    Physical Surface(1) = {114, 120, 130, 132, 142};
    Physical Surface(2) = {124, 126, 128, 122};
    Physical Surface(5) = {138, 140, 134, 136};
    Physical Surface(3) = {116};
    Physical Surface(4) = {118};
EndIf

