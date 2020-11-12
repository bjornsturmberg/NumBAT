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

slab_width = d_in_nm-10;
slab_height = 10;
slab_w = slab_width/d_in_nm;
slab_h = slab_height/d_in_nm;
slab_w_full = 0;
If(slab_w == 1)
    slab_w_full = 1;
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
Point(252) = {0, -hy+slab_h, 0, lc_refine_4};
Point(253) = {d, -hy+slab_h, 0, lc_refine_4};

// Inclusion
Point(5) = {-radius1-b-hx+d/2, -hy+slab_h+radius1y, 0, lc_refine_2};
Point(6) = {-radius1-b-hx+d/2, -hy+slab_h+2*radius1y, 0, lc_refine_1};
Point(7) = {-radius1-b-hx+d/2-radius1, -hy+slab_h+radius1y, 0, lc_refine_1};
Point(8) = {-radius1-b-hx+d/2, -hy+slab_h, 0, lc_refine_1};
Point(9) = {-radius1-b-hx+d/2+radius1, -hy+slab_h+radius1y, 0, lc_refine_1};

Point(10) = {-2*radius1-b-hx+d/2, 0, 0, lc_refine_3};
Point(11) = {-b-hx+d/2, 0, 0, lc_refine_3};
Point(12) = {b-hx+d/2, 0, 0, lc_refine_3};
Point(13) = {2*radius2+b-hx+d/2, 0, 0, lc_refine_3};

Point(18) = {radius2+b-hx+d/2, -hy+slab_h+radius2y, 0, lc_refine_2};
Point(19) = {radius2+b-hx+d/2, -hy+slab_h+2*radius2y, 0, lc_refine_1};
Point(20) = {radius2+b-hx+d/2-radius2, -hy+slab_h+radius2y, 0, lc_refine_1};
Point(21) = {radius2+b-hx+d/2, -hy+slab_h, 0, lc_refine_1};
Point(22) = {radius2+b-hx+d/2+radius2, -hy+slab_h+radius2y, 0, lc_refine_1};

Line(1) = {1, 10};
Line(2) = {10, 11};
Line(3) = {11, 12};
Line(4) = {12, 13};
Line(5) = {13, 4};
Line(6) = {1, 252};
Line(7) = {252, 2};
Line(9) = {3, 253};
Line(10) = {253, 4};
Physical Line(11) = {6, 7};
Physical Line(13) = {9, 10};
Physical Line(14) = {5, 4, 3, 2, 1};
Line(15) = {6, 5};
Line(16) = {5, 8};
Line(17) = {5, 7};
Line(18) = {5, 9};
Line(19) = {19, 18};
Line(20) = {18, 21};
Line(21) = {20, 18};
Line(22) = {18, 22};

If(slab_w_full == 0)
    Point(250) = {d/2-slab_w/2, -hy+slab_h, 0, lc_refine_4};
    Point(251) = {d/2+slab_w/2, -hy+slab_h, 0, lc_refine_4};
    Point(254) = {d/2-slab_w/2, -hy, 0, lc_refine_4};
    Point(255) = {d/2+slab_w/2, -hy, 0, lc_refine_4};
    Line(23) = {2, 254};
    Line(24) = {254, 255};
    Line(25) = {255, 3};
    Physical Line(26) = {23, 24, 25};
EndIf
If(slab_w_full == 1)
    Line(23) = {2, 3};
    Physical Line(24) = {23};
EndIf

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

    If(slab_w_full == 1)
        Line(39) = {252, 8};
        Line(40) = {8, 21};
        Line(41) = {21, 253};
        Line Loop(42) = {6, 39, -29, -35, -1};
        Plane Surface(43) = {42};
        Line Loop(44) = {39, 40, 41, -9, -23, -7};
        Plane Surface(45) = {44};
        Line Loop(46) = {41, 10, -5, 38, -34};
        Plane Surface(47) = {46};
    EndIf
    If(slab_w_full == 0)
        Line(39) = {252, 250};
        Line(40) = {250, 254};
        Line(41) = {250, 8};
        Line(42) = {8, 21};
        Line(43) = {21, 251};
        Line(44) = {251, 253};
        Line(45) = {251, 255};
        Line Loop(46) = {39, 40, -23, -7};
        Plane Surface(47) = {46};
        Line Loop(48) = {6, 39, 41, -29, -35, -1};
        Plane Surface(49) = {48};
        Line Loop(50) = {41, 42, 43, 45, -24, -40};
        Plane Surface(51) = {50};
        Line Loop(52) = {44, -9, -25, -45};
        Plane Surface(53) = {52};
        Line Loop(54) = {44, 10, -5, 38, -34, 43};
        Plane Surface(55) = {54};
    EndIf

    Line Loop(56) = {35, -28, -27, -36, -2};
    Plane Surface(57) = {56};
    Line Loop(60) = {37, -32, -31, -38, -4};
    Plane Surface(61) = {60};
    Line Loop(62) = {32, 21, -19};
    Plane Surface(63) = {62};
    Line Loop(64) = {19, 22, 31};
    Plane Surface(65) = {64};
    Line Loop(66) = {22, -34, -20};
    Plane Surface(67) = {66};
    Line Loop(68) = {20, -33, 21};
    Plane Surface(69) = {68};
    Line Loop(70) = {30, -18, 16};
    Plane Surface(71) = {70};
    Line Loop(72) = {16, -29, -17};
    Plane Surface(73) = {72};
    Line Loop(74) = {17, -28, 15};
    Plane Surface(75) = {74};
    Line Loop(76) = {15, 18, 27};
    Plane Surface(77) = {76};
    Physical Surface(2) = {75, 77, 71, 73};
    Physical Surface(3) = {63, 65, 67, 69};

    If(slab_w_full == 1)
        Line Loop(78) = {36, -30, 40, -33, -37, -3};
        Plane Surface(79) = {78};
        Physical Surface(1) = {43, 57, 79, 61, 47};
        Physical Surface(4) = {45};
    EndIf
    If(slab_w_full == 0)
        Line Loop(78) = {36, -30, 42, -33, -37, -3};
        Plane Surface(79) = {78};
        Physical Surface(1) = {49, 57, 79, 61, 55};
        Physical Surface(5) = {53, 47};
        Physical Surface(4) = {51};
    EndIf
EndIf


If(rect == 1)
    Point(150) = {-radius1-b-hx+d/2+radius1, -hy+slab_h+2*radius1y, 0,lc_refine_2};
    Point(151) = {-radius1-b-hx+d/2-radius1, -hy+slab_h+2*radius1y, 0,lc_refine_2};
    Point(152) = {-radius1-b-hx+d/2+radius1, -hy+slab_h, 0,lc_refine_2};
    Point(153) = {-radius1-b-hx+d/2-radius1, -hy+slab_h, 0,lc_refine_2};

    Point(154) = {radius2+b-hx+d/2+radius2, -hy+slab_h+2*radius2y, 0,lc_refine_2};
    Point(155) = {radius2+b-hx+d/2-radius2, -hy+slab_h+2*radius2y, 0,lc_refine_2};
    Point(156) = {radius2+b-hx+d/2+radius2, -hy+slab_h, 0,lc_refine_2};
    Point(157) = {radius2+b-hx+d/2-radius2, -hy+slab_h, 0,lc_refine_2};

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

    If(slab_w_full == 1)
        Line(35) = {252, 153};
        Line(42) = {156, 253};
        Line Loop(49) = {6, 35, 33, 34, -48, -1};
        Plane Surface(50) = {49};
        Line Loop(51) = {35, -32, -31, -47, -44, -43, 42, -9, -23, -7};
        Plane Surface(52) = {51};
        Line Loop(53) = {42, 10, -5, -39, 40, 41};
        Plane Surface(54) = {53};
    EndIf
    If(slab_w_full == 0)
        Line(49) = {252, 250};
        Line(50) = {250, 153};
        Line(51) = {250, 254};
        Line(52) = {156, 251};
        Line(53) = {251, 253};
        Line(54) = {251, 255};
        Line Loop(55) = {6, 49, 50, 33, 34, -48, -1};
        Plane Surface(56) = {55};
        Line Loop(57) = {49, 51, -23, -7};
        Plane Surface(58) = {57};
        Line Loop(59) = {39, 5, -10, -53, -52, -41, -40};
        Plane Surface(60) = {59};
        Line Loop(61) = {53, -9, -25, -54};
        Plane Surface(62) = {61};
        Line Loop(63) = {50, -32, -31, -47, -44, -43, 52, 54, -24, -51};
        Plane Surface(64) = {63};
    EndIf

    Line Loop(85) = {2, -28, -27, -26, -48};
    Plane Surface(86) = {85};
    Line Loop(87) = {26, 15, 17, 34};
    Plane Surface(88) = {87};
    Line Loop(89) = {15, 18, -29, -27};
    Plane Surface(80) = {89};
    Line Loop(81) = {18, -30, 31, -16};
    Plane Surface(82) = {81};
    Line Loop(83) = {16, 32, 33, -17};
    Plane Surface(84) = {83};
    Line Loop(65) = {28, 3, 36, -46, -45, 47, 30, -29};
    Plane Surface(66) = {65};
    Line Loop(67) = {36, 37, 38, 39, -4};
    Plane Surface(68) = {67};
    Line Loop(69) = {37, 19, -21, 46};
    Plane Surface(70) = {69};
    Line Loop(71) = {19, 22, -40, -38};
    Plane Surface(72) = {71};
    Line Loop(73) = {22, 41, 43, -20};
    Plane Surface(74) = {73};
    Line Loop(75) = {20, 44, 45, 21};
    Plane Surface(76) = {75};
    
    Physical Surface(2) = {88, 80, 82, 84};
    Physical Surface(5) = {70, 72, 76, 74};

    If(slab_w_full == 1)
        Physical Surface(1) = {50, 86, 66, 68, 54};
        Physical Surface(3) = {52};
    EndIf
    If(slab_w_full == 0)
        Physical Surface(1) = {56, 86, 66, 68, 60};
        Physical Surface(3) = {64};
        Physical Surface(4) = {58, 62};
    EndIf
EndIf
