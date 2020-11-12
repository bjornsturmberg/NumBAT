// Template mesh geometry file for a single suspended inclusion.
// Inclusion can be circular/elliptical (default), or square/rectangular.

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
a2 = 10;
a2y = 20;
radius2 = (a2/(2*d_in_nm))*d;
radius2y = (a2y/(2*d_in_nm))*d;
sep = 10;
b = sep/(2*d_in_nm);
yoff = -5;
yoffset = yoff/d_in_nm;
rect = 1;
lc = 0; // 0.501 0.201 0.0701;
lc_refine_1 = lc/1; // on cylinder surfaces
lc_refine_2 = lc/1; // cylinder centres

hy = dy; // Thickness: square profile => hy=d
hx = 0.;


Point(1) = {0, 0, 0, lc};
Point(2) = {-hx, -hy, 0, lc};
Point(3) = {-hx+d, -hy, 0, lc};
Point(4) = {d, 0, 0,lc};

// Vertices
Point(5) = {-radius1-b-hx+d/2, -hy/2, 0, lc_refine_2};
Point(6) = {-radius1-b-hx+d/2, -hy/2+radius1y, 0, lc_refine_1};
Point(7) = {-radius1-b-hx+d/2-radius1, -hy/2, 0, lc_refine_1};
Point(8) = {-radius1-b-hx+d/2, -hy/2-radius1y, 0, lc_refine_1};
Point(9) = {-radius1-b-hx+d/2+radius1, -hy/2, 0, lc_refine_1};

Point(10) = {-radius1-b-hx+d/2, 0, 0, lc_refine_1};
Point(11) = {0,-hy/2, 0, lc};
Point(12) = {-radius1-b-hx+d/2, -hy, 0, lc_refine_1};
Point(13) = {d, -hy/2+yoffset, 0, lc};

Point(14) = {-hx+d/2, -hy/2, 0, lc_refine_1};
Point(15) = {-hx+d/2, 0, 0, lc_refine_1};
Point(16) = {-hx+d/2, -hy, 0, lc_refine_1};
If(yoffset < 0)
    Point(17) = {-hx+d/2, -hy/2+yoffset, 0, lc_refine_1};
EndIf
If(yoffset > 0)
    Point(17) = {-hx+d/2, -hy/2+yoffset, 0, lc_refine_1};
EndIf

Point(18) = {radius2+b-hx+d/2, -hy/2+yoffset, 0, lc_refine_2};
Point(19) = {radius2+b-hx+d/2, -hy/2+radius2y+yoffset, 0, lc_refine_1};
Point(20) = {radius2+b-hx+d/2-radius2, -hy/2+yoffset, 0, lc_refine_1};
Point(21) = {radius2+b-hx+d/2, -hy/2-radius2y+yoffset, 0, lc_refine_1};
Point(22) = {radius2+b-hx+d/2+radius2, -hy/2+yoffset, 0, lc_refine_1};
Point(23) = {radius2+b-hx+d/2, 0, 0, lc_refine_1};
Point(24) = {radius2+b-hx+d/2, -hy, 0, lc_refine_1};

Line(1) = {1,10};
Line(2) = {10,15};
Line(3) = {2,12};
Line(4) = {12,16};
Line(5) = {1,11};
Line(6) = {11,2};
Line(7) = {4,13};
Line(8) = {13,3};
Line(9) = {11,7};
Line(10) = {7,5};
Line(11) = {5,9};
Line(12) = {9,14};
Line(13) = {10,6};
Line(14) = {6,5};
Line(15) = {5,8};
Line(16) = {8,12};
Line(21) = {15, 23};
Line(22) = {23, 4};
Line(23) = {3, 24};
Line(24) = {24, 16};
Line(25) = {21, 24};
Line(29) = {23, 19};
Line(30) = {13, 22};
Line(32) = {20, 18};
Line(33) = {18, 22};
Line(34) = {19, 18};
Line(35) = {21, 18};
If(yoffset < 0)
    Line(26) = {14, 17};
    Line(27) = {17, 16};
    Line(28) = {14, 15};
    Line(31) = {20, 17};
EndIf
If(yoffset > 0)
    Line(26) = {14, 17};
    Line(27) = {14, 16};
    Line(28) = {17, 15};
    Line(31) = {20, 17};
EndIf
If(yoffset == 0)
    Line(26) = {14, 15};
    Line(27) = {14, 16};
    Line(31) = {20, 14};
EndIf

Physical Line(80) = {1, 2, 21, 22};
Physical Line(81) = {7, 8};
Physical Line(82) = {23, 24, 4, 3};
Physical Line(83) = {6, 5};


If(rect == 0)
    Ellipsis(17) = {9,5,6,6};
    Ellipsis(18) = {6,5,7,7};
    Ellipsis(19) = {7,5,8,8};
    Ellipsis(20) = {8,5,9,9};

    Ellipsis(36) = {22, 18, 19, 19};
    Ellipsis(37) = {22, 18, 21, 21};
    Ellipsis(38) = {20, 18, 19, 19};
    Ellipsis(39) = {20, 18, 21, 21};

    Line Loop(40) = {5, 9, -18, -13, -1};
    Plane Surface(41) = {40};
    Line Loop(42) = {9, 19, 16, -3, -6};
    Plane Surface(43) = {42};
    Line Loop(52) = {25, -23, -8, 30, 37};
    Plane Surface(53) = {52};
    Line Loop(54) = {36, -29, 22, 7, 30};
    Plane Surface(55) = {54};
    Line Loop(56) = {36, 34, 33};
    Plane Surface(57) = {56};
    Line Loop(58) = {33, 37, 35};
    Plane Surface(59) = {58};
    Line Loop(60) = {35, -32, 39};
    Plane Surface(61) = {60};
    Line Loop(62) = {32, -34, -38};
    Plane Surface(63) = {62};
    Line Loop(64) = {17, 14, 11};
    Plane Surface(65) = {64};
    Line Loop(66) = {14, -10, -18};
    Plane Surface(67) = {66};
    Line Loop(68) = {10, 15, -19};
    Plane Surface(69) = {68};
    Line Loop(70) = {15, 20, -11};
    Plane Surface(71) = {70};
    If(yoffset < 0)
        Line Loop(72) = {16, 4, -27, -26, -12, -20};
        Plane Surface(73) = {44};
        Line Loop(74) = {12, 28, -2, 13, -17};
        Plane Surface(75) = {46};
        Line Loop(76) = {28, 21, 29, -38, 31, -26};
        Plane Surface(77) = {48};
        Line Loop(78) = {31, 27, -24, -25, -39};
        Plane Surface(79) = {50};
    EndIf
    If(yoffset > 0)
        Line Loop(72) = {13, -17, 12, 26, 28, -2};
        Plane Surface(73) = {72};
        Line Loop(74) = {12, 27, -4, -16, 20};
        Plane Surface(75) = {74};
        Line Loop(76) = {27, -24, -25, -39, 31, -26};
        Plane Surface(77) = {76};
        Line Loop(78) = {28, 21, 29, -38, 31};
        Plane Surface(79) = {78};
    EndIf
    If(yoffset == 0)
        Line Loop(72) = {13, -17, 12, 26, -2};
        Plane Surface(73) = {84};
        Line Loop(74) = {12, 27, -4, -16, 20};
        Plane Surface(75) = {86};
        Line Loop(76) = {27, -24, -25, -39, 31};
        Plane Surface(77) = {88};
        Line Loop(78) = {31, 26, 21, 29, -38};
        Plane Surface(79) = {90};
    EndIf

    Physical Surface(1) = {41, 73, 79, 55, 53, 77, 75, 43};
    Physical Surface(2) = {67, 65, 71, 69};
    Physical Surface(3) = {63, 57, 59, 61};
EndIf

If(rect == 1)
    Point(150) = {-radius1-b-hx+d/2+radius1, -hy/2+radius1y, 0,lc_refine_2};
    Point(151) = {-radius1-b-hx+d/2-radius1, -hy/2+radius1y, 0,lc_refine_2};
    Point(152) = {-radius1-b-hx+d/2+radius1, -hy/2-radius1y, 0,lc_refine_2};
    Point(153) = {-radius1-b-hx+d/2-radius1, -hy/2-radius1y, 0,lc_refine_2};

    Point(154) = {radius2+b-hx+d/2+radius2, -hy/2+radius2y+yoffset, 0,lc_refine_2};
    Point(155) = {radius2+b-hx+d/2-radius2, -hy/2+radius2y+yoffset, 0,lc_refine_2};
    Point(156) = {radius2+b-hx+d/2+radius2, -hy/2-radius2y+yoffset, 0,lc_refine_2};
    Point(157) = {radius2+b-hx+d/2-radius2, -hy/2-radius2y+yoffset, 0,lc_refine_2};

    Line(17) = {151, 6};
    Line(18) = {6, 150};
    Line(19) = {150, 9};
    Line(20) = {9, 152};
    Line(41) = {152, 8};
    Line(42) = {8, 153};
    Line(43) = {153, 7};
    Line(44) = {7, 151};
    Line(84) = {20, 155};
    Line(85) = {155, 19};
    Line(86) = {19, 154};
    Line(87) = {154, 22};
    Line(88) = {156, 22};
    Line(89) = {20, 157};
    Line(90) = {157, 21};
    Line(91) = {21, 156};
    Line Loop(92) = {5, 9, 44, 17, -13, -1};
    Plane Surface(93) = {92};
    Line Loop(94) = {9, -43, -42, 16, -3, -6};
    Plane Surface(95) = {94};
    Line Loop(96) = {29, 86, 87, -30, -7, -22};
    Plane Surface(97) = {96};
    Line Loop(98) = {30, -88, -91, 25, -23, -8};
    Plane Surface(99) = {98};
    Line Loop(100) = {17, 14, -10, 44};
    Plane Surface(101) = {100};
    Line Loop(102) = {14, 11, -19, -18};
    Plane Surface(103) = {102};
    Line Loop(104) = {11, 20, 41, -15};
    Plane Surface(105) = {104};
    Line Loop(106) = {15, 42, 43, 10};
    Plane Surface(107) = {106};
    Line Loop(108) = {84, 85, 34, -32};
    Plane Surface(109) = {108};
    Line Loop(110) = {32, -35, -90, -89};
    Plane Surface(111) = {110};
    Line Loop(112) = {35, 33, -88, -91};
    Plane Surface(113) = {112};
    Line Loop(114) = {87, -33, -34, 86};
    Plane Surface(115) = {114};

    If(yoffset < 0)
        Line Loop(116) = {2, -28, -12, -19, -18, -13};
        Plane Surface(117) = {116};
        Line Loop(118) = {12, 26, 27, -4, -16, -41, -20};
        Plane Surface(119) = {118};
        Line Loop(120) = {31, 27, -24, -25, -90, -89};
        Plane Surface(121) = {120};
        Line Loop(122) = {26, -31, 84, 85, -29, -21, -28};
        Plane Surface(123) = {122};
    EndIf
    If(yoffset > 0)
        Line Loop(116) = {13, 18, 19, 12, 26, 28, -2};
        Plane Surface(117) = {116};
        Line Loop(118) = {12, 27, -4, -16, -41, -20};
        Plane Surface(119) = {118};
        Line Loop(120) = {27, -24, -25, -90, -89, 31, -26};
        Plane Surface(121) = {120};
        Line Loop(122) = {28, 21, 29, -85, -84, 31};
        Plane Surface(123) = {122};
    EndIf
    If(yoffset == 0)
        Line Loop(116) = {13, 18, 19, 12, 26, -2};
        Plane Surface(117) = {116};
        Line Loop(118) = {12, 27, -4, -16, -41, -20};
        Plane Surface(119) = {118};
        Line Loop(120) = {27, -24, -25, -90, -89, 31};
        Plane Surface(121) = {120};
        Line Loop(122) = {31, 26, 21, 29, -85, -84};
        Plane Surface(123) = {122};
    EndIf

    Physical Surface(1) = {93, 117, 123, 97, 99, 121, 119, 95};
    Physical Surface(2) = {101, 103, 105, 107};
    Physical Surface(3) = {109, 115, 113, 111};
EndIf
