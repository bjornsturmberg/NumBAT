// Template mesh geometry file for a single inclusion on two slabs.
// Inclusion can be circular, elliptical square, or rectangular.

d = 1; // grating period
ff = 0;
d_in_nm = 100;
dy_in_nm = 50;
dy = dy_in_nm/d_in_nm;
a1 = 20;
a1y = 10;
radius1 = (a1/(2*d_in_nm))*d;
radius1y = (a1y/(2*d_in_nm))*d;

rect = 1;

slab_width = 80;
slab_height = 5;
slab_w = slab_width/d_in_nm;
slab_h = slab_height/d_in_nm;
slab_w_full = 0;
If(slab_w == 1)
    slab_w_full = 1;
EndIf

slab2_width = 100;
slab2_height = 10;
slab2_w = slab2_width/d_in_nm;
slab2_h = slab2_height/d_in_nm;
slab2_w_full = 0;
If(slab2_w == 1)
    slab2_w_full = 1;
EndIf

lc = 0; // 0.501 0.201 0.0701;
lc2 = lc/1; // on cylinder surfaces
lc3 = lc/1; // cylinder1 centres
lc4 = lc/1; // centres of top and bottom
lc5 = lc/1; // slab

hy = dy; // Thickness: square profile => hy=d
hx = 0.;


Point(1) = {0, 0, 0, lc};
Point(2) = {-hx, -hy, 0, lc};
Point(3) = {-hx+d, -hy, 0, lc};
Point(4) = {d, 0, 0,lc};

// Slab
Point(250) = {d/2-slab_w/2, -hy+slab_h+slab2_h, 0, lc5};
Point(251) = {d/2+slab_w/2, -hy+slab_h+slab2_h, 0, lc5};
Point(350) = {d/2-slab2_w/2, -hy+slab2_h, 0, lc5};
Point(351) = {d/2+slab2_w/2, -hy+slab2_h, 0, lc5};

// Inclusion
Point(5) = {-hx+d/2, -hy+radius1y+slab_h+slab2_h, 0, lc3};
Point(6) = {-hx+d/2, -hy+2*radius1y+slab_h+slab2_h, 0, lc2};
Point(7) = {-hx+d/2-radius1, -hy+radius1y+slab_h+slab2_h, 0, lc2};
Point(8) = {-hx+d/2, -hy+slab_h+slab2_h, 0, lc2};
Point(9) = {-hx+d/2+radius1, -hy+radius1y+slab_h+slab2_h, 0, lc2};
Point(38) = {-hx+d/2, -hy+slab2_h, 0, lc2};

Point(10) = {-hx+d/2, 0, 0, lc4};
Point(11) = {0,-hy+radius1y+slab_h+slab2_h, 0, lc};
Point(12) = {-hx+d/2, -hy, 0, lc4};
Point(13) = {d, -hy+radius1y+slab_h+slab2_h, 0, lc};
Line(1) = {1,10};
Line(2) = {10,4};
Line(5) = {1,11};
Line(7) = {4,13};
Line(9) = {11,7};
Line(10) = {7,5};
Line(11) = {5,9};
Line(12) = {9,13};
Line(13) = {10,6};
Line(14) = {6,5};
Line(15) = {5,8};


If(rect == 0)
    Ellipsis(17) = {9,5,6,6};
    Ellipsis(18) = {6,5,7,7};
    Ellipsis(19) = {7,5,8,8};
    Ellipsis(20) = {8,5,9,9};

    If(slab_w_full == 1)
        Line(6) = {11,250};
        Line(8) = {13,251};
        Line(21) = {2, 250};
        Line(22) = {250, 8};
        Line(23) = {8, 251};
        Line(24) = {251, 3};
        Line(25) = {3, 12};
        Line(26) = {12, 2};

    EndIf

    If(slab_w_full == 0)
        Line(21) = {252, 250};
        Line(22) = {250, 8};
        Line(23) = {8, 251};
        Line(24) = {251, 253};
        Line(25) = {253, 12};
        Line(26) = {12, 252};


    EndIf
EndIf


If(rect == 1)
    Point(150) = {-hx+d/2+radius1, -hy+slab_h+slab2_h+2*radius1y, 0,lc3};
    Point(151) = {-hx+d/2-radius1, -hy+slab_h+slab2_h+2*radius1y, 0,lc3};
    Point(152) = {-hx+d/2+radius1, -hy+slab_h+slab2_h, 0,lc3};
    Point(153) = {-hx+d/2-radius1, -hy+slab_h+slab2_h, 0,lc3};

    If(slab_w_full == 1)
        Line(6) = {11,250};
        Line(8) = {13,251};
        Line(27) = {153, 8};
        Line(28) = {8, 152};
        Line(29) = {153, 7};
        Line(30) = {7, 151};
        Line(31) = {151, 6};
        Line(32) = {6, 150};
        Line(33) = {150, 9};
        Line(34) = {9, 152};

        If(slab2_w_full == 1)
            Line(22) = {250, 153};
            Line(23) = {152, 251};     
            Line(25) = {3, 12};
            Line(26) = {12, 2};   
            Line(35) = {350, 2};
            Line(36) = {350, 250};
            Line(37) = {8, 38};
            Line(38) = {38, 12};
            Line(39) = {351, 3};
            Line(40) = {351, 251};
            Line(41) = {350, 38};
            Line(42) = {38, 351};

            Line Loop(43) = {5, 9, 30, 31, -13, -1};
            Plane Surface(44) = {43};
            Line Loop(45) = {13, 32, 33, 12, -7, -2};
            Plane Surface(46) = {45};
            Line Loop(47) = {12, 8, -23, -34};
            Plane Surface(48) = {47};
            Line Loop(49) = {29, -9, 6, 22};
            Plane Surface(50) = {49};
            Line Loop(51) = {29, 10, 15, -27};
            Plane Surface(52) = {51};
            Line Loop(53) = {15, 28, -34, -11};
            Plane Surface(54) = {53};
            Line Loop(55) = {11, -33, -32, 14};
            Plane Surface(56) = {55};
            Line Loop(57) = {31, 14, -10, 30};
            Plane Surface(58) = {57};
            Line Loop(59) = {22, 27, 37, -41, 36};
            Line Loop(60) = {37, 42, 40, -23, -28};
            Plane Surface(61) = {60};
            Plane Surface(62) = {59};
            Line Loop(63) = {41, 38, 26, -35};
            Plane Surface(64) = {63};
            Line Loop(65) = {38, -25, -39, -42};
            Plane Surface(66) = {65};

            Physical Line(67) = {5, 6, 36, 35};
            Physical Line(68) = {26, 25};
            Physical Line(69) = {39, 40, 8, 7};
            Physical Line(70) = {2, 1};

            Physical Surface(1) = {44, 46, 48, 50};
            Physical Surface(2) = {58, 56, 54, 52};
            Physical Surface(3) = {62, 61};
            Physical Surface(4) = {66, 64};
        EndIf

        If(slab2_w_full == 0)
            Point(352) = {d/2-slab2_w/2, -hy+slab_h+slab2_h, 0, lc5};
            Point(353) = {d/2+slab2_w/2, -hy+slab_h+slab2_h, 0, lc5};
            Point(354) = {d/2-slab2_w/2, -hy, 0, lc5};
            Point(355) = {d/2+slab2_w/2, -hy, 0, lc5};
            Point(356) = {d/2-slab_w/2, -hy+slab2_h, 0, lc5};
            Point(357) = {d/2+slab_w/2, -hy+slab2_h, 0, lc5};

            Line(37) = {8, 38};
            Line(38) = {38, 12};
            Line(41) = {350, 38};
            Line(42) = {38, 351};
            Line(43) = {356, 250};
            Line(44) = {356, 2};
            Line(45) = {2, 354};
            Line(46) = {354, 12};
            Line(47) = {12, 355};
            Line(48) = {355, 3};
            Line(49) = {3, 357};
            Line(50) = {357, 251};
            Line(51) = {351, 355};
            Line(52) = {351, 357};
            Line(53) = {353, 251};
            Line(54) = {353, 152};
            Line(55) = {153, 352};
            Line(56) = {352, 250};
            Line(57) = {356, 350};
            Line(58) = {350, 354};
            
            Line Loop(43) = {5, 9, 30, 31, -13, -1};
            Plane Surface(44) = {43};
            Line Loop(45) = {13, 32, 33, 12, -7, -2};
            Plane Surface(46) = {45};
            Line Loop(59) = {12, 8, -53, 54, -34};
            Plane Surface(60) = {59};
            Line Loop(61) = {29, -9, 6, -56, -55};
            Plane Surface(62) = {61};
            Line Loop(63) = {29, 10, 15, -27};
            Plane Surface(64) = {63};
            Line Loop(65) = {15, 28, -34, -11};
            Plane Surface(66) = {65};
            Line Loop(67) = {33, -11, -14, 32};
            Plane Surface(68) = {67};
            Line Loop(69) = {14, -10, 30, 31};
            Plane Surface(70) = {69};
            Line Loop(71) = {28, -54, 53, -50, -52, -42, -37};
            Plane Surface(72) = {71};
            Line Loop(73) = {37, -41, -57, 43, -56, -55, 27};
            Plane Surface(74) = {73};
            Line Loop(75) = {41, 38, -46, -58};
            Plane Surface(76) = {75};
            Line Loop(77) = {38, 47, -51, -42};
            Plane Surface(78) = {77};
            Line Loop(79) = {51, 48, 49, -52};
            Plane Surface(80) = {79};
            Line Loop(81) = {58, -45, -44, 57};
            Plane Surface(82) = {81};

            Physical Line(83) = {5, 6, 43, 44};
            Physical Line(84) = {45, 46, 47, 48};
            Physical Line(85) = {49, 50, 8, 7};
            Physical Line(86) = {2, 1};

            Physical Surface(1) = {44, 62, 60, 46};
            Physical Surface(2) = {68, 70, 64, 66};
            Physical Surface(3) = {72, 74};
            Physical Surface(4) = {76, 78};
            Physical Surface(5) = {80, 82};
        EndIf
    EndIf

    If(slab_w_full == 0)
        Point(352) = {d/2-slab2_w/2, -hy+slab_h+slab2_h, 0, lc5};
        Point(353) = {d/2+slab2_w/2, -hy+slab_h+slab2_h, 0, lc5};
        Point(354) = {d/2-slab2_w/2, -hy, 0, lc5};
        Point(355) = {d/2+slab2_w/2, -hy, 0, lc5};
        Point(356) = {d/2-slab_w/2, -hy+slab2_h, 0, lc5};
        Point(357) = {d/2+slab_w/2, -hy+slab2_h, 0, lc5};

        Line(17) = {151, 6};
        Line(18) = {6, 150};
        Line(19) = {150, 9};
        Line(20) = {9, 152};
        Line(21) = {152, 8};
        Line(22) = {8, 153};
        Line(23) = {153, 7};
        Line(24) = {7, 151};
        Line(26) = {153, 250};
        Line(30) = {251, 152};

        If(slab2_w_full == 1)
            Line(31) = {11, 352};
            Line(32) = {352, 250};
            Line(33) = {352, 350};
            Line(34) = {350, 356};
            Line(35) = {350, 2};
            Line(36) = {2, 12};
            Line(37) = {12, 3};
            Line(38) = {3, 351};
            Line(39) = {351, 353};
            Line(40) = {353, 13};
            Line(41) = {353, 251};
            Line(42) = {251, 357};
            Line(43) = {357, 351};
            Line(44) = {357, 38};
            Line(45) = {38, 356};
            Line(46) = {250, 356};
            Line(47) = {38, 12};
            Line(48) = {8, 38};
            
            Line Loop(49) = {5, 9, 24, 17, -13, -1};
            Plane Surface(50) = {49};
            Line Loop(51) = {13, 18, 19, 12, -7, -2};
            Plane Surface(52) = {51};
            Line Loop(53) = {12, -40, 41, 30, -20};
            Plane Surface(54) = {53};
            Line Loop(55) = {23, -9, 31, 32, -26};
            Plane Surface(56) = {55};
            Line Loop(57) = {17, 14, -10, 24};
            Plane Surface(58) = {57};
            Line Loop(59) = {18, 19, -11, -14};
            Line Loop(60) = {23, 10, 15, 22};
            Plane Surface(61) = {60};
            Line Loop(62) = {15, -21, -20, -11};
            Plane Surface(63) = {62};
            Plane Surface(64) = {59};
            Line Loop(65) = {30, 21, 48, -44, -42};
            Plane Surface(66) = {65};
            Line Loop(67) = {48, 45, -46, -26, -22};
            Plane Surface(68) = {67};
            Line Loop(69) = {41, 42, 43, 39};
            Plane Surface(70) = {69};
            Line Loop(71) = {46, -34, -33, 32};
            Plane Surface(72) = {71};
            Line Loop(73) = {34, -45, 47, -36, -35};
            Plane Surface(74) = {73};
            Line Loop(75) = {47, 37, 38, -43, 44};
            Plane Surface(76) = {75};

            Physical Line(77) = {5, 31, 33, 35};
            Physical Line(78) = {36, 37};
            Physical Line(79) = {7, 40, 39, 38};
            Physical Line(80) = {2, 1};

            Physical Surface(1) = {50, 52, 54, 56};
            Physical Surface(2) = {58, 64, 63, 61};
            Physical Surface(3) = {66, 68};
            Physical Surface(4) = {72, 70};
            Physical Surface(5) = {76, 74};
        EndIf
    EndIf
EndIf



