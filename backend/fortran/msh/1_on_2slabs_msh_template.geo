// Template mesh geometry file for a single inclusion on two slabs.
// Inclusion can be circular/elliptical (default), or square/rectangular.

d = 1; // grating period
d_in_nm = 100;
dy_in_nm = 50;
dy = dy_in_nm/d_in_nm;
a1 = 20;
a1y = 10;
radius1 = (a1/(2*d_in_nm))*d;
radius1y = (a1y/(2*d_in_nm))*d;

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
Point(38) = {-hx+d/2, -hy+slab2_h, 0, lc_refine_4};
Point(250) = {d/2-slab_w/2, -hy+slab_h+slab2_h, 0, lc_refine_4};
Point(251) = {d/2+slab_w/2, -hy+slab_h+slab2_h, 0, lc_refine_4};
Point(350) = {d/2-slab2_w/2, -hy+slab2_h, 0, lc_refine_4};
Point(351) = {d/2+slab2_w/2, -hy+slab2_h, 0, lc_refine_4};

// Inclusion
Point(5) = {-hx+d/2, -hy+radius1y+slab_h+slab2_h, 0, lc_refine_2};
Point(6) = {-hx+d/2, -hy+2*radius1y+slab_h+slab2_h, 0, lc_refine_1};
Point(7) = {-hx+d/2-radius1, -hy+radius1y+slab_h+slab2_h, 0, lc_refine_1};
Point(8) = {-hx+d/2, -hy+slab_h+slab2_h, 0, lc_refine_1};
Point(9) = {-hx+d/2+radius1, -hy+radius1y+slab_h+slab2_h, 0, lc_refine_1};

Point(10) = {-hx+d/2, 0, 0, lc_refine_3};
Point(11) = {0,-hy+radius1y+slab_h+slab2_h, 0, lc};
Point(12) = {-hx+d/2, -hy, 0, lc_refine_3};
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
    Line Loop(45) = {5, 9, -18, -13, -1};
    Plane Surface(46) = {45};
    Line Loop(47) = {13, -17, 12, -7, -2};
    Plane Surface(48) = {47};
    Line Loop(49) = {17, 14, 11};
    Plane Surface(50) = {49};
    Line Loop(51) = {20, -11, 15};
    Plane Surface(52) = {51};
    Line Loop(53) = {15, -19, 10};
    Plane Surface(54) = {53};
    Line Loop(55) = {10, -14, 18};
    Plane Surface(56) = {55};
    Physical Line(57) = {1, 2};

    Physical Surface(2) = {56, 50, 52, 54};

    If(slab_w_full == 1)
        Line(6) = {11,250};
        Line(8) = {13,251};
        Line(43) = {250, 8};
        Line(44) = {8, 251};

        If(slab2_w_full == 1)    
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

            Line Loop(58) = {9, 19, -43, -6};
            Plane Surface(59) = {58};
            Line Loop(60) = {20, 12, 8, -44};
            Plane Surface(61) = {60};
            Line Loop(62) = {44, -40, -42, -37};
            Plane Surface(63) = {62};
            Line Loop(64) = {37, -41, 36, 43};
            Plane Surface(65) = {64};
            Line Loop(66) = {41, 38, 26, -35};
            Plane Surface(67) = {66};
            Line Loop(68) = {42, 39, 25, -38};
            Plane Surface(69) = {68};

            Physical Line(70) = {5, 6, 36, 35};
            Physical Line(71) = {26, 25};
            Physical Line(72) = {7, 8, 40, 39};

            Physical Surface(1) = {46, 48, 61, 59};
            Physical Surface(3) = {65, 63};
            Physical Surface(4) = {67, 69};
        EndIf

        If(slab2_w_full == 0)  
            Point(352) = {0, -hy+slab2_h, 0, lc_refine_4};
            Point(353) = {d, -hy+slab2_h, 0, lc_refine_4};  
            Point(354) = {d/2-slab2_w/2, -hy, 0, lc_refine_4};
            Point(355) = {d/2+slab2_w/2, -hy, 0, lc_refine_4};  

            Line(37) = {8, 38};
            Line(38) = {38, 12};
            Line(41) = {350, 38};
            Line(42) = {38, 351};
            Line(58) = {250, 352};
            Line(59) = {352, 350};
            Line(60) = {352, 2};
            Line(61) = {2, 354};
            Line(62) = {354, 12};
            Line(63) = {12, 355};
            Line(64) = {355, 3};
            Line(65) = {3, 353};
            Line(66) = {353, 251};
            Line(67) = {351, 353};
            Line(68) = {351, 355};
            Line(69) = {350, 354};

            Line Loop(70) = {19, -43, -6, 9};
            Plane Surface(71) = {70};
            Line Loop(72) = {20, 12, 8, -44};
            Plane Surface(73) = {72};
            Line Loop(74) = {44, -66, -67, -42, -37};
            Plane Surface(75) = {74};
            Line Loop(76) = {43, 37, -41, -59, -58};
            Plane Surface(77) = {76};
            Line Loop(78) = {59, 69, -61, -60};
            Plane Surface(79) = {78};
            Line Loop(80) = {69, 62, -38, -41};
            Plane Surface(81) = {80};
            Line Loop(82) = {38, 63, -68, -42};
            Plane Surface(83) = {82};
            Line Loop(84) = {68, 64, 65, -67};
            Plane Surface(85) = {84};

            Physical Line(86) = {5, 6, 58, 60};
            Physical Line(87) = {61, 62, 63, 64};
            Physical Line(88) = {7, 8, 66, 65};

            Physical Surface(1) = {71, 46, 48, 73};
            Physical Surface(3) = {77, 75};
            Physical Surface(4) = {81, 83};
            Physical Surface(5) = {79, 85};
        EndIf
    EndIf

    If(slab_w_full == 0)
        Point(252) = {0, -hy+slab_h+slab2_h, 0, lc_refine_4};
        Point(253) = {d, -hy+slab_h+slab2_h, 0, lc_refine_4};
        Point(352) = {d/2-slab_w/2, -hy+slab2_h, 0, lc_refine_4};
        Point(353) = {d/2+slab_w/2, -hy+slab2_h, 0, lc_refine_4};

        Line(22) = {250, 8};
        Line(23) = {8, 251};

        If(slab2_w_full == 1)  
            Line(58) = {11, 252};
            Line(59) = {252, 250};
            Line(60) = {251, 253};
            Line(61) = {13, 253};
            Line(62) = {253, 351};
            Line(63) = {351, 3};
            Line(64) = {3, 12};
            Line(65) = {12, 2};
            Line(66) = {2, 350};
            Line(67) = {350, 252};
            Line(68) = {251, 353};
            Line(69) = {353, 351};
            Line(70) = {38, 353};
            Line(71) = {8, 38};
            Line(72) = {38, 352};
            Line(73) = {352, 250};
            Line(74) = {352, 350};
            Line(75) = {38, 12};

            Line Loop(76) = {9, 19, -22, -59, -58};
            Plane Surface(77) = {76};
            Line Loop(78) = {20, 12, 61, -60, -23};
            Plane Surface(79) = {78};
            Line Loop(80) = {68, 69, -62, -60};
            Plane Surface(81) = {80};
            Line Loop(82) = {68, -70, -71, 23};
            Plane Surface(83) = {82};
            Line Loop(84) = {71, 72, 73, 22};
            Plane Surface(85) = {84};
            Line Loop(86) = {59, -73, 74, 67};
            Plane Surface(87) = {86};
            Line Loop(88) = {74, -66, -65, -75, 72};
            Plane Surface(89) = {88};
            Line Loop(90) = {75, -64, -63, -69, -70};
            Plane Surface(91) = {90};

            Physical Line(92) = {5, 58, 67, 66};
            Physical Line(93) = {65, 64};
            Physical Line(94) = {63, 62, 61, 7};

            Physical Surface(1) = {77, 46, 48, 79};
            Physical Surface(3) = {85, 83};
            Physical Surface(4) = {87, 81};
            Physical Surface(5) = {91, 89};
        EndIf

        If(slab2_w_full == 0)
            Point(354) = {d/2-slab2_w/2, -hy, 0, lc_refine_4};
            Point(355) = {d/2+slab2_w/2, -hy, 0, lc_refine_4}; 
            Point(356) = {0, -hy+slab2_h, 0, lc_refine_4};
            Point(357) = {d, -hy+slab2_h, 0, lc_refine_4};   

            If(slab2_w > slab_w)
                Line(58) = {11, 252};
                Line(59) = {252, 250};
                Line(60) = {250, 352};
                Line(61) = {352, 350};
                Line(62) = {350, 356};
                Line(63) = {252, 356};
                Line(64) = {356, 2};
                Line(65) = {2, 354};
                Line(66) = {354, 350};
                Line(67) = {354, 12};
                Line(68) = {8, 38};
                Line(69) = {38, 12};
                Line(70) = {12, 355};
                Line(71) = {355, 351};
                Line(72) = {351, 357};
                Line(73) = {357, 3};
                Line(74) = {3, 355};
                Line(75) = {353, 351};
                Line(76) = {251, 253};
                Line(77) = {13, 253};
                Line(78) = {253, 357};
                Line(79) = {251, 353};
                Line(80) = {353, 38};
                Line(81) = {352, 38};

                Line Loop(82) = {9, 19, -22, -59, -58};
                Plane Surface(83) = {82};
                Line Loop(84) = {20, 12, 77, -76, -23};
                Plane Surface(85) = {84};
                Line Loop(86) = {23, 79, 80, -68};
                Plane Surface(87) = {86};
                Line Loop(88) = {79, 75, 72, -78, -76};
                Plane Surface(89) = {88};
                Line Loop(90) = {71, 72, 73, 74};
                Plane Surface(91) = {90};
                Line Loop(92) = {80, 69, 70, 71, -75};
                Plane Surface(93) = {92};
                Line Loop(94) = {81, -68, -22, 60};
                Plane Surface(95) = {94};
                Line Loop(96) = {60, 61, 62, -63, 59};
                Plane Surface(97) = {96};
                Line Loop(98) = {62, 64, 65, 66};
                Plane Surface(99) = {98};
                Line Loop(100) = {66, -61, 81, 69, -67};
                Plane Surface(101) = {100};

                Physical Line(102) = {5, 58, 63, 64};
                Physical Line(103) = {65, 67, 70, 74};
                Physical Line(104) = {7, 77, 78, 73};

                Physical Surface(1) = {46, 83, 48, 85};
                Physical Surface(3) = {87, 95};
                Physical Surface(4) = {97, 89};
                Physical Surface(5) = {101, 93};
                Physical Surface(6) = {99, 91};
            EndIf

            If(slab2_w < slab_w)
                Line(58) = {11, 252};
                Line(59) = {252, 250};
                Line(60) = {250, 352};
                Line(61) = {352, 356};
                Line(62) = {356, 252};
                Line(63) = {356, 2};
                Line(64) = {2, 354};
                Line(65) = {354, 350};
                Line(66) = {350, 38};
                Line(67) = {38, 8};
                Line(68) = {352, 350};
                Line(69) = {354, 12};
                Line(70) = {12, 355};
                Line(71) = {355, 3};
                Line(72) = {253, 357};
                Line(73) = {357, 3};
                Line(74) = {253, 13};
                Line(75) = {253, 251};
                Line(76) = {353, 357};
                Line(77) = {353, 251};
                Line(78) = {353, 351};
                Line(79) = {351, 38};
                Line(80) = {351, 355};
                Line(81) = {38, 12};

                Line Loop(82) = {9, 19, -22, -59, -58};
                Plane Surface(83) = {82};
                Line Loop(84) = {20, 12, -74, 75, -23};
                Plane Surface(85) = {84};
                Line Loop(86) = {75, -77, 76, -72};
                Plane Surface(87) = {86};
                Line Loop(88) = {77, -23, -67, -79, -78};
                Plane Surface(89) = {88};
                Line Loop(90) = {80, -70, -81, -79};
                Plane Surface(91) = {90};
                Line Loop(92) = {78, 80, 71, -73, -76};
                Plane Surface(93) = {92};
                Line Loop(94) = {81, -69, 65, 66};
                Plane Surface(95) = {94};
                Line Loop(96) = {67, -22, 60, 68, 66};
                Plane Surface(97) = {96};
                Line Loop(98) = {60, 61, 62, 59};
                Plane Surface(99) = {98};
                Line Loop(100) = {61, 63, 64, 65, -68};
                Plane Surface(101) = {100};

                Physical Line(102) = {5, 58, 62, 63};
                Physical Line(103) = {64, 69, 70, 71};
                Physical Line(104) = {7, 74, 72, 73};

                Physical Surface(1) = {46, 48, 85, 83};
                Physical Surface(3) = {97, 89};
                Physical Surface(4) = {99, 87};
                Physical Surface(5) = {91, 95};
                Physical Surface(6) = {101, 93};
            EndIf

            If(slab2_w == slab_w)
                Line(58) = {11, 252};
                Line(59) = {252, 250};
                Line(60) = {252, 356};
                Line(61) = {356, 2};
                Line(62) = {2, 354};
                Line(63) = {354, 12};
                Line(64) = {12, 355};
                Line(65) = {355, 3};
                Line(66) = {3, 357};
                Line(67) = {357, 253};
                Line(68) = {253, 13};
                Line(69) = {251, 253};
                Line(70) = {251, 351};
                Line(71) = {351, 357};
                Line(72) = {351, 38};
                Line(73) = {8, 38};
                Line(74) = {38, 12};
                Line(75) = {351, 355};
                Line(76) = {250, 350};
                Line(77) = {350, 354};
                Line(78) = {356, 350};
                Line(79) = {350, 38};

                Line Loop(80) = {9, 19, -22, -59, -58};
                Plane Surface(81) = {80};
                Line Loop(82) = {59, 76, -78, -60};
                Plane Surface(83) = {82};
                Line Loop(84) = {78, 77, -62, -61};
                Plane Surface(85) = {84};
                Line Loop(86) = {77, 63, -74, -79};
                Plane Surface(87) = {86};
                Line Loop(88) = {74, 64, -75, 72};
                Plane Surface(89) = {88};
                Line Loop(90) = {75, 65, 66, -71};
                Plane Surface(91) = {90};
                Line Loop(92) = {69, -67, -71, -70};
                Plane Surface(93) = {92};
                Line Loop(94) = {23, 70, 72, -73};
                Plane Surface(95) = {94};
                Line Loop(96) = {12, -68, -69, -23, 20};
                Plane Surface(97) = {96};
                Line Loop(98) = {22, 73, -79, -76};
                Plane Surface(99) = {98};

                Physical Line(100) = {5, 58, 60, 61};
                Physical Line(101) = {62, 63, 64, 65};
                Physical Line(102) = {66, 67, 68, 7};

                Physical Surface(1) = {81, 46, 48, 97};
                Physical Surface(3) = {99, 95};
                Physical Surface(4) = {93, 83};
                Physical Surface(5) = {87, 89};
                Physical Surface(6) = {85, 91};
            EndIf
        EndIf
    EndIf
EndIf


If(rect == 1)
    Point(150) = {-hx+d/2+radius1, -hy+slab_h+slab2_h+2*radius1y, 0,lc_refine_1};
    Point(151) = {-hx+d/2-radius1, -hy+slab_h+slab2_h+2*radius1y, 0,lc_refine_1};
    Point(152) = {-hx+d/2+radius1, -hy+slab_h+slab2_h, 0,lc_refine_1};
    Point(153) = {-hx+d/2-radius1, -hy+slab_h+slab2_h, 0,lc_refine_1};

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
            Point(352) = {d/2-slab2_w/2, -hy+slab_h+slab2_h, 0, lc_refine_4};
            Point(353) = {d/2+slab2_w/2, -hy+slab_h+slab2_h, 0, lc_refine_4};
            Point(354) = {d/2-slab2_w/2, -hy, 0, lc_refine_4};
            Point(355) = {d/2+slab2_w/2, -hy, 0, lc_refine_4};
            Point(356) = {d/2-slab_w/2, -hy+slab2_h, 0, lc_refine_4};
            Point(357) = {d/2+slab_w/2, -hy+slab2_h, 0, lc_refine_4};

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
        Point(354) = {d/2-slab2_w/2, -hy, 0, lc_refine_4};
        Point(355) = {d/2+slab2_w/2, -hy, 0, lc_refine_4};
        Point(356) = {d/2-slab_w/2, -hy+slab2_h, 0, lc_refine_4};
        Point(357) = {d/2+slab_w/2, -hy+slab2_h, 0, lc_refine_4};

        Line(17) = {151, 6};
        Line(18) = {6, 150};
        Line(19) = {150, 9};
        Line(20) = {9, 152};
        Line(23) = {153, 7};
        Line(24) = {7, 151};
        Line(26) = {153, 250};
        Line(30) = {251, 152};

        If(2*radius1 < slab_w)
            If(slab2_w_full == 1)
                Point(352) = {d/2-slab2_w/2, -hy+slab_h+slab2_h, 0, lc_refine_4};
                Point(353) = {d/2+slab2_w/2, -hy+slab_h+slab2_h, 0, lc_refine_4};
                Line(21) = {152, 8};
                Line(22) = {8, 153};
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

            If(slab2_w_full == 0)
                Point(358) = {d/2-slab_w/2, -hy, 0, lc_refine_4};
                Point(359) = {d/2+slab_w/2, -hy, 0, lc_refine_4};
                Point(360) = {0, -hy+slab2_h, 0, lc_refine_4};
                Point(361) = {d, -hy+slab2_h, 0, lc_refine_4};
                Point(362) = {0, -hy+slab_h+slab2_h, 0, lc_refine_4};
                Point(363) = {d, -hy+slab_h+slab2_h, 0, lc_refine_4};
                Line(31) = {11, 362};
                Line(32) = {362, 360};
                Line(33) = {360, 2};
                Line(36) = {354, 12};
                Line(38) = {355, 359};
                Line(40) = {3, 361};
                Line(41) = {361, 363};
                Line(42) = {363, 13};
                Line(45) = {356, 350};
                Line(46) = {351, 357};
                Line(47) = {8, 38};
                Line(48) = {38, 12};

                If(slab2_w > slab_w)
                    Line(34) = {2, 354};
                    Line(35) = {358, 354};
                    Line(37) = {12, 359};
                    Line(39) = {355, 3};
                    Line(53) = {153, 8};
                    Line(54) = {8, 152};
                    Line(55) = {251, 363};
                    Line(56) = {251, 357};
                    Line(57) = {351, 361};
                    Line(58) = {351, 355};
                    Line(59) = {357, 38};
                    Line(60) = {356, 38};
                    Line(61) = {350, 354};
                    Line(62) = {362, 250};
                    Line(63) = {250, 356};
                    Line(64) = {360, 350};
                                    
                    Line Loop(65) = {5, 9, 24, 17, -13, -1};
                    Plane Surface(66) = {65};
                    Line Loop(67) = {13, 18, 19, 12, -7, -2};
                    Plane Surface(68) = {67};
                    Line Loop(69) = {12, -42, -55, 30, -20};
                    Plane Surface(70) = {69};
                    Line Loop(71) = {23, -9, 31, 62, -26};
                    Plane Surface(72) = {71};
                    Line Loop(73) = {23, 10, 15, -53};
                    Plane Surface(74) = {73};
                    Line Loop(75) = {15, 54, -20, -11};
                    Plane Surface(76) = {75};
                    Line Loop(77) = {11, -19, -18, 14};
                    Plane Surface(78) = {77};
                    Line Loop(79) = {14, -10, 24, 17};
                    Plane Surface(80) = {79};
                    Line Loop(81) = {62, 63, 45, -64, -32};
                    Plane Surface(82) = {81};
                    Line Loop(83) = {63, 60, -47, -53, 26};
                    Plane Surface(84) = {83};
                    Line Loop(85) = {47, -59, -56, 30, -54};
                    Plane Surface(86) = {85};
                    Line Loop(87) = {56, -46, 57, 41, -55};
                    Plane Surface(88) = {87};
                    Line Loop(93) = {64, 61, -34, -33};
                    Plane Surface(94) = {93};
                    Line Loop(95) = {61, 36, -48, -60, 45};
                    Plane Surface(96) = {95};
                    Line Loop(97) = {59, 48, 37, -38, -58, 46};
                    Plane Surface(98) = {97};
                    Line Loop(105) = {58, 39, 40, -57};
                    Plane Surface(106) = {105};

                    Physical Line(49) = {1, 2};
                    Physical Line(50) = {7, 42, 41, 40};
                    Physical Line(52) = {33, 32, 31, 5};
                    Physical Line(107) = {34, 35, 36, 37, 38, 39};

                    Physical Surface(1) = {66, 68, 70, 72};
                    Physical Surface(2) = {80, 78, 76, 74};
                    Physical Surface(3) = {84, 86};
                    Physical Surface(4) = {88, 82};
                    Physical Surface(5) = {96, 98};
                    Physical Surface(6) = {106, 94};
                EndIf

                If(slab2_w < slab_w)
                    Line(49) = {362, 250};
                    Line(50) = {356, 360};
                    Line(51) = {350, 38};
                    Line(52) = {350, 354};
                    Line(53) = {351, 355};
                    Line(54) = {351, 38};
                    Line(55) = {8, 152};
                    Line(56) = {8, 153};
                    Line(57) = {251, 363};
                    Line(58) = {361, 357};
                    Line(59) = {251, 357};
                    Line(60) = {250, 356};
                    Line(61) = {358, 354};
                    Line(62) = {358, 2};
                    Line(63) = {12, 355};
                    Line(64) = {359, 3};

                    Line Loop(65) = {5, 9, 24, 17, -13, -1};
                    Plane Surface(66) = {65};
                    Line Loop(67) = {13, 18, 19, 12, -7, -2};
                    Plane Surface(68) = {67};
                    Line Loop(69) = {12, -42, -57, 30, -20};
                    Plane Surface(70) = {69};
                    Line Loop(71) = {23, -9, 31, 49, -26};
                    Plane Surface(72) = {71};
                    Line Loop(73) = {23, 10, 15, 56};
                    Plane Surface(74) = {73};
                    Line Loop(75) = {10, -14, -17, -24};
                    Plane Surface(76) = {75};
                    Line Loop(77) = {14, 11, -19, -18};
                    Plane Surface(78) = {77};
                    Line Loop(79) = {11, 20, -55, -15};
                    Plane Surface(80) = {79};
                    Line Loop(81) = {30, -55, 47, -54, 46, -59};
                    Plane Surface(82) = {81};
                    Line Loop(83) = {47, -51, -45, -60, -26, -56};
                    Plane Surface(84) = {83};
                    Line Loop(85) = {60, 50, -32, 49};
                    Plane Surface(86) = {85};
                    Line Loop(87) = {59, -58, 41, -57};
                    Plane Surface(88) = {87};
                    Line Loop(89) = {58, -46, 53, 38, 64, 40};
                    Plane Surface(90) = {89};
                    Line Loop(91) = {50, 33, -62, 61, -52, -45};
                    Plane Surface(92) = {91};
                    Line Loop(93) = {52, 36, -48, -51};
                    Plane Surface(94) = {93};
                    Line Loop(95) = {48, 63, -53, 54};
                    Plane Surface(96) = {95};

                    Physical Line(97) = {62, 61, 36, 63, 38, 64};
                    Physical Line(98) = {40, 41, 42, 7};
                    Physical Line(99) = {1, 2};
                    Physical Line(100) = {5, 31, 32, 33};

                    Physical Surface(1) = {66, 68, 70, 72};
                    Physical Surface(2) = {76, 78, 80, 74};
                    Physical Surface(3) = {84, 82};
                    Physical Surface(4) = {86, 88};
                    Physical Surface(5) = {96, 94};
                    Physical Surface(6) = {92, 90};
                EndIf
                
                If(slab2_w == slab_w)
                    Line(49) = {250, 350};
                    Line(50) = {350, 38};
                    Line(51) = {8, 153};
                    Line(52) = {8, 152};
                    Line(53) = {251, 363};
                    Line(54) = {351, 361};
                    Line(55) = {3, 355};
                    Line(56) = {351, 251};
                    Line(57) = {351, 355};
                    Line(58) = {351, 38};
                    Line(59) = {12, 355};
                    Line(60) = {354, 2};
                    Line(61) = {360, 350};
                    Line(62) = {350, 354};
                    Line(63) = {250, 362};
                    Line Loop(64) = {5, 9, 24, 17, -13, -1};
                    Plane Surface(65) = {64};
                    Line Loop(66) = {13, 18, 19, 12, -7, -2};
                    Plane Surface(67) = {66};
                    Line Loop(68) = {12, -42, -53, 30, -20};
                    Plane Surface(69) = {68};
                    Line Loop(70) = {20, -52, -15, 11};
                    Plane Surface(71) = {70};
                    Line Loop(72) = {19, -11, -14, 18};
                    Plane Surface(73) = {72};
                    Line Loop(74) = {14, -10, 24, 17};
                    Line Loop(75) = {9, -23, 26, 63, -31};
                    Plane Surface(76) = {75};
                    Line Loop(77) = {23, 10, 15, 51};
                    Plane Surface(78) = {77};
                    Plane Surface(79) = {74};
                    Line Loop(80) = {30, -52, 47, -58, 56};
                    Plane Surface(81) = {80};
                    Line Loop(82) = {56, 53, -41, -54};
                    Plane Surface(83) = {82};
                    Line Loop(84) = {47, -50, -49, -26, -51};
                    Plane Surface(85) = {84};
                    Line Loop(86) = {49, -61, -32, -63};
                    Plane Surface(87) = {86};
                    Line Loop(88) = {61, 62, 60, -33};
                    Plane Surface(89) = {88};
                    Line Loop(90) = {36, -48, -50, 62};
                    Plane Surface(91) = {90};
                    Line Loop(92) = {48, 59, -57, 58};
                    Plane Surface(93) = {92};
                    Line Loop(94) = {57, -55, 40, -54};
                    Plane Surface(95) = {94};

                    Physical Line(96) = {1, 2};
                    Physical Line(97) = {7, 42, 41, 40};
                    Physical Line(98) = {55, 59, 36, 60};
                    Physical Line(99) = {33, 32, 31, 5};

                    Physical Surface(1) = {65, 67, 69, 76};
                    Physical Surface(2) = {79, 73, 71, 78};
                    Physical Surface(3) = {85, 81};
                    Physical Surface(4) = {83, 87};
                    Physical Surface(5) = {91, 93};
                    Physical Surface(6) = {89, 95};
                EndIf
            EndIf
        EndIf

        If(2*radius1 > slab_w)
            If(slab2_w_full == 1)
                Point(352) = {d/2-slab2_w/2, -hy+slab_h+slab2_h, 0, lc_refine_1};
                Point(353) = {d/2+slab2_w/2, -hy+slab_h+slab2_h, 0, lc_refine_1};
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
                Line(81) = {250, 8};
                Line(82) = {8, 251};

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
                Plane Surface(64) = {59};
                Line Loop(69) = {41, 42, 43, 39};
                Plane Surface(70) = {69};
                Line Loop(71) = {46, -34, -33, 32};
                Plane Surface(72) = {71};
                Line Loop(73) = {34, -45, 47, -36, -35};
                Plane Surface(74) = {73};
                Line Loop(75) = {47, 37, 38, -43, 44};
                Plane Surface(76) = {75};
                Line Loop(83) = {23, 10, 15, -81, -26};
                Plane Surface(84) = {83};
                Line Loop(85) = {15, 82, 30, -20, -11};
                Plane Surface(86) = {85};
                Line Loop(87) = {82, 42, 44, -48};
                Plane Surface(88) = {87};
                Line Loop(89) = {48, 45, -46, 81};
                Plane Surface(90) = {89};

                Physical Line(77) = {5, 31, 33, 35};
                Physical Line(78) = {36, 37};
                Physical Line(79) = {7, 40, 39, 38};
                Physical Line(80) = {2, 1};

                Physical Surface(1) = {50, 52, 54, 56};
                Physical Surface(2) = {58, 64, 84, 86};
                Physical Surface(3) = {90, 88};
                Physical Surface(4) = {72, 70};
                Physical Surface(5) = {76, 74};
            EndIf

            If(slab2_w_full == 0)
                Point(358) = {d/2-slab_w/2, -hy, 0, lc_refine_4};
                Point(359) = {d/2+slab_w/2, -hy, 0, lc_refine_4};
                Point(360) = {0, -hy+slab2_h, 0, lc_refine_4};
                Point(361) = {d, -hy+slab2_h, 0, lc_refine_4};
                Point(362) = {0, -hy+slab_h+slab2_h, 0, lc_refine_4};
                Point(363) = {d, -hy+slab_h+slab2_h, 0, lc_refine_4};
                Line(31) = {11, 362};
                Line(32) = {362, 360};
                Line(33) = {360, 2};
                Line(36) = {354, 358};
                Line(38) = {355, 359};
                Line(40) = {3, 361};
                Line(41) = {361, 363};
                Line(42) = {363, 13};
                Line(45) = {356, 350};
                Line(46) = {351, 357};
                Line(47) = {8, 38};
                Line(48) = {38, 12};

                If(slab2_w > slab_w)
                    Line(49) = {362, 153};
                    Line(50) = {250, 8};
                    Line(51) = {8, 251};
                    Line(52) = {152, 363};
                    Line(53) = {361, 351};
                    Line(54) = {351, 355};
                    Line(55) = {355, 3};
                    Line(56) = {357, 251};
                    Line(57) = {357, 38};
                    Line(58) = {12, 359};
                    Line(59) = {12, 358};
                    Line(60) = {350, 354};
                    Line(61) = {356, 250};
                    Line(62) = {38, 356};
                    Line(63) = {350, 360};
                    Line(64) = {2, 354};

                    Line Loop(65) = {5, 9, 24, 17, -13, -1};
                    Plane Surface(66) = {65};
                    Line Loop(67) = {13, 18, 19, 12, -7, -2};
                    Plane Surface(68) = {67};
                    Line Loop(69) = {12, -42, -52, -20};
                    Plane Surface(70) = {69};
                    Line Loop(71) = {20, -30, -51, -15, 11};
                    Plane Surface(72) = {71};
                    Line Loop(73) = {11, -19, -18, 14};
                    Plane Surface(74) = {73};
                    Line Loop(75) = {14, -10, 24, 17};
                    Plane Surface(76) = {75};
                    Line Loop(77) = {10, 15, -50, -26, 23};
                    Plane Surface(78) = {77};
                    Line Loop(79) = {23, -9, 31, 49};
                    Plane Surface(80) = {79};
                    Line Loop(81) = {49, 26, -61, 45, 63, -32};
                    Plane Surface(82) = {81};
                    Line Loop(83) = {30, 52, -41, 53, 46, 56};
                    Plane Surface(84) = {83};
                    Line Loop(85) = {56, -51, 47, -57};
                    Plane Surface(86) = {85};
                    Line Loop(87) = {47, 62, 61, 50};
                    Plane Surface(88) = {87};
                    Line Loop(89) = {62, 45, 60, 36, -59, -48};
                    Plane Surface(90) = {89};
                    Line Loop(91) = {48, 58, -38, -54, 46, 57};
                    Plane Surface(92) = {91};
                    Line Loop(93) = {54, 55, 40, 53};
                    Plane Surface(94) = {93};
                    Line Loop(95) = {60, -64, -33, -63};
                    Plane Surface(96) = {95};

                    Physical Line(97) = {64, 36, 59, 58, 38, 55};
                    Physical Line(98) = {5, 31, 32, 33};
                    Physical Line(99) = {1, 2};
                    Physical Line(100) = {7, 42, 41, 40};

                    Physical Surface(1) = {66, 68, 70, 80};
                    Physical Surface(2) = {76, 74, 72, 78};
                    Physical Surface(3) = {88, 86};
                    Physical Surface(4) = {84, 82};
                    Physical Surface(5) = {90, 92};
                    Physical Surface(6) = {96, 94};
                EndIf

                If(slab2_w < slab_w)
                    Line(49) = {350, 354};
                    Line(50) = {354, 12};
                    Line(51) = {12, 355};
                    Line(52) = {351, 355};
                    Line(53) = {351, 38};
                    Line(54) = {38, 350};
                    Line(55) = {356, 360};
                    Line(56) = {2, 358};
                    Line(57) = {359, 3};
                    Line(58) = {361, 357};
                    Line(59) = {357, 251};
                    Line(60) = {251, 8};
                    Line(61) = {8, 250};
                    Line(62) = {250, 356};
                    Line(63) = {153, 362};
                    Line(64) = {152, 363};

                    Line Loop(65) = {5, 9, 24, 17, -13, -1};
                    Plane Surface(66) = {65};
                    Line Loop(67) = {13, 18, 19, 12, -7, -2};
                    Plane Surface(68) = {67};
                    Line Loop(69) = {12, -42, -64, -20};
                    Plane Surface(70) = {69};
                    Line Loop(71) = {9, -23, 63, -31};
                    Plane Surface(72) = {71};
                    Line Loop(73) = {10, -14, -17, -24};
                    Plane Surface(74) = {73};
                    Line Loop(75) = {14, 11, -19, -18};
                    Plane Surface(76) = {75};
                    Line Loop(77) = {11, 20, -30, 60, -15};
                    Plane Surface(78) = {77};
                    Line Loop(79) = {15, 61, -26, 23, 10};
                    Plane Surface(80) = {79};
                    Line Loop(81) = {47, -53, 46, 59, 60};
                    Plane Surface(82) = {81};
                    Line Loop(83) = {47, 54, -45, -62, -61};
                    Plane Surface(84) = {83};
                    Line Loop(85) = {62, 55, -32, -63, 26};
                    Plane Surface(86) = {85};
                    Line Loop(87) = {55, 33, 56, -36, -49, -45};
                    Plane Surface(88) = {87};
                    Line Loop(89) = {49, 50, -48, 54};
                    Plane Surface(90) = {89};
                    Line Loop(91) = {48, 51, -52, 53};
                    Plane Surface(92) = {91};
                    Line Loop(93) = {52, 38, 57, 40, 58, -46};
                    Line Loop(94) = {30, 64, -41, 58, 59};
                    Plane Surface(95) = {94};
                    Plane Surface(96) = {93};

                    Physical Line(97) = {56, 36, 50, 51, 38, 57};
                    Physical Line(98) = {40, 41, 42, 7};
                    Physical Line(99) = {2, 1};
                    Physical Line(100) = {5, 31, 32, 33};

                    Physical Surface(1) = {66, 68, 70, 72};
                    Physical Surface(2) = {74, 80, 78, 76};
                    Physical Surface(3) = {82, 84};
                    Physical Surface(4) = {86, 95};
                    Physical Surface(5) = {92, 90};
                    Physical Surface(6) = {96, 88};
                EndIf

                If(slab2_w == slab_w)
                    Line(49) = {362, 153};
                    Line(50) = {250, 8};
                    Line(51) = {8, 251};
                    Line(52) = {152, 363};
                    Line(53) = {361, 351};
                    Line(54) = {351, 38};
                    Line(55) = {38, 350};
                    Line(56) = {350, 360};
                    Line(57) = {2, 354};
                    Line(58) = {354, 12};
                    Line(59) = {12, 355};
                    Line(60) = {355, 3};
                    Line(61) = {355, 351};
                    Line(62) = {351, 251};
                    Line(63) = {250, 350};
                    Line(64) = {350, 354};

                    Line Loop(65) = {5, 9, 24, 17, -13, -1};
                    Plane Surface(66) = {65};
                    Line Loop(67) = {13, 18, 19, 12, -7, -2};
                    Plane Surface(68) = {67};
                    Line Loop(69) = {12, -42, -52, -20};
                    Plane Surface(70) = {69};
                    Line Loop(71) = {20, -30, -51, -15, 11};
                    Plane Surface(72) = {71};
                    Line Loop(73) = {11, -19, -18, 14};
                    Plane Surface(74) = {73};
                    Line Loop(75) = {17, 14, -10, 24};
                    Plane Surface(76) = {75};
                    Line Loop(77) = {23, 10, 15, -50, -26};
                    Plane Surface(78) = {77};
                    Line Loop(79) = {9, -23, -49, -31};
                    Plane Surface(80) = {79};
                    Line Loop(81) = {32, -56, -63, -26, -49};
                    Plane Surface(82) = {81};
                    Line Loop(83) = {63, -55, -47, -50};
                    Plane Surface(84) = {83};
                    Line Loop(85) = {47, -54, 62, -51};
                    Plane Surface(86) = {85};
                    Line Loop(87) = {52, -41, 53, 62, 30};
                    Plane Surface(88) = {87};
                    Line Loop(89) = {53, -61, 60, 40};
                    Plane Surface(90) = {89};
                    Line Loop(91) = {61, 54, 48, 59};
                    Plane Surface(92) = {91};
                    Line Loop(93) = {48, -58, -64, -55};
                    Plane Surface(94) = {93};
                    Line Loop(95) = {64, -57, -33, -56};
                    Plane Surface(96) = {95};

                    Physical Line(97) = {1, 2};
                    Physical Line(98) = {7, 42, 41, 40};
                    Physical Line(99) = {60, 59, 58, 57};
                    Physical Line(100) = {33, 32, 31, 5};

                    Physical Surface(1) = {66, 68, 70, 80};
                    Physical Surface(2) = {76, 74, 72, 78};
                    Physical Surface(3) = {86, 84};
                    Physical Surface(4) = {88, 82};
                    Physical Surface(5) = {94, 92};
                    Physical Surface(6) = {90, 96};
                EndIf
            EndIf
        EndIf

        If(2*radius1 == slab_w)
            Point(358) = {d/2-slab_w/2, -hy, 0, lc_refine_4};
            Point(359) = {d/2+slab_w/2, -hy, 0, lc_refine_4};
            Point(362) = {0, -hy+slab_h+slab2_h, 0, lc_refine_4};
            Point(363) = {d, -hy+slab_h+slab2_h, 0, lc_refine_4};

            Line(31) = {11, 362};
            Line(34) = {13, 363};

            If(slab2_w > slab_w)
                If(slab2_w_full == 1)
                    Line(35) = {362, 350};
                    Line(36) = {350, 2};
                    Line(37) = {2, 358};
                    Line(38) = {358, 12};
                    Line(39) = {12, 359};
                    Line(40) = {359, 3};
                    Line(41) = {3, 351};
                    Line(42) = {363, 351};
                    Line(43) = {152, 363};
                    Line(44) = {357, 351};
                    Line(45) = {152, 357};
                    Line(46) = {8, 38};
                    Line(47) = {38, 357};
                    Line(48) = {356, 38};
                    Line(49) = {153, 356};
                    Line(50) = {153, 8};
                    Line(51) = {8, 152};
                    Line(52) = {362, 153};
                    Line(53) = {350, 356};

                    Line Loop(54) = {5, 9, 24, 17, -13, -1};
                    Plane Surface(55) = {54};
                    Line Loop(56) = {13, 18, 19, 12, -7, -2};
                    Plane Surface(57) = {56};
                    Line Loop(58) = {12, 34, -43, -20};
                    Plane Surface(59) = {58};
                    Line Loop(60) = {20, -51, -15, 11};
                    Plane Surface(61) = {60};
                    Line Loop(62) = {18, 19, -11, -14};
                    Plane Surface(63) = {62};
                    Line Loop(64) = {17, 14, -10, 24};
                    Plane Surface(65) = {64};
                    Line Loop(66) = {10, 15, -50, 23};
                    Plane Surface(67) = {66};
                    Line Loop(68) = {23, -9, 31, 52};
                    Plane Surface(69) = {68};
                    Line Loop(70) = {52, 49, -53, -35};
                    Plane Surface(71) = {70};
                    Line Loop(72) = {49, 48, -46, -50};
                    Plane Surface(73) = {72};
                    Line Loop(74) = {51, 45, -47, -46};
                    Plane Surface(75) = {74};
                    Line Loop(76) = {43, 42, -44, -45};
                    Plane Surface(77) = {76};
                    Line Loop(78) = {44, -41, -40, -39, -38, -37, -36, 53, 48, 47};
                    Plane Surface(79) = {78};

                    Physical Line(80) = {5, 31, 35, 36};
                    Physical Line(81) = {37, 38, 39, 40};
                    Physical Line(82) = {1, 2};
                    Physical Line(83) = {7, 34, 42, 41};

                    Physical Surface(1) = {55, 57, 59, 69};
                    Physical Surface(2) = {65, 63, 61, 67};
                    Physical Surface(3) = {73, 75};
                    Physical Surface(4) = {77, 71};
                    Physical Surface(5) = {79};
                EndIf

                If(slab2_w_full == 0)
                    Point(360) = {0, -hy+slab2_h, 0, lc_refine_4};
                    Point(361) = {d, -hy+slab2_h, 0, lc_refine_4};
                    Line(32) = {362, 360};
                    Line(33) = {360, 2};
                    Line(35) = {363, 361};
                    Line(36) = {361, 3};
                    Line(37) = {2, 354};
                    Line(38) = {354, 358};
                    Line(39) = {358, 12};
                    Line(40) = {12, 359};
                    Line(41) = {359, 355};
                    Line(42) = {355, 3};
                    Line(43) = {361, 351};
                    Line(44) = {351, 355};
                    Line(45) = {363, 152};
                    Line(46) = {152, 357};
                    Line(47) = {357, 351};
                    Line(49) = {8, 38};
                    Line(50) = {38, 12};
                    Line(51) = {153, 356};
                    Line(53) = {356, 38};
                    Line(54) = {153, 8};
                    Line(55) = {8, 152};
                    Line(56) = {38, 357};
                    Line(57) = {350, 356};
                    Line(58) = {360, 350};
                    Line(59) = {362, 153};
                    Line(60) = {350, 354};

                    Line Loop(61) = {5, 9, 24, 17, -13, -1};
                    Plane Surface(62) = {61};
                    Line Loop(63) = {13, 18, 19, 12, -7, -2};
                    Plane Surface(64) = {63};
                    Line Loop(65) = {12, 34, 45, -20};
                    Plane Surface(66) = {65};
                    Line Loop(67) = {20, -55, -15, 11};
                    Plane Surface(68) = {67};
                    Line Loop(69) = {18, 19, -11, -14};
                    Plane Surface(70) = {69};
                    Line Loop(71) = {14, -10, 24, 17};
                    Plane Surface(72) = {71};
                    Line Loop(73) = {10, 15, -54, 23};
                    Plane Surface(74) = {73};
                    Line Loop(75) = {23, -9, 31, 59};
                    Plane Surface(76) = {75};
                    Line Loop(77) = {59, 51, -57, -58, -32};
                    Plane Surface(78) = {77};
                    Line Loop(79) = {51, 53, -49, -54};
                    Plane Surface(80) = {79};
                    Line Loop(81) = {49, 56, -46, -55};
                    Plane Surface(82) = {81};
                    Line Loop(83) = {45, 46, 47, -43, -35};
                    Plane Surface(84) = {83};
                    Line Loop(85) = {43, 44, 42, -36};
                    Plane Surface(86) = {85};
                    Line Loop(87) = {58, 60, -37, -33};
                    Plane Surface(88) = {87};
                    Line Loop(89) = {57, 53, 50, -39, -38, -60};
                    Plane Surface(90) = {89};
                    Line Loop(91) = {56, 47, 44, -41, -40, -50};
                    Plane Surface(92) = {91};

                    Physical Line(93) = {37, 38, 39, 40, 41, 42};
                    Physical Line(94) = {36, 35, 34, 7};
                    Physical Line(95) = {2, 1};
                    Physical Line(96) = {5, 31, 32, 33};

                    Physical Surface(1) = {62, 64, 66, 76};
                    Physical Surface(2) = {72, 70, 68, 74};
                    Physical Surface(3) = {82, 80};
                    Physical Surface(4) = {84, 78};
                    Physical Surface(5) = {90, 92};
                    Physical Surface(6) = {88, 86};
                EndIf
            EndIf

            If(slab2_w < slab_w)
                Point(360) = {0, -hy+slab2_h, 0, lc_refine_4};
                Point(361) = {d, -hy+slab2_h, 0, lc_refine_4};
                Line(32) = {362, 360};
                Line(33) = {360, 2};
                Line(35) = {363, 361};
                Line(36) = {361, 3};
                Line(37) = {362, 153};
                Line(38) = {153, 8};
                Line(39) = {8, 152};
                Line(40) = {152, 363};
                Line(41) = {361, 357};
                Line(42) = {152, 357};
                Line(43) = {357, 351};
                Line(44) = {351, 355};
                Line(45) = {359, 355};
                Line(46) = {359, 3};
                Line(47) = {351, 38};
                Line(48) = {8, 38};
                Line(49) = {38, 12};
                Line(50) = {350, 354};
                Line(51) = {354, 12};
                Line(52) = {12, 355};
                Line(53) = {38, 350};
                Line(54) = {356, 350};
                Line(55) = {358, 354};
                Line(56) = {358, 2};
                Line(57) = {360, 356};
                Line(58) = {153, 356};

                Line Loop(59) = {5, 9, 24, 17, -13, -1};
                Plane Surface(60) = {59};
                Line Loop(61) = {13, 18, 19, 12, -7, -2};
                Plane Surface(62) = {61};
                Line Loop(63) = {12, 34, -40, -20};
                Plane Surface(64) = {63};
                Line Loop(65) = {20, -39, -15, 11};
                Plane Surface(66) = {65};
                Line Loop(67) = {11, -19, -18, 14};
                Plane Surface(68) = {67};
                Line Loop(69) = {17, 14, -10, 24};
                Plane Surface(70) = {69};
                Line Loop(71) = {10, 15, -38, 23};
                Plane Surface(72) = {71};
                Line Loop(73) = {23, -9, 31, 37};
                Plane Surface(74) = {73};
                Line Loop(75) = {37, 58, -57, -32};
                Plane Surface(76) = {75};
                Line Loop(77) = {58, 54, -53, -48, -38};
                Plane Surface(78) = {77};
                Line Loop(79) = {48, -47, -43, -42, -39};
                Plane Surface(80) = {79};
                Line Loop(81) = {42, -41, -35, -40};
                Plane Surface(82) = {81};
                Line Loop(83) = {41, 43, 44, -45, 46, -36};
                Plane Surface(84) = {83};
                Line Loop(85) = {44, -52, -49, -47};
                Plane Surface(86) = {85};
                Line Loop(87) = {49, -51, -50, -53};
                Plane Surface(88) = {87};
                Line Loop(89) = {50, -55, 56, -33, 57, 54};
                Plane Surface(90) = {89};

                Physical Line(91) = {1, 2};
                Physical Line(92) = {7, 34, 35, 36};
                Physical Line(93) = {46, 45, 52, 51, 55, 56};
                Physical Line(94) = {33, 32, 31, 5};

                Physical Surface(1) = {60, 62, 64, 74};
                Physical Surface(2) = {70, 68, 66, 72};
                Physical Surface(3) = {78, 80};
                Physical Surface(4) = {76, 82};
                Physical Surface(5) = {86, 88};
                Physical Surface(6) = {90, 84};
            EndIf

            If(slab2_w == slab_w)
                Point(360) = {0, -hy+slab2_h, 0, lc_refine_4};
                Point(361) = {d, -hy+slab2_h, 0, lc_refine_4};
                Line(32) = {362, 360};
                Line(33) = {360, 2};
                Line(35) = {363, 361};
                Line(36) = {361, 3};
                Line(37) = {152, 351};
                Line(38) = {152, 363};
                Line(39) = {152, 8};
                Line(40) = {8, 153};
                Line(41) = {362, 153};
                Line(42) = {153, 350};
                Line(43) = {8, 38};
                Line(44) = {38, 12};
                Line(45) = {351, 355};
                Line(46) = {350, 354};

                Line(47) = {360, 350};
                Line(48) = {2, 354};
                Line(49) = {354, 12};
                Line(50) = {12, 355};
                Line(51) = {351, 38};
                Line(52) = {350, 38};
                Line(53) = {351, 361};
                Line(54) = {355, 3};

                Line Loop(55) = {5, 9, 24, 17, -13, -1};
                Plane Surface(56) = {55};
                Line Loop(57) = {13, 18, 19, 12, -7, -2};
                Plane Surface(58) = {57};
                Line Loop(59) = {19, -11, -14, 18};
                Plane Surface(60) = {59};
                Line Loop(61) = {14, -10, 24, 17};
                Plane Surface(62) = {61};
                Line Loop(63) = {23, -9, 31, 41};
                Plane Surface(64) = {63};
                Line Loop(65) = {23, 10, 15, 40};
                Plane Surface(66) = {65};
                Line Loop(67) = {15, -39, -20, -11};
                Plane Surface(68) = {67};
                Line Loop(69) = {20, 38, -34, -12};
                Plane Surface(70) = {69};
                Line Loop(71) = {38, 35, -53, -37};
                Plane Surface(72) = {71};
                Line Loop(73) = {37, 51, -43, -39};
                Plane Surface(74) = {73};
                Line Loop(75) = {43, -52, -42, -40};
                Plane Surface(76) = {75};
                Line Loop(77) = {42, -47, -32, 41};
                Plane Surface(78) = {77};
                Line Loop(79) = {47, 46, -48, -33};
                Plane Surface(80) = {79};
                Line Loop(81) = {46, 49, -44, -52};
                Plane Surface(82) = {81};
                Line Loop(83) = {44, 50, -45, 51};
                Plane Surface(84) = {83};
                Line Loop(85) = {45, 54, -36, -53};
                Plane Surface(86) = {85};

                Physical Line(87) = {1, 2};
                Physical Line(88) = {5, 31, 32, 33};
                Physical Line(89) = {48, 49, 50, 54};
                Physical Line(90) = {7, 34, 35, 36};

                Physical Surface(1) = {56, 58, 70, 64};
                Physical Surface(2) = {62, 60, 68, 66};
                Physical Surface(3) = {76, 74};
                Physical Surface(4) = {72, 78};
                Physical Surface(5) = {82, 84};
                Physical Surface(6) = {80, 86};
            EndIf
        EndIf
    EndIf
EndIf
