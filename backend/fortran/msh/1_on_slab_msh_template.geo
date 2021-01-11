// Template mesh geometry file for a single inclusion on a slab.
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

rect = 1;

slab_width = d_in_nm;
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

hy = dy/2 + (slab_h/2) + radius1y; // 
hx = 0.;


Point(1) = {0, 0, 0, lc};
Point(2) = {-hx, -hy, 0, lc};
Point(3) = {-hx+d, -hy, 0, lc};
Point(4) = {d, 0, 0,lc};


// Slab
Point(250) = {d/2-slab_w/2, -hy+slab_h, 0, lc};
Point(251) = {d/2+slab_w/2, -hy+slab_h, 0, lc};

// Inclusion
Point(5) = {-hx+d/2, -hy+radius1y+slab_h, 0, lc_refine_2};
Point(6) = {-hx+d/2, -hy+2*radius1y+slab_h, 0, lc_refine_1};
Point(7) = {-hx+d/2-radius1, -hy+radius1y+slab_h, 0, lc_refine_1};
Point(8) = {-hx+d/2, -hy+slab_h, 0, lc_refine_1};
Point(9) = {-hx+d/2+radius1, -hy+radius1y+slab_h, 0, lc_refine_1};

Point(10) = {-hx+d/2, 0, 0, lc};
Point(12) = {-hx+d/2, -hy, 0, lc};

Point(990) = {-hx, -dy, 0, lc};
Point(991) = {-hx+d, -dy, 0, lc};
Point(992) = {-hx+d/2, -dy, 0, lc};

Line(1) = {1,10};
Line(2) = {10,4};
Line(10) = {7,5};
Line(11) = {5,9};
Line(13) = {10,6};
Line(14) = {6,5};
Line(15) = {5,8};
Line(16) = {8,12};


If(slab_w_full == 0)
    Point(252) = {d/2-slab_w/2, -hy, 0, lc_refine_4};
    Point(253) = {d/2+slab_w/2, -hy, 0, lc_refine_4};
    Line(3) = {2,252};
    Line(4) = {253,3};
EndIf

If(slab_w_full == 1)
    Line(63) = {2, 990};
    Line(64) = {990, 992};
    Line(65) = {992, 991};
    Line(66) = {991, 3};
    Line(67) = {12, 992};
EndIf

If(rect == 0)
    Ellipsis(17) = {9,5,6,6};
    Ellipsis(18) = {6,5,7,7};
    Ellipsis(19) = {7,5,8,8};
    Ellipsis(20) = {8,5,9,9};

    If(slab_w_full == 1)
        Line(6) = {1,250};
        Line(8) = {4,251};
        Line(21) = {2, 250};
        Line(22) = {250, 8};
        Line(23) = {8, 251};
        Line(24) = {251, 3};
        Line(25) = {3, 12};
        Line(26) = {12, 2};

        Line Loop(39) = {17, 14, 11};
        Plane Surface(40) = {39};
        Line Loop(41) = {14, -10, -18};
        Plane Surface(42) = {41};
        Line Loop(43) = {10, 15, -19};
        Plane Surface(44) = {43};
        Line Loop(45) = {15, 20, -11};
        Plane Surface(46) = {45};
        Line Loop(47) = {23, 24, 25, -16};
        Plane Surface(48) = {47};
        Line Loop(49) = {16, 26, 21, 22};
        Plane Surface(50) = {49};
        Line Loop(51) = {22, -19, -18, -13, -1, 6};
        Plane Surface(52) = {51};
        Line Loop(53) = {20, 17, -13, 2, 8, -23};
        Plane Surface(54) = {53};

        Physical Line(27) = {5, 21};
        Physical Line(28) = {1, 2};
        Physical Line(29) = {7, 24};
        Physical Line(30) = {25, 26};

        Physical Surface(1) = {52, 54};
        Physical Surface(2) = {42, 40, 46, 44};
        Physical Surface(3) = {50, 48};
    EndIf

    If(slab_w_full == 0)            
        Point(352) = {0, -hy+slab_h, 0, lc};
        Point(353) = {d, -hy+slab_h, 0, lc};
        Line(6) = {1,352};
        Line(8) = {4,353};
        Line(21) = {252, 250};
        Line(22) = {250, 8};
        Line(23) = {8, 251};
        Line(24) = {251, 253};
        Line(25) = {253, 12};
        Line(26) = {12, 252};
        Line(27) = {2, 352};
        Line(28) = {352, 250};
        Line(29) = {251, 353};
        Line(30) = {353, 3};

        Line Loop(31) = {22, -19, -18, -13, -1, 6, 28};
        Plane Surface(32) = {31};
        Line Loop(33) = {23, 29, -8, -2, 13, -17, -20};
        Plane Surface(34) = {33};
        Line Loop(35) = {17, 14, 11};
        Plane Surface(36) = {35};
        Line Loop(37) = {11, -20, -15};
        Plane Surface(38) = {37};
        Line Loop(39) = {15, -19, 10};
        Plane Surface(40) = {39};
        Line Loop(41) = {10, -14, 18};
        Plane Surface(42) = {41};
        Line Loop(43) = {24, 25, -16, 23};
        Plane Surface(44) = {43};
        Line Loop(45) = {26, 21, 22, 16};
        Plane Surface(46) = {45};
        Line Loop(47) = {28, -21, -3, 27};
        Plane Surface(48) = {47};
        Line Loop(49) = {24, 4, -30, -29};
        Plane Surface(50) = {49};
        Physical Line(51) = {3, 26, 25, 4};
        Physical Line(52) = {27, 6};
        Physical Line(53) = {1, 2};
        Physical Line(54) = {8, 30};
        Physical Surface(1) = {34, 32};
        Physical Surface(2) = {42, 36, 38, 40};
        Physical Surface(4) = {50, 48};
        Physical Surface(3) = {46, 44};
    EndIf
EndIf


If(rect == 1)
    Point(150) = {-hx+d/2+radius1, -hy+slab_h+2*radius1y, 0,lc_refine_1};
    Point(151) = {-hx+d/2-radius1, -hy+slab_h+2*radius1y, 0,lc_refine_1};
    Point(152) = {-hx+d/2+radius1, -hy+slab_h, 0,lc_refine_1};
    Point(153) = {-hx+d/2-radius1, -hy+slab_h, 0,lc_refine_1};

    If(slab_w > 2*radius1)
        If(slab_w_full == 1)
            Line(6) = {1,250};
            Line(8) = {4,251};
            Line(21) = {2, 250};
            Line(22) = {250, 153};
            Line(23) = {152, 251};
            Line(24) = {251, 3};
            Line(25) = {3, 12};
            Line(26) = {12, 2};        
            Line(27) = {153, 8};
            Line(28) = {8, 152};
            Line(29) = {153, 7};
            Line(30) = {7, 151};
            Line(31) = {151, 6};
            Line(32) = {6, 150};
            Line(33) = {150, 9};
            Line(34) = {9, 152};

            Line Loop(43) = {34, -28, -15, 11};
            Plane Surface(44) = {43};
            Line Loop(45) = {11, -33, -32, 14};
            Plane Surface(46) = {45};
            Line Loop(47) = {14, -10, 30, 31};
            Plane Surface(48) = {47};
            Line Loop(49) = {10, 15, -27, 29};
            Plane Surface(50) = {49};
            Line Loop(51) = {27, 16, 26, 21, 22};
            Plane Surface(52) = {51};
            Line Loop(53) = {16, -25, -24, -23, -28};
            Plane Surface(54) = {53};
            Line Loop(59) = {6, 22, 29, 30, 31, -13, -1};
            Plane Surface(60) = {59};
            Line Loop(61) = {13, 32, 33, 34, 23, -8, -2};
            Plane Surface(62) = {61};

            Physical Line(55) = {6, 21};
            Physical Line(56) = {26, 25};
            Physical Line(57) = {24, 8};
            Physical Line(58) = {2, 1};

            Physical Surface(1) = {60, 62};
            Physical Surface(2) = {48, 46, 44, 50};
            Physical Surface(3) = {52, 54};
        EndIf

        If(slab_w_full == 0)    
            Point(352) = {0, -hy+slab_h, 0, lc};
            Point(353) = {d, -hy+slab_h, 0, lc};
            Line(6) = {1,352};
            Line(8) = {4,353};
            Line(17) = {151, 6};
            Line(18) = {6, 150};
            Line(19) = {150, 9};
            Line(20) = {9, 152};
            Line(21) = {152, 8};
            Line(22) = {8, 153};
            Line(23) = {153, 7};
            Line(24) = {7, 151};
            Line(25) = {252, 250};
            Line(26) = {153, 250};
            Line(27) = {12, 252};
            Line(28) = {12, 253};
            Line(29) = {253, 251};
            Line(30) = {251, 152};
            Line(31) = {352, 2};
            Line(32) = {352, 250};
            Line(33) = {251, 353};
            Line(34) = {3, 353};

            Line Loop(39) = {24, 17, 14, -10};
            Plane Surface(40) = {39};
            Line Loop(41) = {14, 11, -19, -18};
            Plane Surface(42) = {41};
            Line Loop(47) = {33, -34, -4, 29};
            Plane Surface(48) = {47};
            Line Loop(49) = {20, 21, -15, 11};
            Plane Surface(50) = {49};
            Line Loop(51) = {10, 15, 22, 23};
            Plane Surface(52) = {51};
            Line Loop(53) = {21, 16, 28, 29, 30};
            Plane Surface(54) = {53};
            Line Loop(55) = {16, 27, 25, -26, -22};
            Plane Surface(56) = {55};
            Line Loop(57) = {25, -32, 31, 3};
            Plane Surface(58) = {57};
            Line Loop(63) = {26, -32, -6, 1, 13, -17, -24, -23};
            Plane Surface(64) = {63};
            Line Loop(65) = {13, 18, 19, 20, -30, 33, -8, -2};
            Plane Surface(66) = {65};
            Physical Line(59) = {6, 31};
            Physical Line(60) = {3, 27, 28, 4};
            Physical Line(61) = {34, 8};
            Physical Line(62) = {2, 1};

            Physical Surface(1) = {64, 66};
            Physical Surface(2) = {40, 42, 50, 52};
            Physical Surface(4) = {58, 48};
            Physical Surface(3) = {54, 56};

        EndIf
    EndIf

    If(slab_w < 2*radius1)   
        Point(352) = {0, -hy+slab_h, 0, lc};
        Point(353) = {d, -hy+slab_h, 0, lc};
        Line(6) = {1,352};
        Line(8) = {4,353};
        Line(17) = {151, 6};
        Line(18) = {6, 150};
        Line(19) = {150, 9};
        Line(20) = {9, 152};
        Line(21) = {152, 251};
        Line(22) = {251, 253};
        Line(23) = {8, 251};
        Line(24) = {253, 12};
        Line(25) = {12, 252};
        Line(26) = {252, 250};
        Line(27) = {250, 153};
        Line(28) = {153, 7};
        Line(29) = {7, 151};
        Line(30) = {250, 8};
        Line(31) = {352, 153};
        Line(32) = {352, 2};
        Line(33) = {152, 353};
        Line(34) = {353, 3};
        Line Loop(39) = {28, 10, 15, -30, 27};
        Plane Surface(40) = {39};
        Line Loop(41) = {10, -14, -17, -29};
        Plane Surface(42) = {41};
        Line Loop(45) = {19, -11, -14, 18};
        Plane Surface(46) = {45};
        Line Loop(47) = {11, 20, 21, -23, -15};
        Plane Surface(48) = {47};
        Line Loop(49) = {21, 22, 4, -34, -33};
        Plane Surface(50) = {49};
        Line Loop(53) = {23, 22, 24, -16};
        Plane Surface(54) = {53};
        Line Loop(55) = {16, 25, 26, 30};
        Plane Surface(56) = {55};
        Line Loop(57) = {26, 27, -31, 32, 3};
        Plane Surface(58) = {57};
        Line Loop(63) = {6, 31, 28, 29, 17, -13, -1};
        Plane Surface(64) = {63};
        Line Loop(65) = {13, 18, 19, 20, 33, -8, -2};
        Plane Surface(66) = {65};

        Physical Line(59) = {1, 2};
        Physical Line(60) = {8, 34};
        Physical Line(61) = {4, 24, 25, 3};
        Physical Line(62) = {32, 6};
        Physical Surface(1) = {64, 66};
        Physical Surface(2) = {42, 46, 48, 40};
        Physical Surface(4) = {50, 58};
        Physical Surface(3) = {56, 54};

    EndIf
EndIf
