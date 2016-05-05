// Template mesh geometry file for a single inclusion on a slab.
// Inclusion can be circular, elliptical square, or rectangular.

d = 1; // grating period
ff = 0;
d_in_nm = 100;
dy_in_nm = 50;
dy = dy_in_nm/d_in_nm;
a1 = 30;
a1y = 30;
radius1 = (a1/(2*d_in_nm))*d;
radius1y = (a1y/(2*d_in_nm))*d;

rect = 1;

slab_width = 100;
slab_height = 10;
slab_w = slab_width/d_in_nm;
slab_h = slab_height/d_in_nm;
slab_w_full = 0;
If(slab_w == 1)
    slab_w_full = 1;
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
Point(250) = {d/2-slab_w/2, -hy+slab_h, 0, lc5};
Point(251) = {d/2+slab_w/2, -hy+slab_h, 0, lc5};

// Inclusion
Point(5) = {-hx+d/2, -hy+radius1y+slab_h, 0, lc3};
Point(6) = {-hx+d/2, -hy+2*radius1y+slab_h, 0, lc2};
Point(7) = {-hx+d/2-radius1, -hy+radius1y+slab_h, 0, lc2};
Point(8) = {-hx+d/2, -hy+slab_h, 0, lc2};
Point(9) = {-hx+d/2+radius1, -hy+radius1y+slab_h, 0, lc2};

Point(10) = {-hx+d/2, 0, 0, lc4};
Point(11) = {0,-hy+radius1y+slab_h, 0, lc};
Point(12) = {-hx+d/2, -hy, 0, lc4};
Point(13) = {d, -hy+radius1y+slab_h, 0, lc};
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
Line(16) = {8,12};


If(slab_w_full == 0)
    Point(252) = {d/2-slab_w/2, -hy, 0, lc5};
    Point(253) = {d/2+slab_w/2, -hy, 0, lc5};

    Line(3) = {2,252};
    Line(4) = {253,3};
    Line(6) = {11,2};
    Line(8) = {13,3};
EndIf


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

        Line Loop(31) = {5, 9, -18, -13, -1};
        Plane Surface(32) = {31};
        Line Loop(33) = {13, -17, 12, -7, -2};
        Plane Surface(34) = {33};
        Line Loop(35) = {9, 19, -22, -6};
        Plane Surface(36) = {35};
        Line Loop(37) = {20, 12, 8, -23};
        Plane Surface(38) = {37};
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

        Physical Line(27) = {5, 6, 21};
        Physical Line(28) = {1, 2};
        Physical Line(29) = {7, 8, 24};
        Physical Line(30) = {25, 26};

        Physical Surface(1) = {32, 34, 38, 36};
        Physical Surface(2) = {42, 40, 46, 44};
        Physical Surface(3) = {50, 48};
    EndIf

    If(slab_w_full == 0)
        Line(21) = {252, 250};
        Line(22) = {250, 8};
        Line(23) = {8, 251};
        Line(24) = {251, 253};
        Line(25) = {253, 12};
        Line(26) = {12, 252};

        Line Loop(27) = {5, 9, -18, -13, -1};
        Plane Surface(28) = {27};
        Line Loop(29) = {13, -17, 12, -7, -2};
        Plane Surface(30) = {29};
        Line Loop(31) = {12, 8, -4, -24, -23, 20};
        Plane Surface(32) = {31};
        Line Loop(33) = {19, -22, -21, -3, -6, 9};
        Plane Surface(34) = {33};
        Line Loop(35) = {19, -15, -10};
        Plane Surface(36) = {35};
        Line Loop(37) = {15, 20, -11};
        Plane Surface(38) = {37};
        Line Loop(39) = {11, 17, 14};
        Plane Surface(40) = {39};
        Line Loop(41) = {18, 10, -14};
        Plane Surface(42) = {41};
        Line Loop(43) = {22, 16, 26, 21};
        Plane Surface(44) = {43};
        Line Loop(45) = {16, -25, -24, -23};
        Plane Surface(46) = {45};

        Physical Line(47) = {5, 6};
        Physical Line(48) = {3, 26, 25, 4};
        Physical Line(49) = {8, 7};
        Physical Line(50) = {2, 1};

        Physical Surface(1) = {28, 30, 32, 34};
        Physical Surface(2) = {42, 40, 38, 36};
        Physical Surface(3) = {46, 44};
    EndIf
EndIf


If(rect == 1)
    Point(150) = {-hx+d/2+radius1, -hy+slab_h+2*radius1y, 0,lc3};
    Point(151) = {-hx+d/2-radius1, -hy+slab_h+2*radius1y, 0,lc3};
    Point(152) = {-hx+d/2+radius1, -hy+slab_h, 0,lc3};
    Point(153) = {-hx+d/2-radius1, -hy+slab_h, 0,lc3};

    If(slab_w_full == 1)
        Line(6) = {11,250};
        Line(8) = {13,251};
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

        Line Loop(35) = {5, 9, 30, 31, -13, -1};
        Plane Surface(36) = {35};
        Line Loop(37) = {9, -29, -22, -6};
        Plane Surface(38) = {37};
        Line Loop(39) = {13, 32, 33, 12, -7, -2};
        Plane Surface(40) = {39};
        Line Loop(41) = {12, 8, -23, -34};
        Plane Surface(42) = {41};
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

        Physical Line(55) = {5, 6, 21};
        Physical Line(56) = {26, 25};
        Physical Line(57) = {24, 8, 7};
        Physical Line(58) = {2, 1};

        Physical Surface(1) = {36, 40, 42, 38};
        Physical Surface(2) = {48, 46, 44, 50};
        Physical Surface(3) = {54, 52};
    EndIf

    If(slab_w_full == 0)
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

        Line Loop(31) = {5, 9, 24, 17, -13, -1};
        Plane Surface(32) = {31};
        Line Loop(33) = {13, 18, 19, 12, -7, -2};
        Plane Surface(34) = {33};
        Line Loop(35) = {12, 8, -4, 29, 30, -20};
        Plane Surface(36) = {35};
        Line Loop(37) = {23, -9, 6, 3, 25, -26};
        Plane Surface(38) = {37};
        Line Loop(39) = {26, -25, -27, -16, 22};
        Plane Surface(40) = {39};
        Line Loop(41) = {21, 16, 28, 29, 30};
        Plane Surface(42) = {41};
        Line Loop(43) = {15, -21, -20, -11};
        Plane Surface(44) = {43};
        Line Loop(45) = {11, -19, -18, 14};
        Plane Surface(46) = {45};
        Line Loop(47) = {14, -10, 24, 17};
        Plane Surface(48) = {47};
        Line Loop(49) = {10, 15, 22, 23};
        Plane Surface(50) = {49};

        Physical Line(51) = {5, 6};
        Physical Line(52) = {3, 27, 28, 4};
        Physical Line(53) = {8, 7};
        Physical Line(54) = {2, 1};

        Physical Surface(1) = {32, 34, 36, 38};
        Physical Surface(2) = {48, 46, 44, 50};
        Physical Surface(3) = {40, 42};
    EndIf
EndIf

