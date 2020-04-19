// Template mesh geometry file for a single suspended inclusion.
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
lc = 0; // 0.501 0.201 0.0701;
lc_refine_1 = lc/1; // on cylinder surfaces
lc_refine_2 = lc/1; // cylinder1 centres

hy = dy; // Thickness: square profile => hy=d
hx = 0.;


Point(1) = {0, 0, 0, lc};
Point(2) = {-hx, -hy, 0, lc};
Point(3) = {-hx+d, -hy, 0, lc};
Point(4) = {d, 0, 0,lc};

// Vertices
Point(5) = {-hx+d/2, -hy/2, 0, lc_refine_2};
Point(6) = {-hx+d/2, -hy/2+radius1y, 0, lc_refine_1};
Point(7) = {-hx+d/2-radius1, -hy/2, 0, lc_refine_1};
Point(8) = {-hx+d/2, -hy/2-radius1y, 0, lc_refine_1};
Point(9) = {-hx+d/2+radius1, -hy/2, 0, lc_refine_1};

Point(10) = {-hx+d/2, 0, 0, lc};
Point(11) = {0,-hy/2, 0, lc};
Point(12) = {-hx+d/2, -hy, 0, lc};
Point(13) = {d, -hy/2, 0, lc};
Line(1) = {1,10};
Line(2) = {10,4};
Line(3) = {2,12};
Line(4) = {12,3};
Line(5) = {1,11};
Line(6) = {11,2};
Line(7) = {4,13};
Line(8) = {13,3};
Line(9) = {11,7};
Line(10) = {7,5};
Line(11) = {5,9};
Line(12) = {9,13};
Line(13) = {10,6};
Line(14) = {6,5};
Line(15) = {5,8};
Line(16) = {8,12};


If(rect == 0)
    Ellipsis(17) = {9,5,6,6};
    Ellipsis(18) = {6,5,7,7};
    Ellipsis(19) = {7,5,8,8};
    Ellipsis(20) = {8,5,9,9};

    Line Loop(21) = {7,-12,17,-13,2};
    Plane Surface(22) = {21};
    Line Loop(23) = {1,13,18,-9,-5};
    Plane Surface(24) = {23};
    Line Loop(25) = {6,3,-16,-19,-9};
    Plane Surface(26) = {25};
    Line Loop(27) = {4,-8,-12,-20,16};
    Plane Surface(28) = {27};
    Line Loop(29) = {17,14,11};
    Plane Surface(30) = {29};
    Line Loop(31) = {18,10,-14};
    Plane Surface(32) = {31};
    Line Loop(33) = {19,-15,-10};
    Plane Surface(34) = {33};
    Line Loop(35) = {20,-11,15};
    Plane Surface(36) = {35};

    Physical Line(1) = {1,2};
    Physical Line(2) = {3,4};
    Physical Line(3) = {5,6};
    Physical Line(4) = {7,8};

    Physical Surface(1) = {24,26,28,22};
    Physical Surface(2) = {30,32,34,36};
EndIf

If(rect == 1)
    Point(150) = {-hx+d/2+radius1, -hy/2+radius1y, 0,lc_refine_1};
    Point(151) = {-hx+d/2-radius1, -hy/2+radius1y, 0,lc_refine_1};
    Point(152) = {-hx+d/2+radius1, -hy/2-radius1y, 0,lc_refine_1};
    Point(153) = {-hx+d/2-radius1, -hy/2-radius1y, 0,lc_refine_1};

    Line(17) = {151, 6};
    Line(18) = {6, 150};
    Line(19) = {150, 9};
    Line(20) = {9, 152};
    Line(21) = {152, 8};
    Line(22) = {8, 153};
    Line(23) = {153, 7};
    Line(24) = {7, 151};

    Line Loop(25) = {5, 9, 24, 17, -13, -1};
    Plane Surface(26) = {25};
    Line Loop(27) = {6, 3, -16, 22, 23, -9};
    Plane Surface(28) = {27};
    Line Loop(29) = {16, 4, -8, -12, 20, 21};
    Plane Surface(30) = {29};
    Line Loop(31) = {7, -12, -19, -18, -13, 2};
    Plane Surface(32) = {31};
    Line Loop(33) = {24, 17, 14, -10};
    Plane Surface(34) = {33};
    Line Loop(35) = {14, 11, -19, -18};
    Plane Surface(36) = {35};
    Line Loop(37) = {11, 20, 21, -15};
    Plane Surface(38) = {37};
    Line Loop(39) = {22, 23, 10, 15};
    Plane Surface(40) = {39};

    Physical Line(1) = {1,2};
    Physical Line(2) = {3,4};
    Physical Line(3) = {5,6};
    Physical Line(4) = {7,8};

    Physical Surface(1) = {26, 32, 30, 28};
    Physical Surface(2) = {34, 36, 38, 40};
EndIf
