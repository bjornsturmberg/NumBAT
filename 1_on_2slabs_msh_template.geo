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

slab_width = 100;
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

// Inclusion
Point(5) = {-hx+d/2, -hy+radius1y+slab_h+slab2_h, 0, lc3};
Point(6) = {-hx+d/2, -hy+2*radius1y+slab_h+slab2_h, 0, lc2};
Point(7) = {-hx+d/2-radius1, -hy+radius1y+slab_h+slab2_h, 0, lc2};
Point(8) = {-hx+d/2, -hy+slab_h+slab2_h, 0, lc2};
Point(9) = {-hx+d/2+radius1, -hy+radius1y+slab_h+slab2_h, 0, lc2};

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
        Line(22) = {250, 153};
        Line(23) = {152, 251};
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

    EndIf
EndIf

