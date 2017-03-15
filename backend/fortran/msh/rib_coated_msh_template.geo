// Coated rib waveguide mesh template.

// Meshing parameters
lc = 0; // base level meshing parameter
lc2 = lc/1; // inclusion surfaces
lc3 = lc/1; // inclusion centres

// Unit cell
d_in_nm = 1000;
dy_in_nm = 600;
d = 1; // normalised unitcell limits x
h = dy_in_nm/d_in_nm; // unitcell limit y

Point(1) = {0, 0, 0, lc};
Point(2) = {0, -h, 0, lc};
Point(3) = {d, -h, 0, lc};
Point(4) = {d, 0, 0,lc};

// Bottom slab
slab1 = 100;
s1 = slab1/d_in_nm;
Point(5) = {0, -h+s1, 0, lc};
Point(6) = {d, -h+s1, 0, lc};

// Waveguide slab
slab2 = 50;
s2 = slab2/d_in_nm;
Point(7) = {0, -h+s1+s2, 0, lc};
Point(8) = {d, -h+s1+s2, 0, lc};

// Rib
ribx = 200;
riby = 30;
rx = ribx/d_in_nm;
ry = riby/d_in_nm;
Point(9) = {d/2-rx/2, -h+s1+s2, 0, lc};
Point(10) = {d/2+rx/2, -h+s1+s2, 0, lc};
Point(11) = {d/2-rx/2, -h+s1+s2+ry, 0, lc};
Point(12) = {d/2+rx/2, -h+s1+s2+ry, 0, lc};

// Coating
coatx = 20; // around rib
coaty = 20;
cx = coatx/d_in_nm;
cy = coaty/d_in_nm;
Point(13) = {0, -h+s1+s2+cy, 0, lc};
Point(14) = {d, -h+s1+s2+cy, 0, lc};
Point(15) = {d/2-rx/2-cx, -h+s1+s2+cy, 0, lc};
Point(16) = {d/2+rx/2+cx, -h+s1+s2+cy, 0, lc};
Point(17) = {d/2-rx/2-cx, -h+s1+s2+2*cy, 0, lc};
Point(18) = {d/2+rx/2+cx, -h+s1+s2+2*cy, 0, lc};