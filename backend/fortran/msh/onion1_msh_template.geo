// Template mesh geometry file for a single suspended inclusion.
// Inclusion can be circular/elliptical (default), or square/rectangular.

d = 1; // grating period
d_in_nm = 1000;
dy_in_nm = 1000;
dy = dy_in_nm/d_in_nm;
a1 = 100;
a2 = 200;
a3 = 300;
rad1 = (a1/(2*d_in_nm))*d;
rad2 = ((a1+a2)/(2*d_in_nm))*d;
rad3 = ((a1+a2+a3)/(2*d_in_nm))*d;
lc = 0; // 0.501 0.201 0.0701;
lc2 = lc/1; // on cylinder surfaces
lc3 = lc/1; // cylinder1 centres

hy = dy; // Thickness: square profile => hy=d
hx = 0.;


//Box part a 
Point(1) = {0, 0, 0, lc};        // NW
Point(2) = {-hx, -hy, 0, lc};    // SW
Point(3) = {-hx+d, -hy, 0, lc};  // SE
Point(4) = {d, 0, 0,lc};         // E

//Box part b 
Point(6) = {-hx+d/2, 0, 0, lc};       //N
Point(7) = {0,-hy/2, 0, lc};          //W
Point(8) = {-hx+d/2, -hy, 0, lc};     //S
Point(9) = {d, -hy/2, 0, lc};         //E 

// Circle 1 vertices
Point(10) = {-hx+d/2, -hy/2, 0, lc3};        //Center
Point(11) = {-hx+d/2, -hy/2+rad1, 0, lc2};   //N
Point(12) = {-hx+d/2-rad1, -hy/2, 0, lc2};   //W
Point(13) = {-hx+d/2, -hy/2-rad1, 0, lc2};   //S
Point(14) = {-hx+d/2+rad1, -hy/2, 0, lc2};   //E


// Outer Box
Line(1) = {1,6};  // top left
Line(2) = {6,4};  // top right
Line(3) = {2,8};  // bottom left
Line(4) = {8,3};  // bottom right
Line(5) = {1,7};  // left upper
Line(6) = {7,2};  // left lower
Line(7) = {4,9};  // right upper
Line(8) = {9,3};  // right lower

//Circle 1 radii 
Line(11) = {11,10};  // N
Line(12) = {12,10};  // W
Line(13) = {10,13};  // S
Line(14) = {10,14};  // E

// Box to outer circle
Line(2001) = {6,11}; //N
Line(2002) = {7,12}; //W
Line(2003) = {13,8}; //S
Line(2004) = {14,9}; //E

//Circle 1
Ellipsis(211) = {14,10,11,11};  //NE
Ellipsis(212) = {11,10,12,12};  //NW
Ellipsis(213) = {12,10,13,13};  //SW
Ellipsis(214) = {13,10,14,14};  //SE

//
//Surfaces 
//

// Sectors circle 1
Line Loop(511) = {11, 14, 211};  //NE
Plane Surface(516) = {511};
Line Loop(512) = {212, 12, -11}; //NW 
Plane Surface(517) = {512};
Line Loop(513) = {13, -213, 12}; //SW
Plane Surface(518) = {513};
Line Loop(514) = {14, -214, -13}; //SE
Plane Surface(519) = {514};

// Outer circle to boundary surfaces
Line Loop(501) = {2001, -211, 2004, -7, -2}; //NE
Plane Surface(506) = {501};
Line Loop(502) = {212, -2002, -5, 1, 2001};  //NW
Plane Surface(507) = {502};
Line Loop(503) = {2003, -3, -6, 2002, 213};   //SW
Plane Surface(508) = {503};
Line Loop(504) = {2004, 8, -4, -2003, 214};  //SE
Plane Surface(509) = {504};

//
//Physical features
//

// Outer boundary
Physical Line(1) = {1,2};  // Top
Physical Line(2) = {3,4};  // Bottom
Physical Line(3) = {5,6};  // Left
Physical Line(4) = {7,8};  // Right

Physical Surface(1) = {507,506,509,508};     //Outer region
Physical Surface(2) = {517, 516, 519, 518};  //Circle 1



