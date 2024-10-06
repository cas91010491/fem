// Gmsh project created on Fri Sep  3 15:46:47 2021
SetFactory("OpenCASCADE");
//+
// 	GEOMETRY
a = 0.5;			// cube edge
R = 1;				// Sphere Radius
L = 1;				// Indentor length
r = R/Sqrt(3);
r2= R/Sqrt(2);

// MESH
ne  = 3;			// number of elements per edge of pseudo-cube face
nr = 1 ; 			// number of elements along sphere-cylinder radii
alpha = 1.0;   // Simmilarity between longitudinal and external element size
b = ((1-alpha)*a + alpha*Pi*r/2)/2;
nl = Round(ne*L/b);

//+ Points defining Cube
Point(1)  = {a/2, a/2, 0, 1.0};
Point(2)  = {-a/2, a/2, 0, 1.0};
Point(3)  = {-a/2, -a/2, 0, 1.0};
Point(4)  = {a/2, -a/2, 0, 1.0};
Point(5)  = {a/2, a/2, -a/2, 1.0};
Point(6)  = {-a/2, a/2, -a/2, 1.0};
Point(7)  = {-a/2, -a/2, -a/2, 1.0};
Point(8)  = {a/2, -a/2, -a/2, 1.0};
//+ Points defining Sphere
Point(9)  = {r2, r2, 0, 1.0};
Point(10) = {-r2, r2, 0, 1.0};
Point(11) = {-r2, -r2, 0, 1.0};
Point(12) = {r2, -r2, 0, 1.0};
Point(13) = {r, r, -r, 1.0};
Point(14) = {-r, r, -r, 1.0};
Point(15) = {-r, -r, -r, 1.0};
Point(16) = {r, -r, -r, 1.0};
//+ Auxiliary point for arcs
Point(17) = {0.0, 0.0, 0.0};

//+ Cube's edges
Line(1)  = {1, 2};
Line(2)  = {2, 3};
Line(3)  = {3, 4};
Line(4)  = {4, 1};
Line(5)  = {5, 6};
Line(6)  = {6, 7};
Line(7)  = {7, 8};
Line(8)  = {8, 5};
Line(9)  = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};
//+ Cube-Sphere connectors
Line(13) = {1, 9};
Line(14) = {2, 10};
Line(15) = {3, 11};
Line(16) = {4, 12};
Line(17) = {5, 13};
Line(18) = {6, 14};
Line(19) = {7, 15};
Line(20) = {8, 16};

//+ Sphere's arcs
Circle(21) = { 9, 17, 10};
Circle(22) = {10, 17, 11};
Circle(23) = {11, 17, 12};
Circle(24) = {12, 17,  9};
Circle(25) = {13, 17, 14};
Circle(26) = {14, 17, 15};
Circle(27) = {15, 17, 16};
Circle(28) = {16, 17, 13};
Circle(29) = { 9, 17, 13};
Circle(30) = {10, 17, 14};
Circle(31) = {11, 17, 15};
Circle(32) = {12, 17, 16};
Delete{ Point{17}; }  // Not needed anymore

//+ Plane Surfaces

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
//+
Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};
//+
Curve Loop(3) = {1, 10, -5, -9};
Plane Surface(3) = {3};
//+
Curve Loop(4) = {2, 11, -6, -10};
Plane Surface(4) = {4};
//+
Curve Loop(5) = {3, 12, -7, -11};
Plane Surface(5) = {5};
//+
Curve Loop(6) = {4, 9, -8, -12};
Plane Surface(6) = {6};
//+
Curve Loop(7) = {1, 14, -21, -13};
Plane Surface(7) = {7};
//+
Curve Loop(8) = {2, 15, -22, -14};
Plane Surface(8) = {8};
//+
Curve Loop(9) = {3, 16, -23, -15};
Plane Surface(9) = {9};
//+
Curve Loop(10) = {4, 13, -24, -16};
Plane Surface(10) = {10};
//+
Curve Loop(11) = {5, 18, -25, -17};
Plane Surface(11) = {11};
//+
Curve Loop(12) = {6, 19, -26, -18};
Plane Surface(12) = {12};
//+
Curve Loop(13) = {7, 20, -27, -19};
Plane Surface(13) = {13};
//+
Curve Loop(14) = {8, 17, -28, -20};
Plane Surface(14) = {14};
//+
Curve Loop(15) = {9, 17, -29, -13};
Plane Surface(15) = {15};
//+
Curve Loop(16) = {10, 18, -30, -14};
Plane Surface(16) = {16};
//+
Curve Loop(17) = {11, 19, -31, -15};
Plane Surface(17) = {17};
//+
Curve Loop(18) = {12, 20, -32, -16};
Plane Surface(18) = {18};


//+ Spheric Surfaces

//+Auxiliary sphere to create lateral spheric surfaces
Sphere(1) = {-0, -0, 0, R, -Pi/2, Pi/2, 2*Pi};
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/4} { Volume{1}; }
//+
BooleanFragments{ Curve{21}; Curve{22}; Curve{23}; Curve{24}; Curve{25}; Curve{26}; Curve{27}; Curve{28}; Curve{29}; Curve{30}; Curve{31}; Curve{32};} {Surface{19}; Curve{34}; Delete;}
Recursive Delete { Volume{1}; Surface{20}; Surface{23}; }

//+Auxiliary sphere to create bottom spheric surface
Sphere(1) = {-0, -0, 0, R, -Pi/2, Pi/2, 2*Pi};
Rotate {{0, 1, 0}, {0, 0, 0}, -Pi/2} { Volume{1}; }
//+
BooleanFragments{ Curve{25}; Curve{26}; Curve{27}; Curve{28}; } {Surface{26}; Delete;}
Recursive Delete { Volume{1}; Surface{27}; }


//+ Volumes in the Sphere

Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};
//+
Surface Loop(2) = {3, 7, 21, 11, 15, 16};
Volume(2) = {2};
//+
Surface Loop(3) = {4, 8, 24, 12, 16, 17};
Volume(3) = {3};
//+
Surface Loop(4) = {5, 9, 25, 13, 17, 18};
Volume(4) = {4};
//+
Surface Loop(5) = {6, 10, 22, 14, 18, 15};
Volume(5) = {5};
//+
Surface Loop(6) = {2, 11, 12, 13, 14, 28};
Volume(6) = {6};


//+ Creating Cylinder   (Points, Lines, Surfaces, Volumes)

Extrude {0, 0, 1} { Surface{1}; Surface{7}; Surface{8}; Surface{9}; Surface{10}; }


//+ MESHING

//+ wide meshing
Transfinite Line {1,2,3,4,5,6,7,8, 21,22,23,24,25,26,27,28, 35,37,39,40, 44,48,51,52} = ne+1;
Transfinite Line {9,10,11,12, 29,30,31,32} = Round(ne/2)+1;
//+ radial meshing
Transfinite Line {13,14,15,16,17,18,19,20, 42,45,47,50} = nr+1;
//+ longitudinal meshing
Transfinite Line {33,34,36,38,41,43,46,49} = nl+1;

//+
Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";


