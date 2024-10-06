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
ne  = 4;			// number of elements per edge of pseudo-cube face
nr = 3 ; 			// number of elements along sphere-cylinder radii
avgA = (a + Pi*r/2)/2;
nl = Round(ne*L/avgA);


Point(1) = {0, 0, 0, 0.0};
Point(2) = {a/2, a/2, 0, 1.0};
Point(3) = {-a/2, a/2, 0, 1.0};
Point(4) = {-a/2, -a/2, 0, 1.0};
Point(5) = {a/2, -a/2, 0, 1.0};
Point(6) = {a/2, a/2, -a/2, 1.0};
Point(7) = {-a/2, a/2, -a/2, 1.0};
Point(8) = {-a/2, -a/2, -a/2, 1.0};
Point(9) = {a/2, -a/2, -a/2, 1.0};

Point(10) = {r2, r2, 0, 1.0};
Point(11) = {-r2, r2, 0, 1.0};
Point(12) = {-r2, -r2, 0, 1.0};
Point(13) = {r2, -r2, 0, 1.0};
Point(14) = {r, r, -r, 1.0};
Point(15) = {-r, r, -r, 1.0};
Point(16) = {-r, -r, -r, 1.0};
Point(17) = {r, -r, -r, 1.0};

//+
Line(1) = {2, 3};
//+
Line(2) = {3, 7};
//+
Line(3) = {7, 6};
//+
Line(4) = {6, 2};
//+
Line(5) = {2, 5};
//+
Line(6) = {3, 4};
//+
Line(7) = {7, 8};
//+
Line(8) = {6, 9};
//+
Line(9) = {5, 4};
//+
Line(10) = {4, 8};
//+
Line(11) = {8, 9};
//+
Line(12) = {9, 5};
//+
Line(13) = {2, 10};
//+
Line(14) = {3, 11};
//+
Line(15) = {4, 12};
//+
Line(16) = {5, 13};
//+
Line(17) = {6, 14};
//+
Line(18) = {7, 15};
//+
Line(19) = {8, 16};
//+
Line(20) = {9, 17};


//+
Circle(21) = {11, 1, 10};
//+
Circle(22) = {10, 1, 14};
//+
Circle(23) = {14, 1, 15};
//+
Circle(24) = {11, 1, 15};
//+
Circle(25) = {10, 1, 13};
//+
Circle(26) = {13, 1, 17};
//+
Circle(27) = {17, 1, 14};
//+
Circle(28) = {13, 1, 12};
//+
Circle(29) = {12, 1, 16};
//+
Circle(30) = {16, 1, 17};
//+
Circle(31) = {12, 1, 11};
//+
Circle(32) = {15, 1, 16};
//+
Curve Loop(1) = {1, 14, 21, -13};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {13, 25, -16, -5};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {9, 15, -28, -16};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {6, 15, 31, -14};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {1, 6, -9, -5};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {1, 2, 3, 4};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {5, -12, -8, 4};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {9, 10, 11, 12};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {6, 10, -7, -2};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {3, 8, -11, -7};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {3, 17, 23, -18};
//+
Plane Surface(11) = {11};
//+
Curve Loop(12) = {8, 20, 27, -17};
//+
Plane Surface(12) = {12};
//+
Curve Loop(13) = {11, 20, -30, -19};
//+
Plane Surface(13) = {13};
//+
Curve Loop(14) = {7, 19, -32, -18};
//+
Plane Surface(14) = {14};
//+
Curve Loop(15) = {14, 24, -18, -2};
//+
Plane Surface(15) = {15};
//+
Curve Loop(16) = {4, 13, 22, -17};
//+
Plane Surface(16) = {16};
//+
Curve Loop(17) = {12, 16, 26, -20};
//+
Plane Surface(17) = {17};
//+
Curve Loop(18) = {10, 19, -29, -15};
//+
Plane Surface(18) = {18};
//+

//+
Sphere(1) = {-0, -0, 0, R, -Pi/2, Pi/2, 2*Pi};
//+
//+
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/4} { Volume{1}; }


//+
Sphere(2) = {-0, -0, 0, R, -Pi/2, Pi/2, 2*Pi};
Rotate {{0, 1, 0}, {0, 0, 0}, -Pi/2} { Volume{2}; }


//+
BooleanFragments{ Curve{21}; Curve{22}; Curve{23}; Curve{24}; Curve{25}; Curve{26}; Curve{27}; Curve{28}; Curve{29}; Curve{30}; Curve{31}; Curve{32};} {Surface{19}; Curve{34}; Delete;}



//+
Recursive Delete {
  Volume{1}; 
  Surface{20};
  Point{1};
  Surface{24};
}


BooleanFragments{ Curve{23}; Curve{27}; Curve{30}; Curve{32}; } {Surface{20}; Delete;}



//+
Recursive Delete {
  Volume{2};
  Surface{27};
  Surface{21};
}




//+
Surface Loop(2) = {5, 8, 10, 7, 9, 6};
//+
Volume(1) = {2};
//+
Surface Loop(3) = {8, 13, 17, 26, 3, 13};
//+
Volume(2) = {3};
//+
Surface Loop(4) = {17, 12, 23, 2, 16, 7};
//+
Volume(3) = {4};
//+
Surface Loop(5) = {16, 11, 22, 1, 15, 6};
//+
Volume(4) = {5};
//+
Surface Loop(6) = {14, 4, 25, 15, 18, 9};
//+
Volume(5) = {6};
//+
Surface Loop(7) = {10, 28, 14, 11, 12, 13};
//+
Volume(6) = {7};


//+
Extrude {0, 0, 1} {
  Surface{4}; Surface{1}; Surface{2}; Surface{3}; Surface{5}; Curve{9}; Curve{6}; Curve{1}; Curve{5}; Curve{15}; Curve{14}; Curve{13}; Curve{16}; Curve{28}; Curve{31}; Curve{21}; Curve{25}; 
}

//+
Transfinite Line {1,3,5,6,7,8,9,11, 21,23,25,27,28,30,31,32, 35,42,50,51, 39,44,47,52} = ne+1;
//+
Transfinite Line {2,4,10,12, 22,24,26,29} = Round((ne+1)/2);
//+
Transfinite Line {13,14,15,16,17,18,19,20, 37,40,45,49} = nr+1;
Transfinite Line {33,34,36,38,41,43,46,48} = nl+1;

//+
Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";









