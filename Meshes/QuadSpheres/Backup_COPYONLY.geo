// Gmsh project created on Fri Sep  3 15:46:47 2021
SetFactory("OpenCASCADE");
//+

a = 0.5;
R = 1;
r = R/Sqrt(3);

n = 6;


Point(1) = {0, 0, 0, 0.0};
Point(2) = {a/2, a/2, a/2, 1.0};
Point(3) = {-a/2, a/2, a/2, 1.0};
Point(4) = {-a/2, -a/2, a/2, 1.0};
Point(5) = {a/2, -a/2, a/2, 1.0};
Point(6) = {a/2, a/2, -a/2, 1.0};
Point(7) = {-a/2, a/2, -a/2, 1.0};
Point(8) = {-a/2, -a/2, -a/2, 1.0};
Point(9) = {a/2, -a/2, -a/2, 1.0};

Point(10) = {r, r, r, 1.0};
Point(11) = {-r, r, r, 1.0};
Point(12) = {-r, -r, r, 1.0};
Point(13) = {r, -r, r, 1.0};
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
Circle(13) = {11, 1, 10};
//+
Circle(14) = {10, 1, 14};
//+
Circle(15) = {11, 1, 15};
//+
Circle(16) = {14, 1, 15};
//+
Circle(17) = {11, 1, 12};
//+
Circle(18) = {12, 1, 16};
//+
Circle(19) = {16, 1, 15};
//+
Circle(20) = {16, 1, 17};
//+
Circle(21) = {17, 1, 14};
//+
Circle(22) = {12, 1, 13};
//+
Circle(23) = {13, 1, 17};
//+
Circle(24) = {10, 1, 13};
//+
Line(25) = {3, 11};
//+
Line(26) = {4, 12};
//+
Line(27) = {7, 15};
//+
Line(28) = {8, 16};
//+
Line(29) = {9, 17};
//+
Line(30) = {6, 14};
//+
Line(31) = {2, 10};
//+
Line(32) = {5, 13};
//+
Curve Loop(1) = {2, 7, -10, -6};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {4, 5, -12, -8};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {1, 2, 3, 4};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {9, 10, 11, 12};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {5, 9, -6, -1};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {3, 8, -11, -7};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {1, 25, 13, -31};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {2, 27, -15, -25};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {3, 30, 16, -27};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {4, 31, 14, -30};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {9, 26, 22, -32};
//+
Plane Surface(11) = {11};
//+
Curve Loop(12) = {10, 28, -18, -26};
//+
Plane Surface(12) = {12};
//+
Curve Loop(13) = {11, 29, -20, -28};
//+
Plane Surface(13) = {13};
//+
Curve Loop(14) = {12, 32, 23, -29};
//+
Plane Surface(14) = {14};
//+
Curve Loop(15) = {5, 32, -24, -31};
//+
Plane Surface(15) = {15};
//+
Curve Loop(16) = {6, 26, -17, -25};
//+
Plane Surface(16) = {16};
//+
Curve Loop(17) = {7, 28, 19, -27};
//+
Plane Surface(17) = {17};
//+
Curve Loop(18) = {8, 29, 21, -30};
//+
Plane Surface(18) = {18};
//+



//+
Sphere(1) = {-0, -0, 0, 1.0, -Pi/2, Pi/2, 2*Pi};
//+
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/4} {
  Volume{1}; 
}
//+
BooleanFragments{ Curve{15}; Curve{14}; Curve{16}; Curve{13}; Curve{22}; Curve{18}; Curve{23}; Curve{20}; Curve{21}; Curve{24}; Curve{17}; Curve{19}; Surface{19}; Delete; }{ }

//+
Delete {
  Volume{1}; 
}
//+
Delete {
  Curve{34}; 
}
//+
Surface Loop(2) = {4, 3, 2, 5, 1, 6};
//+
Volume(1) = {2};
//+
Surface Loop(3) = {3, 8, 9, 10, 7, 21};
//+
Volume(2) = {3};
//+
Surface Loop(4) = {2, 10, 14, 15, 18, 22};
//+
Volume(3) = {4};

//+
Surface Loop(6) = {4, 14, 12, 13, 11, 25};
//+
Volume(5) = {6};
//+
Surface Loop(7) = {1, 12, 17, 16, 8, 24};
//+
Volume(6) = {7};

//+
Delete {
  Surface{19}; 
}
//+
Delete {
  Curve{37}; Curve{14}; Curve{36}; 
}
//+
Delete {
  Curve{14}; Curve{34}; Curve{36}; Curve{37}; 
}
//+
Delete {
  Curve{14}; Curve{37}; Curve{36}; 
}


Sphere(8) = {-0, -0, 0, 1.0, -Pi/2, Pi/2, 2*Pi};




//+
Rotate {{0, 1, 0}, {0, 0, 0}, Pi/2} {
  Volume{8}; 
}
Rotate {{1, 0, 0}, {0, 0, 0}, Pi/4} {
  Volume{8}; 
}

BooleanFragments{ Curve{15}; Curve{14}; Curve{16}; Curve{13}; Curve{22}; Curve{18}; Curve{23}; Curve{20}; Curve{21}; Curve{24}; Curve{17}; Curve{19}; Surface{26}; Delete; }{ }
//+
Delete {
  Surface{26}; 
  Surface{28};
  Surface{32};
  Surface{20};
  Surface{23};
  Surface{27};
  Surface{30};
}
Delete {
  Volume{8}; 
}

Delete {
  Point{19}; Point{18}; Point{20}; Point{21}; 
}





//+
Surface Loop(9) = {6, 9, 18, 17, 13, 29};
//+
Volume(4) = {9};
//+
Surface Loop(10) = {5, 15, 7, 11, 16, 31};
//+
Volume(7) = {10};




//+
//+
Delete {
  Curve{14}; Curve{36}; Curve{37};  Curve{41}; Curve{39}; Curve{42}; Surface{26};
}


//+
Transfinite Line {13, 14, 15, 16, 20, 18, 23, 22, 21, 24, 19, 17, 30, 31, 4, 3, 1, 2, 27, 25, 6, 7, 10, 28, 26, 9, 11, 5, 8, 12, 32, 29} = n+1;

//+



Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";

Delete {
  Curve{39};
}
//+

Delete {
  Point{1}; Point{18}; Point{19}; Point{20}; Point{21};
}
//+
Delete {
  Point{19}; Point{18}; Point{20}; Point{21}; 
}
//+

