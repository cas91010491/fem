// Gmsh project created on Mon Aug 30 16:40:18 2021
SetFactory("OpenCASCADE");

L = 1.0;	// Length of the cylinder
a = 3.0;	// side of the inner square
r = 5.0;	// Radius of the cylinder
n = 3;		// number of quads in the middle square


na = n;
avgA = (a + Pi*r/2)/2;
//nL = Round(n*L/avgA);
nL = 1;
//nr = Round(n*(r-a/2/Sqrt(2))/avgA);
nr = 1;
h = r/Sqrt(2);



Point(1) = {-L/2, -0,  0, 1.0};
Point(2) = {-L/2,  a/2, -a/2, 1.0};
Point(3) = {-L/2,  a/2,  a/2, 1.0};
Point(4) = {-L/2, -a/2,  a/2, 1.0};
Point(5) = {-L/2, -a/2, -a/2, 1.0};
Point(6) = {-L/2,  h, -h, 1.0};
Point(7) = {-L/2,  h,  h, 1.0};
Point(8) = {-L/2, -h,  h, 1.0};
Point(9) = {-L/2, -h, -h, 1.0};


Line(1) = {2, 3};
Line(2) = {3, 4};
Line(3) = {4, 5};
Line(4) = {5, 2};
Line(5) = {2, 6};
Line(6) = {3, 7};
Line(7) = {4, 8};
Line(8) = {5, 9};
Circle(9) = {6, 1, 7};
Circle(10) = {7, 1, 8};
Circle(11) = {8, 1, 9};
Circle(12) = {9, 1, 6};
Delete {
  Point{1}; 
}


Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Curve Loop(2) = {1, 6, -9, -5};
Plane Surface(2) = {2};
Curve Loop(3) = {2, 7, -10, -6};
Plane Surface(3) = {3};
Curve Loop(4) = {7, 11, -8, -3};
Plane Surface(4) = {4};
Curve Loop(5) = {4, 5, -12, -8};
Plane Surface(5) = {5};


Extrude {L, 0, 0} {
  Point{2}; Point{3}; Point{4}; Point{5}; Point{6}; Point{7}; Point{8}; Point{9}; Curve{1}; Curve{2}; Curve{3}; Curve{4}; Curve{5}; Curve{6}; Curve{7}; Curve{8}; Curve{9}; Curve{10}; Curve{11}; Curve{12}; Surface{1}; Surface{2}; Surface{3}; Surface{4}; Surface{5}; 
}
//+
Transfinite Curve {1, 2, 3, 4, 9, 10, 11, 12, 15, 17, 19, 20, 24, 28, 30, 32} = na+1 Using Progression 1;
//+
//Transfinite Curve {5, 6, 7, 8, 25, 22, 27, 31} = nr+1 Using Progression 1;
Transfinite Curve {5, 6, 7, 8, 25, 22, 27, 31} = nr+1 Using Progression 1;
//+
Transfinite Curve {13, 14, 16, 18, 23, 21, 26, 29} = nL+1 Using Progression 1;


Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";

