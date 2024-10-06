// Gmsh project created on Wed Sep 22 11:21:01 2021
SetFactory("OpenCASCADE");
//+
Box(1) = {-4, -2, -1.5, 8, 4, 3};
//+

Transfinite Curve {1, 3, 5, 7} = 7 Using Progression 1;
Transfinite Curve {2, 4, 6, 8} = 9 Using Progression 1;
Transfinite Curve {9, 10, 11, 12} = 17 Using Progression 1;

Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";
//+

