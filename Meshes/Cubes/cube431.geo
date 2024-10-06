// Gmsh project created on Wed Sep 22 11:21:01 2021
SetFactory("OpenCASCADE");
//+
Box(1) = {-2, -1.5, -0.5, 4, 3, 1};
//+

Transfinite Curve {1, 3, 5, 7} = 2 Using Progression 1;
Transfinite Curve {2, 4, 6, 8} = 4 Using Progression 1;
Transfinite Curve {9, 10, 11, 12} = 5 Using Progression 1;

Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";
//+

