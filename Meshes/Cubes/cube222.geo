// Gmsh project created on Wed Sep 22 11:21:01 2021
SetFactory("OpenCASCADE");
//+
Box(1) = {-1, -1, -1, 2, 2, 2};
//+
//Transfinite Curve {12, 11, 9, 10, 6, 2, 4, 8, 7, 5, 1, 3} = 3 Using Progression 1;
Transfinite Curve "*" = 3 Using Progression 1;


Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";
