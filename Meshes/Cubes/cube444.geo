// Gmsh project created on Wed Sep 22 11:21:01 2021
SetFactory("OpenCASCADE");
//+
Box(1) = {-2, -2, -2, 4, 4, 4};
//+

Transfinite Curve "*" = 5 Using Progression 1;
Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";
