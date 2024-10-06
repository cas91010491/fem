// Gmsh project created on Wed Sep 22 11:21:01 2021
SetFactory("OpenCASCADE");
//+
Box(1) = {-1.5, -1.5, -1.5, 3, 3, 3};
//+

Transfinite Curve "*" = 4 Using Progression 1;
Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";
