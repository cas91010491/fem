// Gmsh project created on Wed Sep 22 11:21:01 2021
SetFactory("OpenCASCADE");
//+
Box(1) = {-0.5, -0.5, -0.5, 1, 1, 1};
//+

Transfinite Curve "*" = 2 Using Progression 1;
Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";
