// Gmsh project created on Wed Sep 22 11:21:01 2021
SetFactory("OpenCASCADE");
//+
Box(1) = {-2.0, -2.0, -2.0, 4, 4, 4};
//+

n = 9;




Transfinite Curve "*" = n+1 Using Progression 1;
Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";
