//+
SetFactory("OpenCASCADE");


dx = 15;
dy = 15;
dz = 5;

nx = 16;
ny = 16;
nz = 8;




Box(1) = {-dx/2, -dy/2, -dz/2, dx, dy, dz};
//+
Transfinite Curve {9, 11, 12, 10} = nx+1 Using Progression 1;
//+
Transfinite Curve {4, 8, 6, 2} = ny+1 Using Progression 1;
//+
Transfinite Curve {1, 5, 7, 3} = nz+1 Using Progression 1;

Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";
//+
Field[1] = Ball;
//+
Field[1].Radius = 3;
//+
Field[1].ZCenter = 2.5;
//+
Field[1].VIn = 0.1;
//+
Field[1].VOut = 0.5;
//+
Field[1].Thickness = 0.05;
//+
Characteristic Length {1} = 0.1;
