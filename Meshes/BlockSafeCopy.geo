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
