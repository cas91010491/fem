//+
SetFactory("OpenCASCADE");


dx = 2;
dy = 1;
dz = 2;

nx = 10;
ny = 1;
nz = 10;




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
