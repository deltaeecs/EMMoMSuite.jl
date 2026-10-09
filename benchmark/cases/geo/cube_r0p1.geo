// uniform dielectric cube, side = 0.2 m (unit: m)
SetFactory("OpenCASCADE");
Box(1) = {-0.1, -0.1, -0.1, 0.2, 0.2, 0.2};

// physical volume so tetrahedra carry the region label (volume cases)
Physical Volume("body") = {1};
