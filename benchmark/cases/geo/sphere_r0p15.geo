// PEC/dielectric sphere, r = 0.15 m (unit: m)
SetFactory("OpenCASCADE");
Sphere(1) = {0, 0, 0, 0.15};

// physical group so mesh tags carry the surface label
Physical Surface("body") = Boundary{ Volume{1}; };

// physical volume so tetrahedra carry the region label (volume cases)
Physical Volume("body") = {1};
