// PEC sphere, r = 0.5 m (unit: m)
SetFactory("OpenCASCADE");
Sphere(1) = {0, 0, 0, 0.5};

// physical group so mesh tags carry the surface label
Physical Surface("body") = Boundary{ Volume{1}; };

// physical volume so tetrahedra carry the region label (volume cases)
Physical Volume("body") = {1};
