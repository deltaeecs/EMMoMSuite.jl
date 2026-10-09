// PEC sphere, r = 0.1 m (unit: m)
SetFactory("OpenCASCADE");
Sphere(1) = {0, 0, 0, 0.1};

// physical group so mesh tags carry the surface label
Physical Surface("body") = Boundary{ Volume{1}; };
