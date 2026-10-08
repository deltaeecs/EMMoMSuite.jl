// PEC cylinder, r = 0.25 m, h = 1 m (unit: m)
SetFactory("OpenCASCADE");
Cylinder(1) = {0, 0, 0, 0, 0, 1, 0.25};

// physical group so mesh tags carry the surface label
Physical Surface("body") = Boundary{ Volume{1}; };
