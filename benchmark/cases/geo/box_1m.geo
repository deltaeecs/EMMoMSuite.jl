// PEC box, 1 m x 1 m x 1 m (unit: m)
SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 1, 1, 1};

// physical group so mesh tags carry the surface label
Physical Surface("body") = Boundary{ Volume{1}; };
