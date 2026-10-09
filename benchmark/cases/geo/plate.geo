// PEC thin plate, 1 m x 1 m x 0.02 m (unit: m)
SetFactory("OpenCASCADE");
Box(1) = {-0.5, -0.5, -0.01, 1, 1, 0.02};

// physical group so mesh tags carry the surface label
Physical Surface("body") = Boundary{ Volume{1}; };
