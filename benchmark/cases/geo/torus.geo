// PEC torus, ring radius 0.4 m, tube radius 0.1 m
SetFactory("OpenCASCADE");
Torus(1) = {0, 0, 0, 0.4, 0.1};

// physical group so mesh tags carry the surface label
Physical Surface("body") = Boundary{ Volume{1}; };
