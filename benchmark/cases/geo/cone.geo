// PEC truncated cone, r1 = 0.4 m -> r2 = 0.05 m, h = 1 m
SetFactory("OpenCASCADE");
Cone(1) = {0, 0, 0, 0, 0, 1, 0.4, 0.05};

// physical group so mesh tags carry the surface label
Physical Surface("body") = Boundary{ Volume{1}; };
