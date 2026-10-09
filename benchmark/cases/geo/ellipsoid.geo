// PEC ellipsoid, semi-axes 0.5/0.35/0.25 m
// (gmsh OCC geo has no Ellipsoid keyword: build from a unit-radius sphere, dilated to the semi-axes)
SetFactory("OpenCASCADE");
Sphere(1) = {0, 0, 0, 0.5};
Dilate {{0, 0, 0}, {1, 0.7, 0.5}}{ Volume{1}; }

// physical group so mesh tags carry the surface label
Physical Surface("body") = Boundary{ Volume{1}; };
