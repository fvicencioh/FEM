mdc = 0.1;
mdi = 0.0025;
mdr = 0.0005;

L = 1;

r = 0.25;
a = L/10;
b = r*a;

Point(1)={L, 0, 0, mdc};
Point(2)={L, 2*L, 0, mdc};
Point(3)={0, 2*L, 0, mdc};
Point(4)={0, b, 0, mdr};
Point(5)={0, 0, 0, mdi};
Point(6)={a, 0, 0, mdr};

Line(1)={1, 2};
Line(2)={2, 3};
Line(3)={3, 4};
Ellipse(4)={4, 5, 6, 6};
Line(5)={6, 1};
Line(6)={6, 5};
Line(7)={5, 4};

Curve Loop(1) = {1, 2, 3, 4, 5};
Curve Loop(2) = {-4, -6, -7};

Plane Surface(1) = {1};
Plane Surface(2) = {2};

Physical Line("right") = {1};
Physical Line("top") = {2};
Physical Line("left") = {3};
Physical Line("boundary") = {4};
Physical Line("bottom") = {5};
Physical Line("ellipse_bottom") = {6};
Physical Line("ellipse_left") = {7};

Physical Surface("Steel") = {1};
Physical Surface("Ellipse") = {2};

Mesh 2;
Mesh.ElementOrder = 2;
Mesh.SurfaceFaces = 1;
Mesh.points = 1;
Mesh.saveAll = 1;
Mesh.SaveGroupsOfElements = 1;
