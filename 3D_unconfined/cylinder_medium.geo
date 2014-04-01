MS=0.25;

Point(1) = {0,0,0,MS};
Point(2) = {1,0,0,MS};
Point(3) = {0,1,0,MS};
Point(4) = {-1,0,0,MS};
Point(5) = {0,-1,0,MS};

Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};

Line Loop(5) = {1,2,3,4};
Plane Surface(6) = {5};

Extrude {0,0,1} {
  Surface{6};
}
