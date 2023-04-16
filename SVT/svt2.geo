//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {-1, -1, 0, 1.0};
//+
Point(3) = {0, -1, 0, 1.0};
//+
Point(4) = {0, 0, 0, 1.0};
//+
Point(5) = {-1, -0, 0, 1.0};
//+
Line(1) = {2, 5};
//+
Line(2) = {5, 1};
//+
Line(3) = {1, 3};
//+
Line(4) = {3, 2};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Transfinite Curve {1} = 10 Using Progression 1;
//+
Transfinite Curve {1} = 10 Using Progression 1;
//+
Physical Surface(1) -= {1};
//+
Line(5) = {3, 5};
//+
Transfinite Curve {5} = 10 Using Progression 1;
//+
Transfinite Curve {2, 1} = 10 Using Progression 1.2;
//+
Transfinite Curve {5, 1} = 10 Using Progression -1.2;
//+
Line(6) = {5, 1};
//+
Line(7) = {5, 2};
//+
Line(8) = {3, 2};
//+
Line(9) = {3, 1};
//+
Line(10) = {5, 3};
//+
Transfinite Curve {7, 6, 9, 8, 10} = 10 Using Progression 1.2;
//+
Curve Loop(2) = {7, -8, -10};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {6, -9, -10};
//+
Plane Surface(3) = {3};
//+
Transfinite Curve {7, 8, 9, 6} = 10 Using Progression 1.2;
//+
Transfinite Curve {10} = 10 Using Progression 1;
