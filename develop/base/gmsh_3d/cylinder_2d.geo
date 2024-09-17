density = DefineNumber[ 4, Name "Parameters/density" ];
density_block = DefineNumber[ 10, Name "Parameters/density_block" ];
thickness = DefineNumber[ 0.1, Name "Parameters/thickness" ];
length = DefineNumber[ 2, Name "Parameters/length" ];
nz = DefineNumber[ 3, Name "Parameters/nz" ];

radius = (2^0.5)/40;
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {radius, radius, 0, 1.0};
//+
Point(3) = {radius, -radius, 0, 1.0};
//+
Point(4) = {-radius, -radius, 0, 1.0};
//+
Point(5) = {-radius, radius, 0, 1.0};
//+



//+
Point(6) = {-0.2, 0.21, 0, 1.0};
//+
Point(7) = {-0.2, -0.2, 0, 1.0};
//+
Point(8) = {0.2, 0.21, 0, 1.0};
//+
Point(9) = {0.2, -0.2, 0, 1.0};
//+

//+
Line(1) = {7, 6};
//+
Line(2) = {6, 8};
//+
Line(3) = {8, 9};
//+
Line(4) = {9, 7};

//+
Line(5) = {7, 4};
//+
Line(6) = {6, 5};
//+
Line(7) = {8, 2};
//+
Line(8) = {9, 3};
//+
Circle(9) = {4, 1, 5};
//+
Circle(10) = {5, 1, 2};
//+
Circle(11) = {2, 1, 3};
//+
Circle(12) = {3, 1, 4};

//+
Extrude {0, 0, thickness} {
  Point{2,3,4,5,6,7,8,9}; Curve{1,2,3,4,5,6,7,8,9,10,11,12};  
}
//+
Curve Loop(1) = {1, 6, -9, -5};
//+
Plane Surface(69) = {1};
//+
Curve Loop(2) = {6, 10, -7, -2};
//+
Plane Surface(70) = {2};
//+
Curve Loop(3) = {7, 11, -8, -3};
//+
Plane Surface(71) = {3};
//+
Curve Loop(4) = {12, -5, -4, 8};
//+
Plane Surface(72) = {4};
//+
Curve Loop(5) = {21, 41, -53, -37};
//+
Plane Surface(73) = {5};
//+
Curve Loop(6) = {41, 57, -45, -25};
//+
Plane Surface(74) = {6};
//+
Curve Loop(7) = {45, 61, -49, -29};
//+
Plane Surface(75) = {7};
//+
Curve Loop(8) = {65, -37, -33, 49};
//+
Plane Surface(76) = {8};
//+
Surface Loop(1) = {69, 24, 73, 56, 40, 44};
//+
Volume(1) = {1};
//+
Surface Loop(2) = {70, 60, 74, 28, 44, 48};
//+
Volume(2) = {2};
//+
Surface Loop(3) = {71, 64, 75, 32, 48, 52};
//+
Volume(3) = {3};
//+
Surface Loop(4) = {76, 68, 72, 36, 52, 40};
//+
Volume(4) = {4};


// right geometry block
//+
Point(20) = {length, -0.2, 0, 1.0};
//+
Point(21) = {length, -0.2, thickness, 1.0};
//+
Point(22) = {length, 0.21, 0, 1.0};
//+
Point(23) = {length, 0.21, thickness, 1.0};

//+
Line(66) = {8, 22};
//+
Line(67) = {22, 23};
//+
Line(68) = {23, 16};
//+
Line(69) = {22, 20};
//+
Line(70) = {23, 21};
//+
Line(71) = {20, 21};
//+
Line(72) = {21, 17};
//+
Line(73) = {20, 9};
//+
//+
Curve Loop(9) = {68, -19, 66, 67};
//+
Plane Surface(77) = {9};
//+
Curve Loop(10) = {67, 70, -71, -69};
//+
Plane Surface(78) = {10};
//+
Curve Loop(11) = {72, -20, -73, 71};
//+
Plane Surface(79) = {11};
//+
Curve Loop(12) = {3, -73, -69, -66};
//+
Plane Surface(80) = {12};
//+
Curve Loop(13) = {29, -72, -70, 68};
//+
Plane Surface(81) = {13};
//+
Surface Loop(5) = {81, 79, 80, 78, 77, 32};
//+
Volume(5) = {5};






//+
Transfinite Curve {13, 14, 15, 16, 17, 18, 19, 20, 67, 71} = nz Using Progression 1.0;
//+
Transfinite Curve {6, 41, 5, 37, 8, 49, 7, 45} = density+1 Using Progression 0.8;
//+
Transfinite Curve {25, 2, 10, 57, 29, 3, 61, 11, 33, 4, 65, 12, 21, 1, 9, 53, 69, 70, 66, 68, 72, 73} = density Using Progression 1;
//+
Transfinite Curve {66, 68, 72, 73} = density_block Using Progression 1;

//+
Transfinite Surface {24,28,32,36,40,44,48,52,69,73,70,74,71,75,72,76,56,60,64,68,77,78,79,80,81};
//+
Transfinite Volume{1,2,3,4,5};


//+
Recombine Surface {24,28,32,36,40,44,48,52,69,73,70,74,71,75,72,76,56,60,64,68,77,78,79,80,81};
//+
Recombine Volume{1,2,3,4,5};

Periodic Surface {70, 69, 72, 71, 80} = {74, 73, 76, 75, 81} Translate {0, 0, -thickness};
//+


//+
Physical Surface("inlet", 77) = {24};
//+
Physical Surface("outlet", 78) = {78};
//+
Physical Surface("top", 79) = {28, 77};
//+
Physical Surface("bottom", 80) = {36, 79};
//+
Physical Surface("left_wall", 81) = {70, 69, 72, 71, 80};
//+
Physical Surface("right_wall", 82) = {74, 73, 76, 75, 81};
//+
Physical Surface("cylinder", 83) = {56, 60, 64, 68};
//+
Physical Volume("fluid", 84) = {1, 2, 3, 4, 5};
//+
