//+
radius = DefineNumber[ 1, Name "Parameters/radius" ];
//+
length = DefineNumber[ 2, Name "Parameters/length" ];
//+
density = DefineNumber[ 2, Name "Parameters/density" ];
//+
nEz = DefineNumber[ 2, Name "Parameters/nEz" ];


//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, 0, -length, 1.0};
//+
Point(3) = {0, radius, 0, 1.0};
//+
Point(4) = {0, radius, -length, 1.0};
//+
Point(5) = {0, -radius, 0, 1.0};
//+
Point(6) = {0, -radius, -length, 1.0};

//+
Point(7) = {radius, 0, 0, 1.0};
//+
Point(8) = {radius, 0, -length, 1.0};
//+
Point(9) = {-radius, 0, 0, 1.0};
//+
Point(10) = {-radius, 0, -length, 1.0};

//+
Point(11) = {0, radius/2, 0, 1.0};
//+
Point(12) = {0, -radius/2, 0, 1.0};
//+
Point(13) = {radius/2, 0, 0, 1.0};
//+
Point(14) = {-radius/2, 0, 0, 1.0};

//+
Point(15) = {0, radius/2, -length, 1.0};
//+
Point(16) = {0, -radius/2, -length, 1.0};
//+
Point(17) = {radius/2, 0, -length, 1.0};
//+
Point(18) = {-radius/2, 0, -length, 1.0};




//+
Circle(1) = {3, 1, 7};
//+
Circle(2) = {7, 1, 5};
//+
Circle(3) = {5, 1, 9};
//+
Circle(4) = {9, 1, 3};

//+
Line(5) = {3, 4};
//+
Line(6) = {5, 6};
//+
Line(7) = {7, 8};
//+
Line(8) = {9, 10};

//+
Circle(9) = {4, 2, 8};
//+
Circle(10) = {8, 2, 6};
//+
Circle(11) = {6, 2, 10};
//+
Circle(12) = {10, 2, 4};





//+
Line(13) = {11, 15};
//+
Line(14) = {13, 17};
//+
Line(15) = {12, 16};
//+
Line(16) = {14, 18};
//+
Line(17) = {11, 3};
//+
Line(18) = {15, 4};
//+
Line(19) = {13, 7};
//+
Line(20) = {17, 8};
//+
Line(21) = {12, 5};
//+
Line(22) = {16, 6};
//+
Line(23) = {14, 9};
//+
Line(24) = {18, 10};
//+

//Circle(25) = {11, 1, 13};
Line(25) = {11, 13};
//+
//Circle(26) = {13, 1, 12};
Line(26) = {13, 12};
//+
//Circle(27) = {12, 1, 14};
Line(27) = {12, 14};
//+
//Circle(28) = {14, 1, 11};
Line(28) = {14, 11};
//+
//Circle(29) = {15, 2, 17};
Line(29) = {15, 17};
//+
//Circle(30) = {17, 2, 16};
Line(30) = {17, 16};
//+
//Circle(31) = {16, 2, 18};
Line(31) = {16, 18};
//+
//Circle(32) = {18, 2, 15};
Line(32) = {18, 15};


// inner cylinder

//+
Curve Loop(1) = {25, 26, 27, 28};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {29, 30, 31, 32};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {13, 29, -14, -25};
//+
Surface(3) = {3};
//+
Curve Loop(4) = {14, 30, -15, -26};
//+
Surface(4) = {4};
//+
Curve Loop(5) = {15, 31, -16, -27};
//+
Surface(5) = {5};
//+
Curve Loop(6) = {16, 32, -13, -28};
//+
Surface(6) = {6};
//+
Surface Loop(1) = {3, 6, 5, 4, 2, 1};
//+
Volume(1) = {1};





//+
Transfinite Curve {25, 28, 27, 26, 29, 32, 31, 30} = density+1 Using Progression 1;
//+
Transfinite Curve {13, 16, 15, 14} = nEz+1 Using Progression 1;


//+
Transfinite Surface {1, 2, 3, 4, 5, 6};

//+
Transfinite Volume{1};

//+
Recombine Surface {1, 2, 3, 4, 5, 6};






// outer cylinder

//+
Curve Loop(7) = {5, -18, -13, 17};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {14, 20, -7, -19};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {15, 22, -6, -21};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {16, 24, -8, -23};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {1, -19, -25, 17};
//+
Plane Surface(11) = {11};
//+
Curve Loop(12) = {2, -21, -26, 19};
//+
Plane Surface(12) = {12};
//+
Curve Loop(13) = {3, -23, -27, 21};
//+
Plane Surface(13) = {13};
//+
Curve Loop(14) = {4, -17, -28, 23};
//+
Plane Surface(14) = {14};
//+
Curve Loop(15) = {9, -20, -29, 18};
//+
Plane Surface(15) = {15};
//+
Curve Loop(16) = {10, -22, -30, 20};
//+
Plane Surface(16) = {16};
//+
Curve Loop(17) = {11, -24, -31, 22};
//+
Plane Surface(17) = {17};
//+
Curve Loop(18) = {12, -18, -32, 24};
//+
Plane Surface(18) = {18};
//+
Curve Loop(19) = {5, 9, -7, -1};
//+
Surface(19) = {19};
//+
Curve Loop(20) = {7, 10, -6, -2};
//+
Surface(20) = {20};
//+
Curve Loop(21) = {6, 11, -8, -3};
//+
Surface(21) = {21};
//+
Curve Loop(22) = {8, 12, -5, -4};
//+
Surface(22) = {22};
//+
Surface Loop(2) = {11, 19, 15, 7, 3, 8};
//+
Volume(2) = {2};
//+
Surface Loop(3) = {20, 16, 12, 4, 8, 9};
//+
Volume(3) = {3};
//+
Surface Loop(4) = {21, 17, 13, 5, 10, 9};
//+
Volume(4) = {4};
//+
Surface Loop(5) = {22, 18, 14, 6, 10, 7};
//+
Volume(5) = {5};




//+
Transfinite Curve {1, 2, 3, 4, 9, 10, 11, 12} = density+1 Using Progression 1;
//+
Transfinite Curve {17, 19, 21, 23, 20, 22, 24, 18} = density+1 Using Progression 1;
//+
Transfinite Curve {5, 7, 6, 8} = nEz+1 Using Progression 1;
//+

Transfinite Surface {7, 8, 9, 10, 11, 12, 13, 14, 15};
//+
Transfinite Surface {16, 17, 18, 19, 20, 21, 22};

//+
Transfinite Volume{2, 3, 4, 5};

//+
Recombine Surface {7, 8, 9, 10, 11, 12, 13, 14, 15};
//+
Recombine Surface {16, 17, 18, 19, 20, 21, 22};

Recombine Volume {1, 2, 3, 4, 5};

Periodic Surface {2, 15, 16, 17, 18} = {1, 11, 12, 13, 14} Translate {0, 0, -length};
//+


//+
Physical Surface("wall", 33) = {19, 20, 21, 22};
//+
Physical Surface("inlet", 34) = {1, 11, 12, 13, 14};
//+
Physical Surface("outlet", 35) = {2, 15, 16, 17, 18};
//+
Physical Volume("fluid", 36) = {1, 2, 3, 4, 5};
