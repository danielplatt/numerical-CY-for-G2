Q:=RationalField();
P5<x,y,z,u,v> := ProjectiveSpace(Q,4);
epsilon:=1/10;
cubic:=Random(0,1)^3*u^3 + Random(0,1)^3*u^2*v + Random(0,1)^3*u*v^2 + Random(0,1)^3*v^3 + 
 Random(0,1)^3*u^2*x + Random(0,1)^3*u*v*x + Random(0,1)^3*v^2*x + Random(0,1)^3*u*x^2 + 
 Random(0,1)^3*v*x^2 + Random(0,1)^3*x^3 + Random(0,1)^3*u^2*y + Random(0,1)^3*u*v*y + 
 Random(0,1)^3*v^2*y + Random(0,1)^3*u*x*y + Random(0,1)^3*v*x*y + Random(0,1)^3*x^2*y + 
 Random(0,1)^3*u*y^2 + Random(0,1)^3*v*y^2 + Random(0,1)^3*x*y^2 + Random(0,1)^3*y^3 + 
 Random(0,1)^3*u^2*z + Random(0,1)^3*u*v*z + Random(0,1)^3*v^2*z + Random(0,1)^3*u*x*z + 
 Random(0,1)^3*v*x*z + Random(0,1)^3*x^2*z + Random(0,1)^3*u*y*z + Random(0,1)^3*v*y*z + 
 Random(0,1)^3*x*y*z + Random(0,1)^3*y^2*z + Random(0,1)^3*u*z^2 + Random(0,1)^3*v*z^2 + 
 Random(0,1)^3*x*z^2 + Random(0,1)^3*y*z^2 + Random(0,1)^3*z^3;
V:=Scheme(P5,cubic);
cubic;
IsNonsingular(V);