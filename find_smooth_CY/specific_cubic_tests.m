Q:=RationalField();
P5<x,y,z,u,v> := ProjectiveSpace(Q,4);
cubic:=x*(y^2+z^2+u^2+v^2)-(y^3+z^3+u^3-(1/2)*v^3);
V:=Scheme(P5,cubic);
cubic;
print "The variety is smooth: ";
IsNonsingular(V);

smoothing:=cubic+(1/4)*x^3;
W:=Scheme(P5,smoothing);
smoothing;
print "The variety is smooth: ";
IsNonsingular(W);