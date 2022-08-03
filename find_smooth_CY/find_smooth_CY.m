Q:=RationalField();
P5<x,y,z,u,v,w> := ProjectiveSpace(Q,5);
circ:=x^2+y^2-w^2;
epsilon:=1/2;
circperturb:=epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u^2 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u*v + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*v^2 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u*w +
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*v*w + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*w^2 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u*x + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*v*x +
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*w*x + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*x^2 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u*y + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*v*y +
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*w*y + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*x*y + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*y^2 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u*z +
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*v*z + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*w*z + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*x*z + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*y*z +
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*z^2;
quart:=z^4+u^4+v^4-w^4;
quartperturb:=epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u^4 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u^3*v + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u^2*v^2 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u*v^3 + 
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*v^4 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u^3*w + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u^2*v*w + 
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u*v^2*w + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*v^3*w + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u^2*w^2 + 
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u*v*w^2 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*v^2*w^2 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u*w^3 + 
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*v*w^3 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*w^4 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u^3*x + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u^2*v*x + 
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u*v^2*x + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*v^3*x + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u^2*w*x + 
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u*v*w*x + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*v^2*w*x + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u*w^2*x + 
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*v*w^2*x + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*w^3*x + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u^2*x^2 + 
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u*v*x^2 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*v^2*x^2 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u*w*x^2 + 
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*v*w*x^2 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*w^2*x^2 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u*x^3 + 
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*v*x^3 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*w*x^3 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*x^4 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u^3*y + 
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u^2*v*y + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u*v^2*y + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*v^3*y + 
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u^2*w*y + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u*v*w*y + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*v^2*w*y + 
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u*w^2*y + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*v*w^2*y + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*w^3*y + 
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u^2*x*y + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u*v*x*y + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*v^2*x*y + 
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u*w*x*y + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*v*w*x*y + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*w^2*x*y + 
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u*x^2*y + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*v*x^2*y + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*w*x^2*y + 
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*x^3*y + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u^2*y^2 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u*v*y^2 + 
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*v^2*y^2 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u*w*y^2 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*v*w*y^2 + 
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*w^2*y^2 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u*x*y^2 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*v*x*y^2 + 
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*w*x*y^2 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*x^2*y^2 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u*y^3 + 
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*v*y^3 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*w*y^3 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*x*y^3 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*y^4 + 
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u^3*z + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u^2*v*z + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u*v^2*z + 
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*v^3*z + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u^2*w*z + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u*v*w*z + 
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*v^2*w*z + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u*w^2*z + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*v*w^2*z + 
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*w^3*z + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u^2*x*z + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u*v*x*z + 
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*v^2*x*z + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u*w*x*z + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*v*w*x*z + 
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*w^2*x*z + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u*x^2*z + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*v*x^2*z + 
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*w*x^2*z + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*x^3*z + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u^2*y*z + 
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u*v*y*z + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*v^2*y*z + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u*w*y*z + 
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*v*w*y*z + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*w^2*y*z + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u*x*y*z + 
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*v*x*y*z + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*w*x*y*z + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*x^2*y*z + 
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u*y^2*z + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*v*y^2*z + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*w*y^2*z + 
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*x*y^2*z + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*y^3*z + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u^2*z^2 + 
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u*v*z^2 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*v^2*z^2 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u*w*z^2 + 
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*v*w*z^2 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*w^2*z^2 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u*x*z^2 + 
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*v*x*z^2 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*w*x*z^2 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*x^2*z^2 + 
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u*y*z^2 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*v*y*z^2 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*w*y*z^2 + 
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*x*y*z^2 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*y^2*z^2 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*u*z^3 + 
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*v*z^3 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*w*z^3 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*x*z^3 + epsilon*Random(0,1)*Random(0,1)*Random(0,1)*y*z^3 + 
 epsilon*Random(0,1)*Random(0,1)*Random(0,1)*z^4;
V:=Scheme(P5,[circ+circperturb,quart+quartperturb]);
[circ+circperturb,quart+quartperturb];
IsNonsingular(V);