Q := Rationals();
P4<x0,x1,x2,x3,x4> := ProjectiveSpace(Q,4);
fstar := x0*(x1^2+x2^2+x3^2-x4^2)-(x1^3+x2^3+x3^3-1/2*x4^3);
f := fstar - 1/4*x0^3;
D := Scheme(P4, f);
time IsNonsingular(D); // X=Z(f) is smooth

gstar := (x0^2+x1^2+x2^2+x3^2+x4^2)*f;
Xstar := Scheme(P4, gstar);
time IsNonsingular(Xstar); // X=Z(gstar) is singular

p := 9/1000*(x1^5+x2^5+x3^5);
fPerturbed := gstar + p;
fPerturbed;
X := Scheme(P4, fPerturbed);
time IsNonsingular(X); //