loadPackage "SheafStability";

R = QQ[x,y,z,w,MonomialOrder=>GRevLex]
I = ideal(x^3+y^3+w^2*z,w^3+z^3+w^2*x)
S = R/I
M = matrix{{ x^3, w^2*y, z^2*x }}

x = computeSemistability(M)