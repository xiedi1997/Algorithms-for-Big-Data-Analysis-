% solve the following problem for variable Z
% Z = argmin  L(X,Z,Y), L is augmented Lagrangian function 


function Z = operator_Z(X, Y, rho, beta)
Z = shrink_matrix(X + beta * Y, beta * rho);