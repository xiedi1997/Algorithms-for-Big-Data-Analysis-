% use the ADMM to solve following problem
% min -log(det(X)) + Tr(SX) + rho||X||_1     s.t.    X>=0
% X0 ,Y0, Z0 are the intital values of algorithm
% S,  rho are the given values
% beta,gamma are parameters of augmented lagrangian function
% eps is the precision of augmented lagrangian method

function X1 = ADMM_X(S, X0, Y0, Z0, rho, beta, gamma, eps)
X1 = operator_X(S, Y0, Z0, beta);
Z1 = operator_Z(X1, Y0, rho, beta);
Y1 = Y0 + (gamma / beta) * (X1 - Z1);
while  abs(obj_X(S, X1, rho) - obj_X(S, X0, rho)) > eps
    X0 = X1;
    Z0 = Z1;
    Y0 = Y1;
    X1 = operator_X(S, Y0, Z0, beta);
    Z1 = operator_Z(X1, Y0, rho, beta);
    Y1 = Y0 + (gamma / beta) * (X1 - Z1);
    obj_X(S, X1, rho)
end
