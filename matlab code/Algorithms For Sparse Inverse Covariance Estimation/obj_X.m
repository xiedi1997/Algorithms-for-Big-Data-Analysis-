% objective function
% -log(det(X)) + Tr(SX) + rho||X||_1
% S, rho ,mu are the given values
% X is a iterative value

function f = obj_X(S, X, rho)
n = length(S);
A = ones(n, n);
f = -log(det(X)) + trace(S * X) + rho * trace(A * abs(X));