% use the CVX to solve following CS problem
% max log(det(X)) - Tr(SX) - rho||X||_1     s.t.    X>=0
% S is a given  covariance matrix
% rho is a given parameter
% X is the optimal solution
% out is a struct which saves all other output information.

function [X,out] = CS_cvx(S, rho)
n = size(S,1);
A = ones(n, n);
cvx_begin
    cvx_solver mosek;
    variable X(n,n) symmetric;
    minimize( -log_det(X) + trace(S * X) + rho * trace(A * abs(X)));
    subject to
        X >= 0;
cvx_end
out.optval = cvx_optval;
out.cputime = cvx_cputime;