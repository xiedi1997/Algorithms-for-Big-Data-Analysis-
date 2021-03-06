% use the gurobi solver of CVX to solve following BP problem
% min mu*||x||_1 + ||Ax-b||_2
% x0 is a given input initial solution(It is not used in this function)
% opts is a struct which stores the options of the algorithm(It is not used in this function)
% x is the optimal solution
% out is a struct which saves all other output information.

function [x, out] = BP_cvx_gurobi(x0, A, b, mu, opts)
n = size(A, 2);
cvx_begin
    cvx_solver gurobi
    cvx_save_prefs
    variable x(n)
    minimize(norm(A*x-b) + mu*norm(x,1))
cvx_end
out.optval = cvx_optval;
out.cputime = cvx_cputime;