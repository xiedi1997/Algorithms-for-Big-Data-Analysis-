% use the gurobi to solve following BP problem
% min mu*||x||_1 + ||Ax-b||_2
% x0 is a given input initial solution(It is not used in this function)
% opts is a struct which stores the options of the algorithm(It is not used in this function)
% x is the optimal solution
% out is a struct which saves all other output information.

function [x, out] = BP_gurobi(x0, A, b, mu, opts)

m = size(A, 1);
n = size(A, 2);

% Set objective : (c^T) * x
model . obj = [ zeros(1, n), zeros(1, n + m), ones(1, n + 1) ];  % c
model . lb = [-inf * ones(1, n + n + m), zeros(1, n + 1)];  % Lower bound of x
model . ub = inf*ones(1, n + n + m + n + 1); % upper bound of x
model . modelsense = 'min';  % minimize the objective function

% Add constraint : Ax-b=0
model .A = sparse([mu * eye(n), -eye(n), zeros(n, m), zeros(n, n + 1); A, zeros(m, n), -eye(m), zeros(m, n + 1)]);
model .rhs = [zeros(1, n),b'];
model . sense = '=';

% Add second - order cone : (x^T)Qx + (q^T)x + c >= 0
% Add second - order cone : \bar{v}_i^2 <= v_{i0}^2     i = 1,...,n
for i = 1:n
    Q = zeros(n + n + m + n + 1, n + n + m + n + 1);
    Q(2 * n + m + i, 2 * n + m + i) = -1;
    Q(n + i, n + i) = 1;
    model . quadcon (i). Qc = sparse(Q);
    model . quadcon (i). q = zeros(n + n + m + n + 1, 1);
    model . quadcon (i). rhs = 0;
    i
end
% Add second - order cone : \sum_i^n{\bar{v}_i^2} <= v_{n+1 0}^2
Q = zeros(n + n + m + n + 1, n + n + m + n + 1);
Q(2 * n + m + n + 1, 2 * n + m + n + 1) = -1;
Q(2 * n + 1:2 * n + m, 2 * n + 1:2 * n + m) = eye(m);

model . quadcon (n + 1). Qc = sparse(Q);
model . quadcon (n + 1). q = zeros(n + n + m + n + 1, 1);  
model . quadcon (n + 1). rhs = 0;

% solve the problem
gurobi_write (model , 'qcp .lp');



out = gurobi ( model );
x = out.x(1:n);
