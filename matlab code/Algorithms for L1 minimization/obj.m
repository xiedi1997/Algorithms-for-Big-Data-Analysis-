% objective function
% mu * ||x||_1 + ||Ax-b||_2
% A, b ,mu are the given values
% x is a iterative value

function f = obj(A, x, b, mu)
f = mu * norm(x, 1) + norm(A * x- b, 2);