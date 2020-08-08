% augmented Lagrangian function
% mu * ||x||_1 + ||z||_2 + y'(Ax-b-z) + (||Ax-b-z||_2^2)/2t
% A, b ,mu are the given values
% t is parameter of augmented lagrangian function
% x ,y, z are the iterative value

function f = L(A, b, x, y, z, t, mu)
f = mu * norm(x, 1) + norm(z, 2) + (norm(A * x - b - z, 2)^2) / (2 * t) + y' * (A * x - b - z);