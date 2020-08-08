% The gradient of f with respect to x
% A1 = A * A'
% A, b are the given values
% t is parameter of augmented lagrangian function
% x ,y, z are the iterative value

function g = gradient_1(A, A1, b, x, y, z, t)
% g = (A' * (A * x - b - z)) / t + A' * y;
g = (A1 * x) / t  - A' * ((b + z) / t - y); 