% The gradient of f with respect to z
% A, b are the given values
% t is parameter of augmented lagrangian function
% x ,y, z are the iterative value

function g = gradient_2(A, b, x, y, z, t)
g = (A * x - b - z) / (-t) - y;