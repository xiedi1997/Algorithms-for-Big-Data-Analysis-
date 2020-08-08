% Differentiable part of augmented Lagrangian function
% L(x,z,y) = r(x,z,y) + f(x,z,y)
% A, b are the given values
% t is parameter of augmented lagrangian function
% x ,y, z are the iterative value

function f = f(A, b, x, y, z, t)
f = (norm(A * x - b - z, 2)^2) / (2 * t) + y' * (A * x - b - z);