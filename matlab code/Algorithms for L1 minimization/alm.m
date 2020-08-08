% use the Augmented Lagrangian method to solve following BP problem
% min mu*||x||_1 + ||Ax-b||_2
% A1 = A' * A
% x0 ,y0, z0 are the intital values of algorithm
% A, b ,mu are the given values
% t is parameter of augmented lagrangian function
% eps2 is the precision of augmented lagrangian method
% tau is parameter of proximal gradient method
% eps1 is the precision of  proximal gradient method
% x1 is the optimal solution

function x1 = alm(A, A1, b, x0, y0, z0, t, mu, tau, eps1, eps2)
[x1, z1] = prox_method_12(A, A1, b, x0, y0, z0, t, mu, tau, eps1);
y1 = y0 + (A * x1- b - z1) / t;
while abs(obj(A, x1, b, mu) - obj(A, x0, b, mu)) > eps2
    y0 = y1;
    x0 = x1;
    z0 = z1;
    [x1, z1] = prox_method_12(A, A1, b, x0, y0, z0, t, mu, tau, eps1);
    y1 = y0 + (A * x1- b - z1) / t;
end