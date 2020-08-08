% proximal gradient method for min mu*|x|_1 + f(x,z,y)
% x0 ,y0, z0 are the intital values of algorithm
% A, b ,mu are the given values
% t is parameter of augmented lagrangian function
% tau is parameter of proximal gradient method
% eps is the precision of  proximal gradient method

function x1 = prox_method_1(A, A1, b, x0, y0, z0, t, mu, tau, eps)
g1 = gradient_1(A, A1, b, x0, y0, z0, t);
x1 = shrink(x0 - tau * g1, mu * tau);
while abs(L(A, b, x1, y0, z0, t, mu) - L(A, b, x0, y0, z0, t, mu)) > eps
    x0 = x1;
    g1 = gradient_1(A, A1, b, x0, y0, z0, t);
    x1 = shrink(x0 - tau * g1, mu * tau);
end