% proximal gradient method for min |z|_2 + f(x,z,y)
% x0 ,y0, z0 are the intital values of algorithm
% A, b ,mu are the given values
% t is parameter of augmented lagrangian function
% tau is parameter of proximal gradient method
% eps is the precision of  proximal gradient method

function z1 = prox_method_2(A, b, x0, y0, z0, t, tau, eps)
g2 = gradient_2(A, b, x0, y0, z0, t);
z1 = proximal_2(z0 - tau * g2, tau);
while norm(z1 - z0, 2) > eps
    z0 = z1;
    g2 = gradient_2(A, b, x0, y0, z0, t);
    z1 = proximal_2(z0 - tau * g2, tau);
end