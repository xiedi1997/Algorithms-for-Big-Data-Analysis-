% solve the following problem for variable X
% X = argmin  L(X,Z,Y), L is augmented Lagrangian function 

function X = operator_X(S, Y, Z, beta)
n = length(S);
lambda = zeros(n, 1);
A = S + Y';
[U, D] = eig(Z);
mu = diag(D);
B = U' * A * U;
for i = 1:n
    lambda(i) = (mu(i) - B(i,i) * beta + ((B(i,i) * beta - mu(i))^2 + 4 * beta)^0.5) / 2;
end
X = U * diag(lambda) * U';