% generate covariance matrix S using model 1
% p is the dimensions of the generated data
% n is the number of sampling

function S = model_2(p, n)
B = binornd(1, 0.1, p, p);
B = tril(B, -1) + tril(B, -1)';
B = B / 2;
lambda_B = eig(B);
lambda1 = max(lambda_B);
lambda2 = min(lambda_B);
delta = (lambda1 - lambda2 * p)/(p - 1);
Omega = B + delta * eye(p);
Omega_0 = Omega / delta;
sigma_0 = inv(Omega_0);
z = mvrnnd(zeros(1, p), sigma_0, n);
S = cov(z);