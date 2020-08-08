% generate covariance matrix S using model 1
% p is the dimensions of the generated data
% n is the number of sampling

function S = model_1(p,n)
Omega_0 = zeros(p);
for i = 1:p
    for j = 1:p
        Omega_0(i, j) = 0.6^(abs(i-j));
    end
end
sigma_0 = inv(Omega_0);
z = mvnrnd(zeros(1, p), sigma_0, n);
S = cov(z);