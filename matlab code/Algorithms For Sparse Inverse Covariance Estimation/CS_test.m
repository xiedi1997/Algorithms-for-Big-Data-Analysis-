% function Test_BP

% max log(det(X)) - Tr(SX) - rho||X||_1     s.t.    X>=0

% generate data
clear;clc
p = 30;   % dimensions of the generated data
n = 100;   % number of sampling
S = model_1(p,n);   % generate covariance matrix S using model 1

rho = [10; 0.1; 0.01];  
beta = [0.1; 10; 1000];
gamma = 1;

eps = 10^-9;

objfun = @(X, rho)(-log(det(X)) + trace(S * X) + rho * trace(ones(p, p) * abs(X)));

time1 = zeros(3, 10);
time2 = zeros(3, 10);
obj_CVX = zeros(3, 10);
obj_ADMM = zeros(3, 10);
for i = 1:3
    for j = 1:10
        U = orth(randn(p, p));  % random orthogonal matrix
        X0 = U * diag(rand(p, 1)) * U';  % random positive definite matrix
        Y0 = randn(p, p);     % random matrix
        U = orth(randn(p, p));    % random orthogonal matrix
        Z0 = U * diag(rand(p, 1)) * U';    % random positive definite matrix
        
        % cvx calling mosek
        [X1,out1] = CS_cvx(S ,rho(i));
        time1(i,j) = out1.cputime;
        obj_CVX(i,j) = objfun(X1, rho(i));

        % AMDD
        start = cputime;
        X2 = ADMM_X(S, X0, Y0, Z0, rho(i), beta(i), gamma, eps);
        time2(i,j) = cputime - start;
        obj_ADMM(i,j) = objfun(X2, rho(i));
    end
end
for i = 1:3
    for j = 1:10
        fprintf('cvx:        cpu: %5.2f, obj: %5.2e\t ADMM:       cpu: %5.2f, obj: %5.2e\n', time1(i,j), obj_CVX(i,j), time2(i,j), real(obj_ADMM(i,j)));
    end
     fprintf('\n\n\n');
end
