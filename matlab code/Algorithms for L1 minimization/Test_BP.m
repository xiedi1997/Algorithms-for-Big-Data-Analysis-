% function Test_BP

% min mu*||x||_1 + ||Ax-b||_2

% generate data
clear;clc
n = 1024;
m = 512;

A = randn(m,n);
u = sprandn(n,1,0.1);
b = A*u;
mu = 1e-2;
A1 = A'*A;

x0 = randn(n, 1);
y0 = randn(m, 1);
z0 = randn(m, 1);

t = 10000;
tau = 5;
eps1 = 0.1;
eps2 = 10^-9;

errfun = @(x1, x2) norm(x1-x2)/(1+norm(x1));
resfun = @(x) norm(A*x-b);
nrm1fun = @(x) norm(x,1);
objfun = @(x)(mu * norm(x, 1) + norm(A * x- b, 2));

% cvx directly
opts1 = []; 
[x1, out1] = BP_cvx(x0, A, b, mu, opts1);
time1 = out1.cputime;

% cvx calling mosek
opts2 = []; 
[x2, out2] = BP_cvx_mosek(x0, A, b, mu, opts2);
time2 = out2.cputime;

% cvx calling gurobi
opts3 = []; 
[x3, out3] = BP_cvx_gurobi(x0, A, b, mu, opts3);
time3 = out3.cputime;

% call mosek directly
opts4 = []; 
start = cputime;
[x4, out4] = BP_mosek(x0, A, b, mu, opts4);
time4 = cputime - start;

% call gurobi directly
opts5 = []; 
[x5, out5] = BP_gurobi(x0, A, b, mu, opts5);
time5 = out5.runtime;

% Augmented Lagrangian method
start = cputime;
x6 = alm(A, A1, b, x0, y0, z0, t, mu, tau, eps1, eps2);
time6 = cputime - start;

% AMDD
start = cputime;
x7 = ADMM(A, A1, b, x0, y0, z0, t, mu, tau, eps1, eps2);
time7 = cputime - start;

% print comparison results with u
fprintf('cvx:        nrm1: %3.2e, res: %3.2e, cpu: %5.2f, err-to-u: %3.2e, obj: %5.2e\n', ... 
        resfun(x1), nrm1fun(x1), time1, errfun(u, x1), objfun(x1));
fprintf('cvx_mosek:  nrm1: %3.2e, res: %3.2e, cpu: %5.2f, err-to-u: %3.2e, obj: %5.2e\n', ...
        resfun(x2), nrm1fun(x2), time2, errfun(u, x2), objfun(x2));
fprintf('cvx_gurobi: nrm1: %3.2e, res: %3.2e, cpu: %5.2f, err-to-u: %3.2e, obj: %5.2e\n', ...
        resfun(x3), nrm1fun(x3), time3, errfun(u, x3), objfun(x3));
fprintf('mosek:      nrm1: %3.2e, res: %3.2e, cpu: %5.2f, err-to-u: %3.2e, obj: %5.2e\n', ...
        resfun(x4), nrm1fun(x4), time4, errfun(u, x4), objfun(x4));
fprintf('gurobi:     nrm1: %3.2e, res: %3.2e, cpu: %5.2f, err-to-u: %3.2e, obj: %5.2e\n', ...
        resfun(x5), nrm1fun(x5), time5, errfun(u, x5), objfun(x5));
fprintf('AML:        nrm1: %3.2e, res: %3.2e, cpu: %5.2f, err-to-u: %3.2e, obj: %5.2e\n', ...
        resfun(x6), nrm1fun(x6), time6, errfun(u, x6), objfun(x6));
fprintf('ADMM:       nrm1: %3.2e, res: %3.2e, cpu: %5.2f, err-to-u: %3.2e, obj: %5.2e\n', ...
        resfun(x7), nrm1fun(x7), time7, errfun(u, x7), objfun(x7));