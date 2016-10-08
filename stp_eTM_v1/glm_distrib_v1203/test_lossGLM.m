
addpath(genpath('C:\Users\ian\Documents\MATLAB\general\DERIVESTsuite'))

global X
global y

b = randn(10,1);
X = randn(100,10);

% loss_fun = @lossGLM_poiss_exp;
% y = poissrnd(exp(X*b));

loss_fun = @lossGLM_poiss_log;
y = poissrnd(log(1+exp(X*b)));

% loss_fun = @lossGLM_binom_logistic;
% y = rand(100,1)>(1./(1+exp(-X*b)));

% loss_fun = @lossGLM_gauss_iden;
% y = X*b+randn(size(X,1),1);


%%%% Double check gradient and hessian

% Analytic
[llhd,dx,H] = feval(loss_fun,b);

% Numerical Estimates
[grad,e]=gradest(loss_fun,b);
[hess,e]=hessian(loss_fun,b);

[grad' dx]
grad_rmse = sqrt(mean((grad'-dx).^2))

[hess(:) H(:)]
hess_rmse = sqrt(mean((hess(:)-H(:)).^2))