% Loss function for GLM: Binomial noise, logisitc nonlinearity
% opts=optimset('Gradobj','on','Hessian','on');
% fminunc(@lossGLM_binom_logistic,b0,opts,X,y,nu)
%
% Note that likelihood is per bin, prior is per parameter
% (This makes it easier to use the same hyperparam for different datasets)

function [llhd, dx, H] = lossGLM_binom_logistic(b,X,y,nu,offset,penalty,penalty_type)

% % For testing...
% global X
% global y

%%%% Defaults
if nargin<4, nu=0; end                                  % No Regularization
if nargin<5, offset=0; end                              % No offset
if nargin<6 || isempty(penalty), penalty=b*0+1; end     % Equal penalty on everything
if nargin<7, penalty_type='none'; end                   % MLE

% % For testing
% nu=1;
% penalty_type = 'l1';

[n,p] = size(X);
lam = 1./(1+exp(-(X*b+offset)));
ymax = max(y);

% MLE
llhd = -sum(y.*log(lam+(lam==0)) + (ymax-y).*log(1-lam))/n;
dx = - (X'*(y.*(1-lam) - (ymax-y).*lam))/n;
Htmp = -y.*lam.*(1-lam) - (ymax-y).*lam.*(1-lam);
H = -X'*bsxfun(@times,Htmp,X)/n;

if strcmp(penalty_type,'l2')        % MAP - L2
    llhd  = llhd + nu*sum((b.*penalty).^2)/p;
    dx = dx + nu*2*(b.*penalty)/p;
    H = H + diag(2*nu*penalty)/p;
elseif strcmp(penalty_type,'l1')    % MAP - L1
    llhd  = llhd + nu*sum(abs(b.*penalty))/p;
    dx = dx + nu*sign(b.*penalty)/p;
end