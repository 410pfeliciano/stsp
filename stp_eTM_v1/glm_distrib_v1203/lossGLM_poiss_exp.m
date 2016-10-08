% Loss function for GLM: Poisson noise, exp nonlinearity
% opts=optimset('Gradobj','on','Hessian','on');
% fminunc(@lossGLM_poiss_exp,b0,opts,X,y,nu)
%
% Note that likelihood is per bin, prior is per parameter
% (This makes it easier to use the same hyperparam for different datasets)

function [llhd, dx, H] = lossGLM_poiss_exp(b,X,y,nu,offset,penalty,penalty_type)

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
% penalty_type = 'l2';

[n,p] = size(X);
lam = exp(X*b+offset);
idx = isfinite(lam);

% MLE
llhd  = sum(lam(idx) - log(lam(idx)+(lam(idx)==0)).*y(idx))/n;
dx = X(idx,:)'*(lam(idx) - y(idx))/n;
H = X(idx,:)'*bsxfun(@times,lam(idx),X(idx,:))/n;

if strcmp(penalty_type,'l2')        % MAP - L2
    llhd  = llhd + nu*sum((b.*penalty).^2)/p;
    dx = dx + nu*2*(b.*penalty)/p;
    H = H + diag(2*nu*penalty)/p;
elseif strcmp(penalty_type,'l1')    % MAP - L1
    llhd  = llhd + nu*sum(abs(b.*penalty))/p;
    dx = dx + nu*sign(b.*penalty)/p;
end