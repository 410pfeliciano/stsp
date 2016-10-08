% Loss function for GLM: Poisson noise, log(1+exp(Xb)) nonlinearity
% opts=optimset('Gradobj','on','Hessian','on');
% fminunc(@lossGLM_poiss_log,b0,opts,X,y,nu)
%
% Note that likelihood is per bin, prior is per parameter
% (This makes it easier to use the same hyperparam for different datasets)

function [llhd, dx, H] = lossGLM_poiss_log(b,X,y,nu,offset,penalty,penalty_type)

% % For testing...
% global X
% global y

%%%% Defaults
if nargin<4, nu=0; end                                  % No Regularization
if nargin<5, offset=0*y; end                            % No offset
if nargin<6 || isempty(penalty), penalty=b*0+1; end     % Equal penalty on everything
if nargin<7, penalty_type='none'; end                   % MLE

% % For testing
% nu=1;
% penalty_type = 'l2';

[n,p] = size(X);
eXb = exp(X*b+offset);
lam = log(1+eXb);

idx = ~isfinite(lam);
lam(idx) = X(idx,:)*b+offset(idx);
tmp = eXb./(1+eXb);
tmp(idx) = 1;

% MLE
llhd  = sum(lam - log(lam+(lam==0)).*y)/n;
dxtmp = tmp.*(1 - y./lam);
idx = isfinite(dxtmp)&(lam>0)&isfinite(eXb);
dx = X(idx,:)'*dxtmp(idx)/n;
Htmp = ((y./(lam.^2+(lam==0))).*eXb+(1-y./(lam+(lam==0)))).*tmp./(1+eXb);
H = X(idx,:)'*bsxfun(@times,Htmp(idx),X(idx,:))/n;

if strcmp(penalty_type,'l2')        % MAP - L2
    llhd  = llhd + nu*sum((b.*penalty).^2)/p;
    dx = dx + nu*2*(b.*penalty)/p;
    H = H + diag(2*nu*penalty)/p;
elseif strcmp(penalty_type,'l1')    % MAP - L1
    llhd  = llhd + nu*sum(abs(b.*penalty))/p;
    dx = dx + nu*sign(b.*penalty)/p;
end

% Stop for errors...
if any(~isfinite(dx)) || any(~isfinite(b)) || ~isfinite(llhd) || any(~isfinite(H(:)))
    keyboard
    llhd = Inf;
    dx = b*0;
end