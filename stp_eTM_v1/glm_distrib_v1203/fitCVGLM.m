% Fit a regularized GLM (input X, output y) using n-fold Cross Validation
% Returns a structure m, with params and log-likelihood information
% m has fields...
%   b       [p x nfoldcv x length(nu_vec)] matrix of parameters
%   breg    [p x nfoldcv]                  MAP estimates at the optimized hyperparam
%   llhdr       training-set negative log-likelihoods (per bin)
%   llhdt       test-set negative log-likelihoods (per bin)
%   llhdt0      test-set negative log-likelihoods (per bin) for homogeneous process
%   llrt        test-set log likelihood ratios relative to homogeneous process
%   numspks     test-set spikes to get llrt/spk
%   numbins     test-set bins to get llrt/s
%   exitflag    output from fminunc to determine convergence

function m = fitCVGLM(X,y,nfoldcv,nu_vec,loss_fun)

% Use poisson noise with log(1+Xb) output nonlinearity by default
if nargin<5,
    loss_fun = 'lossGLM_poiss_log';
end

opts=optimset('Gradobj','on','Hessian','on');

% Don't penalize the baseline parameter...
penalty = ones(size(X,2)+1,1);
penalty(1) = 0;

% Cross-Validation
k=floor(size(y,1)/nfoldcv);
for i=1:nfoldcv
    fprintf('Fold %03i/%03i...',i,nfoldcv)
    tic
    idx_tr = [1:(k*(i-1)) k*i:size(y,1)];
    idx_ts = [(k*(i-1)+1):(k*i-1)];
    
    % Fit a homogeneous model to compare others against...
    [bmu, f, exitflag] = fminunc(loss_fun,0,opts,[X(idx_tr,1)*0+1],y(idx_tr));
    m.llhdt0(i) = feval(loss_fun,bmu,ones(length(idx_ts),1),y(idx_ts));

    % Fit MLE to use as a starting point...
    b0 = zeros(size(X,2)+1,1); b0(1) = bmu;
    b0 = fminunc(loss_fun,b0,opts,[X(idx_tr,1)*0+1 X(idx_tr,:)],y(idx_tr),0);

    % Loop over regularization hyper-parameters...
    for j=1:length(nu_vec)
        [m.b(:,i,j), f, m.exitflag(i,j)]   = fminunc(loss_fun,b0,opts,[X(idx_tr,1)*0+1 X(idx_tr,:)],y(idx_tr),nu_vec(j),0*y(idx_tr),penalty,'l1');
        m.llhdr(i,j) = feval(loss_fun,m.b(:,i,j),[ones(length(idx_tr),1) X(idx_tr,:)],y(idx_tr));
        m.llhdt(i,j) = feval(loss_fun,m.b(:,i,j),[ones(length(idx_ts),1) X(idx_ts,:)],y(idx_ts));
        m.llrt(i,j)  = (-m.llhdt(i,j)+m.llhdt0(i))/log(2);
        toc
    end
    [m.llrtmax(i),id]=max(m.llrt(i,:));
    m.breg(:,i) = m.b(:,i,id);
    m.numspks(i) = sum(y(idx_ts));
    m.numbins(i) = length(idx_ts);
end