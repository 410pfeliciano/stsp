function [interv,theta] = stpGBLM_v2(s_pre,s_post,ifreg,ifcv,iflambd,n_exp)
if nargin<3; ifreg=0;end % choose whether regularize parameters or not
if nargin<4; ifcv=0;end
if nargin<5; iflambd=0;else Lambd=[5e-06];end
if nargin<6; n_exp=7;end % timeconstant selection of exp in w(t)

L = length(s_pre);
if L~=length(s_post);warning('s_pre and s_post should have the same length');end
dt   = 0.001; % binsize (don't change it) 
T    = L*dt;   % seconds

% basis functions
delay=[150 200]/(dt*1000);
nf = [5 5];
Fc = orth(getBasis('rcos',nf(1),delay(1),20,0)')';
Fh = orth(getBasis('rcos',nf(2),delay(2),20,0)')';

Xc = getX(s_pre,Fc,0,1,0)'; 
Xh = getX(s_post,Fh,0,1,0)';

% finding the nonweighted stp plasticity trace of each interval
NN = 50; % number of intervals suggested values 8 - 60
interv= logspace(log10(1/100 ),log10(1),NN); % [1 400]hz
Tlist = find(s_pre==1)*dt;
for i = 1: NN
    if i == 1
        index = find(diff(Tlist)<interv(i));
    else
        index = find(diff(Tlist)<interv(i) & diff(Tlist)>interv(i-1));
    end
    Tplas{i} = Tlist(index+1);
end

S1= double(getSpkMat(Tplas,dt,T,1));

mprops1.nfilt = 40;
mprops1.delay = 1000/(dt*1000);
mprops1.basis = getBasis('exp',mprops1.nfilt,mprops1.delay,20,0);
e = getX(S1,mprops1.basis(n_exp,:)/max(mprops1.basis(n_exp,:)),0,1,0)'; % 

% spline expansion
Nq=8; % 8-20 suggested but the more baseline = overfitting
q_base = getBasis('rcos',Nq,2000,20,0);
% q_base = [ones(1,2000);linspace(0,1,2000)];
c_gam = q_base(:,round(interv*1000));

nu = c_gam*e';
wt = ones(T/dt,1)'; % check if the T/dt is an integer
dev=0;
tic
for i = 2:20 % could change the maximum number of iteration 
    Xc_ = repmat(wt',1,nf(1)).*Xc;
    bta = glmfit([Xh Xc Xc_],s_post,'poisson');
    a_offset = bta(1) + Xh*bta(2:6) + Xc*bta(7:11) ;%repmat(ones(T/dt,1),1,nf(1)).*
    W_ = repmat(Xc*bta(12:end),1,Nq).*nu';
    if ifreg==1 % regularized gblm
        if ifcv==1 % cross-validated
            if iflambd==1
                [alph,stats] = lassoglm( W_ ,s_post,'poisson',...
                    'Alpha',1,'Lambda',Lambd,'Offset',a_offset,'CV',4);
            else
                [alph,stats] = lassoglm( W_ ,s_post,'poisson',...
                    'Alpha',1,'NumLambda',8,'Offset',a_offset,'CV',4);
            end
        else
            if iflambd==1
                [alph,stats] = lassoglm( W_ ,s_post,'poisson',...
                    'Alpha',1,'Lambda',Lambd,'Offset',a_offset);
            else
                [alph,stats] = lassoglm( W_ ,s_post,'poisson',...
                    'Alpha',1,'NumLambda',8,'Offset',a_offset);
            end
        end
        [dev(i),b] = min(stats.Deviance);
%         wt = 1+alph(:,b)' * nu;
        wt = alph(:,b)' * nu;
    else
        [alph,dev(i),stats] = glmfit( W_ ,s_post,'poisson',...
            'offset',a_offset,'constant','off');
        wt = alph' * nu;
    end
    fprintf('Dev difference: %04.01f in %02.01f \n',abs(dev(i)-dev(i-1)),toc);
    if abs(dev(i)-dev(i-1))<.1;break;end % could change the convergence limit
end
% keyboard
%%
% [yfit,dlo,dhi] = glmval(alph,X_temp,'log',stats,.95,[],a_offset,'off');
%%
if ifreg==1
    theta.bw = 1+alph(:,b)'*q_base(:,100:1001);
    theta.nb = b;
    theta.fh = bta(2:6)'*Fh;
    theta.b0 = bta(1);
    theta.offset=a_offset;
    theta.X=[ones(size(W_,1),1) W_];
    theta.b=alph;

else
    theta.bw = 1+alph'*q_base(:,100:1001);
    theta.stats = stats.covb;
    % theta.Xc_ = repmat(wt',1,nf(1)).*Xc;
    theta.fh = bta(2:6)'*Fh;
    theta.b0 = bta(1);
    theta.offset=a_offset;
    theta.X=W_;
    theta.b=alph;
end