function [interv,theta] = stpGBLM1(s_pre,s_post,varargin)
if (~isempty(varargin))
    bb = varargin{1}; %% determines the shape of the exponential in w(t)
else
    bb = 7;
end 
L = length(s_pre);
if L~=length(s_post);warning('s_pre and s_post should have the same length');end
dt   = 0.001; % binsize (don't change it) 
T    = L*dt;   % seconds

% basis functions
delay=[100 100]/(dt*1000);
nf = [5 5];
Fc = orth(getBasis('rcos',nf(1),delay(1),20,0)')';
Fh = orth(getBasis('rcos',nf(2),delay(2),20,0)')';

Xc = getX(s_pre,Fc,0,1,0)'; 
Xh = getX(s_post,Fh,0,1,0)';

% finding the nonweighted stp plasticity trace of each interval
NN = 35; % number of intervals suggested values 8 - 60
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
e = getX(S1,mprops1.basis(bb,:)/max(mprops1.basis(bb,:)),0,1,0)'; % 

% spline expansion
Nq=8; % 8-20 suggested but the more baseline = overfitting
q_base = getBasis('rcos',Nq,2000,20,0);
c_gam = q_base(:,round(interv*1000));

nu = c_gam*e';
wt = ones(T/dt,1)'; % check if the T/dt is an integer
dev=0;
tic
for i = 2:40 % could change the maximum number of iteration 
    Xc_ = repmat(wt',1,nf(1)).*Xc;
    bta = glmfit([Xh Xc_],s_post,'poisson');
    a_offset = repmat(ones(T/dt,1),1,nf(1)).*Xc*bta(7:end)+bta(1)+Xh*bta(2:6);
    W_ = repmat(Xc*bta(7:end),1,Nq).*nu';
    [alph,dev(i),stats] = glmfit( W_ ,s_post,'poisson','offset',a_offset,'constant','off');
    wt = 1+alph' * nu;
    fprintf('Dev difference: %04.01f in %02.01f \n',abs(dev(i)-dev(i-1)),toc);
    if abs(dev(i)-dev(i-1))<.1;break;end % could change the convergence limit
end
%%
% [yfit,dlo,dhi] = glmval(alph,X_temp,'log',stats,.95,[],a_offset,'off');
%%
theta.bw = 1+alph'*q_base(:,100:1001);
theta.stats = stats.covb;
% theta.Xc_ = repmat(wt',1,nf(1)).*Xc;
theta.fh = bta(2:6)'*Fh;
theta.b0 = bta(1);
theta.offset=a_offset;
theta.X=W_;
theta.b=alph;
