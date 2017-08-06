function theta = stp_gblm(s_pre,s_post,varargin)% 1- loading the spiking activity

dt   = 0.001; % binsize (don't change it) 
L = length(s_pre);
isi = diff(find(s_pre>0)*dt); % interspike intervals

% defaults
delay=[100 100]/(dt*1000);
nfilt = [5 5];
NumLambda = 20;
numCV = 8;
toleranceValue=.1;

if (~isempty(varargin))
    c = 1 ;
    % user defined
    while c <= length(varargin)
        switch varargin{c}
            case {'delay'}
                delay = varargin{c+1};
            case {'nfilt'}
                nfilt = varargin{c+1};
            case {'NumLambda'}
                NumLambda = varargin{c+1};
            case {'numCV'}
                numCV = varargin{c+1};
            case {'toleranceValue'}
                toleranceValue = varargin{c+1};
        end % switch
        c = c + 2;
    end % for
end % if

% 2- create basis functions for coupling and post-spike history filter
Fc = orth(getBasis('rcos',nfilt(1),delay(1),20,0)')';
Fh = orth(getBasis('rcos',nfilt(2),delay(2),20,0)')';
Xc = getX(s_pre,Fc,0,1,0)'; 
Xh = getX(s_post,Fh,0,1,0)';

% 3- spline expansion
Nq=8; % 8-20 suggested but the more baseline = overfitting
Bm = getBasis('rcos',Nq,500,20,0);
Bm = padarray(Bm,[0 5000]);
Bm=Bm(:,5001:end);
Bm_dt = Bm(:,round(isi*1000));
s=zeros(Nq,L);
for i =1: Nq
    s(i,s_pre>0)=[0 Bm_dt(i,:)];
end

% 4- transient effect of stp carried by exp function
% from splined spikes to exponential ( from "s" [Nq x T] --> "e" [Nq x T] )
tau= .2; % timeconstant for decay
x0 = linspace(0,1,1000);
kern_stp = exp(-x0/tau);
e = filter(kern_stp,1,s');

% 5- gblm
wt = ones(L,1); % check if the T/dt is an integer
dev=0;
tic
for i = 2:40 % could change the maximum number of iteration 
    Xc_ = repmat(wt,1,nfilt(1)).*Xc; % modified coupling term
%     bta = glmfit([Xh Xc_],s_post,'poisson'); % glm1 -> estimtes coupling and history filter based on modified coupling term
    if i==2
        % initial glm
        bta = glmfit([Xh Xc_],s_post,'poisson');
    else
        % warm-start (using previous bta)
        bta = glmfit([Xh Xc_],s_post,'poisson','B0',bta);
    end
    a_offset = Xc*bta(7:end)+bta(1)+Xh*bta(2:6); % static part of the lambda
    W_ = repmat(Xc*bta(7:end),1,Nq).*e; % plastic part of lambda
    % glm2 -> estimated the modification term (affected by the plasticity)
%     [alph,stats] = lassoglm( W_ ,s_post,'poisson',...
%                     'Alpha',1,'NumLambda',NumLambda,'Offset',a_offset,'CV',numCV,'RelTol',1e-3);
    % for faster response ignore NumLambda and provide a single one
    [alph,stats] = lassoglm( W_ ,s_post,'poisson',...
                    'Alpha',1,'Lambda',1e-5,'Offset',a_offset,'CV',numCV,'RelTol',1e-3);
    [dev(i),b] = min(stats.Deviance);
    wt = 1 + e*alph(:,b);
    fprintf('Dev difference: %04.01f in %02.01f \n',(dev(i)-dev(i-1)),toc);
    if abs(dev(i)-dev(i-1))<toleranceValue && i>10;break;end % could change the convergence limit
end
for b = 1:size(alph,2)
    theta.wt(:,b) = 1+e*alph(:,b);
end
theta.alph = alph;
theta.bta = bta;
theta.modif_fxn = alph'*Bm;
theta.stats = stats;
theta.couplingFilt = Fc;
theta.histFilt = Fh;