function [param_bta, param, optimres]= stp_tmglm(population,T,varargin)
%-INPUTS-------------------------------------------------------------------
% population: a cell array [2x1]
%           - first cell: presynaptic neuron's time-stamps of the APs
%           - second cell: postsynaptic neuron's time-stamps of the APs
% T: duration of the recorded pre and postsynaptic activity [scalar]

%-OUTPUT-------------------------------------------------------------------
% param: estimated parameters [D,F,U,f]

% default
beta_h=[];
isConstrained = false;
delay =[150 150];
dt = .001;
hist.nfilt = 5;
coupling.nfilt = 5;
noise=[];
nrestart=1;

% change default
if (~isempty(varargin))
    c = 1 ;
    
    % user defined
    while c <= length(varargin)
        switch varargin{c}
            case {'beta'} % adds the fixed beta and doesn't compute the beta each time
                beta_h = varargin{c+1};
                
            case {'constrained'}
                isConstrained = true;
                upper = varargin{c+1};
                
            case {'delay'}
                delay = varargin{c+1};
                
            case {'nfilt'}
                hist.nfilt = varargin{c+1}(2);
                coupling.nfilt = varargin{c+1}(1);
                
            case ('dt')
                dt = varargin{c+1};
                
            case ('noise') % this could be LFP recording to help the estimation process
                noise = varargin{c+1};
            
            case ('nrestart') % this could be LFP recording to help the estimation process
                nrestart = varargin{c+1};
                
        end % switch
        c = c + 2;
    end % for
end % if

% Basis functions...
%  history
hist.delay = delay(2)/(dt*1000);
hist.basis = getBasis('rcos',hist.nfilt,hist.delay,20,0);
hist.basis = orth(hist.basis')';

%  coupling
coupling.delay = delay(1)/(dt*1000);
coupling.basis = getBasis('rcos',coupling.nfilt,coupling.delay,20,0);

S = double(getSpkMat(population,dt,T,1));
population{1}=find(S(1,:))*dt;

XY = getX(S(2,:),hist.basis,0,1,0)';    % Generate covariates using the basis functions...

optimres=[]; fmin=Inf;
tic
for rand_res=1:nrestart
    x([1:2 5]) = abs(randn(3,1));
    x(3:4) = rand(2,1);
    x(5) = x(5)*5*sign(randn(1)); % what is this?
    
    ff=[]; xx=[];  llhd_val=[]; bb=[]; 
    xx(1,:)=x;
    
    for rep=1:20
        
        % Fit GLM...
        [PSP,~,~] = eTM_modified(population{1},[1 x(3)],x);
        mP=mean(PSP);
        PSP=PSP/abs(mP);
        PSP_full = S(1,:);
        PSP_full(PSP_full~=0) = PSP;
        XX = getX(PSP_full,coupling.basis,0,1,0)';    % Generate covariates using the basis functions...
        X=[XX XY noise];
        
        if rep==1 && rand_res==1
            % initial glm
            [bta,~] = glmfit(X,S(2,:)','poisson');
        else
            % warm-start (using previous bta)
            [bta,~] = glmfit(X,S(2,:)','poisson','B0',bta);
        end
        coupling_filter = bta(2:(size(coupling.basis,1)+1))'*coupling.basis;
        cfnorm = sqrt(sum(coupling_filter.^2));
        if rep==1 && rand_res==1, x(end)=cfnorm*mP; xx(1,:)=x; end % use estimated amplitude as starting point
        
        coupling_filter = coupling_filter/cfnorm;
        offset = [XY noise]*bta((size(coupling.basis,1)+2):end)+bta(1);
        lambda = exp(X*bta(2:end)+bta(1));
        llhd_val(rep+1) = -(S(2,:)*log(lambda+(lambda==0))-sum(lambda));
        
        % Fit eTM params...        
        
        funObj = @(x) fun_opt_TM(x,population,T,S,dt,coupling_filter,offset);
        LB = [zeros(4,1)-10;-1000];
        UB = [log([5;5;1;1]-10e-3);1000];
        PLB = [zeros(4,1)-10;-100];
        PUB = [log([2;2;1;1]-10e-3);100];
    
        [x,f] = bads(funObj,[log(x(1:4)) x(5)],LB',UB',PLB',PUB');
        x(1:4)=exp(x(1:4));

        % Store path...
        ff(rep+1) = f;
        xx(rep+1,:)=x';
        bb(rep+1,:)=bta';
        
        fprintf('.')
        if rep>2 && abs(ff(end)-ff(end-1))<10e-6, break; end
    end
    if ff(end)<fmin
        param=x';
        param_bta=bta;
        fmin=ff(end);
    end
    optimres(rand_res).param_path=xx;
    optimres(rand_res).bta_path=bb;
    optimres(rand_res).param_ff=ff;
    optimres(rand_res).bta_ff=llhd_val;
    optimres(rand_res).coupling = coupling;
    optimres(rand_res).hist = hist;

    toc
end