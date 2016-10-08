function llhd = fun_opt(x,population,T,mpropsC,XY,S,str,isConstrained,beta_h,dt,noise)
% this is the function that designed to have the variables to be searched
% in the liklihood space 

%INPUTS
%x: the four parameters that we are going to search for in ther 4D space
%population: the pre and post synaptic timestamps of the APs
%mprops: contains the filters
%   mprops.delay
%   mprops.nfilt
%   mrops.basis

%OUTPUT
%llhd_opt: the scalar value of the optimization function which is the
%negative of the liklihood function since we are minimizing the function

param.D = x(1);               % depression time constant
param.F = x(2);               % facilitation time const
param.U = x(3);               % baseline probability of vesicle release
param.f = x(4);               % speed coefficent for decaying to the baseline 
param.A = 1;               % Amplitude factor related to the number of relewase sites 
param.dt = dt;
param.sim_time = T;

initial.R0 = 1;
initial.u0 = param.U;
[PSP,~,~] = eTM_modified(population{1},initial,param);
PSP = PSP/mean(PSP);
PSP_full = getSpkMat(population{1},dt,T,1);
PSP_full = double(PSP_full);
PSP_full(PSP_full==1) = PSP;
XX = getX(PSP_full,mpropsC.basis,0,1,0)';    % Generate covariates using the basis functions...

X=[XX XY noise];

if (~isempty(beta_h))
    bta = beta_h;
else
    if rand<.2
        bta = glmfit(X,S(2,:),'poisson');
        save(['beta',str,'.mat'], 'bta')
    else
        if exist(['beta',str,'.mat'],'file')
            load(['beta',str,'.mat'], 'bta');
        else
            bta = glmfit(X,S(2,:),'poisson');
            save(['beta',str,'.mat'], 'bta')
        end
    end
end

lambda = exp(X*bta(2:end)+bta(1));

if isConstrained
    llhd = -(S(2,:)*log(lambda)-sum(lambda));
else
    llhd = -(S(2,:)*log(lambda)-sum(lambda) ...
        +(log(cauchypdf(param.D,0,2.5))+log(cauchypdf(param.F,0,2.5))...
        +log(betapdf(param.U,2,2))+log(betapdf(param.f,.9,2))));
end