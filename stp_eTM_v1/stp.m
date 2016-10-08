function [bta, param]= stp(population,T,varargin)

%-INPUTS-------------------------------------------------------------------
% population: a cell array [2x1]
%           - first cell: presynaptic neuron's time-stamps of the APs
%           - second cell: postsynaptic neuron's time-stamps of the APs
% T: duration of the recorded pre and postsynaptic activity [scalar]

%-OUTPUT-------------------------------------------------------------------
% param: estimated parameters [D,F,U,f]

%-Description--------------------------------------------------------------
%/ here we first define an interval to search for the parameters that
%gives the maximum liklihood and then using the simplex method search for
%the values in the 4D space

% default
beta_h=[];
theta0 = [2*rand 2*rand rand rand];
isConstrained = false;
delay =[150 150];
dt = .001;
mprops.nfilt = 5;
mpropsC.nfilt = 5;
noise=[]; 
if (~isempty(varargin))
    c = 1 ;

    % user defined
    while c <= length(varargin)
        switch varargin{c}
            case {'theta0'}
                theta0 = varargin{c+1};
                c = c +1;
            case {'beta'} % adds the fixed beta and doesn't compute the beta each time
                beta_h = varargin{c+1};
                c = c +1;
            case {'constrained'}
                isConstrained = true;
                upper = varargin{c+1};
                c = c +1;            
            case {'delay'}
                delay = varargin{c+1};
                c = c +1;
            case {'nfilt'}
                mprops.nfilt = varargin{c+1}(2);
                mpropsC.nfilt = varargin{c+1}(1);
                c = c +1;
            case ('dt')
                dt = varargin{c+1};
                c=c+1;
            case ('noise') % this could be LFP recording to help the estimation process
                noise = varargin{c+1};
                c=c+1;
                
        end % switch
        c = c + 1;
    end % for
end % if


% Basis functions...
% for history
mprops.delay = delay(2)/(dt*1000);
mprops.basis = getBasis('rcos',mprops.nfilt,mprops.delay,20,0);
mprops.basis = orth(mprops.basis')';
% for coupling
mpropsC.delay = delay(1)/(dt*1000);
mpropsC.basis = getBasis('rcos',mpropsC.nfilt,mpropsC.delay,20,0);
mpropsC.basis = orth(mpropsC.basis')';

% coupling for common input
if length(delay)>2
    mpropsC.delay2 = delay(3)/(dt*1000);
    mpropsC.basis2 = getBasis('rcos',mprops.nfilt,mpropsC.delay2,20,0);
    mpropsC.basis2 = orth(mpropsC.basis2')';
    mpropsC.basis = padarray(mpropsC.basis,[0,length(mpropsC.basis2)-length(mpropsC.basis)],'post');
    mpropsC.basis = [mpropsC.basis;mpropsC.basis2];
end

S = getSpkMat(population,dt,T,1);
S= double(S);

XY = getX(S(2,:),mprops.basis,0,1,0)';    % Generate covariates using the basis functions...
str=num2str(floor(rand*1000000)); % to save and load the beta in optimization step
%%
options = [];
nVars = 4;
LB = zeros(nVars,1);
UB = [2;2;1;1];

funObj = @(x) fun_opt(x,population,T,mpropsC,XY,S,str,isConstrained,beta_h,dt,noise);
options.numDiff = 1;
% options.verbose = 0;
options.maxIter = 100;
param = minConf_TMP(funObj,theta0',LB,UB,options);
%%
if isempty(beta_h)
    load(['beta',str,'.mat'], 'bta');
    delete(['beta',str,'.mat']);
else
    bta = beta_h;
end