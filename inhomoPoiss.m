% random rate generator

% gets :
% T : time
% nknots : number of knots
% isCirc : is it circular
% fr : average firing rate

function [tspike,lambda] = inhomoPoiss(T,dt,nknots,fr)
% nknots = 50;
% T=50;
% fr=10;
% dt=.0001;
y=randn(1,nknots);
xi=linspace(0,T,T/dt);
x=linspace(0,T,nknots);
yi=spline(x,y,xi);
lambda = exp(yi)/mean(exp(yi))*fr;
tspike=spike_poiss2(T,dt,lambda);

%% plot
% stem(tspike,'.')
% hold on
% plot(lambda/max(lambda),'r','LineWidth',2)
