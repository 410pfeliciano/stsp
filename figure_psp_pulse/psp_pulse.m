clc;close all;
% this file creates the PSP amplitudes of a Markram Stimuli for a train
% spikes with a fixed interval
freq = [5 30]; % Hz
dt = 1./freq; % \Delta_t
Ns = [5 5]; % number of spikes
par = [ 1.7 .02 .7  .05;
        .5  .05 .5  .05;
        .2  .2  .25 .3;
        .05 .5  .15 .15;
        .02 1.7 .1  .11];
%         .05 1 .01  .1];
%         0   0   1   0];
C=colormap(hot(10));
for k=1:2
PSP=[];
for j =1:5
cArray = num2cell(par(j,:));
[D,F,U,f] = deal(cArray{:});
spk = [0 20:dt(k):20+Ns(k)*dt(k)];
R(1) = 1;
u(1) = U;
for i = 1:Ns(k)
    R(i+1) = 1 - ( 1 - R(i) * (1-u(i)) ) * exp(- (spk(i+1)-spk(i)) / D ) ;
    u(i+1) = U + ( u(i) + f * ( 1 - u(i)) - U  )* exp(- (spk(i+1)-spk(i)) / F ) ;
    PSP(i+1) = R(i+1) * u(i+1) ;
end
subplot(1,2,k)
hold on
plot((spk(3:end)-spk(3))*1000,PSP(2:end)/PSP(2),'-*','Color',C(6-j,:))
end
ylim([0 4])
end