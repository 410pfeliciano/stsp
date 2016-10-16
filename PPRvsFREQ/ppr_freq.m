% creates the ground truth of freq vs ppr for short term synaptic
% plasticity
% assumed that there is a long interaval between the two spikes
clc;clear all;close all
par = [ 1.7 .02 .7  .05;
    .5  .05 .5  .05;
    .2  .2  .25 .3;
    .05 .5  .15 .15;
    .02 1.7 .1  .11];
C=colormap(hot(10));
for j=1:5
    for freq=1:50

    PSP=[];
    cArray = num2cell(par(j,:));
    [D,F,U,f] = deal(cArray{:});
    dt=1/freq;
    spk = [0 20 20+dt];
    R(1) = 1;
    u(1) = U;
    for i = 1:2
        R(i+1) = 1 - ( 1 - R(i) * (1-u(i)) ) * exp(- (spk(i+1)-spk(i)) / D ) ;
        u(i+1) = U + ( u(i) + f * ( 1 - u(i)) - U  )* exp(- (spk(i+1)-spk(i)) / F ) ;
        PSP(i+1) = R(i+1) * u(i+1) ;
    end
    ppr(freq)=PSP(end)/PSP(end-1);
    end
    hold on
    plot(ppr,'Color',C(6-j,:))
end
