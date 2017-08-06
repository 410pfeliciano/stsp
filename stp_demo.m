clc;clear all;close all;
addpath(genpath('stp_glm'))

% simulating pre and postsynaptic responses ...
T=50;

% TM parameters [D F U f]

% true_params = [ 1.7 .02 .7  .05]; % strong depression
% true_params = [.5  .05 .5  .05];%  depression
% true_params = [.2  .2  .25 .3];% facilitation/depression
% true_params = [.05 .5  .15 .15];% facilitation
true_params = [.05 1 .01  .1];% strong facilitation

[Tpre, Tpost] = LIFoutput(T,20,50,true_params,1);
population{1}=Tpre;
population{2}=Tpost;

%% tm-glm
delay = [70 200];
nfilt = [6 5];
dt=.001;
[param_bta, est_params, optimres]= ...
    stp_tmglm(population,T,'dt',dt,'delay',delay,'nfilt',nfilt,'nrestart',1);

%% gblm
S = double(getSpkMat(population,dt,T,1));
theta = ...
    stp_gblm(S(1,:),S(2,:),'delay',delay,'nfilt',nfilt,'numCV',4,'toleranceValue',5);

%% plot
figure,
subplot(3,3,1)
plot((est_params(5)*optimres(1).bta_path(end,2:nfilt(1)+1)*optimres(1).coupling.basis))
title('coupling filter');xlim([0 100]);
subplot(3,3,2)
plot((optimres(1).bta_path(end,nfilt(1)+2:end)*optimres(1).hist.basis))
title('history filter')

[corr,~] = corr_fast(Tpre, Tpost,-.01,.1,110);
subplot(3,3,4)
bar(linspace(-.01,.1,110),corr,'k');xlim([-.01 .1]);
title('cross-correlation pre/post')

[corr,~] = corr_fast(Tpost, Tpost,-.1,.1,200);
subplot(3,3,5)
bar(linspace(-.1,.1,200),corr,'k');xlim([-.1 .1]);
title('auto-correlation - post')

[V_est,current_est,t] = markram_response(est_params(1:4)');
[V,current] = markram_response(true_params);
% [V,current,t] = markram_response([.05 1 .01  .1]);
subplot(3,3,3);hold on
plot(t,current_est)
plot(t,current,'r')
title('markram response')
legend('estimated','true')

subplot(3,3,6);hold on
plot(t,V_est(1:length(t)))
plot(t,V(1:length(t)),'r')
title('markram response')

[corr_sim_short,corr_sim_long] = short_hist(Tpre,Tpost,-.005,.05,150);
t=linspace(-5,50,150);
subplot(3,3,7)
bar(t,corr_sim_long,'k');xlim([0 50]);
yl1 = ylim;
title('splitxcorr - recovered')

subplot(3,3,8)
bar(t,corr_sim_short,'k');xlim([0 50]);ylim(yl1)
title('splitxcorr - burst')
subplot(3,3,9)
plot((theta.modif_fxn'));xlim([0 1000])
title('modification function')
