function [V,current,time] = markram_response(params)
%% make input current

T = 3; % seconds
dt = .001; % should be less than the tau in LIF function
% total number of steps
N = round(T ./ dt);
% time - for plotting
time = linspace(0, T, N);
freq = 20;
s_pre = MarkramStimuli(T,dt,freq);

initial.R0 = 1;
initial.u0 = params(3);%(STP_class,3);
Tpre= find(s_pre~=0)*dt;
% numReleaseSites = 64;
[PSP,~,~] = eTM_modified(Tpre,initial,[params 1]);
PSP = PSP/mean(PSP);
t_spike = s_pre;
t_spike(t_spike~=0)=PSP;

% generate current
tau = .01;
T_kernel = .25;
x0 = linspace(0,T_kernel,T_kernel/dt);
kernel_psc = exp(-x0/tau)-exp(-x0/(tau/10));
current = filter(kernel_psc,1,t_spike);

%% generate spike with the generated current
% 
I = .5e-9*current;
[~,V] = LIF_SFA(I);
