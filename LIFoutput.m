function [Tpre, Tpost] = LIFoutput(T,rate,nknots,params,Adap_class)

%  poisson input - could be changed to the inhomogeneous poisson spiking
% rate = 20; % hz
% T = 500; % seconds
dt = .0001; % should be less than the tau in LIF function
% total number of steps
N = round(T ./ dt);
% time - for plotting
% time = linspace(0, T, N);

% generate poisson spike time for presynaptic neuron
% nknots =50;
[s_pre,~] = inhomoPoiss(T,dt,nknots,rate);

% generate weighted spikes using TM model

initial.R0 = 1;
initial.u0 = params(3);
Tpre= find(s_pre~=0)*dt;
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

% LIF output
I = 1e-9*current;
I0 = 1e-9*ones(1,N);
randI = 1e-9 .* random('Normal', 0, 1, [1, N]);
% keyboard
if Adap_class==1 % no adaptation
	[ s_post_poiss , ~ ] = LIF_SFA(.8*I+.8*I0+5*randI,0);
else % with adaptation
	[ s_post_poiss , ~ ] = LIF_SFA(1.5*I+1*I0+10*randI,100e-9);
end

% Tpost_poiss_noAdap = find(s_post_poiss_noAdap~=0)*dt;
Tpost = find(s_post_poiss~=0)*dt;