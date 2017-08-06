function [spikes,V] = LIF_SFA(I,varargin)
prog=0; % report progress on loop

E=-70e-3;
Rm=10e6;
tau_m=10e-3;
Vth=-54e-3;
Vreset=-80e-3;
Vspike = 20e-3;
E_K = -80e-3;       % reversal potential for spike-rate adaptation (SRA) and 
                    % refractory current

% Iapp = 4e-9;        
dt=0.0001;
% tmax=0.7;
% tpulse=0.1;
% lengthpulse=0.5;
% Nt=tmax/dt;
Nt = length(I);
tmax = Nt*dt;
T=0:dt:tmax;
% I=zeros(size(T));
V=NaN(size(T));
gsra=zeros(size(T));    % spike-rate adaptaion conductance
tau_sra = 0.1;          % time for spike-rate adaptation to decay
% increase in spike-rate conductance per spike
if length(varargin)==1
    delta_gsra = varargin{1};
else
    delta_gsra = 100e-9;     % increase in refractory conductance per spike
end
gref=zeros(size(T));    % spike refractory conductance
tau_ref = 0.002;        % time for refractory conductance to decay
delta_gref = 200e-9;     % increase in refractory conductance per spike

V(1)=E;                 % begin simulation with membrane potential at leak level
gsra(1) = 0;            % initially no spike-rate adaptation
gref(1) = 0;            % initially no refractory current

spikes=zeros(size(T));  % vector to record spike times
% for i = tpulse/dt : (tpulse+lengthpulse)/dt
%     I(i) = Iapp;        % applied current during a pulse
% end
for i = 2:Nt
    if mod(100*i/Nt,5)==0 && prog==1
        fprintf('%2.2f \n',100*i/Nt)
    end
    % Next line solves tau_sra*(d gsra/dt) = -gsra
    gsra(i) = gsra(i-1)*(1-dt/tau_sra);
    % Next line solves tau_ref*(d gref/dt) = -gref
    gref(i) = gref(i-1)*(1-dt/tau_ref);
    
    if ( V(i-1) == Vspike ) 
        V(i) = Vreset;      % reset after spike
    else
        V(i) = V(i-1) + ...     %otherwise integrate 
            dt/tau_m*(E - V(i-1) +Rm*I(i-1) - ...
        (V(i-1)-E_K)*Rm*(gsra(i-1)+gref(i-1)) ); % current for SRA is (V-E_K)*gsra
    
    end
    
    if V(i) > Vth       % if voltage is greater than threshold spike
        V(i) = Vspike;  % change voltage so we see the spike
        spikes(i) = 1;  % record the time of the spike
        gsra(i) = gsra(i-1) + delta_gsra;   % increase the spike-rate-daptation conductance
        gref(i) = gref(i-1) + delta_gref;   % increase the refractory conductance
    end
end