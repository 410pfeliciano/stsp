function [PSP,R,u] = eTM_modified(spk,initial,param)
% gets the plasticity parameteres and gives back the PSPs for those
 
 
% by: Abed Ghanbari 
%     abed.ghanbari@uconn.edu
% -------------------------------------------------------------------------
 
% Costa, R. P., Sjöström, P. J., & van Rossum, M. C. W. (2013). 
% Probabilistic inference of short-term synaptic plasticity in neocortical 
% microcircuits. 
% Frontiers in Computational Neuroscience, 7(June), 75. 
% doi:10.3389/fncom.2013.00075
 
% -------------------------------------------------------------------------
 
% spk : spike's time vector
% initial.R0, initial.u0 : initial values
% psp :  post synaptic potential
 
% param.D : depression time constant
% param.F : facilitation time constant
% param.f : change in number of vesicles after facilitation process
% param.U : change in number of vesicles after depression process
% param.A : is an amplitude factor that includes the number of release 
% sites, the properties and number of postsynaptic receptors, and cable 
% filtering.
 
% -------------------------------------------------------------------------
 
D = param.D ;
F = param.F ;
f = param.f ;
U = param.U ;
A = param.A ;
 
L = length(spk);
 
R = zeros(1,L);
u = zeros(1,L);
 
R(1) = initial.R0;
u(1) = initial.u0;
 
PSP = zeros(1,L);
PSP(1)=A*R(1)*u(1);
for i = 1 : L-1
    
    R(i+1) = 1 - ( 1 - R(i) * (1-u(i)) ) * exp(- (spk(i+1)-spk(i)) / D ) ;
    
    u(i+1) = U + ( u(i) + f * ( 1 - u(i)) - U  )* exp(- (spk(i+1)-spk(i)) / F ) ;
    
    PSP(i+1) = A * R(i+1) * u(i+1) ;
end

