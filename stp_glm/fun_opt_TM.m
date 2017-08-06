function llhd = fun_opt_TM(x,population,T,S,dt,coupling_filter,offset)

if any(~isfinite(x))
    keyboard
end

x(1:4) = exp(x(1:4));
param.D = x(1);               % depression time constant
param.F = x(2);               % facilitation time const
param.U = x(3);               % baseline probability of vesicle release
param.f = x(4);               % speed coefficent for decaying to the baseline 
param.A = x(5);               % Amplitude factor related to the number of relewase sites 
param.dt = dt;
param.sim_time = T;

initial.R0 = 1;
initial.u0 = param.U;
[PSP,~,~] = eTM_modified(population{1},initial,param);

PSP_full = S(1,:);
PSP_full(PSP_full==1) = PSP;

Xc = getX(PSP_full,coupling_filter,0,1,0)';
lambda = exp(Xc+offset);
lambda(lambda>2000)=2000;

llhd = -(S(2,:)*log(lambda+(lambda==0))-sum(lambda));

    if param.f>1 || param.f<10e-10 || param.U>1 || param.U<10e-10
        llhd = 10e6;
    end


    if ~isfinite(llhd)
        keyboard
    end
end