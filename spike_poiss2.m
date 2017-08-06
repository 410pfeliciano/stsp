function spike=spike_poiss2(T,dt,rr)
% just for simulation of a spike train
N = round(T/dt);
spike=zeros(1,N);
t=linspace(0,T,N);

for i=6:N
    if rand<1-exp(-rr(i)*dt) && sum(spike(i-5:i-1))==0;
       spike(i)=1;
    end
end