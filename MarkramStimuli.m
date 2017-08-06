function s_blk = MarkramStimuli(T,dt,freq)
% dt = .0001;
N = T/dt ;
time = linspace(0,T,N);
spk = zeros(N,1);
indx = zeros(N,1);
for i = 1 : N
    if mod(i,round(1/dt/freq))==0
        spk(i) = 1;
    end
    indx(time>.5 & time<1)=1;
    indx(end-1/dt/2)=1;
    
end
s_blk=indx.*spk;