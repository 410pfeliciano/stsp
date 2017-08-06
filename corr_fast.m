function [d,Ctot] = corr_fast(Tpre,Tpost,Ta,Tb,bin)
if length(Tpre)>1e5 || length(Tpost)>1e5;keyboard;end
% Tlist_pre : presynaptic spiking times
% Tlist_post : postsynaptic spiking times
% Ta : starting point of the histogram before the presynaptic onset
% Tb : ending interval 
% bin : number of bins
% plot : 0/1 plot/don't plot
k=1;
Ctot = [];
Tchunk = 1000;
while (k-1)*Tchunk<max(Tpre)
    Tlist_pre = [];
    Tlist_post = [];
    Tlist_pre = Tpre(Tpre>(k-1)*Tchunk & Tpre<k*Tchunk);
    Tlist_post = Tpost(Tpost>(k-1)*Tchunk & Tpost<k*Tchunk);

    if size(Tlist_post,1)==1
        A = repmat(Tlist_post,length(Tlist_pre),1);
    else
        A = repmat(Tlist_post',length(Tlist_pre),1);
    end

    if size(Tlist_pre,1)==1
        B = repmat(Tlist_pre',1,length(Tlist_post));
    else
        B = repmat(Tlist_pre,1,length(Tlist_post));
    end

    C = A - B ;

    C(C>Tb | C<Ta | C==0) = NaN;
    C = C(:);
    Ctot = [Ctot C'];
    k = k+1;

end
d = hist(Ctot,bin);