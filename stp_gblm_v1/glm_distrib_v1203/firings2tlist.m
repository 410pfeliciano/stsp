
function Tlist = firings2tlist(firings)

un = unique(firings(:,2));
for neuron=1:max(un)
    Tlist{neuron} = firings(firings(:,2)==neuron,1);
end
