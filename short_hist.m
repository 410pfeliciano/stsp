function [corr_short,corr_long,C1,C2] = short_hist(T_pre,T_post,ta,tb,bin_corr)

% calculates the split cross-correlograms of the two neurons.
% T_pre: spike times of the presynaptic neuron
% T_post: spike times of the postsynaptic neuron
% ta: time(seconds) where the correlation window starts (ex: -.05)
% tb: time (seconds) where the window ends (ex: .15)
% bin_corr:number of bins for our correlation window

dt=.001;      
isi = diff(T_pre);

% spikes with short and long intervals defined by the quartiles
Ts = prctile(isi,25);
Tl = prctile(isi,75);
Tshort = T_pre(find(isi<Ts)+1);
Tlong = T_pre(find(isi>Tl)+1);

% spliding window for sequnces with short and long intervals
[corr_short,C1] = corr_fast(Tshort, T_post,ta,tb,bin_corr); 
[corr_long,C2] = corr_fast(Tlong, T_post,ta,tb,bin_corr);

% keyboard
corr_short=corr_short/length(Tshort)/dt; % gives firing rates in Hz
corr_long=corr_long/length(Tlong)/dt;