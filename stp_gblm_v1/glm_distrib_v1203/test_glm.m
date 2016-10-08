
%% Simulate some glm neurons...

num_neurons = 3;
sim_time    = 100;   % seconds
params.dt   = 0.001; % binsize

% Basis functions...
mprops.nfilt = 5;
mprops.delay = 500/(params.dt*1000);
mprops.basis = getBasis('rcos',mprops.nfilt,mprops.delay,20,0);

figure(1); clf
subplot(2,3,1)
plot(mprops.basis')
xlabel('Time [bins]')
title('Basis Functions')

%% Random sparse connectivity...
C=[];
for i=1:num_neurons
    alph=[];
    for j=1:num_neurons
        if i==j % Self-connections
            alph = [alph (randn(1,mprops.nfilt)-1)*mprops.basis];
        elseif rand(1)>0.5 % Sparsify
            alph = [alph zeros(1,mprops.nfilt)*mprops.basis];
        else % All other connections
            alph = [alph randn(1,mprops.nfilt)*mprops.basis];
        end
    end
    C=[C;alph];
end
C = [rand(num_neurons,1)*0 C];    % Baseline firing rates

subplot(2,3,2:3)
plot(C(:,2:end)'/5+repmat([1:num_neurons],size(C,2)-1,1))
axis tight
title('Simulated Connectivity')

firings = simLNP_trsX(C,params.dt,sim_time,[],[]);
Tlist = firings2tlist(firings);
subplot(2,3,4:6)
drawRaster(Tlist,1,1,[0 120])
title('Simulated Spikes')
ylabel('Neuron')
xlabel('Time [s]')

%% Get connectivity Covariates
S = getSpkMat(Tlist,params.dt,sim_time,1);

mprops.nfilt = 5;
mprops.delay = 500/(params.dt*1000);
mprops.basis = getBasis('rcos',mprops.nfilt,mprops.delay,20,0);
mprops.obasis = orth(mprops.basis')';
X = getX(S,mprops.obasis,0,1,0)';    % Generate covariates using the basis functions...

%% Fit..

% Regularization hyper-parameters...
% nu_vec = [0 logspace(-3,1,10)];
nu_vec = 0;

% For each neuron
for neuron=1:num_neurons
    % L1-regularized GLM with cross-validation...
    fprintf('\n\nNeuron %i/%i...\n',neuron,num_neurons)
    m(neuron) = fitCVGLM(X,S(neuron,:)',10,nu_vec,'lossGLM_poiss_exp');
end

%% Plot results...
Chat=[];
for i=1:num_neurons
    for j=1:num_neurons
        idx = ((j-1)*size(mprops.basis,1)+1):(j*size(mprops.basis,1));
        odx = ((j-1)*size(mprops.basis,2)+1):(j*size(mprops.basis,2));
        Chat(i,odx+1) = mprops.obasis'*(m(i).breg(idx+1,1));
    end
    Chat(i,1) = (m(i).breg(1,1));
end

figure(2); clf
plot(C(:,2:end)'/5+repmat([1:num_neurons],size(C,2)-1,1))
axis tight
hold on
plot(Chat(:,2:end)'/5+repmat([1:num_neurons],size(C,2)-1,1),':')
title('True and Estimated Connectivity')

%%
figure(3); clf
for neuron=1:num_neurons
%     plot(mean(m(neuron).llrt,1))
    plot(m(neuron).llrt)
    hold on
end
hold off
axis tight
title('Hyper-Parameter Optimization')
ylabel('Cross-Validated Log Likelihood Ratio')
xlabel('Hyperparameter Index')

%%
figure(4); clf
% for neuron=1:num_neurons
    plot(squeeze(m(1).b(2:end,1,:))')
%     hold on
% end
% hold off
% axis tight