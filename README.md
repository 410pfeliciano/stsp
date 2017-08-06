# Estimating short-term synaptic plasticity from pre- and postsynaptic spiking

for convinience all codes related to estimating STP parameters are in ```stp_glm``` folder

```addpath(genpath('stp_glm'))```

here we generated pre and postsynaptic spikes from a LIF neuron - replace that with your own data

```[Tpre, Tpost] = LIFoutput(T,20,50,true_params,1);```


a sample of estimated parameters for a facilitating neuron with only few hundreds of spikes and T=50

![Marginals](https://raw.githubusercontent.com/abedghanbari2/stsp/master/facilitation_screenshot.png)


# References

- Acerbi, L. & Ma, W. J. (2017) Bayesian Adaptive Direct Search (BADS) optimization algorithm for model fitting in MATLAB (https://github.com/lacerbi/bads)

