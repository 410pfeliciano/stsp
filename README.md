# Estimating short-term synaptic plasticity from pre- and postsynaptic spiking

for convinience all codes related to estimating STP parameters are in ```stp_glm``` folder

```addpath(genpath('stp_glm'))```

- data for the in vitro experiment is available in ```in vitro data```

# Demo
```stp_demo.m``` generates the pre and postsynaptic spiking activity in a connection with TM parameters in ```true_params```
and estimates the parameters using both TM-GLM and GBLM descibed in:

Ghanbari, A. & Malyshev, A. & Volgushev, M. & Stevenson, I. (2017)
Estimating short-term synaptic plasticity from pre- and postsynaptic spiking (http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005738)

here we generated pre and postsynaptic spikes from an LIF neuron - replace that with your own data

```[Tpre, Tpost] = LIFoutput(T,20,50,true_params,1);```


a sample of estimated parameters for a facilitating neuron with only few hundreds of spikes and T=100 sec

![Marginals](https://raw.githubusercontent.com/abedghanbari2/stsp/master/facilitation_screenshot.png)


# References

- Acerbi, L. & Ma, W. J. (2017) Bayesian Adaptive Direct Search (BADS) optimization algorithm for model fitting in MATLAB (https://github.com/lacerbi/bads)

