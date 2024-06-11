The code simulates the neuronal activity of an excitatory and inhibitory network following topological details presented in Menesse and Torres (2024). 

The script EEG_LIF_model.jl contains functions for creating the network and simulate its dynamic.

The script tau_mu_map.jl execute the functions to run the simulation given the values of two parameters of the model, the noise level mu and the synaptic resources recovery time scale tau_{rec}.

The output of the script gives:
  - File with average membrane potential time series of 5 groups of excitatory and inhibitory neurons.
  - File with binned time series of spiking state of N neurons of the network.

Reference:
  Menesse G, Torres J. 2023. Information dynamics efficiently discriminates high γ-rhythms in EEG brain waves. DOI: https://arxiv.org/abs/2311.13977 

The Information Dynamics analysis software used in the Reference is available in https://github.com/GuEMM/PhiID_Tools.git
