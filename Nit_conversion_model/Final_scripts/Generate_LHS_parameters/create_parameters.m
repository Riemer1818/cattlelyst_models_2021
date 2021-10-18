%% create parameters is a function to generate random parameter sets for simulations of HNAD for
%Pseudomonas stutzeri YZN_001, XL2
%N = amount of parameters to be created

function create_parameters(N)
%specify parameter bounds for all parameters make sure vectors are of
%similar length
min_values = [-3 -3 -3 -3 -3 -3 -3 -3 -3 -3 -3 -3 -3 -3 -3 -3 -3 -3 -3 -3 -3 -3 -3 -3 -3 -3 -3 -3 -3 -3 -3 -3 -4 -4 -3 -3];
max_values = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 -2 -2 1 1];

% Sampling parameters using Latin hypercube
power_p_sets = lhsdesign_modified(N,min_values,max_values);
p_sets = 10.^(power_p_sets);

% save parameter sets
save('parameters_$$date_.mat','p_sets')
end