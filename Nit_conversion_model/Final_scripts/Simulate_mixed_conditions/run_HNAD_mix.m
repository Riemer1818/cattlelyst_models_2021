%% run_HNAD_mix is a function to call the simulation of nitrogen dynamics in
%mixed conditions
%n = amount of initial conditions to test

function run_HNAD_mix(n)
%load model parameters (1*36 double)
load('p_sets_best.mat','p_sets_SA')
% convert to linear scale
p_sets = 10.^(p_sets_SA);

%specify which system to simulate

%call sim_HNAD_mix to simulate system for ZN1
[output] = sim_HNAD_mix_ZN1(p_sets(1,:),n);

%call sim_HNAD_mix_SDU10 to simulate system for SDU10
%[output] = sim_HNAD_mix_SDU10(p_sets(1,:),n);

% save output into file specify k if multiple files are generated
% sequentially
k = 1;

%specify file name
save(sprintf("Output_mixed_condition_%i.mat",k),'output');
end

