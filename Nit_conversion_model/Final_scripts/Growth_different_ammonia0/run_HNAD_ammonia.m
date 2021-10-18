%% run_HNAD_ammonia is a function to run model with multiple ammonia concentrations
%n = amount of initial conditions to test, now 7 are specified in
%sim_HNAD_am

function run_HNAD_ammonia(n)
%load model parameters (1*36 double)
load('p_sets_best.mat','p_sets_SA')
% convert to linear scale
p_sets = 10.^(p_sets_SA);

%call function to simulate the system
[output] = sim_HNAD_am(p_sets(1,:),n);

% save output into file specify k if multiple files are generated
% sequentially
k = 1;

%specify name
save(sprintf("Output_ammonia_range_%i.mat",k),'output');
end

