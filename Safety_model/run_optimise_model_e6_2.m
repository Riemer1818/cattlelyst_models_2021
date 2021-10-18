%% Goal: Run model extension 6.1 on server
% function used to run model extension 6.1 on the server with the generated
% parameter sets. The parameter sets can be split in chunks to allow for
% running on multiple screens. 
% N = number of parameter sets per chunk (total amount of parameter sets = 10,000)
% q = number of the chunk (1-10)
%% run_optimise_model_e6_2 function
function run_optimise_model_e6_2(N,q) 
    % load parameter sets model extension 6.1 (same as parameters extension 6)
    load('e6_parameters_o.mat','p_sets')
    
    % make chunk of the parameter sets
    params = p_sets(1+N*(q-1):N*q,:);
    
    % call optimise model extension 6.1 function
    [sim] = optimise_model_e6_2(params,N,q);
    
    % save simulations
    save(sprintf("e6_2_output_o_%i.mat",q),'params','sim','-v7.3');
end