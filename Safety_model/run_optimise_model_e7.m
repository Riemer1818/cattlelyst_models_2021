%% Goal: Run model extension 7 on server
% function used to run model extension 7 on the server with the generated
% parameter sets. The parameter sets can be split in chunks to allow for
% running on multiple screens. 
% N = number of parameter sets per chunk (total amount of parameter sets = 10,000)
% q = number of the chunk (1-10)
%% run_optimise_model_e6 function
function run_optimise_model_e7(N,q)
    % load parameter sets model extension 7
    load('e7_parameters_o.mat','p_sets')
    
    % make chunk of the parameter sets
    params = p_sets(1+N*(q-1):N*q,:);
    
    % call optimise model extension 7 function
    [sim] = optimise_model_e7(params,N,q);
    
    % save simulations
    save(sprintf("e7_output_o_%i.mat",q),'params','sim','-v7.3');
end