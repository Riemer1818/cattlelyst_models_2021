%% Goal: Run full model on server
% function used to run the full model on the server with the generated
% parameter sets. The parameter sets can be split in chunks to allow for
% running on multiple screens. (used later on as well to simulate the
% parameter sets that gave a good hok/sok ratio after 12 hours for 48
% hours)
% N = number of parameter sets per chunk (total amount of parameter sets = 100,000)
% q = number of the chunk (1-10)
%% run_optimise_model function
function run_optimise_model(N,q) 
    % load parameter sets full model
    load('model_parameters.mat','p_sets')
    %load('median_low_hoksokratio2','params_GhighF_ratio')
    
    % make chunk of the parameter sets
    params = p_sets(1+N*(q-1):N*q,:);
    %params = params_GhighF_ratio(1+N*(q-1):N*q,:);
    
    % call optimise model function
    [sim] = optimise_model(params,N,q);
    
    % save simulations
    save(sprintf("FM_output_%i.mat",q),'params','sim');
    %save(sprintf("median_low_hoksokratio2_48hsim",q),'params','sim','-v7.3');
end