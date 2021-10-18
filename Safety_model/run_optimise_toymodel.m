%% Goal: Run toy model on server
% function used to run the toy model on the server with the generated
% parameter sets. The parameter sets can be split in chunks to allow for
% running on multiple screens. 
% N = number of parameter sets per chunk (total amount of parameter sets = 100,000)
% q = number of the chunk (1-10)
%% run_optimise_toymodel function
function run_optimise_toymodel(N,q)
    % load parameter sets toy model
    load('TM_parameters1.mat')

    % load data hok/sok ratio toy model. Note: both data and scores are not
    % used in analysis.
    load('data.mat')

    % make chunk of the parameter sets
    params = p_sets(1+N*(q-1):N*q,:);

    % call optimise toymodel function
    [score_ij,score,sim] = optimise_toymodel(data,params,N,q);

    % save simulations
    save(sprintf("TM_output1_%i.mat",q),'score_ij','score','params','sim');
end