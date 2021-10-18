%% Goal: Run sensitivity analysis full model on server
% function used to run the sensitivity analysis of the full model on the
% server with the good parameter sets. The SA can be split in chunks to allow for
% running on multiple screens. 
% N = number of parameter sets per chunk (total amount of parameter sets = 4482)
% q = chuck of parameters. (parameters are divided in groups of two, so
% using q=1 simulates sensitivity analysis of parameter 1 and 2)
%% run FM_sensitivity_analysis function
function FM_sensitivity_analysis(N,q)
% load good parameter sets
load('median_low_hoksokratio2','params_GhighF_ratio')
for i=q*2-1:q*2
     p_new = params_GhighF_ratio(:,1:22);
     %divide new parameters one by one by ten
     p_new(:,i) = params_GhighF_ratio(:,i)./10;
     %simulate with new parameter sets
     [sim] = optimise_model(p_new,N,i);
     save(sprintf("outputFM_SAGlowRatio_div10_%i.mat",i),'p_new','sim','-v7.3');
end
for i=q*2-1:q*2
     p_new = params_GhighF_ratio(:,1:22);
     %multiply new parameters one by one by ten
     p_new(:,i) = params_GhighF_ratio(:,i).*10;
     %simulate with new parameter sets
     [sim] = optimise_model(p_new,N,i);
     save(sprintf("outputFM_SAGlowRatio_times10_%i.mat",i),'p_new','sim','-v7.3');
end
end
