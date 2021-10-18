%% Goal: create random parameter sets for model extension 7 using latin hypercube sampling
%Latin hypercube creates random numbers between 0-1. 
%N = number of parameter samples to simulate
%% Function to create parameter sets
function create_parameters_e7(N)
    %load good parameter sets from the full model
    load('median_low_hoksokratio2.mat','params_GhighF_ratio','median_lowhoksok');
    
    %calculate average values for certain parameter types
    average_a_l = (median_lowhoksok(1,2)+median_lowhoksok(1,4))/2;
    average_a_m = (median_lowhoksok(1,1)+median_lowhoksok(1,3)+median_lowhoksok(1,5))/3;
    average_b_m = (median_lowhoksok(1,6)+median_lowhoksok(1,8)+median_lowhoksok(1,10)+median_lowhoksok(1,11)+median_lowhoksok(1,12))/5;
    average_b_p = (median_lowhoksok(1,7)+median_lowhoksok(1,9))/2;
    average_g = (median_lowhoksok(1,13)+median_lowhoksok(1,14))/2;
    average_k_h = (median_lowhoksok(1,17)+median_lowhoksok(1,20))/2;
    
    % parameter nr 1056 is outlier in sensitivity. parameter nr 2241 is the
    % parameter set giving the median hoksok ratio 
    nr = 1056;
    
    p_sets = lhsdesign(N,29);

    % Parameter a1 alpha T7
    p_sets(:,1) = params_GhighF_ratio(nr,1);
    % Parameter a2 alpha frm leakage
    p_sets(:,2) = params_GhighF_ratio(nr,2);
    % Parameter a3 alpha frm max
    p_sets(:,3) = params_GhighF_ratio(nr,3);
    % Parameter a4 lac promoter leakage
    p_sets(:,4) = params_GhighF_ratio(nr,4);
    % Parameter a5 lac max transcription rate
    p_sets(:,5) = params_GhighF_ratio(nr,5);
    % Parameter b1 degradation rate mFrmR
    p_sets(:,6) = params_GhighF_ratio(nr,6);
    % Parameter b2 degradation rate FrmR
    p_sets(:,7) = params_GhighF_ratio(nr,7);
    % Parameter b3 degradation rate mLacI
    p_sets(:,8) = params_GhighF_ratio(nr,8);
    % Parameter b4 degradation rate LacI
    p_sets(:,9) = params_GhighF_ratio(nr,9);
    % Parameter b5 degradation rate sok
    p_sets(:,10) = params_GhighF_ratio(nr,10);
    % Parameter b6 degradation rate hoksok
    p_sets(:,11) = params_GhighF_ratio(nr,11);
    % Parameter b7 degradation rate hok
    p_sets(:,12) = params_GhighF_ratio(nr,12);
    % Parameter g1 translation rate mfrmR
    p_sets(:,13) = params_GhighF_ratio(nr,13);
    % Parameter g2 translation rate mlacI
    p_sets(:,14) = params_GhighF_ratio(nr,14);
    % Parameter k1 unbinding formaldehyde-frmR complex
    p_sets(:,15) = params_GhighF_ratio(nr,15);
    % Parameter k2 binding formaldehyde and frmR
    p_sets(:,16) = params_GhighF_ratio(nr,16);
    % Parameter k3 binding frmR in hill equation
    p_sets(:,17) = params_GhighF_ratio(nr,17);
    % Parameter k4 binding hok and sok
    p_sets(:,18) = params_GhighF_ratio(nr,18);
    % Parameter k5 unbinding hok and sok
    p_sets(:,19) = params_GhighF_ratio(nr,19);
    % Parameter k6 binding LacI in hill equation
    p_sets(:,20) = params_GhighF_ratio(nr,20);
    % Parameter n1 hill equation frmR
    p_sets(:,21) = params_GhighF_ratio(nr,21);
    % Parameter n2 hill equation lacI
    p_sets(:,22) = params_GhighF_ratio(nr,22);
    % New parameters
    % Parameter a6
    p_sets(:,23) = average_a_l*10.^((2*p_sets(:,23))-1);
    % Parameter a7
    p_sets(:,24) = average_a_m*10.^((2*p_sets(:,24))-1);
    % Parameter b8
    p_sets(:,25) = average_b_m*10.^((2*p_sets(:,25))-1);
    % Parameter b9
    p_sets(:,26) = average_b_p*10.^((2*p_sets(:,26))-1);
    % Parameter g3
    p_sets(:,27) = average_g*10.^((2*p_sets(:,27))-1);
    % Parameter k7
    p_sets(:,28) = average_k_h*10.^((2*p_sets(:,28))-1);
    % Parameter n3
    p_sets(:,29) = round(2*p_sets(:,29)+1);

    %save p-sets
    save('e7_parameters_o.mat','p_sets')

end