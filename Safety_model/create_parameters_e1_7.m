%% Goal: create parameter sets for model extension 1+7 
%N = number of parameter samples to simulate (1)
%% Function to create parameter sets
function create_parameters_e1_7(N)
    %load good parameter sets from the full model
    load('median_low_hoksokratio2.mat','params_GhighF_ratio','median_lowhoksok');
    
    %load parameter sets from the extensions and give them different names
    load('e1_parameters_o.mat','p_sets');
    p_sets_e1 = p_sets;
    load('e7_parameters_o.mat','p_sets');
    p_sets_e7 = p_sets;
    
    % parameter numbers used from the full model and the extentions.
    % parameter nr 1056 is outlier in sensitivity from the full model.
    nr = 1056;
    nr1 = 4331; %p-set used from extension 1
    %nr7 = 3082; %max hoksok from extension 7
    nr7 = 2425; %hoksok ratio 6000 from extension 7
    
    p_sets = lhsdesign(N,35);

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
    % New parameters 1
    % Parameter b8
    p_sets(:,23) = p_sets_e1(nr1,23);
    % Parameter b9
    p_sets(:,24) = p_sets_e1(nr1,24);
    % Parameter b10
    p_sets(:,25) = p_sets_e1(nr1,25);
    % Parameter g3
    p_sets(:,26) = p_sets_e1(nr1,26);
    % Parameter k7
    p_sets(:,27) = p_sets_e1(nr1,27);
    % Parameter k8
    p_sets(:,28) = p_sets_e1(nr1,28);
    % new parameters e7
    % parameter a6
    p_sets(:,29) = p_sets_e7(nr7,23);
    % parameter a7 
    p_sets(:,30) = p_sets_e7(nr7,24);
    % parameter b11
    p_sets(:,31) = p_sets_e7(nr7,25);
    % parameter b12
    p_sets(:,32) = p_sets_e7(nr7,26);
    % parameter g4 
    p_sets(:,33) = p_sets_e7(nr7,27);
    % parameter k9
    p_sets(:,34) = p_sets_e7(nr7,28);
    % parameter n3
    p_sets(:,35) = p_sets_e7(nr7,29);

    %save p-sets
    save('e1_7_parameters_o.mat','p_sets')

end