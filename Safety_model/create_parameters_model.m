%% Goal: create random parameter sets for the full model using latin hypercube sampling
%Latin hypercube creates random numbers between 0-1. 
%N = number of parameter samples to simulate
%% Function to create parameter sets
function create_parameters_model(N)
    %load 95% confidence intervals from the toy model
    load('TM_CI_pwr.mat','TM_pwr_all');
    
    p_sets = lhsdesign(N,22);

    % Parameter a1 alpha constitutive promoter
    p_sets(:,1) = 10.^((4*p_sets(:,1))-12);
    % Parameter a2 alpha frm leakage
    p_sets(:,2) = 10.^((2*p_sets(:,2))-14);
    % Parameter a3 alpha frm max
    p_sets(:,3) = 10.^((4*p_sets(:,3))-12);
    % Parameter a4 lac promoter leakage
    p_sets(:,4) = 10.^((TM_pwr_all(1,4)*p_sets(:,4))+TM_pwr_all(1,2));
    % Parameter a5 lac max transcription rate
    p_sets(:,5) = 10.^((TM_pwr_all(2,4)*p_sets(:,5))+TM_pwr_all(2,2));
    % Parameter b1 degradation rate mFrmR
    p_sets(:,6) = 10.^((6*p_sets(:,6))-6);
    % Parameter b2 degradation rate FrmR
    p_sets(:,7) = 10.^((6*p_sets(:,7))-6);
    % Parameter b3 degradation rate mLacI
    p_sets(:,8) = 10.^((6*p_sets(:,8))-6);
    % Parameter b4 degradation rate LacI
    p_sets(:,9) = 2.3*10.^((2*p_sets(:,9))-5);
    % Parameter b5 degradation rate sok
    p_sets(:,10) = 10.^((TM_pwr_all(4,4)*p_sets(:,10))+TM_pwr_all(4,2));
    % Parameter b6 degradation rate hoksok
    p_sets(:,11) = 10.^((TM_pwr_all(5,4)*p_sets(:,11))+TM_pwr_all(5,2));
    % Parameter b7 degradation rate hok
    p_sets(:,12) = 10.^((TM_pwr_all(6,4)*p_sets(:,12))+TM_pwr_all(6,2));
    % Parameter g1 translation rate mfrmR
    p_sets(:,13) = 10.^((5*p_sets(:,13))-2);
    % Parameter g2 translation rate mlacI
    p_sets(:,14) = 10.^((5*p_sets(:,14))-2);
    % Parameter k1 unbinding formaldehyde-frmR complex
    p_sets(:,15) = 10.^((6*p_sets(:,15))-1);
    % Parameter k2 binding formaldehyde and frmR
    p_sets(:,16) = 10.^((5*p_sets(:,16))-2);
    % Parameter k3 binding frmR in hill equation
    p_sets(:,17) = 26.3*10.^((2*p_sets(:,17))-7);
    % Parameter k4 binding hok and sok
    p_sets(:,18) = 10.^((TM_pwr_all(7,4)*p_sets(:,18))+TM_pwr_all(7,2));
    % Parameter k5 unbinding hok and sok
    p_sets(:,19) = 10.^((TM_pwr_all(8,4)*p_sets(:,19))+TM_pwr_all(8,2));
    % Parameter k6 binding LacI in hill equation
    p_sets(:,20) = 10.^((3*p_sets(:,20))-12);
    % Parameter n1 hill equation frmR
    p_sets(:,21) = round(1*p_sets(:,21)+1);
    % Parameter n2 hill equation lacI
    p_sets(:,22) = round(2*p_sets(:,22)+1);

    %save p-sets
    save('model_parameters.mat','p_sets')

end