%% Goal: create random parameter sets for the toy model using latin hypercube sampling
%Latin hypercube creates random numbers between 0-1. 
%N = number of parameter samples to simulate
%data from literature is saved (file is not used in analysis)
%% Function to create parameter sets
function create_parameters_toymodel(N)
    p_sets = lhsdesign(N,10);

    % Parameter a4 hok leakage transcription rate
    p_sets(:,1) = 10.^((2*p_sets(:,1))-14);
    % Parameter a5 hok max transcription rate
    p_sets(:,2) = 10.^((4*p_sets(:,2))-12);
    % Parameter a6 sok transcription rate
    p_sets(:,3) = 10.^((4*p_sets(:,2))-12);
    % Parameter b5 degradation rate sok
    p_sets(:,4) = 2.3*10.^((2*p_sets(:,4))-3);
    % Parameter b6 degradation rate hoksok
    p_sets(:,5) = 10.^((6*p_sets(:,5))-6);
    % Parameter b7 degradation rate hok
    p_sets(:,6) = 5.8*10.^((2*p_sets(:,6))-5);
    % Parameter k5 binding hok and sok
    p_sets(:,7) = 3*10.^((4*p_sets(:,7))+2);
    % Parameter k6 unbinding hok and sok
    p_sets(:,8) = 10.^((5*p_sets(:,8))-2);
    % Parameter k7 binding IPTG in hill equation
    p_sets(:,9) = 10.^((7*p_sets(:,9))-12);
    % Parameter n3
    p_sets(:,10) = round(4*p_sets(:,10)+1);

    %save p-sets
    save('TM_parameters1.mat','p_sets')

    %save data
    data (:,:,1) = [0 18; 3600 200];
    data (:,:,2) = [0 18; 3600 6000];
    save('data.mat','data');
end