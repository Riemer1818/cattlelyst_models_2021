%% file used to analyse the output of the simulations of the toy model
% contains:
% 1. Select good parameter sets
% 2. Histograms parameter values
% 3. Scatter plots (parameter value vs parameter value)
% 4. Median values + 95% confidence intervals
% 5. time series plot toy model of selected parameter set

% input files: 
%     TM_output1_1.mat
% output files:
%     TM_gp.mat
%     TM_CI_pwr.mat
%% 1. Select good parameter sets

% load output toy model 
    load('TM_output1_1.mat')
% get all simulation values at t=1h (101) 
    no_IPTG = zeros(100000, 4); 
    with_IPTG= zeros(100000, 4);
    for i=1:length(sim)
         no_IPTG(i,1:length(sim{i}{1}(101,:))) = sim{i}{1}(101,:);
         with_IPTG(i,1:length(sim{i}{2}(101,:))) = sim{i}{2}(101,:);            
    end
% make a matrix in which the parameter sets are numbered
    params_all = [[1:100000]' score score_ij params];

% make a matrix with the numbers of the parameter sets (column 1), time
% (2), sok concentration (3), hok-sok concentration (4), hok concentration
% (5) and hok/sok ratio (6)
    A = [[1:100000]' with_IPTG (with_IPTG(:,3)+with_IPTG(:,4))./(with_IPTG(:,3)+with_IPTG(:,2))];
    A_no = [[1:100000]' no_IPTG (no_IPTG(:,3)+no_IPTG(:,4))./(no_IPTG(:,3)+no_IPTG(:,2))];
% Select parameters with a hoksok ratio without IPTG between 10-30
    B = sortrows(A,6);
    gp_A = B(58112:61859,1);
    gp_B = A_no(gp_A,:);
% Select parameters with a hoksok ratio with IPTG between 5000-7000
    C = sortrows(gp_B,6);
    gp = C(2667:3173,1);
    gp_params = params_all(gp,:);
    GP_A = A(gp,:);
    GP_A_no =A_no(gp,:);

% Save good parameter sets
    save('TM_gp.mat','gp_params','GP_A','GP_A_no');

%% 2. Histograms parameter values
% Used to check the distribution of the values of the parameters of the good parameter sets

% load good parameter sets
    load('TM_gp.mat');

% histograms
    figure(1)
    subplot(3,1,1);   
    [~,edges] = histcounts(log10(gp_params(:,5)));
    histogram(gp_params(:,5),10.^edges);
    set(gca, 'XScale', 'log');
    title('hok promoter leakage');
    subplot(3,1,2);
    [~,edges] = histcounts(log10(gp_params(:,6)));
    histogram(gp_params(:,6),10.^edges);
    set(gca, 'XScale', 'log');
    title('hok max transcription rate');
    subplot(3,1,3);
    [~,edges] = histcounts(log10(gp_params(:,7)));
    histogram(gp_params(:,7),10.^edges);
    set(gca, 'XScale', 'log');
    title('sok promoter leakage');

    figure(2)
    subplot(3,1,1);
    edges = [2.3*10^-3 2.3*10^-2.75 2.3*10^-2.5 2.3*10^-2.25 2.3*10^-2 2.3*10^-1.75 2.3*10^-1.5 2.3*10^-1.25 2.3*10^-1];
    histogram(gp_params(:,8),edges);
    set(gca, 'XScale', 'log');
    title('degradation rate sok');
    subplot(3,1,2);
    [~,edges] = histcounts(log10(gp_params(:,9)));
    histogram(gp_params(:,9),10.^edges);
    set(gca, 'XScale', 'log');
    title('degradation rate hoksok');
    subplot(3,1,3);
    edges = [5.8*10^-5 5.8*10^-4.75 5.8*10^-4.5 5.8*10^-4.25 5.8*10^-4 5.8*10^-3.75 5.8*10^-3.5 5.8*10^-3.25 5.8*10^-3];
    N = histcounts(log10(gp_params(:,10)),edges);
    histogram(gp_params(:,10),edges);
    set(gca, 'XScale', 'log');
    title('degradation rate hok');

    figure(3)
    subplot(2,1,1);
    [~,edges] = histcounts(log10(gp_params(:,11)));
    histogram(gp_params(:,11),10.^edges);
    set(gca, 'XScale', 'log');
    title('binding hok and sok');
    subplot(2,1,2);
    [~,edges] = histcounts(log10(gp_params(:,12)));
    histogram(gp_params(:,12),10.^edges);
    set(gca, 'XScale', 'log');
    title('unbinding hok and sok');

    figure(4)
    subplot(2,1,1);
    [~,edges] = histcounts(log10(gp_params(:,13)));
    histogram(gp_params(:,13),10.^edges);
    set(gca, 'XScale', 'log');
    title('binding IPTG in hill equation');
    subplot(2,1,2);
    [~,edges] = histcounts((gp_params(:,14)));
    histogram(gp_params(:,14),edges);
    title('n3 hill equation');


%% 3. Scatter plots (parameter value vs parameter value)
% from the good parameter sets, scatter plots are made to check if there is
% a relation between the values of the good parameter sets. Value of one
% parameter is on the x-axis and the values of the other parameters per
% subplot on the y-axis. 

% load good parameter sets
    load('TM_gp.mat');
    
% scatter plots
    figure(5);
    subplot(3,3,2);
    scatter(gp_params(:,5),gp_params(:,6));
    set(gca,'xscale','log');
    set(gca,'yscale','log','FontSize',13);
    xlabel('a2');
    ylabel('a3');
    title('max. transcription rate hok');
    subplot(3,3,1);
    scatter(gp_params(:,5),gp_params(:,7));
    set(gca,'xscale','log');
    set(gca,'yscale','log','FontSize',13);
    xlabel('a2');
    ylabel('a1');
    title('max. transcription rate sok');
    subplot(3,3,3);
    scatter(gp_params(:,5),gp_params(:,8));
    set(gca,'xscale','log');
    set(gca,'yscale','log','FontSize',13);
    xlabel('a2');
    ylabel('b1');
    title('degradation rate sok');
    subplot(3,3,4);
    scatter(gp_params(:,5),gp_params(:,9));
    set(gca,'xscale','log');
    set(gca,'yscale','log','FontSize',13);
    xlabel('a2');
    ylabel('b2');
    title('degradation rate hok-sok');
    subplot(3,3,5);
    scatter(gp_params(:,5),gp_params(:,10));
    set(gca,'xscale','log');
    set(gca,'yscale','log','FontSize',13);
    xlabel('a2');
    ylabel('b3');
    title('degradation rate hok');
    subplot(3,3,6);
    scatter(gp_params(:,5),gp_params(:,11));
    set(gca,'xscale','log');
    set(gca,'yscale','log','FontSize',13);
    xlabel('a2');
    ylabel('k1');
    title('binding rate hok and sok');
    subplot(3,3,7);
    scatter(gp_params(:,5),gp_params(:,12));
    set(gca,'xscale','log');
    set(gca,'yscale','log','FontSize',13);
    xlabel('a2');
    ylabel('k2');
    title('unbinding rate hok-sok');
    subplot(3,3,8);
    scatter(gp_params(:,5),gp_params(:,13));
    set(gca,'xscale','log');
    set(gca,'yscale','log','FontSize',13);
    xlabel('a2');
    ylabel('k3');
    title('binding IPTG in Hill function');
    subplot(3,3,9);
    scatter(gp_params(:,5),gp_params(:,14));
    set(gca,'xscale','log','FontSize',13);
    xlabel('a2');
    ylabel('n');
    title('Hill coefficient');
    sgtitle('background hok production (a2) vs values other parameters','Fontsize',22);

    figure(6);
    subplot(3,3,2);
    scatter(gp_params(:,6),gp_params(:,5));
    set(gca,'xscale','log');
    set(gca,'yscale','log','FontSize',13);
    xlabel('a3');
    ylabel('a2');
    title('background hok production');
    subplot(3,3,1);
    scatter(gp_params(:,6),gp_params(:,7));
    set(gca,'xscale','log');
    set(gca,'yscale','log','FontSize',13);
    xlabel('a3');
    ylabel('a1');
    title('max. transcription rate sok');
    subplot(3,3,3);
    scatter(gp_params(:,6),gp_params(:,8));
    set(gca,'xscale','log');
    set(gca,'yscale','log','FontSize',13);
    xlabel('a3');
    ylabel('b1');
    title('degradation rate sok');
    subplot(3,3,4);
    scatter(gp_params(:,6),gp_params(:,9));
    set(gca,'xscale','log');
    set(gca,'yscale','log','FontSize',13);
    xlabel('a3');
    ylabel('b2');
    title('degradation rate hok-sok');
    subplot(3,3,5);
    scatter(gp_params(:,6),gp_params(:,10));
    set(gca,'xscale','log');
    set(gca,'yscale','log','FontSize',13);
    xlabel('a3');
    ylabel('b3');
    title('degradation rate hok');
    subplot(3,3,6);
    scatter(gp_params(:,6),gp_params(:,11));
    set(gca,'xscale','log');
    set(gca,'yscale','log','FontSize',13);
    xlabel('a3');
    ylabel('k1');
    title('binding rate hok and sok');
    subplot(3,3,7);
    scatter(gp_params(:,6),gp_params(:,12));
    set(gca,'xscale','log');
    set(gca,'yscale','log','FontSize',13);
    xlabel('a3');
    ylabel('k2');
    title('unbinding rate hok-sok');
    subplot(3,3,8);
    scatter(gp_params(:,6),gp_params(:,13));
    set(gca,'xscale','log');
    set(gca,'yscale','log','FontSize',13);
    xlabel('a3');
    ylabel('k3');
    title('binding IPTG in Hill function');
    subplot(3,3,9);
    scatter(gp_params(:,6),gp_params(:,14));
    set(gca,'xscale','log','FontSize',13);
    xlabel('a3');
    ylabel('n');
    title('Hill coefficient');
    sgtitle('max. transcription rate hok (a3) vs values other parameters','Fontsize',22);

    figure(7);
    subplot(3,3,1);
    scatter(gp_params(:,7),gp_params(:,5));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('a1');
    ylabel('a2');
    title('background hok production');
    subplot(3,3,2);
    scatter(gp_params(:,7),gp_params(:,6));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('a1');
    ylabel('a3');
    title('max. transcription rate hok');
    subplot(3,3,3);
    scatter(gp_params(:,7),gp_params(:,8));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('a1');
    ylabel('b1');
    title('degradation rate sok');
    subplot(3,3,4);
    scatter(gp_params(:,7),gp_params(:,9));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('a1');
    ylabel('b2');
    title('degradation rate hok-sok');
    subplot(3,3,5);
    scatter(gp_params(:,7),gp_params(:,10));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('a1');
    ylabel('b3');
    title('degradation rate hok');
    subplot(3,3,6);
    scatter(gp_params(:,7),gp_params(:,11));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('a1');
    ylabel('k1');
    title('binding rate hok and sok');
    subplot(3,3,7);
    scatter(gp_params(:,7),gp_params(:,12));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('a1');
    ylabel('k2');
    title('unbinding rate hok-sok');
    subplot(3,3,8);
    scatter(gp_params(:,7),gp_params(:,13));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('a1');
    ylabel('k3');
    title('binding IPTG in Hill function');
    subplot(3,3,9);
    scatter(gp_params(:,7),gp_params(:,14));
    set(gca,'xscale','log','FontSize',13);
    xlabel('a1');
    ylabel('n');
    title('Hill coefficient');
    sgtitle('max. transcription rate sok (a1) vs values other parameters','Fontsize',22);

    figure(8);
    subplot(3,3,2);
    scatter(gp_params(:,8),gp_params(:,5));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('b1');
    ylabel('a2');
    title('background hok production');
    subplot(3,3,3);
    scatter(gp_params(:,8),gp_params(:,6));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('b1');
    ylabel('a3');
    title('max. transcription rate hok');
    subplot(3,3,1);
    scatter(gp_params(:,8),gp_params(:,7));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('b1');
    ylabel('a1');
    title('max. transcription rate sok');
    subplot(3,3,4);
    scatter(gp_params(:,8),gp_params(:,9));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('b1');
    ylabel('b2');
    title('degradation rate hok-sok');
    subplot(3,3,5);
    scatter(gp_params(:,8),gp_params(:,10));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('b1');
    ylabel('b3');
    title('degradation rate hok');
    subplot(3,3,6);
    scatter(gp_params(:,8),gp_params(:,11));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('b1');
    ylabel('k1');
    title('binding rate hok and sok');
    subplot(3,3,7);
    scatter(gp_params(:,8),gp_params(:,12));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('b1');
    ylabel('k2');
    title('unbinding rate hok-sok');
    subplot(3,3,8);
    scatter(gp_params(:,8),gp_params(:,13));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('b1');
    ylabel('k3');
    title('binding IPTG in Hill function');
    subplot(3,3,9);
    scatter(gp_params(:,8),gp_params(:,14));
    set(gca,'xscale','log','FontSize',13);
    xlabel('b1');
    ylabel('n');
    title('Hill coefficient');
    sgtitle('degradation rate sok (b1) vs values other parameters','Fontsize',22);

    figure(9);
    subplot(3,3,2);
    scatter(gp_params(:,9),gp_params(:,5));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('b2');
    ylabel('a2');
    title('background hok production');
    subplot(3,3,3);
    scatter(gp_params(:,9),gp_params(:,6));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('b2');
    ylabel('a3');
    title('max. transcription rate hok');
    subplot(3,3,1);
    scatter(gp_params(:,9),gp_params(:,7));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('b2');
    ylabel('a1');
    title('max. transcription rate sok');
    subplot(3,3,4);
    scatter(gp_params(:,9),gp_params(:,8));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('b2');
    ylabel('b1');
    title('degradation rate sok');
    subplot(3,3,5);
    scatter(gp_params(:,9),gp_params(:,10));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('b2');
    ylabel('b3');
    title('degradation rate hok');
    subplot(3,3,6);
    scatter(gp_params(:,9),gp_params(:,11));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('b2');
    ylabel('k1');
    title('binding rate hok and sok');
    subplot(3,3,7);
    scatter(gp_params(:,9),gp_params(:,12));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('b2');
    ylabel('k2');
    title('unbinding rate hok-sok');
    subplot(3,3,8);
    scatter(gp_params(:,9),gp_params(:,13));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('b2');
    ylabel('k3');
    title('binding IPTG in Hill function');
    subplot(3,3,9);
    scatter(gp_params(:,9),gp_params(:,14));
    set(gca,'xscale','log','FontSize',13);
    xlabel('b2');
    ylabel('n');
    title('Hill coefficient');
    sgtitle('degradation rate hok-sok (b2) vs values other parameters','Fontsize',22);

    figure(10);
    subplot(3,3,2);
    scatter(gp_params(:,10),gp_params(:,5));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('b3');
    ylabel('a2');
    title('background hok production');
    subplot(3,3,3);
    scatter(gp_params(:,10),gp_params(:,6));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('b3');
    ylabel('a3');
    title('max. transcription rate hok');
    subplot(3,3,1);
    scatter(gp_params(:,10),gp_params(:,7));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('b3');
    ylabel('a1');
    title('max. transcription rate sok');
    subplot(3,3,4);
    scatter(gp_params(:,10),gp_params(:,8));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('b3');
    ylabel('b1');
    title('degradation rate sok');
    subplot(3,3,5);
    scatter(gp_params(:,10),gp_params(:,9));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('b3');
    ylabel('b2');
    title('degradation rate hok-sok');
    subplot(3,3,6);
    scatter(gp_params(:,10),gp_params(:,11));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('b3');
    ylabel('k1');
    title('binding rate hok and sok');
    subplot(3,3,7);
    scatter(gp_params(:,10),gp_params(:,12));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('b3');
    ylabel('k2');
    title('unbinding rate hok-sok');
    subplot(3,3,8);
    scatter(gp_params(:,10),gp_params(:,13));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('b3');
    ylabel('k3');
    title('binding IPTG in Hill function');
    subplot(3,3,9);
    scatter(gp_params(:,10),gp_params(:,14));
    set(gca,'xscale','log','FontSize',13);
    xlabel('b3');
    ylabel('n');
    title('Hill coefficient');
    sgtitle('degradation rate hok (b3) vs values other parameters','Fontsize',22);

    figure(11);
    subplot(3,3,2);
    scatter(gp_params(:,11),gp_params(:,5));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('k1');
    ylabel('a2');
    title('background hok production');
    subplot(3,3,3);
    scatter(gp_params(:,11),gp_params(:,6));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('k1');
    ylabel('a3');
    title('max. transcription rate hok');
    subplot(3,3,1);
    scatter(gp_params(:,11),gp_params(:,7));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('k1');
    ylabel('a1');
    title('max. transcription rate sok');
    subplot(3,3,4);
    scatter(gp_params(:,11),gp_params(:,8));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('k1');
    ylabel('b1');
    title('degradation rate sok');
    subplot(3,3,5);
    scatter(gp_params(:,11),gp_params(:,9));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('k1');
    ylabel('b2');
    title('degradation rate hok-sok');
    subplot(3,3,6);
    scatter(gp_params(:,11),gp_params(:,10));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('k1');
    ylabel('b3');
    title('degradation rate hok');
    subplot(3,3,7);
    scatter(gp_params(:,11),gp_params(:,12));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('k1');
    ylabel('k2');
    title('unbinding rate hok-sok');
    subplot(3,3,8);
    scatter(gp_params(:,11),gp_params(:,13));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('k1');
    ylabel('k3');
    title('binding IPTG in Hill function');
    subplot(3,3,9);
    scatter(gp_params(:,11),gp_params(:,14));
    set(gca,'xscale','log','FontSize',13);
    xlabel('k1');
    ylabel('n');
    title('Hill coefficient');
    sgtitle('binding rate hok and sok (k1) vs values other parameters','Fontsize',22);

    figure(12);
    subplot(3,3,2);
    scatter(gp_params(:,12),gp_params(:,5));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('k2');
    ylabel('a2');
    title('background hok production');
    subplot(3,3,3);
    scatter(gp_params(:,12),gp_params(:,6));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('k2');
    ylabel('a3');
    title('max. transcription rate hok');
    subplot(3,3,1);
    scatter(gp_params(:,12),gp_params(:,7));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('k2');
    ylabel('a1');
    title('max. transcription rate sok');
    subplot(3,3,4);
    scatter(gp_params(:,12),gp_params(:,8));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('k2');
    ylabel('b1');
    title('degradation rate sok');
    subplot(3,3,5);
    scatter(gp_params(:,12),gp_params(:,9));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('k2');
    ylabel('b2');
    title('degradation rate hok-sok');
    subplot(3,3,6);
    scatter(gp_params(:,12),gp_params(:,10));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('k2');
    ylabel('b3');
    title('degradation rate hok');
    subplot(3,3,7);
    scatter(gp_params(:,12),gp_params(:,11));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('k2');
    ylabel('k1');
    title('binding rate hok and sok');
    subplot(3,3,8);
    scatter(gp_params(:,12),gp_params(:,13));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('k2');
    ylabel('k3');
    title('binding IPTG in Hill function');
    subplot(3,3,9);
    scatter(gp_params(:,12),gp_params(:,14));
    set(gca,'xscale','log','FontSize',13);
    xlabel('k2');
    ylabel('n');
    title('Hill coefficient');
    sgtitle('unbinding rate hok-sok (k2) vs values other parameters','Fontsize',22);

    figure(13);
    subplot(3,3,2);
    scatter(gp_params(:,13),gp_params(:,5));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('k3');
    ylabel('a2');
    title('background hok production');
    subplot(3,3,3);
    scatter(gp_params(:,13),gp_params(:,6));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('k3');
    ylabel('a3');
    title('max. transcription rate hok');
    subplot(3,3,1);
    scatter(gp_params(:,13),gp_params(:,7));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('k3');
    ylabel('a1');
    title('max. transcription rate sok');
    subplot(3,3,4);
    scatter(gp_params(:,13),gp_params(:,8));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('k3');
    ylabel('b1');
    title('degradation rate sok');
    subplot(3,3,5);
    scatter(gp_params(:,13),gp_params(:,9));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('k3');
    ylabel('b2');
    title('degradation rate hok-sok');
    subplot(3,3,6);
    scatter(gp_params(:,13),gp_params(:,10));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('k3');
    ylabel('b3');
    title('degradation rate hok');
    subplot(3,3,7);
    scatter(gp_params(:,13),gp_params(:,11));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('k3');
    ylabel('k1');
    title('binding rate hok and sok');
    subplot(3,3,8);
    scatter(gp_params(:,13),gp_params(:,12));
    set(gca,'xscale','log','FontSize',13);
    set(gca,'yscale','log');
    xlabel('k3');
    ylabel('k2');
    title('unbinding rate hok-sok');
    subplot(3,3,9);
    scatter(gp_params(:,13),gp_params(:,14));
    set(gca,'xscale','log','FontSize',13);
    xlabel('k3');
    ylabel('n');
    title('Hill coefficient');
    sgtitle('binding IPTG Hill function (k3) vs values other parameters','Fontsize',22);

%% 4. Median values + 95% confidence intervals
% The CI_lower and CI_higher values are calculated using the formulas
% described in this source; https://www-users.york.ac.uk/~mb55/intro/cicent.htm

%load good parameter sets from the toy model
    load('TM_gp.mat');

%calculate with parameter number denotes the 95% confidence intervals. 
    CI_lower = round(507*0.5-1.96*(507*0.5*(1-0.5))^(0.5));
    CI_higher = round(507*0.5+1.96*(507*0.5*(1-0.5))^(0.5));

%sort the parameter values of the good parameter sets
    gp_params_a4 = sortrows(gp_params,5);
    gp_params_a5 = sortrows(gp_params,6);
    gp_params_a6 = sortrows(gp_params,7);
    gp_params_b5 = sortrows(gp_params,8);
    gp_params_b6 = sortrows(gp_params,9);
    gp_params_b7 = sortrows(gp_params,10);
    gp_params_k5 = sortrows(gp_params,11);
    gp_params_k6 = sortrows(gp_params,12);
    gp_params_k7 = sortrows(gp_params,13);
    gp_params_n3 = sortrows(gp_params,14);

% create a matrix (M_all) containin in column 1 the median value of the parameters,
% in column 2 the lower confidence interval and in column 3 the upper
% confidence interval. 
    M_a4 = [median(gp_params_a4(:,5)) gp_params_a4(CI_lower,5) gp_params_a4(CI_higher,5)];
    M_a5 = [median(gp_params_a5(:,6)) gp_params_a5(CI_lower,6) gp_params_a5(CI_higher,6)];
    M_a6 = [median(gp_params_a6(:,7)) gp_params_a6(CI_lower,7) gp_params_a6(CI_higher,7)];
    M_b5 = [median(gp_params_b5(:,8)) gp_params_b5(CI_lower,8) gp_params_b5(CI_higher,8)];
    M_b6 = [median(gp_params_b7(:,9)) gp_params_b7(CI_lower,9) gp_params_b7(CI_higher,9)];
    M_b7 = [median(gp_params_b7(:,10)) gp_params_b7(CI_lower,10) gp_params_b7(CI_higher,10)];
    M_k5 = [median(gp_params_k5(:,11)) gp_params_k5(CI_lower,11) gp_params_k5(CI_higher,11)];
    M_k6 = [median(gp_params_k6(:,12)) gp_params_k6(CI_lower,12) gp_params_k6(CI_higher,12)];
    M_k7 = [median(gp_params_k7(:,13)) gp_params_k7(CI_lower,13) gp_params_k7(CI_higher,13)];
    M_n3 = [median(gp_params_n3(:,14)) gp_params_n3(CI_lower,14) gp_params_n3(CI_higher,14)];
    M_all = [M_a4; M_a5; M_a6; M_b5; M_b6; M_b7; M_k5; M_k6; M_k7; M_n3];

%This part calculates the x if you convert all median and CI values to 10^x. 
%First column; median. Second column; CI_lower. Third column; CI_higher. 
%The last column is the difference in x between CI_higher and CI_lower. 
    TM_pwr_all = [log10(M_all(:,1)) log10(M_all(:,2)) log10(M_all(:,3)) log10(M_all(:,3))-log10(M_all(:,2))];

%save the TM_pwr_all matrix
    save('TM_CI_pwr.mat','TM_pwr_all');

%% 5. time series plot toy model of selected parameter set
%load simulation files toy model
    load('TM_output1_1.mat')

%select start time (min value = 2) + which parameter set to simulate
    tstart = 2;
    nr = 9561;

%calculate hoksok ratio (column 2) 
    hoksok_ratio_TS_noIPTG = [[tstart:length(sim{1,nr}{2}(:,1))]' (sim{1,nr}{1}(tstart:end,4)+sim{1,nr}{2}(tstart:end,3))./(sim{1,nr}{2}(tstart:end,3)+sim{1,nr}{2}(tstart:end,2))];
    hoksok_ratio_TS_withIPTG= [[tstart:length(sim{1,nr}{2}(:,1))]' (sim{1,nr}{2}(tstart:end,4)+sim{1,nr}{1}(tstart:end,3))./(sim{1,nr}{1}(tstart:end,3)+sim{1,nr}{1}(tstart:end,2))];

%time series plots
    figure(14);
    subplot(2,2,1);
    plot(sim{1,nr}{1,1}(tstart:end,1),sim{1,1}{1,1}(tstart:end,2),'-o','MarkerIndices',1:10:length(sim{1,nr}{1,2}(tstart:end,2)),'LineWidth',1);
    hold on
    plot(sim{1,nr}{1,2}(tstart:end,1),sim{1,1}{1,2}(tstart:end,2),'LineWidth',1);
    hold off
    set(gca,'FontSize',14);
    xlabel('time (h)','Fontsize',16);
    ylabel('[sok] (nM)','Fontsize',16);
    legend('no IPTG','1mM IPTG');
    title('sok','Fontsize',18);
    subplot(2,2,2);
    plot(sim{1,nr}{1,1}(tstart:end,1),sim{1,1}{1,1}(tstart:end,3),'LineWidth',1);
    hold on
    plot(sim{1,nr}{1,2}(tstart:end,1),sim{1,1}{1,2}(tstart:end,3),'LineWidth',1);
    hold off
    set(gca,'FontSize',14);
    xlabel('time (h)','Fontsize',16);
    ylabel('[hok-sok] (nM)','Fontsize',16);
    set(gca,'yscale','log');
    legend('no IPTG','1mM IPTG');
    title('hok-sok complex','Fontsize',18);
    subplot(2,2,3);
    plot(sim{1,nr}{1,1}(tstart:end,1),sim{1,1}{1,1}(tstart:end,4),'LineWidth',1);
    hold on
    plot(sim{1,nr}{1,2}(tstart:end,1),sim{1,1}{1,2}(tstart:end,4),'LineWidth',1);
    hold off
    set(gca,'FontSize',14);
    xlabel('time (h)','Fontsize',16);
    ylabel('[hok] (nM)','Fontsize',16);
    set(gca,'yscale','log');
    legend('no IPTG','1mM IPTG');
    title('hok','Fontsize',18);
    subplot(2,2,4);
    plot(sim{1,nr}{1,1}(tstart:end,1),hoksok_ratio_TS_noIPTG(:,2),'LineWidth',1);
    hold on
    plot(sim{1,nr}{1,2}(tstart:end,1),hoksok_ratio_TS_withIPTG(:,2),'LineWidth',1);
    hold off
    set(gca,'FontSize',14);
    xlabel('time (h)','Fontsize',16);
    ylabel('hok/sok ratio','Fontsize',16);
    set(gca,'yscale','log');
    legend('no IPTG','1mM IPTG');
    title('hok/sok ratio','Fontsize',18);
    sgtitle('Time series plots - example good parameter set toy model','Fontsize',22);