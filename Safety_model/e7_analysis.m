%% file used to analyse the output of the simulations of the model extension 7
% contains:
% 1. create a matrix with end concentrations for the parameter sets that have the median p-set of the full model as first 22 parameter values
% 2. create a matrix with end concentrations for the parameter sets that have the sensitivity outlier p-set of the full model as first 22 parameter values
% 3. find the parameter sets with a the max sensitivity + corresponding timepoint
% 4. scatterplot hoksok ratio vs parameter value
% 5. Time series plot e7

% input files: 
%     e7_output_m_1.mat
%     e7_output_o_1.mat
%     e7_parameters_o.mat
%     median_low_hoksokratio2.mat
% output files:
%     e7_max_sensitivity_analysis.mat
%% 1. create a matrix with end concentrations for the parameter sets that have the median p-set of the full model as first 22 parameter values
% load output file median parameter set
    load('e7_output_m_1.mat','sim')
% creating a matrix with end concentrations for the simulations of each
% p-set.
    low_F_e = zeros(10000, 11); 
    high_F_e= zeros(10000, 11);
    for i=1:length(sim)
         low_F_e(i,1:length(sim{i}{1}(end,:))) = sim{i}{1}(end,:);
         high_F_e(i,1:length(sim{i}{2}(end,:))) = sim{i}{2}(end,:);            
    end
% Matrix A contains a new numbering of the simulations (first column) + end time (column 2) + end
% concentrations of the compounds (column 3-12) + hoksok ratio (column 13).
% Matrix B contains matrix A (column 1-13) + hoksok ratio [F]low - hoksok
% ratio [F]high + hoksok ratio [F]low/hoksok ratio [F]high.
    A_high_F_e = [[1:1000]' high_F_e (high_F_e(:,9)+high_F_e(:,8))./(high_F_e(:,8)+high_F_e(:,7))];
    A_low_F_e = [[1:1000]' low_F_e (low_F_e(:,9)+low_F_e(:,8))./(low_F_e(:,8)+low_F_e(:,7))];
    B_high_F_e = [A_high_F_e A_low_F_e(:,13)-A_high_F_e(:,13) A_low_F_e(:,13)./A_high_F_e(:,13)];
    B_low_F_e = [A_low_F_e A_low_F_e(:,13)-A_high_F_e(:,13) A_low_F_e(:,13)./A_high_F_e(:,13)];
    sortB = sortrows(B_high_F_e,14);
    C_high_F_e = [B_high_F_e A_low_F_e(:,11)./A_high_F_e(:,11) A_low_F_e(:,7)./A_high_F_e(:,7)];

%% 2. create a matrix with end concentrations for the parameter sets that have the sensitivity outlier p-set of the full model as first 22 parameter values

load('e7_output_o_1.mat','sim')
% creating a matrix with end concentrations for the simulations of each
% p-set.
    low_F_e = zeros(10000, 11); 
    high_F_e= zeros(10000, 11);
    for i=1:length(sim)
         low_F_e(i,1:length(sim{i}{1}(end,:))) = sim{i}{1}(end,:);
         high_F_e(i,1:length(sim{i}{2}(end,:))) = sim{i}{2}(end,:);            
    end
% Matrix A contains a new numbering of the simulations (first column) + end time (column 2) + end
% concentrations of the compounds (column 3-12) + hoksok ratio (column 13).
% Matrix B contains matrix A (column 1-13) + hoksok ratio [F]low - hoksok
% ratio [F]high + hoksok ratio [F]low/hoksok ratio [F]high.
    A_high_F_e = [[1:10000]' high_F_e (high_F_e(:,9)+high_F_e(:,8))./(high_F_e(:,8)+high_F_e(:,7))];
    A_low_F_e = [[1:10000]' low_F_e (low_F_e(:,9)+low_F_e(:,8))./(low_F_e(:,8)+low_F_e(:,7))];
    B_high_F_e = [A_high_F_e A_low_F_e(:,13)-A_high_F_e(:,13) A_low_F_e(:,13)./A_high_F_e(:,13)];
    B_low_F_e = [A_low_F_e A_low_F_e(:,13)-A_high_F_e(:,13) A_low_F_e(:,13)./A_high_F_e(:,13)];
    sortB = sortrows(B_high_F_e,15);
    C_high_F_e = [B_high_F_e A_low_F_e(:,11)./A_high_F_e(:,11) A_low_F_e(:,7)./A_high_F_e(:,7)];
    
%% 3. find the parameter sets with a the max sensitivity + corresponding timepoint
load('e7_output_o_1.mat','sim')

% add hoksok ratio in the sim array (12) and the sensitivity (13)
for i=1:length(sim)
     sim{i}{1}(:,12) = (sim{i}{1}(:,9)+sim{i}{1}(:,8))./(sim{i}{1}(:,8)+sim{i}{1}(:,7));
     sim{i}{2}(:,12) = (sim{i}{2}(:,9)+sim{i}{2}(:,8))./(sim{i}{2}(:,8)+sim{i}{2}(:,7));
     sim{i}{1}(:,13) = (sim{i}{1}(:,12))./(sim{i}{2}(:,12));
     sim{i}{2}(:,13) = (sim{i}{1}(:,12))./(sim{i}{2}(:,12));
end
% find the max sensitivities for each simulation and make a vector with the corresponding timepoints 
    low_F_e = zeros(10000, 11); 
    high_F_e= zeros(10000, 11);
    for i=1:length(sim)
         sim_max(i,1) = max(sim{i}{1}(:,13));
         t(i,1) = find(sim{i}{1}(:,13)==sim_max(i,1))./100;
    end
% sort the max sensitivities, to get the most sensitive p-set
    sim_max_n = [[1:10000]' t sim_max];
    sort_max = sortrows(sim_max_n,3);
save('e7_max_sensitivity_analysis.mat','sort_max')

%% 4. scatterplot hoksok ratio vs parameter value
load('e7_parameters.mat')

figure(4);
subplot(3,3,1);
scatter(A_high_F_e(:,13),p_sets(:,23));
set(gca,'yscale','log');
%set(gca,'xscale','log');
xlabel('ratio hoksok [F]_h_i_g_h');
ylabel('leakage expression promoter X');
title('alpha T7 (a1)');
subplot(3,3,2);
scatter(A_high_F_e(:,13),p_sets(:,24));
set(gca,'yscale','log');
%set(gca,'xscale','log');
xlabel('ratio hoksok [F]_h_i_g_h');
ylabel('max expression promoter X');
title('alpha frm leakage (a2)');
subplot(3,3,3);
scatter(A_high_F_e(:,13),p_sets(:,25));
set(gca,'yscale','log');
%set(gca,'xscale','log');
xlabel('ratio hoksok [F]_h_i_g_h');
ylabel('degradation rate mX');
title('alpha frm max (a3)');
subplot(3,3,4);
scatter(A_high_F_e(:,13),p_sets(:,26));
set(gca,'yscale','log');
%set(gca,'xscale','log');
xlabel('ratio hoksok [F]_h_i_g_h');
ylabel('degradation rate X');
title('lac promoter leakage (a4)');
subplot(3,3,5);
scatter(A_high_F_e(:,13),p_sets(:,27));
set(gca,'yscale','log');
%set(gca,'xscale','log');
xlabel('ratio hoksok [F]_h_i_g_h');
ylabel('translation rate mX');
title('lac max transcription rate (a5)');
subplot(3,3,6);
scatter(A_high_F_e(:,13),p_sets(:,28));
set(gca,'yscale','log');
%set(gca,'xscale','log');
xlabel('ratio hoksok [F]_h_i_g_h');
ylabel('X binding in hill equation');
title('degradation rate mFrmR (b1)');
subplot(3,3,7);
scatter(A_high_F_e(:,13),p_sets(:,29));
set(gca,'yscale','log');
%set(gca,'xscale','log');
xlabel('ratio hoksok [F]_h_i_g_h');
ylabel('n hill equation');
title('degradation rate FrmR (b2)');
sgtitle('e7 - hoksok ratio vs parameter value');

%% 5. Time series plot e7
load('median_low_hoksokratio2','params_GhighF_ratio');
load('e7_parameters_o','p_sets')
N = 1;
q = 1;

% select parameter set to simulate
% p_set = 309; %highest sensitivity
p_set = 3082; %highest hoksok ratio
%p_set = 2425; %hoksok ratio 6000
%p_set = 1561; %highest max sensitivity
p_set_m = 1056; %original biosafety circuit

% make simulation of given paremeter set for extension 7 (new) and original
% biosafety circuit (old)
[sim] = optimise_model_e7(p_sets(p_set,:),N,q);
[sim_old] = optimise_model(params_GhighF_ratio(p_set_m,:),N,q);
for i=1:length(sim_old)
     low_F(i,1:length(sim_old{i}{1}(end,:))) = sim_old{i}{1}(end,:);
     high_F(i,1:length(sim_old{i}{2}(end,:))) = sim_old{i}{2}(end,:);            
end
for i=1:length(sim)
     low_F_e(i,1:length(sim{i}{1}(end,:))) = sim{i}{1}(end,:);
     high_F_e(i,1:length(sim{i}{2}(end,:))) = sim{i}{2}(end,:);            
end

% Matrix A contains a time vector (first column) + time series data + hoksok ratio (last column).
% Matrix B contains matrix A (column 1-13) + hoksok ratio [F]low - hoksok
% ratio [F]high + hoksok ratio [F]low/hoksok ratio [F]high.
    A_high_F_TSnew = [[2:length(sim{1}{2}(:,1))]' sim{1}{2}(2:end,:) (sim{1}{2}(2:end,9)+sim{1}{2}(2:end,8))./(sim{1}{2}(2:end,8)+sim{1}{2}(2:end,7))];
    A_low_F_TSnew = [[2:length(sim{1}{1}(:,1))]' sim{1}{1}(2:end,:) (sim{1}{1}(2:end,9)+sim{1}{1}(2:end,8))./(sim{1}{1}(2:end,8)+sim{1}{1}(2:end,7))];
    B_high_F_TSnew = [A_high_F_TSnew A_low_F_TSnew(:,13)-A_high_F_TSnew(:,13) A_low_F_TSnew(:,13)./A_high_F_TSnew(:,13)];
    B_low_F_TSnew = [A_low_F_TSnew A_low_F_TSnew(:,13)-A_high_F_TSnew(:,13) A_low_F_TSnew(:,13)./A_high_F_TSnew(:,13)];

    A_high_F_TSold = [[2:length(sim_old{1}{2}(:,1))]' sim_old{1}{2}(2:end,:) (sim_old{1}{2}(2:end,9)+sim_old{1}{2}(2:end,8))./(sim_old{1}{2}(2:end,8)+sim_old{1}{2}(2:end,7))];
    A_low_F_TSold = [[2:length(sim_old{1}{1}(:,1))]' sim_old{1}{1}(2:end,:) (sim_old{1}{1}(2:end,9)+sim_old{1}{1}(2:end,8))./(sim_old{1}{1}(2:end,8)+sim_old{1}{1}(2:end,7))];
    B_high_F_TSold = [A_high_F_TSold A_low_F_TSold(:,11)-A_high_F_TSold(:,11) A_low_F_TSold(:,11)./A_high_F_TSold(:,11)];
    B_low_F_TSold = [A_low_F_TSold A_low_F_TSold(:,11)-A_high_F_TSold(:,11) A_low_F_TSold(:,11)./A_high_F_TSold(:,11)];

% time series plots
stp_indices = 1000; %changes frequency of the circles in the lines
tstart = 1;
figure(11);
subplot(2,2,1);
plot(sim{1,1}{1,1}(tstart:end,1),sim{1,1}{1,1}(tstart:end,2),'-o','MarkerIndices',1:stp_indices:length(sim{1,1}{1,1}(tstart:end,2)),'LineWidth',1);
hold on
plot(sim{1,1}{1,2}(tstart:end,1),sim{1,1}{1,2}(tstart:end,2),'LineWidth',1);
hold on
plot(sim_old{1,1}{1,1}(tstart:end,1),sim_old{1,1}{1,1}(tstart:end,2),'-o','MarkerIndices',1:stp_indices:length(sim{1,1}{1,2}(tstart:end,2)),'LineWidth',1);
hold on
plot(sim_old{1,1}{1,2}(tstart:end,1),sim_old{1,1}{1,2}(tstart:end,2),'LineWidth',1);
hold off
set(gca,'FontSize',14);
xlabel('time (h)','FontSize',16);
ylabel('[mFrmR] (nM)','FontSize',16);
legend('e7 [FA]_l_o_w','e7 [FA]_h_i_g_h','FM [FA]_l_o_w','FM [FA]_h_i_g_h');
title('mFrmR','FontSize',18);
subplot(2,2,2);
plot(sim{1,1}{1,1}(tstart:end,1),sim{1,1}{1,1}(tstart:end,3),'-o','MarkerIndices',1:stp_indices:length(sim{1,1}{1,1}(tstart:end,3)),'LineWidth',1);
hold on
plot(sim{1,1}{1,2}(tstart:end,1),sim{1,1}{1,2}(tstart:end,3),'LineWidth',1);
hold on
plot(sim_old{1,1}{1,1}(tstart:end,1),sim_old{1,1}{1,1}(tstart:end,3),'-o','MarkerIndices',1:stp_indices:length(sim{1,1}{1,2}(tstart:end,3)),'LineWidth',1);
hold on
plot(sim_old{1,1}{1,2}(tstart:end,1),sim_old{1,1}{1,2}(tstart:end,3),'LineWidth',1);
hold off
set(gca,'FontSize',14);
xlabel('time (h)','FontSize',16);
ylabel('[FrmR] (nM)','FontSize',16);
legend('e7 [FA]_l_o_w','e7 [FA]_h_i_g_h','FM [FA]_l_o_w','FM [FA]_h_i_g_h');
title('FrmR','FontSize',18);
subplot(2,2,3);
plot(sim{1,1}{1,1}(tstart:end,1),sim{1,1}{1,1}(tstart:end,4),'-o','MarkerIndices',1:stp_indices:length(sim{1,1}{1,1}(tstart:end,4)),'LineWidth',1);
hold on
plot(sim{1,1}{1,2}(tstart:end,1),sim{1,1}{1,2}(tstart:end,4),'LineWidth',1);
hold on
plot(sim_old{1,1}{1,1}(tstart:end,1),sim_old{1,1}{1,1}(tstart:end,4),'-o','MarkerIndices',1:stp_indices:length(sim{1,1}{1,2}(tstart:end,4)),'LineWidth',1);
hold on
plot(sim_old{1,1}{1,2}(tstart:end,1),sim_old{1,1}{1,2}(tstart:end,4),'LineWidth',1);
hold off
set(gca,'FontSize',14);
xlabel('time (h)','FontSize',16);
ylabel('[FA-FrmR] (nM)','FontSize',16);
legend('e7 [FA]_l_o_w','e7 [FA]_h_i_g_h','FM [FA]_l_o_w','FM [FA]_h_i_g_h');
title('FA-FrmR','FontSize',18);
subplot(2,2,4);
plot(sim{1,1}{1,1}(tstart:end,1),sim{1,1}{1,1}(tstart:end,5),'-o','MarkerIndices',1:stp_indices:length(sim{1,1}{1,1}(tstart:end,5)),'LineWidth',1);
hold on
plot(sim{1,1}{1,2}(tstart:end,1),sim{1,1}{1,2}(tstart:end,5),'LineWidth',1);
hold on
plot(sim_old{1,1}{1,1}(tstart:end,1),sim_old{1,1}{1,1}(tstart:end,5),'-o','MarkerIndices',1:stp_indices:length(sim{1,1}{1,2}(tstart:end,5)),'LineWidth',1);
hold on
plot(sim_old{1,1}{1,2}(tstart:end,1),sim_old{1,1}{1,2}(tstart:end,5),'LineWidth',1);
hold off
set(gca,'FontSize',14);
xlabel('time (h)','FontSize',16);
ylabel('[mLacI] (nM)','FontSize',16);
legend('e7 [FA]_l_o_w','e7 [FA]_h_i_g_h','FM [FA]_l_o_w','FM [FA]_h_i_g_h');
title('mLacI','FontSize',18);
sgtitle('Time series extension 7 - p-set with highest hok/sok ratio at t=48h','FontSize',22);

figure(12);
subplot(2,2,1);
plot(sim{1,1}{1,1}(tstart:end,1),sim{1,1}{1,1}(tstart:end,6),'-o','MarkerIndices',1:stp_indices:length(sim{1,1}{1,1}(tstart:end,6)),'LineWidth',1);
hold on
plot(sim{1,1}{1,2}(tstart:end,1),sim{1,1}{1,2}(tstart:end,6),'LineWidth',1);
hold on
plot(sim_old{1,1}{1,1}(tstart:end,1),sim_old{1,1}{1,1}(tstart:end,6),'-o','MarkerIndices',1:stp_indices:length(sim{1,1}{1,2}(tstart:end,6)),'LineWidth',1);
hold on
plot(sim_old{1,1}{1,2}(tstart:end,1),sim_old{1,1}{1,2}(tstart:end,6),'LineWidth',1);
hold off
set(gca,'FontSize',14);
xlabel('time (h)','FontSize',16);
ylabel('[LacI] (nM)','FontSize',16);
legend('e7 [FA]_l_o_w','e7 [FA]_h_i_g_h','FM [FA]_l_o_w','FM [FA]_h_i_g_h');
title('LacI','FontSize',18);
subplot(2,2,2);
plot(sim{1,1}{1,1}(tstart:end,1),sim{1,1}{1,1}(tstart:end,7),'-o','MarkerIndices',1:stp_indices:length(sim{1,1}{1,1}(tstart:end,7)),'LineWidth',1);
hold on
plot(sim{1,1}{1,2}(tstart:end,1),sim{1,1}{1,2}(tstart:end,7),'LineWidth',1);
hold on
plot(sim_old{1,1}{1,1}(tstart:end,1),sim_old{1,1}{1,1}(tstart:end,7),'-o','MarkerIndices',1:stp_indices:length(sim{1,1}{1,2}(tstart:end,7)),'LineWidth',1);
hold on
plot(sim_old{1,1}{1,2}(tstart:end,1),sim_old{1,1}{1,2}(tstart:end,7),'LineWidth',1);
hold off
set(gca,'yscale','log','FontSize',14);
xlabel('time (h)','FontSize',16);
ylabel('[msok] (nM)','FontSize',16);
legend('e7 [FA]_l_o_w','e7 [FA]_h_i_g_h','FM [FA]_l_o_w','FM [FA]_h_i_g_h');
title('msok','FontSize',18);
subplot(2,2,3);
plot(sim{1,1}{1,1}(tstart:end,1),sim{1,1}{1,1}(tstart:end,8),'-o','MarkerIndices',1:stp_indices:length(sim{1,1}{1,1}(tstart:end,8)),'LineWidth',1);
hold on
plot(sim{1,1}{1,2}(tstart:end,1),sim{1,1}{1,2}(tstart:end,8),'LineWidth',1);
hold on
plot(sim_old{1,1}{1,1}(tstart:end,1),sim_old{1,1}{1,1}(tstart:end,8),'-o','MarkerIndices',1:stp_indices:length(sim{1,1}{1,2}(tstart:end,8)),'LineWidth',1);
hold on
plot(sim_old{1,1}{1,2}(tstart:end,1),sim_old{1,1}{1,2}(tstart:end,8),'LineWidth',1);
hold off
set(gca,'FontSize',14);
xlabel('time (h)','FontSize',16);
ylabel('[hok-sok] (nM)','FontSize',16);
legend('e7 [FA]_l_o_w','e7 [FA]_h_i_g_h','FM [FA]_l_o_w','FM [FA]_h_i_g_h');
title('hok-sok','FontSize',18);
subplot(2,2,4);
plot(sim{1,1}{1,1}(tstart:end,1),sim{1,1}{1,1}(tstart:end,9),'-o','MarkerIndices',1:stp_indices:length(sim{1,1}{1,1}(tstart:end,9)),'LineWidth',1);
hold on
plot(sim{1,1}{1,2}(tstart:end,1),sim{1,1}{1,2}(tstart:end,9),'LineWidth',1);
hold on
plot(sim_old{1,1}{1,1}(tstart:end,1),sim_old{1,1}{1,1}(tstart:end,9),'-o','MarkerIndices',1:stp_indices:length(sim{1,1}{1,2}(tstart:end,9)),'LineWidth',1);
hold on
plot(sim_old{1,1}{1,2}(tstart:end,1),sim_old{1,1}{1,2}(tstart:end,9),'LineWidth',1);
hold off
set(gca,'FontSize',14);
xlabel('time (h)','FontSize',16);
ylabel('[mhok] (nM)','FontSize',16);
%set(gca,'yscale','log');
legend('e7 [FA]_l_o_w','e7 [FA]_h_i_g_h','FM [FA]_l_o_w','FM [FA]_h_i_g_h');
title('mhok','FontSize',18);
sgtitle('Time series extension 7 - p-set with highest hok/sok ratio at t=48h','FontSize',22);

figure(13);
subplot(2,2,1);
plot(sim{1,1}{1,1}(tstart:end,1),sim{1,1}{1,1}(tstart:end,10),'-o','MarkerIndices',1:stp_indices:length(sim{1,1}{1,1}(tstart:end,10)),'LineWidth',1);
hold on
plot(sim{1,1}{1,2}(tstart:end,1),sim{1,1}{1,2}(tstart:end,10),'LineWidth',1);
hold off
set(gca,'FontSize',14);
xlabel('time (h)','FontSize',16);
ylabel('[mX] (nM)','FontSize',16);
legend('e7 [FA]_l_o_w','e7 [FA]_h_i_g_h');
title('mX','FontSize',18);
subplot(2,2,2);
plot(sim{1,1}{1,1}(tstart:end,1),sim{1,1}{1,1}(tstart:end,11),'-o','MarkerIndices',1:stp_indices:length(sim{1,1}{1,1}(tstart:end,11)),'LineWidth',1);
hold on
plot(sim{1,1}{1,2}(tstart:end,1),sim{1,1}{1,2}(tstart:end,11),'LineWidth',1);
hold off
set(gca,'FontSize',14);
xlabel('time (h)','FontSize',16);
ylabel('[X] (nM)','FontSize',16);
legend('e7 [FA]_l_o_w','e7 [FA]_h_i_g_h');
title('X','FontSize',18);
subplot(2,2,3);
plot(B_low_F_TSnew(:,2),B_low_F_TSnew(:,13),'-o','MarkerIndices',1:stp_indices:length(B_low_F_TSnew(:,2)),'LineWidth',1);
hold on
plot(B_high_F_TSnew(:,2),B_high_F_TSnew(:,13),'LineWidth',1);
hold on
plot(B_low_F_TSold(:,2),B_low_F_TSold(:,11),'-o','MarkerIndices',1:stp_indices:length(B_high_F_TSnew(:,2)),'LineWidth',1);
hold on
plot(B_high_F_TSold(:,2),B_high_F_TSold(:,11),'LineWidth',1);
hold off
set(gca,'yscale','log','FontSize',14);
xlabel('time (h)','FontSize',16);
ylabel('hok/sok ratio','FontSize',16);
legend('e7 [FA]_l_o_w','e7 [FA]_h_i_g_h','FM [FA]_l_o_w','FM [FA]_h_i_g_h');
title('hok/sok ratio','FontSize',18);
subplot(2,2,4);
plot(B_low_F_TSnew(:,2),B_low_F_TSnew(:,15),'-o','MarkerIndices',1:stp_indices:length(B_low_F_TSnew(:,2)),'LineWidth',1);
hold on
plot(B_low_F_TSold(:,2),B_low_F_TSold(:,13),'LineWidth',1);
hold off
set(gca,'yscale','log','FontSize',14);
xlabel('time (h)','FontSize',16);
ylabel('sensitivity','FontSize',16);
legend('e7','FM');
title('sensitivity','FontSize',18);
%sgtitle('extension 7 - p-set 3082 (outlier hoksok ratio)','FontSize',22);
%sgtitle('extension 7 - p-set 1561 (highest max sensitivity)','FontSize',22);
sgtitle('Time series extension 7 - p-set with highest hok/sok ratio at t=48h','FontSize',22);
%sgtitle('Time series extension 7 - p-set with hok/sok ratio 6000 at t=48h','FontSize',22);

figure(14);
subplot(2,2,1);
plot(B_low_F_TSnew(:,2),B_low_F_TSnew(:,13),'-o','MarkerIndices',1:stp_indices:length(B_low_F_TSnew(:,2)),'LineWidth',1);
hold on
plot(B_high_F_TSnew(:,2),B_high_F_TSnew(:,13),'LineWidth',1);
hold on
plot(B_low_F_TSold(:,2),B_low_F_TSold(:,11),'-o','MarkerIndices',1:stp_indices:length(B_high_F_TSnew(:,2)),'LineWidth',1);
hold on
plot(B_high_F_TSold(:,2),B_high_F_TSold(:,11),'LineWidth',1);
hold off
set(gca,'yscale','log','FontSize',14);
xlabel('time (h)','FontSize',16);
ylabel('hok/sok ratio','FontSize',16);
legend('e7 [FA]_l_o_w','e7 [FA]_h_i_g_h','FM [FA]_l_o_w','FM [FA]_h_i_g_h');
title('hok/sok ratio','FontSize',18);
subplot(2,2,2);
plot(B_low_F_TSnew(:,2),B_low_F_TSnew(:,15),'-o','MarkerIndices',1:stp_indices:length(B_low_F_TSnew(:,2)),'LineWidth',1);
hold on
plot(B_low_F_TSold(:,2),B_low_F_TSold(:,13),'LineWidth',1);
hold off
set(gca,'yscale','log','FontSize',14);
xlabel('time (h)','FontSize',16);
ylabel('sensitivity','FontSize',16);
legend('e7','FM');
title('sensitivity','FontSize',18);
sgtitle('Time series extension 7 - p-set with highest hok/sok ratio at t=48h','FontSize',22);