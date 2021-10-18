%% file used to analyse the output of the simulations of the model extension 1
% contains:
% 1. create a matrix with end concentrations for the parameter sets that have the median p-set of the full model as first 22 parameter values
% 2. create a matrix with end concentrations for the parameter sets that have the sensitivity outlier p-set of the full model as first 22 parameter values
% 3. find the parameter sets with a the max sensitivity + corresponding timepoint
% 4. Time series plot e1

% input files: 
%     e1_output_m_1.mat
%     e1_output_o_1.mat
%     e1_parameters_o.mat
%     median_low_hoksokratio2.mat
% output files:
%     e1_max_sensitivity_analysis.mat
%% 1. create a matrix with end concentrations for the parameter sets that have the median p-set of the full model as first 22 parameter values
% load output file median parameter set
load('e1_output_m_1.mat','sim')
% creating a matrix with end concentrations for the simulations of each
% p-set.
    low_F_e = zeros(10000, 12); 
    high_F_e= zeros(10000, 12);
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
    sortB = sortrows(B_high_F_e,16);
    C_high_F_e = [B_high_F_e A_low_F_e(:,12)./A_high_F_e(:,12) A_low_F_e(:,7)./A_high_F_e(:,7)];

%save('e2_k.mat','B_high_F_e','B_low_F_e','e2_k');
    
%% 2. create a matrix with end concentrations for the parameter sets that have the sensitivity outlier p-set of the full model as first 22 parameter values
load('e1_output_o_1.mat','sim')
% creating a matrix with end concentrations for the simulations of each
% p-set.
    low_F_e = zeros(10000, 12); 
    high_F_e= zeros(10000, 12);
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
    B_high_F_e = [A_high_F_e A_low_F_e(:,14)-A_high_F_e(:,14) A_low_F_e(:,14)./A_high_F_e(:,14)];
    B_low_F_e = [A_low_F_e A_low_F_e(:,14)-A_high_F_e(:,14) A_low_F_e(:,14)./A_high_F_e(:,14)];
    sortB = sortrows(B_high_F_e,16);
    C_high_F_e = [B_high_F_e A_low_F_e(:,12)./A_high_F_e(:,12) A_low_F_e(:,7)./A_high_F_e(:,7)];
%% 3. find the parameter sets with a the max sensitivity + corresponding timepoint
load('e1_output_o_1.mat','sim')

% add hoksok ratio in the sim array (13) and the sensitivity (14)
for i=1:length(sim)
     sim{i}{1}(:,13) = (sim{i}{1}(:,9)+sim{i}{1}(:,8))./(sim{i}{1}(:,8)+sim{i}{1}(:,7));
     sim{i}{2}(:,13) = (sim{i}{2}(:,9)+sim{i}{2}(:,8))./(sim{i}{2}(:,8)+sim{i}{2}(:,7));
     sim{i}{1}(:,14) = (sim{i}{1}(:,13))./(sim{i}{2}(:,13));
     sim{i}{2}(:,14) = (sim{i}{1}(:,13))./(sim{i}{2}(:,13));
end
% find the max sensitivities for each simulation and make a vector with the corresponding timepoints 
    low_F_e = zeros(10000, 12); 
    high_F_e= zeros(10000, 12);
    for i=1:length(sim)
         sim_max(i,1) = max(sim{i}{1}(:,14));
         t(i,1) = find(sim{i}{1}(:,14)==sim_max(i,1))./100;
    end
% sort the max sensitivities, to get the most sensitive p-set
    sim_max_n = [[1:10000]' t sim_max];
    sort_max = sortrows(sim_max_n,3);
save('e1_max_sensitivity_analysis.mat','sim','sort_max')

%% 4. Time series plot e1
load('median_low_hoksokratio2','params_GhighF_ratio');
load('e1_parameters_o','p_sets')
N = 1;
q = 1;

% select parameter set to simulate 
p_set = 3893; %highest sensitivity at t=48h
%p_set = 4331; %highest max sensitivity
p_set_m = 1056; %original biosafety circuit

% make simulation of given paremeter set for extension 1 (new) and original
% biosafety circuit (old)
[sim] = optimise_model_e1(p_sets(p_set,:),N,q);
[sim_old] = optimise_model(params_GhighF_ratio(p_set_m,:),N,q);
for i=1:length(sim)
     low_F(i,1:length(sim{i}{1}(end,:))) = sim{i}{1}(end,:);
     high_F(i,1:length(sim{i}{2}(end,:))) = sim{i}{2}(end,:);            
end
A_high_F = [[1:1]' high_F (high_F(:,9)+high_F(:,8))./(high_F(:,8)+high_F(:,7))];
A_low_F = [[1:1]' low_F (low_F(:,9)+low_F(:,8))./(low_F(:,8)+low_F(:,7))];
B_high_F = [A_high_F A_low_F(:,11)-A_high_F(:,11) A_low_F(:,11)./A_high_F(:,11)];
B_low_F = [A_low_F A_low_F(:,11)-A_high_F(:,11) A_low_F(:,11)./A_high_F(:,11)];

% Matrix A contains a time vector (first column) + time series data + hoksok ratio (last column).
% Matrix B contains matrix A (column 1-13) + hoksok ratio [F]low - hoksok
% ratio [F]high + hoksok ratio [F]low/hoksok ratio [F]high.
    A_high_F_TSnew = [[2:length(sim{1}{2}(:,1))]' sim{1}{2}(2:end,:) (sim{1}{2}(2:end,9)+sim{1}{2}(2:end,8))./(sim{1}{2}(2:end,8)+sim{1}{2}(2:end,7))];
    A_low_F_TSnew = [[2:length(sim{1}{2}(:,1))]' sim{1}{2}(2:end,:) (sim{1}{1}(2:end,9)+sim{1}{1}(2:end,8))./(sim{1}{1}(2:end,8)+sim{1}{1}(2:end,7))];
    B_high_F_TSnew = [A_high_F_TSnew A_low_F_TSnew(:,14)-A_high_F_TSnew(:,14) A_low_F_TSnew(:,14)./A_high_F_TSnew(:,14)];
    B_low_F_TSnew = [A_low_F_TSnew A_low_F_TSnew(:,14)-A_high_F_TSnew(:,14) A_low_F_TSnew(:,14)./A_high_F_TSnew(:,14)];

    A_high_F_TSold = [[2:length(sim_old{1}{2}(:,1))]' sim_old{1}{2}(2:end,:) (sim_old{1}{2}(2:end,9)+sim_old{1}{2}(2:end,8))./(sim_old{1}{2}(2:end,8)+sim_old{1}{2}(2:end,7))];
    A_low_F_TSold = [[2:length(sim_old{1}{2}(:,1))]' sim_old{1}{2}(2:end,:) (sim_old{1}{1}(2:end,9)+sim_old{1}{1}(2:end,8))./(sim_old{1}{1}(2:end,8)+sim_old{1}{1}(2:end,7))];
    B_high_F_TSold = [A_high_F_TSold A_low_F_TSold(:,11)-A_high_F_TSold(:,11) A_low_F_TSold(:,11)./A_high_F_TSold(:,11)];
    B_low_F_TSold = [A_low_F_TSold A_low_F_TSold(:,11)-A_high_F_TSold(:,11) A_low_F_TSold(:,11)./A_high_F_TSold(:,11)];

%time series plots
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
legend('e1 [FA]_l_o_w','e1 [FA]_h_i_g_h','FM [FA]_l_o_w','FM [FA]_h_i_g_h');
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
legend('e1 [FA]_l_o_w','e1 [FA]_h_i_g_h','FM [FA]_l_o_w','FM [FA]_h_i_g_h');
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
legend('e1 [FA]_l_o_w','e1 [FA]_h_i_g_h','FM [FA]_l_o_w','FM [FA]_h_i_g_h');
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
legend('e1 [FA]_l_o_w','e1 [FA]_h_i_g_h','FM [FA]_l_o_w','FM [FA]_h_i_g_h');
title('mLacI','FontSize',18);
%sgtitle('Time series extension 1 - p-set with highest max sensitivity','FontSize',22);
sgtitle('Time series extension 1 - p-set with highest sensitivity at t=48h','FontSize',22);

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
legend('e1 [FA]_l_o_w','e1 [FA]_h_i_g_h','FM [FA]_l_o_w','FM [FA]_h_i_g_h');
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
set(gca,'FontSize',14);
xlabel('time (h)','FontSize',16);
ylabel('[msok] (nM)','FontSize',16);
legend('e1 [FA]_l_o_w','e1 [FA]_h_i_g_h','FM [FA]_l_o_w','FM [FA]_h_i_g_h');
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
legend('e1 [FA]_l_o_w','e1 [FA]_h_i_g_h','FM [FA]_l_o_w','FM [FA]_h_i_g_h');
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
legend('e1 [FA]_l_o_w','e1 [FA]_h_i_g_h','FM [FA]_l_o_w','FM [FA]_h_i_g_h');
title('mhok','FontSize',18);
%sgtitle('extension 1 - p-set 4331 (outlier max sensitivity)','FontSize',22);
%sgtitle('Time series extension 1 - p-set with highest max sensitivity','FontSize',22);
sgtitle('Time series extension 1 - p-set with highest sensitivity at t=48h','FontSize',22);

figure(13);
subplot(2,2,1);
plot(sim{1,1}{1,1}(tstart:end,1),sim{1,1}{1,1}(tstart:end,10),'-o','MarkerIndices',1:stp_indices:length(sim{1,1}{1,1}(tstart:end,10)),'LineWidth',1);
hold on
plot(sim{1,1}{1,2}(tstart:end,1),sim{1,1}{1,2}(tstart:end,10),'LineWidth',1);
hold off
set(gca,'FontSize',14);
xlabel('time (h)','FontSize',16);
ylabel('[mC] (nM)','FontSize',16);
legend('e1 [FA]_l_o_w','e1 [FA]_h_i_g_h');
title('mC','FontSize',18);
subplot(2,2,2);
plot(sim{1,1}{1,1}(tstart:end,1),sim{1,1}{1,1}(tstart:end,11),'-o','MarkerIndices',1:stp_indices:length(sim{1,1}{1,1}(tstart:end,11)),'LineWidth',1);
hold on
plot(sim{1,1}{1,2}(tstart:end,1),sim{1,1}{1,2}(tstart:end,11),'LineWidth',1);
hold off
set(gca,'FontSize',14);
xlabel('time (h)','FontSize',16);
ylabel('[C] (nM)','FontSize',16);
legend('e1 [FA]_l_o_w','e1 [FA]_h_i_g_h');
title('C','FontSize',18);
subplot(2,2,3);
plot(sim{1,1}{1,1}(tstart:end,1),sim{1,1}{1,1}(tstart:end,12),'-o','MarkerIndices',1:stp_indices:length(sim{1,1}{1,1}(tstart:end,12)),'LineWidth',1);
hold on
plot(sim{1,1}{1,2}(tstart:end,1),sim{1,1}{1,2}(tstart:end,12),'LineWidth',1);
hold off
set(gca,'FontSize',14);
xlabel('time (h)','FontSize',16);
ylabel('[C-FrmR] (nM)','FontSize',16);
legend('e1 [FA]_l_o_w','e1 [FA]_h_i_g_h');
title('C-FrmR','FontSize',18);
subplot(2,2,4);
plot(B_low_F_TSnew(:,2),B_low_F_TSnew(:,14),'-o','MarkerIndices',1:stp_indices:length(B_low_F_TSnew(:,2)),'LineWidth',1);
hold on
plot(B_high_F_TSnew(:,2),B_high_F_TSnew(:,14),'LineWidth',1);
hold on
plot(B_low_F_TSold(:,2),B_low_F_TSold(:,11),'-o','MarkerIndices',1:stp_indices:length(B_high_F_TSnew(:,2)),'LineWidth',1);
hold on
plot(B_high_F_TSold(:,2),B_high_F_TSold(:,11),'LineWidth',1);
hold off
set(gca,'FontSize',14);
xlabel('time (h)','FontSize',16);
ylabel('hok/sok ratio','FontSize',16);
%set(gca,'yscale','log');
legend('e1 [FA]_l_o_w','e1 [FA]_h_i_g_h','FM [FA]_l_o_w','FM [FA]_h_i_g_h');
title('hok/sok ratio','FontSize',18);
%sgtitle('Time series extension 1 - p-set with highest max sensitivity','FontSize',22);
sgtitle('Time series extension 1 - p-set with highest sensitivity at t=48h','FontSize',22);

figure(14);
subplot(2,2,1);
plot(B_low_F_TSnew(:,2),B_low_F_TSnew(:,16),'-o','MarkerIndices',1:stp_indices:length(B_low_F_TSnew(:,2)),'LineWidth',1);
hold on
plot(B_low_F_TSold(:,2),B_low_F_TSold(:,13),'LineWidth',1);
hold off
set(gca,'FontSize',14);
xlabel('time (h)','FontSize',16);
ylabel('sensitivity','FontSize',16);
legend('e1','FM');
title('sensitivity','FontSize',18);
%sgtitle('Time series extension 1 - p-set with highest max sensitivity','FontSize',22);
sgtitle('Time series extension 1 - p-set with highest sensitivity at t=48h','FontSize',22);

figure(15);
subplot(2,2,1);
plot(B_low_F_TSnew(:,2),B_low_F_TSnew(:,14),'-o','MarkerIndices',1:stp_indices:length(B_low_F_TSnew(:,2)),'LineWidth',1);
hold on
plot(B_high_F_TSnew(:,2),B_high_F_TSnew(:,14),'LineWidth',1);
hold on
plot(B_low_F_TSold(:,2),B_low_F_TSold(:,11),'-o','MarkerIndices',1:stp_indices:length(B_high_F_TSnew(:,2)),'LineWidth',1);
hold on
plot(B_high_F_TSold(:,2),B_high_F_TSold(:,11),'LineWidth',1);
hold off
set(gca,'FontSize',14);
xlabel('time (h)','FontSize',16);
ylabel('hok/sok ratio','FontSize',16);
%set(gca,'yscale','log');
legend('e1 [FA]_l_o_w','e1 [FA]_h_i_g_h','FM [FA]_l_o_w','FM [FA]_h_i_g_h');
title('hok/sok ratio','FontSize',18);
subplot(2,2,2);
plot(B_low_F_TSnew(:,2),B_low_F_TSnew(:,16),'-o','MarkerIndices',1:stp_indices:length(B_low_F_TSnew(:,2)),'LineWidth',1);
hold on
plot(B_low_F_TSold(:,2),B_low_F_TSold(:,13),'LineWidth',1);
hold off
set(gca,'FontSize',14);
xlabel('time (h)','FontSize',16);
ylabel('sensitivity','FontSize',16);
legend('e1','FM');
title('sensitivity','FontSize',18);
%sgtitle('Time series extension 1 - p-set with highest max sensitivity','FontSize',22);
sgtitle('Time series extension 1 - p-set with highest sensitivity at t=48h','FontSize',22);
