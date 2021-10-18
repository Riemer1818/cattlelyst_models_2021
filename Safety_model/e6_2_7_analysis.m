%% file used to analyse the output of the simulations of the model extension 6.1
% contains:
% 1. Time series plot e6.1+7
% 2. find the max hoksok ratio + corresponding timepoint and max sensitivity + timepoint
% 3. Sensitivity analysis
% 4. Plot sensitivity analysis
% 5. run first cycle with initial conditions 0
% 6. run multiple cycles with initial concentration = half end concentration previous cycle 

% input files: 
%     e6_2_7_parameters_o.mat
%     median_low_hoksokratio2.mat
%     e6_2_7_conc_d24h_2c
%     e6_2_7_conc_d24h_2c_1
%     e6_2_7_conc_d24h_0
%     e6_2_7_output.mat
% output files: 
%     e6_2_7_conc.mat
%     e6_2_7_conc_d24h_2c
%     e6_2_7_conc_d24h_2c_1
%     e6_2_7_conc_d24h_0
%     e6_2_7_output.mat
%     e6_2_7_SA_times1.1.mat
%     e6_2_7_SA_div1.1.mat
%     e6_2_7_SA_times10.mat
%     e6_2_7_SA_div10.mat
%% 1. Time series plot e6.1+7
load('median_low_hoksokratio2','params_GhighF_ratio');
load('e6_2_7_parameters_o','p_sets')
N = 1;
q = 1;

% make simulation of given paremeter set for extension 6.1+7 (new)
%[sim] = optimise_model_e6_2_7(p_sets(:,:),N,q);
load('e6_2_7_conc_d24h_2c_2','sim_all');  
sim = sim_all;

% select parameter set of original biosafety circuit to simulate
p_set_m = 1056; 

% make simulation of original biosafety circuit (old)
[sim_old] = optimise_model(params_GhighF_ratio(p_set_m,:),N,q);

% Matrix A contains a time vector (first column) + time series data + hoksok ratio (last column).
% Matrix B contains matrix A (column 1-13) + hoksok ratio [F]low - hoksok
% ratio [F]high + hoksok ratio [F]low/hoksok ratio [F]high.
    A_high_F_TSnew = [[2:length(sim{1}{2}(:,1))]' sim{1}{2}(2:end,:) (sim{1}{2}(2:end,9)+sim{1}{2}(2:end,8))./(sim{1}{2}(2:end,8)+sim{1}{2}(2:end,7))];
    A_low_F_TSnew = [[2:length(sim{1}{1}(:,1))]' sim{1}{1}(2:end,:) (sim{1}{1}(2:end,9)+sim{1}{1}(2:end,8))./(sim{1}{1}(2:end,8)+sim{1}{1}(2:end,7))];
    B_high_F_TSnew = [A_high_F_TSnew A_low_F_TSnew(:,15)-A_high_F_TSnew(:,15) A_low_F_TSnew(:,15)./A_high_F_TSnew(:,15)];
    B_low_F_TSnew = [A_low_F_TSnew A_low_F_TSnew(:,15)-A_high_F_TSnew(:,15) A_low_F_TSnew(:,15)./A_high_F_TSnew(:,15)];

    A_high_F_TSold = [[2:length(sim_old{1}{2}(:,1))]' sim_old{1}{2}(2:end,:) (sim_old{1}{2}(2:end,9)+sim_old{1}{2}(2:end,8))./(sim_old{1}{2}(2:end,8)+sim_old{1}{2}(2:end,7))];
    A_low_F_TSold = [[2:length(sim_old{1}{1}(:,1))]' sim_old{1}{1}(2:end,:) (sim_old{1}{1}(2:end,9)+sim_old{1}{1}(2:end,8))./(sim_old{1}{1}(2:end,8)+sim_old{1}{1}(2:end,7))];
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
%legend('e6.1+7 [FA]_l_o_w','e6.1+7 [FA]_h_i_g_h[FA]_l_o_w','FM [FA]_l_o_w','FM [FA]_h_i_g_h');
legend('e6.1+7 [FA]_l_o_w','e6.1+7 [FA]_h_i_g_h','FM [FA]_l_o_w','FM [FA]_h_i_g_h');
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
set(gca,'yscale','log','FontSize',14);
xlabel('time (h)','FontSize',16);
ylabel('[FrmR] (nM)','FontSize',16);
%legend('e6.1+7 [FA]_l_o_w','e6.1+7 [FA]_h_i_g_h[FA]_l_o_w','FM [FA]_l_o_w','FM [FA]_h_i_g_h');
legend('e6.1+7 [FA]_l_o_w','e6.1+7 [FA]_h_i_g_h','FM [FA]_l_o_w','FM [FA]_h_i_g_h');
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
set(gca,'yscale','log','FontSize',14);
xlabel('time (h)','FontSize',16);
ylabel('[FA-FrmR] (nM)','FontSize',16);
%legend('e6.1+7 [FA]_l_o_w','e6.1+7 [FA]_h_i_g_h[FA]_l_o_w','FM [FA]_l_o_w','FM [FA]_h_i_g_h');
legend('e6.1+7 [FA]_l_o_w','e6.1+7 [FA]_h_i_g_h','FM [FA]_l_o_w','FM [FA]_h_i_g_h');
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
%legend('e6.1+7 [FA]_l_o_w','e6.1+7 [FA]_h_i_g_h[FA]_l_o_w','FM [FA]_l_o_w','FM [FA]_h_i_g_h');
legend('e6.1+7 [FA]_l_o_w','e6.1+7 [FA]_h_i_g_h','FM [FA]_l_o_w','FM [FA]_h_i_g_h');
title('mLacI','FontSize',18);
%sgtitle('extension 6.2+7 - p-set 2439 +2425 (outlier sensitivity + hoksok ratio 6000)','FontSize',22);
%sgtitle('extension 6.2+7 - p-set 2439 +3082 (outlier sensitivity + max hoksok ratio)','FontSize',22);
%sgtitle('Time series extension 6.1+7 - [F] = 9uM at t=0.21h','FontSize',22);
sgtitle('Time series extension 6.1+7','FontSize',22);

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
%legend('e6.1+7 [FA]_l_o_w','e6.1+7 [FA]_h_i_g_h[FA]_l_o_w','FM [FA]_l_o_w','FM [FA]_h_i_g_h');
legend('e6.1+7 [FA]_l_o_w','e6.1+7 [FA]_h_i_g_h','FM [FA]_l_o_w','FM [FA]_h_i_g_h');
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
%legend('e6.1+7 [FA]_l_o_w','e6.1+7 [FA]_h_i_g_h[FA]_l_o_w','FM [FA]_l_o_w','FM [FA]_h_i_g_h');
legend('e6.1+7 [FA]_l_o_w','e6.1+7 [FA]_h_i_g_h','FM [FA]_l_o_w','FM [FA]_h_i_g_h');
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
%legend('e6.1+7 [FA]_l_o_w','e6.1+7 [FA]_h_i_g_h[FA]_l_o_w','FM [FA]_l_o_w','FM [FA]_h_i_g_h');
legend('e6.1+7 [FA]_l_o_w','e6.1+7 [FA]_h_i_g_h','FM [FA]_l_o_w','FM [FA]_h_i_g_h');
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
set(gca,'yscale','log','FontSize',14);
xlabel('time (h)','FontSize',16);
ylabel('[mhok] (nM)','FontSize',16);
%legend('e6.1+7 [FA]_l_o_w','e6.1+7 [FA]_h_i_g_h[FA]_l_o_w','FM [FA]_l_o_w','FM [FA]_h_i_g_h');
legend('e6.1+7 [FA]_l_o_w','e6.1+7 [FA]_h_i_g_h','FM [FA]_l_o_w','FM [FA]_h_i_g_h');
title('mhok','FontSize',18);
%sgtitle('extension 6.2+7 - p-set 2439 +2425 (outlier sensitivity + hoksok ratio 6000)','FontSize',22);
%sgtitle('extension 6.2+7 - p-set 2439 +3082 (outlier sensitivity + max hoksok ratio)','FontSize',22);
%sgtitle('Time series extension 6.1+7 - [F] = 9uM at t=0.21h','FontSize',22);
sgtitle('Time series extension 6.1+7','FontSize',22);

figure(13);
subplot(2,2,1);
plot(sim{1,1}{1,1}(tstart:end,1),sim{1,1}{1,1}(tstart:end,10),'-o','MarkerIndices',1:stp_indices:length(sim{1,1}{1,1}(tstart:end,10)),'LineWidth',1);
hold on
plot(sim{1,1}{1,2}(tstart:end,1),sim{1,1}{1,2}(tstart:end,10),'LineWidth',1);
hold off
set(gca,'FontSize',14);
xlabel('time (h)','FontSize',16);
ylabel('[mX] (nM)','FontSize',16);
%legend('e6.1+7 [FA]_l_o_w','e6.1+7 [FA]_h_i_g_h[FA]_l_o_w');
legend('e6.1+7 [FA]_l_o_w','e6.1+7 [FA]_h_i_g_h');
title('mX','FontSize',18);
subplot(2,2,2);
plot(sim{1,1}{1,1}(tstart:end,1),sim{1,1}{1,1}(tstart:end,11),'-o','MarkerIndices',1:stp_indices:length(sim{1,1}{1,1}(tstart:end,11)),'LineWidth',1);
hold on
plot(sim{1,1}{1,2}(tstart:end,1),sim{1,1}{1,2}(tstart:end,11),'LineWidth',1);
hold off
set(gca,'FontSize',14);
xlabel('time (h)','FontSize',16);
ylabel('[X] (nM)','FontSize',16);
%legend('e6.1+7 [FA]_l_o_w','e6.1+7 [FA]_h_i_g_h[FA]_l_o_w');
legend('e6.1+7 [FA]_l_o_w','e6.1+7 [FA]_h_i_g_h');
title('X','FontSize',18);
subplot(2,2,3);
plot(sim{1,1}{1,1}(tstart:end,1),sim{1,1}{1,1}(tstart:end,12),'-o','MarkerIndices',1:stp_indices:length(sim{1,1}{1,1}(tstart:end,11)),'LineWidth',1);
hold on
plot(sim{1,1}{1,2}(tstart:end,1),sim{1,1}{1,2}(tstart:end,12),'LineWidth',1);
hold off
set(gca,'FontSize',14);
xlabel('time (h)','FontSize',16);
ylabel('[mY] (nM)','FontSize',16);
%legend('e6.1+7 [FA]_l_o_w','e6.1+7 [FA]_h_i_g_h[FA]_l_o_w');
legend('e6.1+7 [FA]_l_o_w','e6.1+7 [FA]_h_i_g_h');
title('mY','FontSize',18);
subplot(2,2,4);
plot(sim{1,1}{1,1}(tstart:end,1),sim{1,1}{1,1}(tstart:end,13),'-o','MarkerIndices',1:stp_indices:length(sim{1,1}{1,1}(tstart:end,11)),'LineWidth',1);
hold on
plot(sim{1,1}{1,2}(tstart:end,1),sim{1,1}{1,2}(tstart:end,13),'LineWidth',1);
hold off
set(gca,'FontSize',14);
xlabel('time (h)','FontSize',16);
ylabel('[Y] (nM)','FontSize',16);
%legend('e6.1+7 [FA]_l_o_w','e6.1+7 [FA]_h_i_g_h[FA]_l_o_w');
legend('e6.1+7 [FA]_l_o_w','e6.1+7 [FA]_h_i_g_h');
title('Y','FontSize',18);
%sgtitle('extension 6.2+7 - p-set 2439 +2425 (outlier sensitivity + hoksok ratio 6000)','FontSize',22);
%sgtitle('extension 6.2+7 - p-set 2439 +3082 (outlier sensitivity + max hoksok ratio)','FontSize',22);
%sgtitle('Time series extension 6.1+7 - [F] = 9uM at t=0.21h','FontSize',22);
sgtitle('Time series extension 6.1+7','FontSize',22);

figure(14);
subplot(2,2,1);
plot(B_low_F_TSnew(:,2),B_low_F_TSnew(:,15),'-o','MarkerIndices',1:stp_indices:length(B_low_F_TSnew(:,2)),'LineWidth',1);
hold on
plot(B_high_F_TSnew(:,2),B_high_F_TSnew(:,15),'LineWidth',1);
hold on
plot(B_low_F_TSold(:,2),B_low_F_TSold(:,11),'-o','MarkerIndices',1:stp_indices:length(B_high_F_TSnew(:,2)),'LineWidth',1);
hold on
plot(B_high_F_TSold(:,2),B_high_F_TSold(:,11),'LineWidth',1);
hold off
set(gca,'yscale','log','FontSize',14);
xlabel('time (h)','FontSize',16);
ylabel('hok/sok ratio','FontSize',16);
%legend('e6.1+7 [FA]_l_o_w','e6.1+7 [FA]_h_i_g_h[FA]_l_o_w','FM [FA]_l_o_w','FM [FA]_h_i_g_h');
legend('e6.1+7 [FA]_l_o_w','e6.1+7 [FA]_h_i_g_h','FM [FA]_l_o_w','FM [FA]_h_i_g_h');
title('hok/sok ratio','FontSize',18);
subplot(2,2,2);
plot(B_low_F_TSnew(:,2),B_low_F_TSnew(:,17),'-o','MarkerIndices',1:stp_indices:length(B_low_F_TSnew(:,2)),'LineWidth',1);
hold on
plot(B_low_F_TSold(:,2),B_low_F_TSold(:,13),'LineWidth',1);
hold off
set(gca,'yscale','log','FontSize',14);
xlabel('time (h)','FontSize',16);
ylabel('sensitivity','FontSize',16);
legend('e6.1+7','FM');
title('sensitivity','FontSize',18);
subplot(2,2,3);
plot(sim{1,1}{1,1}(tstart:end,1),sim{1,1}{1,1}(tstart:end,11),'-o','MarkerIndices',1:stp_indices:length(sim{1,1}{1,1}(tstart:end,11)),'LineWidth',1);
hold on
plot(sim{1,1}{1,2}(tstart:end,1),sim{1,1}{1,2}(tstart:end,11),'LineWidth',1);
hold off
set(gca,'FontSize',14);
xlabel('time (h)','FontSize',16);
ylabel('[X] (nM)','FontSize',16);
%legend('e6.1+7 [FA]_l_o_w','e6.1+7 [FA]_h_i_g_h[FA]_l_o_w');
legend('e6.1+7 [FA]_l_o_w','e6.1+7 [FA]_h_i_g_h');
title('X','FontSize',18);
subplot(2,2,4);
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
%legend('e6.1+7 [FA]_l_o_w','e6.1+7 [FA]_h_i_g_h[FA]_l_o_w','FM [FA]_l_o_w','FM [FA]_h_i_g_h');
legend('e6.1+7 [FA]_l_o_w','e6.1+7 [FA]_h_i_g_h','FM [FA]_l_o_w','FM [FA]_h_i_g_h');
title('msok','FontSize',18);
%sgtitle('extension 6.2+7 - p-set 2439 +2425 (outlier sensitivity + hoksok ratio 6000)','FontSize',22);
%sgtitle('extension 6.2+7 - p-set 2439 +3082 (outlier sensitivity + max hoksok ratio)','FontSize',22);
%sgtitle('Time series extension 6.1+7 - [F] = 9uM at t=0.21h','FontSize',22);
sgtitle('Time series extension 6.1+7','FontSize',22);

%% 2. find the max hoksok ratio + corresponding timepoint and max sensitivity + timepoint
load('e6_2_7_parameters_o','p_sets')
N = 1;
q = 1;
[sim] = optimise_model_e6_2_7(p_sets(:,:),N,q);

% add hoksok ratio in the sim array (15) and the sensitivity (16)
for i=1:length(sim)
     sim{i}{1}(:,14) = (sim{i}{1}(:,9)+sim{i}{1}(:,8))./(sim{i}{1}(:,8)+sim{i}{1}(:,7));
     sim{i}{2}(:,14) = (sim{i}{2}(:,9)+sim{i}{2}(:,8))./(sim{i}{2}(:,8)+sim{i}{2}(:,7));
     sim{i}{1}(:,15) = (sim{i}{1}(:,14))./(sim{i}{2}(:,14));
     sim{i}{2}(:,15) = (sim{i}{1}(:,14))./(sim{i}{2}(:,14));
end
% find the max hoksok ratio (2) and the corresponding timepoint (3) + time point at
% which 95% of this max hoksok ratio is reached (4). 
     sim_max = max(sim{1}{2}(:,14)); %(2)
     t = find(sim{1}{2}(:,14)==sim_max(1,1))./100; %(3)
     thr(1,1) = 0.95*sim_max(1,1);
     t95(1,1) = (find(sim{1}{2}(:,14)>thr(1,1),1,'first'))./100; %(4)
     sim_hoksok_max = [[1:1]' sim_max t t95];
     
% find the max sensitivity (2) and the corresponding timepoint (3) + first time point at
% which a sensitivity > 100 is reached (4), last timepoint as which the sensitivity is > 100 (5) and the difference between these timepoints (6). 
     sim_max = max(sim{1}{2}(:,15)); %(2)
     t = find(sim{1}{2}(:,15)==sim_max(1,1))./100; %(3)
     t100_f(1,1) = (find(sim{1}{2}(:,15)>100,1,'first'))./100; %(4)
     t100_l(1,1) = (find(sim{1}{2}(:,15)>100,1,'last'))./100; %(5)
     delta_t = [t100_l(1,1) - t100_f(1,1)]; %(6)
     sim_sens_max = [[1:1]' sim_max t t100_f t100_l delta_t];
save('e6_2_7_output.mat','sim','sim_hoksok_max','sim_sens_max');
%% 3. Sensitivity analysis
load('e6_2_7_parameters_o','p_sets')
N=36; %nr of parameters
q=1;
for i=1:(length(p_sets(1,:)))
    p_new(i,:)=p_sets;
    p_new(i,i) = p_sets(:,i).*1.1; %multiply the parameters one by one by 10
end
    [sim] = optimise_model_e6_2_7(p_new,N,q); %simulate with new p-sets
    %add hoksok raio is column 15 and sensitivity in column 16
for i=1:length(sim)
    sim{i}{1}(:,14) = (sim{i}{1}(:,9)+sim{i}{1}(:,8))./(sim{i}{1}(:,8)+sim{i}{1}(:,7));
    sim{i}{2}(:,14) = (sim{i}{2}(:,9)+sim{i}{2}(:,8))./(sim{i}{2}(:,8)+sim{i}{2}(:,7));
    sim{i}{1}(:,15) = (sim{i}{1}(:,14))./(sim{i}{2}(:,14));
    sim{i}{2}(:,15) = (sim{i}{1}(:,14))./(sim{i}{2}(:,14));
end

% find the max hoksok ratio for each simulation and the corresponding
% timepoints + time at which 95% of the max hoksok ratio is reached
for i=1:length(sim)
     sim_max(i,1) = max(sim{i}{2}(:,14));
     t(i,1) = find(sim{i}{2}(:,14)==sim_max(i,1))./100;
     thr(i,1) = 0.95*sim_max(i,1);
     t95(i,1) = (find(sim{i}{2}(:,14)>thr(i,1),1,'first'))./100;
end
% combine the sim_max, t and t95 in one table
    sim_hoksok_max = [[1:36]' sim_max t t95];

% find the max sensitivities for each simulation and the corresponding timepoints 
for i=1:length(sim)
     sim_max(i,1) = max(sim{i}{1}(:,15));
     t(i,1) = find(sim{i}{1}(:,15)==sim_max(i,1))./100;
     thr(i,1) = 0.95*sim_max(i,1);
     t95(i,1) = (find(sim{i}{2}(:,15)>thr(i,1),1,'first'))./100;
end
% combine the sim_max, t and t95 in one table
    sim_sens_max = [[1:36]' sim_max t t95];
    save('e6_2_7_SA_times1.1.mat','p_new','sim','sim_hoksok_max','sim_sens_max');

for i=1:(length(p_sets(1,:)))
    p_new(i,:)=p_sets;
    p_new(i,i) = p_sets(:,i)./1.1; %divide the parameters one by one by 10
end
    [sim] = optimise_model_e6_2_7(p_new,N,q); %simulate with new p-sets
    %add hoksok raio is column 15 and sensitivity in column 16
for i=1:length(sim)
    sim{i}{1}(:,14) = (sim{i}{1}(:,9)+sim{i}{1}(:,8))./(sim{i}{1}(:,8)+sim{i}{1}(:,7));
    sim{i}{2}(:,14) = (sim{i}{2}(:,9)+sim{i}{2}(:,8))./(sim{i}{2}(:,8)+sim{i}{2}(:,7));
    sim{i}{1}(:,15) = (sim{i}{1}(:,14))./(sim{i}{2}(:,14));
    sim{i}{2}(:,15) = (sim{i}{1}(:,14))./(sim{i}{2}(:,14));
end

% find the max hoksok ratio for each simulation and the corresponding
% timepoints + time at which 95% of the max hoksok ratio is reached
for i=1:length(sim)
     sim_max(i,1) = max(sim{i}{2}(:,14));
     t(i,1) = find(sim{i}{2}(:,14)==sim_max(i,1))./100;
     thr(i,1) = 0.95*sim_max(i,1);
     t95(i,1) = (find(sim{i}{2}(:,14)>thr(i,1),1,'first'))./100;
end
% combine the sim_max, t and t95 in one table
    sim_hoksok_max = [[1:36]' sim_max t t95];

% find the max sensitivities for each simulation and the corresponding timepoints 
for i=1:length(sim)
     sim_max(i,1) = max(sim{i}{1}(:,15));
     t(i,1) = find(sim{i}{1}(:,15)==sim_max(i,1))./100;
     thr(i,1) = 0.95*sim_max(i,1);
     t95(i,1) = (find(sim{i}{2}(:,15)>thr(i,1),1,'first'))./100;
end
% combine the sim_max, t and t95 in one table
    sim_sens_max = [[1:36]' sim_max t t95];
    save('e6_2_7_SA_div1.1.mat','p_new','sim','sim_hoksok_max','sim_sens_max');
   
%% 4. Plot sensitivity analysis
load('e6_2_7_SA_times1.1.mat','sim','sim_hoksok_max','sim_sens_max');
sim_sens_max_t = sim_sens_max;
%get end concentrations of the simulations
for i=1:length(sim)
     end_conc_low_F_times10(i,1:length(sim{i}{1}(end,:))) = sim{i}{1}(end,:);
     end_conc_high_F_times10(i,1:length(sim{i}{2}(end,:))) = sim{i}{2}(end,:);            
end

load('e6_2_7_SA_div1.1.mat','sim','sim_hoksok_max','sim_sens_max');
sim_sens_max_d = sim_sens_max;
%get end concentrations of the simulations
for i=1:length(sim)
     end_conc_low_F_div10(i,1:length(sim{i}{1}(end,:))) = sim{i}{1}(end,:);
     end_conc_high_F_div10(i,1:length(sim{i}{2}(end,:))) = sim{i}{2}(end,:);            
end

%make matrix for scatterplot
load('e6_2_7_output.mat','sim','sim_hoksok_max','sim_sens_max');
sim_max = max(sim{1}{2}(:,15));
RP_low  = [ones(1, 36).*sim{1}{1}(end,14); ones(1, 36).*sim{1}{1}(end,15); ones(1, 36).*sim_sens_max(2)]; %end hoksok ratio (:,1), end sensitivity (:,2) and max sensitivity (:,3) of the original parameter set. 
RP_high = [ones(1, 36).*sim{1}{2}(end,14); ones(1, 36).*sim{1}{2}(end,15)]; %end hoksok ratio and end sensitivity of the original parameter set. 
Nr = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36];

%scatterplot comparing end sensitivities of the SA 
figure(5);
scatter(Nr(1,:),RP_high(2,:),50,'g','filled','d');
hold on
scatter(Nr(1,:),end_conc_low_F_div10(:,15),50,'b','filled','o');
hold on
scatter(Nr(1,:),end_conc_low_F_times10(:,15),50,'r','filled','d');
hold off
set(gca,'FontSize',14);
xlabel('parameter nr');
ylabel('sensitivity');
legend('RF','div1.1','times1.1');
title('Extension 6.1+7 - SA - sensitivity','FontSize',22);

%scatterplot comparing end hoksok ratio of the SA (times 10)
figure(6);
scatter(Nr(1,:),RP_high(1,:),50,'g','filled','d');
hold on
scatter(Nr(1,:),RP_low(1,:),50,'y','filled','o');
hold on
scatter(Nr(1,:),end_conc_high_F_times10(:,14),50,'b','filled','d');
hold on
scatter(Nr(1,:),end_conc_low_F_times10(:,14),50,'r','filled','o');
hold off
set(gca,'FontSize',14);
xlabel('parameter nr');
ylabel('hok/sok ratio');
set(gca,'yscale','log');
legend('RF [FA]_h_i_g_h','RF [FA]_l_o_w','times10 [FA]_h_i_g_h','times10 [FA]_l_o_w');
title('Extension 6.1+7 - SA - hok/sok ratio - times 10','FontSize',22);

%scatterplot comparing end hoksok ratio of the SA (div 10)
figure(7);
scatter(Nr(1,:),RP_high(1,:),50,'g','filled','d');
hold on
scatter(Nr(1,:),RP_low(1,:),50,'y','filled','o');
hold on
scatter(Nr(1,:),end_conc_high_F_div10(:,14),50,'b','filled','d');
hold on
scatter(Nr(1,:),end_conc_low_F_div10(:,14),50,'r','filled','o');
hold off
set(gca,'FontSize',14);
xlabel('parameter nr');
ylabel('hok/sok ratio');
set(gca,'yscale','log');
legend('RF [FA]_h_i_g_h','RF [FA]_l_o_w','div10 [FA]_h_i_g_h','div10 [FA]_l_o_w');
title('Extension 6.1+7 - SA - hok/sok ratio - div 10','FontSize',22);

%scatterplot comparing max sensitivities
figure(8);
scatter(Nr(1,:),RP_low(3,:),50,'g','filled','d');
hold on
scatter(Nr(1,:),sim_sens_max_d(:,2),50,'b','filled','o');
hold on
scatter(Nr(1,:),sim_sens_max_t(:,2),50,'r','filled','d');
hold off
set(gca,'FontSize',14);
xlabel('parameter nr');
ylabel('sensitivity');
%set(gca,'yscale','log');
legend('RF','div10','times10');
title('Extension 6.1+7 - SA - max sensitivity','FontSize',22);

%% 5. run first cycle with initial conditions 0
% use this and (9) to obtain simulations where the concentration of
% compounds in the cells is halved every x hours. This is the first cycle
% with initial concentrations 0 (don't forget to edit the
% optimise_model_e1_7 file)
load('e6_2_7_parameters_o','p_sets')
N = 1;
q = 1;
[sim] = optimise_model_e6_2_7(p_sets(:,:),N,q);

    % This part can be used to save end concentrations of the simulation
    c_highF = sim{1}{2}(end,2:end);
    c_lowF = sim{1}{1}(end,2:end);
    save('e6_2_7_conc.mat','c_highF','c_lowF');
    save('e6_2_7_conc_d24h_0.mat','c_highF','c_lowF','sim');
%% 6. run multiple cycles with initial concentration = half end concentration previous cycle 
% use this and (8) to obtain simulations where the concentration of
% compounds in the cells is halved every x hours. This are the other cycles
% (don't forget to edit the optimise_model_e1_7 file)
load('e6_2_7_conc_d24h_0.mat','c_highF','c_lowF','sim');
sim_all = sim;
for i=1:1 %1:(nr of cycles -1)
    [sim] = optimise_model_e6_2_7(p_sets(:,:),N,q);
    c_highF = sim{1}{2}(end,2:end);
    c_lowF = sim{1}{1}(end,2:end);
    %add concentrations to the previous cycle
    sim_all{1}{1}((2+(length(sim{1}{1}(:,1))-1)*i):(length(sim{1}{1}(:,1))+(length(sim{1}{1}(:,1))-1)*i),1) = sim{1}{1}(2:end,1)+(length(sim{1}{1}(:,1))-1)*i/100;
    sim_all{1}{2}((2+(length(sim{1}{1}(:,1))-1)*i):(length(sim{1}{1}(:,1))+(length(sim{1}{1}(:,1))-1)*i),1) = sim{1}{2}(2:end,1)+(length(sim{1}{1}(:,1))-1)*i/100;
    sim_all{1}{1}((2+(length(sim{1}{1}(:,1))-1)*i):(length(sim{1}{1}(:,1))+(length(sim{1}{1}(:,1))-1)*i),2:13) = sim{1}{1}(2:end,2:13);
    sim_all{1}{2}((2+(length(sim{1}{1}(:,1))-1)*i):(length(sim{1}{1}(:,1))+(length(sim{1}{1}(:,1))-1)*i),2:13) = sim{1}{2}(2:end,2:13);
    save('e6_2_7_conc.mat','c_highF','c_lowF');
    save('e6_2_7_conc_d24h_2c_2','c_highF','c_lowF','sim_all');           
end