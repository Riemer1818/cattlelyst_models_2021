%% Goal: make time-series plots of the original biosafety circuit
% contains:
% 1. Time series plots full model
% 2. Time series plot SA 4482 p-sets

% input files;
%    median_low_hoksokratio2
%    outputFM_SAGlowRatio_times10_16
%    outputFM_SAGlowRatio_times10_1
%    outputFM_SAGlowRatio_div10_7

%% 1. Time series plots full model
load('median_low_hoksokratio2','params_GhighF_ratio'); %load good parameter sets
%load('FM_Gsim_high_low.mat','Gsim_high','Gsim_low','Gparams'); %load all parameter sets
N = 1;
q = 1;

% select parameter set
p_set = 3048;
%p_set = 1056;
%p_set = 4293;

%simulate model with the parameter set
[sim] = optimise_model(params_GhighF_ratio(p_set,:),N,q);
for i=1:length(sim)
     low_F(i,1:length(sim{i}{1}(end,:))) = sim{i}{1}(end,:);
     high_F(i,1:length(sim{i}{2}(end,:))) = sim{i}{2}(end,:);            
end
A_high_F = [[1:1]' high_F (high_F(:,9)+high_F(:,8))./(high_F(:,8)+high_F(:,7))];
A_low_F = [[1:1]' low_F (low_F(:,9)+low_F(:,8))./(low_F(:,8)+low_F(:,7))];
B_high_F = [A_high_F A_low_F(:,11)-A_high_F(:,11) A_low_F(:,11)./A_high_F(:,11)];
B_low_F = [A_low_F A_low_F(:,11)-A_high_F(:,11) A_low_F(:,11)./A_high_F(:,11)];

% Matrix A contains a time vector (first column) + time series data + hoksok ratio (last column).
% Matrix B contains matrix A + hoksok ratio [F]low - hoksok ratio [F]high + hoksok ratio [F]low/hoksok ratio [F]high.
    A_high_F_TS = [[2:length(sim{1}{2}(:,1))]' sim{1}{2}(2:end,:) (sim{1}{2}(2:end,9)+sim{1}{2}(2:end,8))./(sim{1}{2}(2:end,8)+sim{1}{2}(2:end,7))];
    A_low_F_TS = [[2:length(sim{1}{2}(:,1))]' sim{1}{2}(2:end,:) (sim{1}{1}(2:end,9)+sim{1}{1}(2:end,8))./(sim{1}{1}(2:end,8)+sim{1}{1}(2:end,7))];
    B_high_F_TS = [A_high_F_TS A_low_F_TS(:,11)-A_high_F_TS(:,11) A_low_F_TS(:,11)./A_high_F_TS(:,11)];
    B_low_F_TS = [A_low_F_TS A_low_F_TS(:,11)-A_high_F_TS(:,11) A_low_F_TS(:,11)./A_high_F_TS(:,11)];

stp_indices = 1000;
tstart = 1;
figure(11);
subplot(2,2,1);
plot(sim{1,1}{1,1}(tstart:end,1),sim{1,1}{1,1}(tstart:end,2),'LineWidth',1);
hold on
plot(sim{1,1}{1,2}(tstart:end,1),sim{1,1}{1,2}(tstart:end,2),'LineWidth',1);
hold off
set(gca,'FontSize',14);
xlabel('time (h)','Fontsize',16);
ylabel('[mFrmR] (nM)','Fontsize',16);
legend('[FA]_l_o_w','[FA]_h_i_g_h');
title('mFrmR','Fontsize',18);
subplot(2,2,2);
plot(sim{1,1}{1,1}(tstart:end,1),sim{1,1}{1,1}(tstart:end,3),'LineWidth',1);
hold on
plot(sim{1,1}{1,2}(tstart:end,1),sim{1,1}{1,2}(tstart:end,3),'LineWidth',1);
hold off
set(gca,'FontSize',14);
xlabel('time (h)','Fontsize',16);
ylabel('[FrmR] (nM)','Fontsize',16);
legend('[FA]_l_o_w','[FA]_h_i_g_h');
title('FrmR','Fontsize',18);
subplot(2,2,3);
plot(sim{1,1}{1,1}(tstart:end,1),sim{1,1}{1,1}(tstart:end,4),'LineWidth',1);
hold on
plot(sim{1,1}{1,2}(tstart:end,1),sim{1,1}{1,2}(tstart:end,4),'LineWidth',1);
hold off
set(gca,'FontSize',14);
xlabel('time (h)','Fontsize',16);
ylabel('[FA-FrmR] (nM)','Fontsize',16);
legend('[FA]_l_o_w','[FA]_h_i_g_h');
title('FA-FrmR','Fontsize',18);
subplot(2,2,4);
plot(sim{1,1}{1,1}(tstart:end,1),sim{1,1}{1,1}(tstart:end,5),'LineWidth',1);
hold on
plot(sim{1,1}{1,2}(tstart:end,1),sim{1,1}{1,2}(tstart:end,5),'LineWidth',1);
hold off
set(gca,'FontSize',14);
xlabel('time (h)','Fontsize',16);
ylabel('[mLacI] (nM)','Fontsize',16);
legend('[FA]_l_o_w','[FA]_h_i_g_h');
title('mLacI','Fontsize',18);
sgtitle('Time series plots original biosafety circuit - p-set with highest sensitivity at t=48h','Fontsize',22);

figure(12);
subplot(2,2,1);
plot(sim{1,1}{1,1}(tstart:end,1),sim{1,1}{1,1}(tstart:end,6),'LineWidth',1);
hold on
plot(sim{1,1}{1,2}(tstart:end,1),sim{1,1}{1,2}(tstart:end,6),'LineWidth',1);
hold off
set(gca,'FontSize',14);
xlabel('time (h)','Fontsize',16);
ylabel('[LacI] (nM)','Fontsize',16);
legend('[FA]_l_o_w','[FA]_h_i_g_h');
title('LacI','Fontsize',18);
subplot(2,2,2);
plot(sim{1,1}{1,1}(tstart:end,1),sim{1,1}{1,1}(tstart:end,7),'LineWidth',1);
hold on
plot(sim{1,1}{1,2}(tstart:end,1),sim{1,1}{1,2}(tstart:end,7),'LineWidth',1);
hold off
set(gca,'FontSize',14);
xlabel('time (h)','Fontsize',16);
ylabel('[msok] (nM)','Fontsize',16);
legend('[FA]_l_o_w','[FA]_h_i_g_h');
title('msok','Fontsize',18);
subplot(2,2,3);
plot(sim{1,1}{1,1}(tstart:end,1),sim{1,1}{1,1}(tstart:end,8),'LineWidth',1);
hold on
plot(sim{1,1}{1,2}(tstart:end,1),sim{1,1}{1,2}(tstart:end,8),'LineWidth',1);
hold off
set(gca,'FontSize',14);
xlabel('time (h)','Fontsize',16);
ylabel('[hok-sok] (nM)','Fontsize',16);
legend('[FA]_l_o_w','[FA]_h_i_g_h');
title('hok-sok','Fontsize',18);
subplot(2,2,4);
plot(sim{1,1}{1,1}(tstart:end,1),sim{1,1}{1,1}(tstart:end,9),'LineWidth',1);
hold on
plot(sim{1,1}{1,2}(tstart:end,1),sim{1,1}{1,2}(tstart:end,9),'LineWidth',1);
hold off
set(gca,'FontSize',14);
xlabel('time (h)','Fontsize',16);
ylabel('[mhok] (nM)','Fontsize',16);
%set(gca,'yscale','log');
legend('[FA]_l_o_w','[FA]_h_i_g_h');
title('mhok','Fontsize',18);
sgtitle('Time series plots original biosafety circuit - p-set with highest sensitivity at t=48h','Fontsize',22);

figure(13);
subplot(2,2,1);
plot(B_low_F_TS(:,2),B_low_F_TS(:,11),'LineWidth',1);
hold on
plot(B_high_F_TS(:,2),B_high_F_TS(:,11),'LineWidth',1);
hold off
set(gca,'FontSize',14);
xlabel('time (h)','Fontsize',16);
ylabel('hok/sok ratio','Fontsize',16);
%set(gca,'yscale','log');
legend('[FA]_l_o_w','[FA]_h_i_g_h');
title('hok/sok ratio','Fontsize',18);
subplot(2,2,2);
plot(B_low_F_TS(:,2),B_low_F_TS(:,13),'LineWidth',1);
hold off
set(gca,'FontSize',14);
xlabel('time (h)','Fontsize',16);
ylabel('sensitivity','Fontsize',16);
title('sensitivity','Fontsize',18);
sgtitle('Time series plots original biosafety circuit - p-set with highest sensitivity at t=48h','Fontsize',22);
%% 2. Time series plot SA 4482 p-sets
load('median_low_hoksokratio2','params_GhighF_ratio');
%load('outputFM_SAGlowRatio_times10_16.mat','p_new')
load('outputFM_SAGlowRatio_times10_1.mat','p_new')
%load('outputFM_SAGlowRatio_div10_7.mat','p_new')
N = 1;
q = 1;

%select parameter set from SA to simulate
%p_set = 784; %max hok/sok ratio increase for parameter 7 divided by 10
p_set = 328; %max hok/sok ratio increase for parameter 1 times 10
%p_set = 1056; %max increase in formaldehyde sensitivity
[sim] = optimise_model(p_new(p_set,:),N,q);
[sim_old] = optimise_model(params_GhighF_ratio(p_set,:),N,q);
%2-4801
for i=1:length(sim)
     low_F(i,1:length(sim{i}{1}(end,:))) = sim{i}{1}(end,:);
     high_F(i,1:length(sim{i}{2}(end,:))) = sim{i}{2}(end,:);            
end
A_high_F = [[1:1]' high_F (high_F(:,9)+high_F(:,8))./(high_F(:,8)+high_F(:,7))];
A_low_F = [[1:1]' low_F (low_F(:,9)+low_F(:,8))./(low_F(:,8)+low_F(:,7))];
B_high_F = [A_high_F A_low_F(:,11)-A_high_F(:,11) A_low_F(:,11)./A_high_F(:,11)];
B_low_F = [A_low_F A_low_F(:,11)-A_high_F(:,11) A_low_F(:,11)./A_high_F(:,11)];

% Matrix A contains a time vector (first column) + time series data + hoksok ratio (last column).
% Matrix B contains matrix A + hoksok ratio [F]low - hoksok ratio [F]high + hoksok ratio [F]low/hoksok ratio [F]high.
% Matrix C contains matrix A + new sensitivity if parameter is made formaldehyde dependent. 
    A_high_F_TSold = [[2:length(sim_old{1}{2}(:,1))]' sim_old{1}{2}(2:end,:) (sim_old{1}{2}(2:end,9)+sim_old{1}{2}(2:end,8))./(sim_old{1}{2}(2:end,8)+sim_old{1}{2}(2:end,7))];
    A_low_F_TSold = [[2:length(sim_old{1}{2}(:,1))]' sim_old{1}{2}(2:end,:) (sim_old{1}{1}(2:end,9)+sim_old{1}{1}(2:end,8))./(sim_old{1}{1}(2:end,8)+sim_old{1}{1}(2:end,7))];
    B_high_F_TSold = [A_high_F_TSold A_low_F_TSold(:,11)-A_high_F_TSold(:,11) A_low_F_TSold(:,11)./A_high_F_TSold(:,11)];
    B_low_F_TSold = [A_low_F_TSold A_low_F_TSold(:,11)-A_high_F_TSold(:,11) A_low_F_TSold(:,11)./A_high_F_TSold(:,11)];

    A_high_F_TSnew = [[2:length(sim{1}{2}(:,1))]' sim{1}{2}(2:end,:) (sim{1}{2}(2:end,9)+sim{1}{2}(2:end,8))./(sim{1}{2}(2:end,8)+sim{1}{2}(2:end,7))];
    A_low_F_TSnew = [[2:length(sim{1}{2}(:,1))]' sim{1}{2}(2:end,:) (sim{1}{1}(2:end,9)+sim{1}{1}(2:end,8))./(sim{1}{1}(2:end,8)+sim{1}{1}(2:end,7))];
    B_high_F_TSnew = [A_high_F_TSnew A_low_F_TSnew(:,11)-A_high_F_TSnew(:,11) A_low_F_TSnew(:,11)./A_high_F_TSnew(:,11)];
    B_low_F_TSnew = [A_low_F_TSnew A_low_F_TSnew(:,11)-A_high_F_TSnew(:,11) A_low_F_TSnew(:,11)./A_high_F_TSnew(:,11)];
    C_high_F_TSnew = [A_high_F_TSnew A_low_F_TSnew(:,11)./A_high_F_TSold(:,11)];
    C_low_F_TSnew = [A_low_F_TSnew A_low_F_TSnew(:,11)./A_high_F_TSold(:,11)];

%time series plots
stp_indices = 1000; %changes frequency of the circles in the lines
tstart = 1000;
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
legend('new p-set [FA]_l_o_w','new p-set [FA]_h_i_g_h','old p-set [FA]_l_o_w','old p-set [FA]_h_i_g_h');
title('mFrmR');
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
legend('new p-set [FA]_l_o_w','new p-set [FA]_h_i_g_h','old p-set [FA]_l_o_w','old p-set [FA]_h_i_g_h');
title('FrmR');
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
legend('new p-set [FA]_l_o_w','new p-set [FA]_h_i_g_h','old p-set [FA]_l_o_w','old p-set [FA]_h_i_g_h');
title('FA-FrmR');
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
legend('new p-set [FA]_l_o_w','new p-set [FA]_h_i_g_h','old p-set [FA]_l_o_w','old p-set [FA]_h_i_g_h');
title('mLacI');
sgtitle('Time series plots - sensitivity analysis - parameter 1 times 10','Fontsize',22);

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
legend('new p-set [FA]_l_o_w','new p-set [FA]_h_i_g_h','old p-set [FA]_l_o_w','old p-set [FA]_h_i_g_h');
title('LacI');
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
legend('new p-set [FA]_l_o_w','new p-set [FA]_h_i_g_h','old p-set [FA]_l_o_w','old p-set [FA]_h_i_g_h');
title('msok');
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
legend('new p-set [FA]_l_o_w','new p-set [FA]_h_i_g_h','old p-set [FA]_l_o_w','old p-set [FA]_h_i_g_h');
title('hok-sok');
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
legend('new p-set [FA]_l_o_w','new p-set [FA]_h_i_g_h','old p-set [FA]_l_o_w','old p-set [FA]_h_i_g_h');
title('mhok');
sgtitle('Time series plots - sensitivity analysis - parameter 1 times 10','Fontsize',22);

figure(13);
subplot(2,2,1);
plot(B_low_F_TSnew(:,2),B_low_F_TSnew(:,11),'-o','MarkerIndices',1:stp_indices:length(B_low_F_TSnew(:,2)),'LineWidth',1);
hold on
plot(B_high_F_TSnew(:,2),B_high_F_TSnew(:,11),'LineWidth',1);
hold on
plot(B_low_F_TSold(:,2),B_low_F_TSold(:,11),'-o','MarkerIndices',1:stp_indices:length(B_high_F_TSnew(:,2)),'LineWidth',1);
hold on
plot(B_high_F_TSold(:,2),B_high_F_TSold(:,11),'LineWidth',1);
hold off
set(gca,'FontSize',14);
xlabel('time (h)','FontSize',16);
ylabel('hok/sok ratio','FontSize',16);
%set(gca,'yscale','log');
legend('new p-set [FA]_l_o_w','new p-set [FA]_h_i_g_h','old p-set [FA]_l_o_w','old p-set [FA]_h_i_g_h');
title('hok/sok ratio');
subplot(2,2,2);
plot(C_low_F_TSnew(:,2),C_low_F_TSnew(:,12),'-o','MarkerIndices',1:stp_indices:length(C_low_F_TSnew(:,2)),'LineWidth',1);
hold on
plot(B_low_F_TSold(:,2),B_low_F_TSold(:,13),'LineWidth',1);
hold off
set(gca,'FontSize',14);
xlabel('time (h)','FontSize',16);
ylabel('sensitivity','FontSize',16);
legend('new p-set','old p-set');
title('sensitivity');
sgtitle('Time series plots - sensitivity analysis - parameter 1 times 10','Fontsize',22);
