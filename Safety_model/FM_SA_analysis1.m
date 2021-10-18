%% file used to analyse the sensitivity analysis of the original biosafety circuit
% contains:
% 1. Select the parameter sets and simulations that give a good hok/sok ratio after 48 hours
% 2. analyse results SA times 10
% 3. analyse results SA div 10
% 4. Plot sensitivity analysis - scatterplot
% 5. Plot sensitivity analysis - boxplot
% 6. Select outliers

% input files: 
%     median_low_hoksokratio2_48hsim.mat
%     outputFM_SAGlowRatio_times10_*.mat
%     outputFM_SAGlowRatio_div10_*.mat
% output files:
%     FM_SA_GlowRatio_*D_2.mat
%     FM_SA_GlowRatio_*T_2.mat
%% 1. Select the parameter sets and simulations that give a good hok/sok ratio after 48 hours
load('median_low_hoksokratio2_48hsim.mat', 'sim')
% creating a matrix with end concentrations for the simulations of each
% p-set.
    low_F_SA = zeros(4482, 9); 
    high_F_SA= zeros(4482, 9);
    for i=1:length(sim)
         low_F_SA(i,1:length(sim{i}{1}(end,:))) = sim{i}{1}(end,:);
         high_F_SA(i,1:length(sim{i}{2}(end,:))) = sim{i}{2}(end,:);            
    end
% Matrix A contains a new numbering of the simulations (first column) + end time (column 2) + end
% concentrations of the compounds (column 3-10) + hoksok ratio (column 11).
% Matrix B contains matrix A (column 1-11) + hoksok ratio [F]low - hoksok
% ratio [F]high + hoksok ratio [F]low/hoksok ratio [F]high.
    A_high_F_SA = [[1:4482]' high_F_SA (high_F_SA(:,9)+high_F_SA(:,8))./(high_F_SA(:,8)+high_F_SA(:,7))];
    A_low_F_SA = [[1:4482]' low_F_SA (low_F_SA(:,9)+low_F_SA(:,8))./(low_F_SA(:,8)+low_F_SA(:,7))];
    sort_A_high_SA = sortrows(A_high_F_SA,11);
    %select parameter sets that give a good hok/sok ratio after 48h of
    %simulation
    GP_48h = sort_A_high_SA(2:3120,1);
    B_high_SA = [A_high_F_SA(GP_48h,:) A_low_F_SA(GP_48h,11)-A_high_F_SA(GP_48h,11) A_low_F_SA(GP_48h,11)./A_high_F_SA(GP_48h,11)];
    B_low_SA = [A_low_F_SA(GP_48h,:) A_low_F_SA(GP_48h,11)-A_high_F_SA(GP_48h,11) A_low_F_SA(GP_48h,11)./A_high_F_SA(GP_48h,11)];
    
%% 2. analyse results SA times 10
load('outputFM_SAGlowRatio_times10_19.mat','sim')
% parameter sets that give the upper or lower 95% confidence interval. 
    CI_lower = round(3119*0.5-1.96*(3119*0.5*(1-0.5))^(0.5));
    CI_upper = round(3119*0.5+1.96*(3119*0.5*(1-0.5))^(0.5));
% creating a matrix with end concentrations for the simulations of each
% p-set.
    low_F_SABP_times10 = zeros(4482, 9); 
    high_F_SABP_times10= zeros(4482, 9);
    for i=1:length(sim)
         low_F_SABP_times10(i,1:length(sim{i}{1}(end,:))) = sim{i}{1}(end,:);
         high_F_SABP_times10(i,1:length(sim{i}{2}(end,:))) = sim{i}{2}(end,:);            
    end
% Matrix A contains a new numbering of the simulations (first column) + end time (column 2) + end
% concentrations of the compounds (column 3-10) + hoksok ratio (column 11).
% Matrix B contains matrix A (column 1-11) + hoksok ratio [F]low - hoksok
% ratio [F]high + hoksok ratio [F]low/hoksok ratio [F]high.
    A_high_F_SABP_times10 = [[1:4482]' high_F_SABP_times10 (high_F_SABP_times10(:,9)+high_F_SABP_times10(:,8))./(high_F_SABP_times10(:,8)+high_F_SABP_times10(:,7))];
    A_low_F_SABP_times10 = [[1:4482]' low_F_SABP_times10 (low_F_SABP_times10(:,9)+low_F_SABP_times10(:,8))./(low_F_SABP_times10(:,8)+low_F_SABP_times10(:,7))];
    B_high_SABP_times10_19T = [A_high_F_SABP_times10(GP_48h,:) A_low_F_SABP_times10(GP_48h,11)-A_high_F_SABP_times10(GP_48h,11) A_low_F_SABP_times10(GP_48h,11)./A_high_F_SABP_times10(GP_48h,11)];
    B_low_SABP_times10_19T = [A_low_F_SABP_times10(GP_48h,:) A_low_F_SABP_times10(GP_48h,11)-A_high_F_SABP_times10(GP_48h,11) A_low_F_SABP_times10(GP_48h,11)./A_high_F_SABP_times10(GP_48h,11)];
% Get the hoksok ratio with new p-set/old p-set
    %load('median_low_hoksokratio2','sim_GhighF_ratio')
    ratio_HFHF_hoksok19T = [B_high_SABP_times10_19T(:,11)./B_high_SA(:,11)];
    ratio_LFHF_hoksok19T = [B_low_SABP_times10_19T(:,11)./B_high_SA(:,11)];
    ratio_sensitivity19T = [B_high_SABP_times10_19T(:,13)./B_high_SA(:,13)];
% Get the median values and 95% confidence intervals (column 1: HFHF,
% column 2 LFHF, column 3 sensitiviy)
    sort_ratio_HFHF_hoksok = sortrows(ratio_HFHF_hoksok19T);
    sort_ratio_LFHF_hoksok = sortrows(ratio_LFHF_hoksok19T);
    sort_ratio_sensitivity = sortrows(ratio_sensitivity19T);
    median19T = [median(sort_ratio_HFHF_hoksok) median(sort_ratio_LFHF_hoksok) median(sort_ratio_sensitivity)];
    CI19T_l = [sort_ratio_HFHF_hoksok(CI_lower) sort_ratio_LFHF_hoksok(CI_lower) sort_ratio_sensitivity(CI_lower)];
    CI19T_u = [sort_ratio_HFHF_hoksok(CI_upper) sort_ratio_LFHF_hoksok(CI_upper) sort_ratio_sensitivity(CI_upper)];

%save('FM_SA_GlowRatio_19T_2.mat','B_high_SABP_times10_19T','B_low_SABP_times10_19T','ratio_HFHF_hoksok19T','ratio_LFHF_hoksok19T','ratio_sensitivity19T','median19T','CI19T_l','CI19T_u');
%% 3. analyse results SA div 10
load('outputFM_SAGlowRatio_div10_20.mat','sim')
% parameter sets that give the upper or lower 95% confidence interval. 
    CI_lower = round(3119*0.5-1.96*(3119*0.5*(1-0.5))^(0.5));
    CI_upper = round(3119*0.5+1.96*(3119*0.5*(1-0.5))^(0.5));
% creating a matrix with end concentrations for the simulations of each
% p-set.
    low_F_SABP_div10 = zeros(4482, 9); 
    high_F_SABP_div10= zeros(4482, 9);
    for i=1:length(sim)
         low_F_SABP_div10(i,1:length(sim{i}{1}(end,:))) = sim{i}{1}(end,:);
         high_F_SABP_div10(i,1:length(sim{i}{2}(end,:))) = sim{i}{2}(end,:);            
    end
% Matrix A contains a new numbering of the simulations (first column) + end time (column 2) + end
% concentrations of the compounds (column 3-10) + hoksok ratio (column 11).
% Matrix B contains matrix A (column 1-11) + hoksok ratio [F]low - hoksok
% ratio [F]high + hoksok ratio [F]low/hoksok ratio [F]high.
    A_high_F_SABP_div10 = [[1:4482]' high_F_SABP_div10 (high_F_SABP_div10(:,9)+high_F_SABP_div10(:,8))./(high_F_SABP_div10(:,8)+high_F_SABP_div10(:,7))];
    A_low_F_SABP_div10 = [[1:4482]' low_F_SABP_div10 (low_F_SABP_div10(:,9)+low_F_SABP_div10(:,8))./(low_F_SABP_div10(:,8)+low_F_SABP_div10(:,7))];
    B_high_SABP_div10_20D = [A_high_F_SABP_div10(GP_48h,:) A_low_F_SABP_div10(GP_48h,11)-A_high_F_SABP_div10(GP_48h,11) A_low_F_SABP_div10(GP_48h,11)./A_high_F_SABP_div10(GP_48h,11)];
    B_low_SABP_div10_20D = [A_low_F_SABP_div10(GP_48h,:) A_low_F_SABP_div10(GP_48h,11)-A_high_F_SABP_div10(GP_48h,11) A_low_F_SABP_div10(GP_48h,11)./A_high_F_SABP_div10(GP_48h,11)];
% Get the hoksok ratio with new p-set/old p-set
    %load('median_low_hoksokratio2','sim_GhighF_ratio')
    ratio_HFHF_hoksok20D = [B_high_SABP_div10_20D(:,11)./B_high_SA(:,11)];
    ratio_LFHF_hoksok20D = [B_low_SABP_div10_20D(:,11)./B_high_SA(:,11)];
    ratio_sensitivity20D = [B_high_SABP_div10_20D(:,13)./B_high_SA(:,13)];
% Get the median values and 95% confidence intervals (column 1: HFHF,
% column 2 LFHF, column 3 sensitiviy)
    sort_ratio_HFHF_hoksok = sortrows(ratio_HFHF_hoksok20D);
    sort_ratio_LFHF_hoksok = sortrows(ratio_LFHF_hoksok20D);
    sort_ratio_sensitivity = sortrows(ratio_sensitivity20D);
    median20D = [median(sort_ratio_HFHF_hoksok) median(sort_ratio_LFHF_hoksok) median(sort_ratio_sensitivity)];
    CI20D_l = [sort_ratio_HFHF_hoksok(CI_lower) sort_ratio_LFHF_hoksok(CI_lower) sort_ratio_sensitivity(CI_lower)];
    CI20D_u = [sort_ratio_HFHF_hoksok(CI_upper) sort_ratio_LFHF_hoksok(CI_upper) sort_ratio_sensitivity(CI_upper)];

save('FM_SA_GlowRatio_20D_2.mat','B_high_SABP_div10_20D','B_low_SABP_div10_20D','ratio_HFHF_hoksok20D','ratio_LFHF_hoksok20D','ratio_sensitivity20D','median20D','CI20D_l','CI20D_u');    
%% 4. Plot sensitivity analysis - scatterplot
% Load all files containing the median and confidence intervals. 
    load('FM_SA_GlowRatio_1T_2.mat','median1T','CI1T_l','CI1T_u'); 
    load('FM_SA_GlowRatio_2T_2.mat','median2T','CI2T_l','CI2T_u'); 
    load('FM_SA_GlowRatio_3T_2.mat','median3T','CI3T_l','CI3T_u'); 
    load('FM_SA_GlowRatio_4T_2.mat','median4T','CI4T_l','CI4T_u'); 
    load('FM_SA_GlowRatio_5T_2.mat','median5T','CI5T_l','CI5T_u'); 
    load('FM_SA_GlowRatio_6T_2.mat','median6T','CI6T_l','CI6T_u'); 
    load('FM_SA_GlowRatio_7T_2.mat','median7T','CI7T_l','CI7T_u'); 
    load('FM_SA_GlowRatio_8T_2.mat','median8T','CI8T_l','CI8T_u'); 
    load('FM_SA_GlowRatio_9T_2.mat','median9T','CI9T_l','CI9T_u'); 
    load('FM_SA_GlowRatio_10T_2.mat','median10T','CI10T_l','CI10T_u'); 
    load('FM_SA_GlowRatio_11T_2.mat','median11T','CI11T_l','CI11T_u'); 
    load('FM_SA_GlowRatio_12T_2.mat','median12T','CI12T_l','CI12T_u'); 
    load('FM_SA_GlowRatio_13T_2.mat','median13T','CI13T_l','CI13T_u'); 
    load('FM_SA_GlowRatio_14T_2.mat','median14T','CI14T_l','CI14T_u'); 
    load('FM_SA_GlowRatio_15T_2.mat','median15T','CI15T_l','CI15T_u'); 
    load('FM_SA_GlowRatio_16T_2.mat','median16T','CI16T_l','CI16T_u'); 
    load('FM_SA_GlowRatio_17T_2.mat','median17T','CI17T_l','CI17T_u'); 
    load('FM_SA_GlowRatio_18T_2.mat','median18T','CI18T_l','CI18T_u'); 
    load('FM_SA_GlowRatio_19T_2.mat','median19T','CI19T_l','CI19T_u'); 
    load('FM_SA_GlowRatio_20T_2.mat','median20T','CI20T_l','CI20T_u'); 

    load('FM_SA_GlowRatio_1D_2.mat','median1D','CI1D_l','CI1D_u'); 
    load('FM_SA_GlowRatio_2D_2.mat','median2D','CI2D_l','CI2D_u'); 
    load('FM_SA_GlowRatio_3D_2.mat','median3D','CI3D_l','CI3D_u'); 
    load('FM_SA_GlowRatio_4D_2.mat','median4D','CI4D_l','CI4D_u'); 
    load('FM_SA_GlowRatio_5D_2.mat','median5D','CI5D_l','CI5D_u'); 
    load('FM_SA_GlowRatio_6D_2.mat','median6D','CI6D_l','CI6D_u'); 
    load('FM_SA_GlowRatio_7D_2.mat','median7D','CI7D_l','CI7D_u'); 
    load('FM_SA_GlowRatio_8D_2.mat','median8D','CI8D_l','CI8D_u'); 
    load('FM_SA_GlowRatio_9D_2.mat','median9D','CI9D_l','CI9D_u'); 
    load('FM_SA_GlowRatio_10D_2.mat','median10D','CI10D_l','CI10D_u'); 
    load('FM_SA_GlowRatio_11D_2.mat','median11D','CI11D_l','CI11D_u'); 
    load('FM_SA_GlowRatio_12D_2.mat','median12D','CI12D_l','CI12D_u'); 
    load('FM_SA_GlowRatio_13D_2.mat','median13D','CI13D_l','CI13D_u'); 
    load('FM_SA_GlowRatio_14D_2.mat','median14D','CI14D_l','CI14D_u'); 
    load('FM_SA_GlowRatio_15D_2.mat','median15D','CI15D_l','CI15D_u'); 
    load('FM_SA_GlowRatio_16D_2.mat','median16D','CI16D_l','CI16D_u'); 
    load('FM_SA_GlowRatio_17D_2.mat','median17D','CI17D_l','CI17D_u'); 
    load('FM_SA_GlowRatio_18D_2.mat','median18D','CI18D_l','CI18D_u'); 
    load('FM_SA_GlowRatio_19D_2.mat','median19D','CI19D_l','CI19D_u'); 
    load('FM_SA_GlowRatio_20D_2.mat','median20D','CI20D_l','CI20D_u');

% Combining the median and CI in one matrix
    Nr = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20];
    median_T = [median1T; median2T; median3T; median4T; median5T; median6T; median7T; median8T; median9T; median10T; median11T; median12T; median13T; median14T; median15T; median16T; median17T; median18T; median19T; median20T];
    CI_T_l = [CI1T_l; CI2T_l; CI3T_l; CI4T_l; CI5T_l; CI6T_l; CI7T_l; CI8T_l; CI9T_l; CI10T_l; CI11T_l; CI12T_l; CI13T_l; CI14T_l; CI15T_l; CI16T_l; CI17T_l; CI18T_l; CI19T_l; CI16T_l];
    CI_T_u = [CI1T_u; CI2T_u; CI3T_u; CI4T_u; CI5T_u; CI6T_u; CI7T_u; CI8T_u; CI9T_u; CI10T_u; CI11T_u; CI12T_u; CI13T_u; CI14T_u; CI15T_u; CI16T_u; CI17T_u; CI18T_u; CI19T_u; CI20T_u];

    median_D = [median1D; median2D; median3D; median4D; median5D; median6D; median7D; median8D; median9D; median10D; median11D; median12D; median13D; median14D; median15D; median16D; median17D; median18D; median19D; median20D];
    CI_D_l = [CI1D_l; CI2D_l; CI3D_l; CI4D_l; CI5D_l; CI6D_l; CI7D_l; CI8D_l; CI9D_l; CI10D_l; CI11D_l; CI12D_l; CI13D_l; CI14D_l; CI15D_l; CI16D_l; CI17D_l; CI18D_l; CI19D_l; CI20D_l];
    CI_D_u = [CI1D_u; CI2D_u; CI3D_u; CI4D_u; CI5D_u; CI6D_u; CI7D_u; CI8D_u; CI9D_u; CI10D_u; CI11D_u; CI12D_u; CI13D_u; CI14D_u; CI15D_u; CI16D_u; CI17D_u; CI18D_u; CI19D_u; CI20D_u];

% Plot senstivity analysis 
    % Choose how to plot the senstivitiy analysis 
        % 1=hoksok ratio [F]high_new/[F]high_old
        % 2=hoksok ratio [F]low_new/[F]high_old 
        % 3=sensitivity_new/sensitivity_old
        c = 3;
    
    figure(1);
    scatter(Nr(1,:),median_T(:,c),40,'r','filled','o');
    hold on
    scatter(Nr(1,:),median_D(:,c),30,'b','filled','d');
    hold on
    yline(1);
    hold on
    er = errorbar(Nr(1,:),median_T(:,c),(median_T(:,c)-CI_T_l(:,c)),(CI_T_u(:,c)-median_T(:,c)));
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';
    hold on
    er2 = errorbar(Nr(1,:),median_D(:,c),(median_D(:,c)-CI_D_l(:,c)),(CI_D_u(:,c)-median_D(:,c)));
    er2.Color = [0 0 0];                            
    er2.LineStyle = 'none';
    hold off
    set(gca,'Fontsize',14);
    xlabel('parameter nr.','Fontsize',16);
    %ylabel('(hok/sok ratio new p-sets ([F]_l_o_w))/(hok/sok ratio old p-sets ([F]_h_i_g_h))','Fontsize',16);
    ylabel('(sensitivity new p-sets)/(sensitivity old p-sets)');
    legend('times 10','div 10');
    title('Sensitivity analysis - effect of parameters on formaldehyde sensitivity','Fontsize',22);
%% 5. Plot sensitivity analysis - boxplot
% load all the files that contain the hoksok ratio change and sensitivity
% ratio change. 
    load('FM_SA_GlowRatio_1T_2.mat','B_high_SABP_times10_1T','B_low_SABP_times10_1T','ratio_HFHF_hoksok1T','ratio_LFHF_hoksok1T','ratio_sensitivity1T'); 
    load('FM_SA_GlowRatio_2T_2.mat','B_high_SABP_times10_2T','B_low_SABP_times10_2T','ratio_HFHF_hoksok2T','ratio_LFHF_hoksok2T','ratio_sensitivity2T'); 
    load('FM_SA_GlowRatio_3T_2.mat','B_high_SABP_times10_3T','B_low_SABP_times10_3T','ratio_HFHF_hoksok3T','ratio_LFHF_hoksok3T','ratio_sensitivity3T'); 
    load('FM_SA_GlowRatio_4T_2.mat','B_high_SABP_times10_4T','B_low_SABP_times10_4T','ratio_HFHF_hoksok4T','ratio_LFHF_hoksok4T','ratio_sensitivity4T'); 
    load('FM_SA_GlowRatio_5T_2.mat','B_high_SABP_times10_5T','B_low_SABP_times10_5T','ratio_HFHF_hoksok5T','ratio_LFHF_hoksok5T','ratio_sensitivity5T'); 
    load('FM_SA_GlowRatio_6T_2.mat','B_high_SABP_times10_6T','B_low_SABP_times10_6T','ratio_HFHF_hoksok6T','ratio_LFHF_hoksok6T','ratio_sensitivity6T'); 
    load('FM_SA_GlowRatio_7T_2.mat','B_high_SABP_times10_7T','B_low_SABP_times10_7T','ratio_HFHF_hoksok7T','ratio_LFHF_hoksok7T','ratio_sensitivity7T'); 
    load('FM_SA_GlowRatio_8T_2.mat','B_high_SABP_times10_8T','B_low_SABP_times10_8T','ratio_HFHF_hoksok8T','ratio_LFHF_hoksok8T','ratio_sensitivity8T'); 
    load('FM_SA_GlowRatio_9T_2.mat','B_high_SABP_times10_9T','B_low_SABP_times10_9T','ratio_HFHF_hoksok9T','ratio_LFHF_hoksok9T','ratio_sensitivity9T'); 
    load('FM_SA_GlowRatio_10T_2.mat','B_high_SABP_times10_10T','B_low_SABP_times10_10T','ratio_HFHF_hoksok10T','ratio_LFHF_hoksok10T','ratio_sensitivity10T'); 
    load('FM_SA_GlowRatio_11T_2.mat','B_high_SABP_times10_11T','B_low_SABP_times10_11T','ratio_HFHF_hoksok11T','ratio_LFHF_hoksok11T','ratio_sensitivity11T'); 
    load('FM_SA_GlowRatio_12T_2.mat','B_high_SABP_times10_12T','B_low_SABP_times10_12T','ratio_HFHF_hoksok12T','ratio_LFHF_hoksok12T','ratio_sensitivity12T');
    load('FM_SA_GlowRatio_13T_2.mat','B_high_SABP_times10_13T','B_low_SABP_times10_13T','ratio_HFHF_hoksok13T','ratio_LFHF_hoksok13T','ratio_sensitivity13T');
    load('FM_SA_GlowRatio_14T_2.mat','B_high_SABP_times10_14T','B_low_SABP_times10_14T','ratio_HFHF_hoksok14T','ratio_LFHF_hoksok14T','ratio_sensitivity14T');
    load('FM_SA_GlowRatio_15T_2.mat','B_high_SABP_times10_15T','B_low_SABP_times10_15T','ratio_HFHF_hoksok15T','ratio_LFHF_hoksok15T','ratio_sensitivity15T'); 
    load('FM_SA_GlowRatio_16T_2.mat','B_high_SABP_times10_16T','B_low_SABP_times10_16T','ratio_HFHF_hoksok16T','ratio_LFHF_hoksok16T','ratio_sensitivity16T');
    load('FM_SA_GlowRatio_17T_2.mat','B_high_SABP_times10_17T','B_low_SABP_times10_17T','ratio_HFHF_hoksok17T','ratio_LFHF_hoksok17T','ratio_sensitivity17T');
    load('FM_SA_GlowRatio_18T_2.mat','B_high_SABP_times10_18T','B_low_SABP_times10_18T','ratio_HFHF_hoksok18T','ratio_LFHF_hoksok18T','ratio_sensitivity18T');
    load('FM_SA_GlowRatio_19T_2.mat','B_high_SABP_times10_19T','B_low_SABP_times10_19T','ratio_HFHF_hoksok19T','ratio_LFHF_hoksok19T','ratio_sensitivity19T');
    load('FM_SA_GlowRatio_20T_2.mat','B_high_SABP_times10_20T','B_low_SABP_times10_20T','ratio_HFHF_hoksok20T','ratio_LFHF_hoksok20T','ratio_sensitivity20T');

    load('FM_SA_GlowRatio_1D_2.mat','B_high_SABP_div10_1D','B_low_SABP_div10_1D','ratio_HFHF_hoksok1D','ratio_LFHF_hoksok1D','ratio_sensitivity1D'); 
    load('FM_SA_GlowRatio_2D_2.mat','B_high_SABP_div10_2D','B_low_SABP_div10_2D','ratio_HFHF_hoksok2D','ratio_LFHF_hoksok2D','ratio_sensitivity2D'); 
    load('FM_SA_GlowRatio_3D_2.mat','B_high_SABP_div10_3D','B_low_SABP_div10_3D','ratio_HFHF_hoksok3D','ratio_LFHF_hoksok3D','ratio_sensitivity3D'); 
    load('FM_SA_GlowRatio_4D_2.mat','B_high_SABP_div10_4D','B_low_SABP_div10_4D','ratio_HFHF_hoksok4D','ratio_LFHF_hoksok4D','ratio_sensitivity4D'); 
    load('FM_SA_GlowRatio_5D_2.mat','B_high_SABP_div10_5D','B_low_SABP_div10_5D','ratio_HFHF_hoksok5D','ratio_LFHF_hoksok5D','ratio_sensitivity5D'); 
    load('FM_SA_GlowRatio_6D_2.mat','B_high_SABP_div10_6D','B_low_SABP_div10_6D','ratio_HFHF_hoksok6D','ratio_LFHF_hoksok6D','ratio_sensitivity6D'); 
    load('FM_SA_GlowRatio_7D_2.mat','B_high_SABP_div10_7D','B_low_SABP_div10_7D','ratio_HFHF_hoksok7D','ratio_LFHF_hoksok7D','ratio_sensitivity7D'); 
    load('FM_SA_GlowRatio_8D_2.mat','B_high_SABP_div10_8D','B_low_SABP_div10_8D','ratio_HFHF_hoksok8D','ratio_LFHF_hoksok8D','ratio_sensitivity8D'); 
    load('FM_SA_GlowRatio_9D_2.mat','B_high_SABP_div10_9D','B_low_SABP_div10_9D','ratio_HFHF_hoksok9D','ratio_LFHF_hoksok9D','ratio_sensitivity9D'); 
    load('FM_SA_GlowRatio_10D_2.mat','B_high_SABP_div10_10D','B_low_SABP_div10_10D','ratio_HFHF_hoksok10D','ratio_LFHF_hoksok10D','ratio_sensitivity10D'); 
    load('FM_SA_GlowRatio_11D_2.mat','B_high_SABP_div10_11D','B_low_SABP_div10_11D','ratio_HFHF_hoksok11D','ratio_LFHF_hoksok11D','ratio_sensitivity11D'); 
    load('FM_SA_GlowRatio_12D_2.mat','B_high_SABP_div10_12D','B_low_SABP_div10_12D','ratio_HFHF_hoksok12D','ratio_LFHF_hoksok12D','ratio_sensitivity12D');
    load('FM_SA_GlowRatio_13D_2.mat','B_high_SABP_div10_13D','B_low_SABP_div10_13D','ratio_HFHF_hoksok13D','ratio_LFHF_hoksok13D','ratio_sensitivity13D');
    load('FM_SA_GlowRatio_14D_2.mat','B_high_SABP_div10_14D','B_low_SABP_div10_14D','ratio_HFHF_hoksok14D','ratio_LFHF_hoksok14D','ratio_sensitivity14D');
    load('FM_SA_GlowRatio_15D_2.mat','B_high_SABP_div10_15D','B_low_SABP_div10_15D','ratio_HFHF_hoksok15D','ratio_LFHF_hoksok15D','ratio_sensitivity15D'); 
    load('FM_SA_GlowRatio_16D_2.mat','B_high_SABP_div10_16D','B_low_SABP_div10_16D','ratio_HFHF_hoksok16D','ratio_LFHF_hoksok16D','ratio_sensitivity16D');
    load('FM_SA_GlowRatio_17D_2.mat','B_high_SABP_div10_17D','B_low_SABP_div10_17D','ratio_HFHF_hoksok17D','ratio_LFHF_hoksok17D','ratio_sensitivity17D');
    load('FM_SA_GlowRatio_18D_2.mat','B_high_SABP_div10_18D','B_low_SABP_div10_18D','ratio_HFHF_hoksok18D','ratio_LFHF_hoksok18D','ratio_sensitivity18D');
    load('FM_SA_GlowRatio_19D_2.mat','B_high_SABP_div10_19D','B_low_SABP_div10_19D','ratio_HFHF_hoksok19D','ratio_LFHF_hoksok19D','ratio_sensitivity19D');
    load('FM_SA_GlowRatio_20D_2.mat','B_high_SABP_div10_20D','B_low_SABP_div10_20D','ratio_HFHF_hoksok20D','ratio_LFHF_hoksok20D','ratio_sensitivity20D');

    Nr = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20];
    x1 = [B_high_SABP_times10_1T(:,13), B_high_SABP_times10_2T(:,13), B_high_SABP_times10_3T(:,13), B_high_SABP_times10_4T(:,13), B_high_SABP_times10_5T(:,13), B_high_SABP_times10_6T(:,13), B_high_SABP_times10_7T(:,13), B_high_SABP_times10_8T(:,13), B_high_SABP_times10_9T(:,13), B_high_SABP_times10_10T(:,13), B_high_SABP_times10_11T(:,13), B_high_SABP_times10_12T(:,13), B_high_SABP_times10_13T(:,13), B_high_SABP_times10_14T(:,13), B_high_SABP_times10_15T(:,13), B_high_SABP_times10_16T(:,13), B_high_SABP_times10_17T(:,13), B_high_SABP_times10_18T(:,13), B_high_SABP_times10_19T(:,13), B_high_SABP_times10_20T(:,13)];
    x2 = [B_high_SABP_times10_1T(:,11), B_high_SABP_times10_2T(:,11), B_high_SABP_times10_3T(:,11), B_high_SABP_times10_4T(:,11), B_high_SABP_times10_5T(:,11), B_high_SABP_times10_6T(:,11), B_high_SABP_times10_7T(:,11), B_high_SABP_times10_8T(:,11), B_high_SABP_times10_9T(:,11), B_high_SABP_times10_10T(:,11), B_high_SABP_times10_11T(:,11), B_high_SABP_times10_12T(:,11), B_high_SABP_times10_13T(:,11), B_high_SABP_times10_14T(:,11), B_high_SABP_times10_15T(:,11), B_high_SABP_times10_16T(:,11), B_high_SABP_times10_17T(:,11), B_high_SABP_times10_18T(:,11), B_high_SABP_times10_19T(:,11), B_high_SABP_times10_20T(:,11)];
    x3 = [B_low_SABP_times10_1T(:,11), B_low_SABP_times10_2T(:,11), B_low_SABP_times10_3T(:,11), B_low_SABP_times10_4T(:,11), B_low_SABP_times10_5T(:,11), B_low_SABP_times10_6T(:,11), B_low_SABP_times10_7T(:,11), B_low_SABP_times10_8T(:,11), B_low_SABP_times10_9T(:,11), B_low_SABP_times10_10T(:,11), B_low_SABP_times10_11T(:,11), B_low_SABP_times10_12T(:,11), B_low_SABP_times10_13T(:,11), B_low_SABP_times10_14T(:,11), B_low_SABP_times10_15T(:,11), B_low_SABP_times10_16T(:,11), B_low_SABP_times10_17T(:,11), B_low_SABP_times10_18T(:,11), B_low_SABP_times10_19T(:,11), B_low_SABP_times10_20T(:,11)];
    x4 = [ratio_HFHF_hoksok1T, ratio_HFHF_hoksok2T, ratio_HFHF_hoksok3T, ratio_HFHF_hoksok4T, ratio_HFHF_hoksok5T, ratio_HFHF_hoksok6T, ratio_HFHF_hoksok7T, ratio_HFHF_hoksok8T, ratio_HFHF_hoksok9T, ratio_HFHF_hoksok10T, ratio_HFHF_hoksok11T, ratio_HFHF_hoksok12T, ratio_HFHF_hoksok13T, ratio_HFHF_hoksok14T, ratio_HFHF_hoksok15T, ratio_HFHF_hoksok16T, ratio_HFHF_hoksok17T, ratio_HFHF_hoksok18T, ratio_HFHF_hoksok19T, ratio_HFHF_hoksok20T];
    x5 = [ratio_LFHF_hoksok1T, ratio_LFHF_hoksok2T, ratio_LFHF_hoksok3T, ratio_LFHF_hoksok4T, ratio_LFHF_hoksok5T, ratio_LFHF_hoksok6T, ratio_LFHF_hoksok7T, ratio_LFHF_hoksok8T, ratio_LFHF_hoksok9T, ratio_LFHF_hoksok10T, ratio_LFHF_hoksok11T, ratio_LFHF_hoksok12T, ratio_LFHF_hoksok13T, ratio_LFHF_hoksok14T, ratio_LFHF_hoksok15T, ratio_LFHF_hoksok16T, ratio_LFHF_hoksok17T, ratio_LFHF_hoksok18T, ratio_LFHF_hoksok19T, ratio_LFHF_hoksok20T]; 
    x6 = [ratio_sensitivity1T, ratio_sensitivity2T, ratio_sensitivity3T, ratio_sensitivity4T, ratio_sensitivity5T, ratio_sensitivity6T, ratio_sensitivity7T, ratio_sensitivity8T, ratio_sensitivity9T, ratio_sensitivity10T, ratio_sensitivity11T, ratio_sensitivity12T, ratio_sensitivity13T, ratio_sensitivity14T, ratio_sensitivity15T, ratio_sensitivity16T, ratio_sensitivity17T, ratio_sensitivity18T, ratio_sensitivity19T, ratio_sensitivity20T];
            
    % Choose how to plot the senstivitiy analysis (times 10)
        % x1=new sensitivity
        % x2=hoksok_ratio_[F]high_new
        % x3=hoksok_ratio_[F]low_new
        % x4=hoksok_ratio_[F]high_new/hoksok_ratio_[F]high_old
        % x5=hoksok_ratio_[F]low_new/hoksok_ratio_[F]high_old 
        % x6=sensitivity_new/sensitivity_old
        c = x1;
        
    figure(2);      
    z = [c];
    %boxplot(z,Nr,'labels',{'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'},'Colors',[0 0.2 0.7410],'Symbol','+b');
    boxplot(z,Nr,'labels',{'α1','α2','α3','α4','α5','β1','β2','β3','β4','β5','β6','β7','γ1','γ2','K-1','K1','K3','K4','K-4','K6'},'Colors',[0 0.2 0.7410],'Symbol','+b');
    set(gca,'yscale','log','Fontsize',14);
    ylabel('(new hok/sok ratio [FA]_l_o_w)/(old hok/sok ratio [FA]_h_i_g_h)','Fontsize',16);
    %ylabel('new sensitivity/old sensitivity','Fontsize',16);
    xlabel('parameter','Fontsize',16);
    title('Sensitivity analysis - times 10','Fontsize',22);
       
    y1 = [B_high_SABP_div10_1D(:,13), B_high_SABP_div10_2D(:,13), B_high_SABP_div10_3D(:,13), B_high_SABP_div10_4D(:,13), B_high_SABP_div10_5D(:,13), B_high_SABP_div10_6D(:,13), B_high_SABP_div10_7D(:,13), B_high_SABP_div10_8D(:,13), B_high_SABP_div10_9D(:,13), B_high_SABP_div10_10D(:,13), B_high_SABP_div10_11D(:,13), B_high_SABP_div10_12D(:,13), B_high_SABP_div10_13D(:,13), B_high_SABP_div10_14D(:,13), B_high_SABP_div10_15D(:,13), B_high_SABP_div10_16D(:,13), B_high_SABP_div10_17D(:,13), B_high_SABP_div10_18D(:,13), B_high_SABP_div10_19D(:,13), B_high_SABP_div10_20D(:,13)];
    y2 = [B_high_SABP_div10_1D(:,11), B_high_SABP_div10_2D(:,11), B_high_SABP_div10_3D(:,11), B_high_SABP_div10_4D(:,11), B_high_SABP_div10_5D(:,11), B_high_SABP_div10_6D(:,11), B_high_SABP_div10_7D(:,11), B_high_SABP_div10_8D(:,11), B_high_SABP_div10_9D(:,11), B_high_SABP_div10_10D(:,11), B_high_SABP_div10_11D(:,11), B_high_SABP_div10_12D(:,11), B_high_SABP_div10_13D(:,11), B_high_SABP_div10_14D(:,11), B_high_SABP_div10_15D(:,11), B_high_SABP_div10_16D(:,11), B_high_SABP_div10_17D(:,11), B_high_SABP_div10_18D(:,11), B_high_SABP_div10_19D(:,11), B_high_SABP_div10_20D(:,11)];
    y3 = [B_low_SABP_div10_1D(:,11), B_low_SABP_div10_2D(:,11), B_low_SABP_div10_3D(:,11), B_low_SABP_div10_4D(:,11), B_low_SABP_div10_5D(:,11), B_low_SABP_div10_6D(:,11), B_low_SABP_div10_7D(:,11), B_low_SABP_div10_8D(:,11), B_low_SABP_div10_9D(:,11), B_low_SABP_div10_10D(:,11), B_low_SABP_div10_11D(:,11), B_low_SABP_div10_12D(:,11), B_low_SABP_div10_13D(:,11), B_low_SABP_div10_14D(:,11), B_low_SABP_div10_15D(:,11), B_low_SABP_div10_16D(:,11), B_low_SABP_div10_17D(:,11), B_low_SABP_div10_18D(:,11), B_low_SABP_div10_19D(:,11), B_low_SABP_div10_20D(:,11)];
    y4 = [ratio_HFHF_hoksok1D, ratio_HFHF_hoksok2D, ratio_HFHF_hoksok3D, ratio_HFHF_hoksok4D, ratio_HFHF_hoksok5D, ratio_HFHF_hoksok6D, ratio_HFHF_hoksok7D, ratio_HFHF_hoksok8D, ratio_HFHF_hoksok9D, ratio_HFHF_hoksok10D, ratio_HFHF_hoksok11D, ratio_HFHF_hoksok12D, ratio_HFHF_hoksok13D, ratio_HFHF_hoksok14D, ratio_HFHF_hoksok15D, ratio_HFHF_hoksok16D, ratio_HFHF_hoksok17D, ratio_HFHF_hoksok18D, ratio_HFHF_hoksok19D, ratio_HFHF_hoksok20D];
    y5 = [ratio_LFHF_hoksok1D, ratio_LFHF_hoksok2D, ratio_LFHF_hoksok3D, ratio_LFHF_hoksok4D, ratio_LFHF_hoksok5D, ratio_LFHF_hoksok6D, ratio_LFHF_hoksok7D, ratio_LFHF_hoksok8D, ratio_LFHF_hoksok9D, ratio_LFHF_hoksok10D, ratio_LFHF_hoksok11D, ratio_LFHF_hoksok12D, ratio_LFHF_hoksok13D, ratio_LFHF_hoksok14D, ratio_LFHF_hoksok15D, ratio_LFHF_hoksok16D, ratio_LFHF_hoksok17D, ratio_LFHF_hoksok18D, ratio_LFHF_hoksok19D, ratio_LFHF_hoksok20D]; 
    y6 = [ratio_sensitivity1D, ratio_sensitivity2D, ratio_sensitivity3D, ratio_sensitivity4D, ratio_sensitivity5D, ratio_sensitivity6D, ratio_sensitivity7D, ratio_sensitivity8D, ratio_sensitivity9D, ratio_sensitivity10D, ratio_sensitivity11D, ratio_sensitivity12D, ratio_sensitivity13D, ratio_sensitivity14D, ratio_sensitivity15D, ratio_sensitivity16D, ratio_sensitivity17D, ratio_sensitivity18D, ratio_sensitivity19D, ratio_sensitivity20D];
    
    % Choose how to plot the senstivitiy analysis (div 10)
        % y1=new sensitivity
        % y2=hoksok_ratio_[F]high_new
        % y3=hoksok_ratio_[F]low_new
        % y4=hoksok_ratio_[F]high_new/hoksok_ratio_[F]high_old
        % y5=hoksok_ratio_[F]low_new/hoksok_ratio_[F]high_old 
        % y6=sensitivity_new/sensitivity_old
        d = y1;
    
    figure(3);
    z = [d];
    %boxplot(z,Nr,'labels',{'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'},'Colors',[0 0.2 0.7410],'Symbol','+b');
    boxplot(z,Nr,'labels',{'α1','α2','α3','α4','α5','β1','β2','β3','β4','β5','β6','β7','γ1','γ2','K-1','K1','K3','K4','K-4','K6'},'Colors',[0 0.2 0.7410],'Symbol','+b');
    set(gca,'yscale','log','Fontsize',14);
    ylabel('(new hok/sok ratio [FA]_l_o_w)/(old hok/sok ratio [FA]_h_i_g_h)','Fontsize',16);
    %ylabel('new sensitivity/old sensitivity','Fontsize',16);
    xlabel('parameter','Fontsize',16);
    title('Sensitivity analysis - div 10','Fontsize',22);
    
%% 6. Select outliers 
load('FM_SA_GlowRatio_7D_2.mat','B_high_SABP_div10_7D','B_low_SABP_div10_7D','ratio_HFHF_hoksok7D','ratio_LFHF_hoksok7D','ratio_sensitivity7D'); 
load('FM_SA_GlowRatio_1T_2.mat','B_high_SABP_times10_1T','B_low_SABP_times10_1T','ratio_HFHF_hoksok1T','ratio_LFHF_hoksok1T','ratio_sensitivity1T'); 
load('FM_SA_GlowRatio_16T_2.mat','B_high_SABP_times10_16T','B_low_SABP_times10_16T','ratio_HFHF_hoksok16T','ratio_LFHF_hoksok16T','ratio_sensitivity16T');
load('outputFM_SAGlowRatio_times10_16.mat','sim')

% give numbers to the parameter sets and sort them to select the set that
% gives the biggest increase in sensitivity
ratio_LFHF_hoksok7D_nr = [GP_48h ratio_LFHF_hoksok7D];
sort_ratio_LFHF_hoksok7D_nr = sortrows(ratio_LFHF_hoksok7D_nr,2);
ratio_LFHF_hoksok1T_nr = [GP_48h ratio_LFHF_hoksok1T];
sort_ratio_LFHF_hoksok1T_nr = sortrows(ratio_LFHF_hoksok1T_nr,2);
ratio_sensitivity16T_nr = [GP_48h ratio_sensitivity16T];
sort_ratio_sensitivity16T_nr = sortrows(ratio_sensitivity16T_nr,2);