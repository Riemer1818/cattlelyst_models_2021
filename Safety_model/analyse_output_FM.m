%% Goal: make time-series plots of the original biosafety circuit
% contains:
% 1. Combine data from simulations
% 2. remove NaN parameter sets
% 3. Select parameters which give a high or low hoksok ratio
% 4. Boxplots parameter values giving low or high hok/sok ratio
% 5. scatter plots parameter value vs FA senstivity (all good p-sets)

% input files;
%    FM_output_*.mat
%    FM_sim_high_low.mat
%    FM_Gsim_high_low.mat
%    median_low_hoksokratio2
% output files;
%    FM_sim_high_low.mat
%    FM_Gsim_high_low.mat
%    median_low_hoksokratio2

%% 1. Combine data from simulations
load('FM_output_1.mat')
O1 = [[1:10000]' params];
S1 = [sim];
load('FM_output_2.mat')
O2 = [[10001:20000]' params];
S2 = [sim];
load('FM_output_3.mat')
O3 = [[20001:30000]' params];
S3 = [sim];
load('FM_output_4.mat')
O4 = [[30001:40000]' params];
S4 = [sim];
load('FM_output_5.mat')
O5 = [[40001:50000]' params];
S5 = [sim];
load('FM_output_6.mat')
O6 = [[50001:60000]' params];
S6 = [sim];
load('FM_output_7.mat')
O7 = [[60001:70000]' params];
S7 = [sim];
load('FM_output_8.mat')
O8 = [[70001:80000]' params];
S8 = [sim];
load('FM_output_9.mat')
O9 = [[80001:90000]' params];
S9 = [sim];
load('FM_output_10.mat')
O10 = [[90001:100000]' params];
S10 = [sim];

%combine files
params_all = [O1; O2; O3; O4; O5; O6; O7; O8; O9; O10];
sim_all = [S1 S2 S3 S4 S5 S6 S7 S8 S9 S10];

% make a matrix with end concentrations.
low_F = zeros(100000, 9); 
high_F= zeros(100000, 9);
for i=1:length(sim_all)
     low_F(i,1:length(sim_all{i}{1}(end,:))) = sim_all{i}{1}(end,:);
     high_F(i,1:length(sim_all{i}{2}(end,:))) = sim_all{i}{2}(end,:);            
end
A_high_F = [[1:100000]' high_F (high_F(:,9)+high_F(:,8))./(high_F(:,8)+high_F(:,7))];
A_low_F = [[1:100000]' low_F (low_F(:,9)+low_F(:,8))./(low_F(:,8)+low_F(:,7))];
save('FM_sim_high_low.mat','A_high_F','A_low_F','params_all');
%% 2. remove NaN parameter sets

load('FM_sim_high_low.mat','A_high_F','A_low_F','params_all')

B_high = [A_high_F A_low_F(:,11)-A_high_F(:,11) A_low_F(:,11)./A_high_F(:,11)];
B_low = [A_low_F A_low_F(:,11)-A_high_F(:,11) A_low_F(:,11)./A_high_F(:,11)];

%removing NaN and negative simulations from B_low, B_high and params_all

sort_B_high = sortrows(B_high,13);
sort_B_low = sortrows(B_low,13);
% NaN_p = sort_B_high(99119:100000,1);
%sort_C_high = sortrows(B_high,8);
%sort_C_low = sortrows(B_low,8);
% N_p_high = sort_C_high(1:86,1);
% N_p_low = sort_C_low(1:430,1);
% BP_all = [NaN_p; N_p_high; N_p_low];
% BP = unique(BP_all);
Gsim_high = sortrows(B_high,1);
Gsim_low = sortrows(B_low,1);
% Gsim_high(BP,:)=[];
% Gsim_low(BP,:)=[];
Gparams = params_all;
% Gparams(BP,:)=[]; 
sort_Gsim_high = sortrows(Gsim_high,13);
sort_Gsim_low = sortrows(Gsim_low,13);
save('FM_Gsim_high_low.mat','Gsim_high','Gsim_low','Gparams');

%% 3. Select parameters which give a high or low hoksok ratio

load('FM_Gsim_high_low.mat','Gsim_high','Gsim_low','Gparams');
%first select the parameters that give a low hoksok ratio with high
%formaldehyde concentrations (cut off ratio 15-25)
sort_Gsim_high = sortrows(Gsim_high,11);
p_sets_GhighF_ratio = sort_Gsim_high(51742:56223,1);
Gparams_ratio = [Gparams Gsim_high(:,13)];
params_GhighF_ratio = Gparams_ratio(p_sets_GhighF_ratio,2:24);
sim_GhighF_ratio = Gsim_high(p_sets_GhighF_ratio,:);
%select parameters that give a high hoksok ratio with low formaldehyde
%concentrations (cut off ratio 5000-7000)
sort_Gsim_low = sortrows(Gsim_low,11);
p_sets_GlowF_ratio = sort_Gsim_low(97423:98913,1);
params_GlowF_ratio = Gparams_ratio(p_sets_GlowF_ratio,2:24);

%select median parameter values
median_lowhoksok = median(params_GhighF_ratio);
median_highhoksok = median(params_GlowF_ratio);
N = 1;
q = 1;
[sim] = optimise_model(median_lowhoksok(1:22),N,q);
for i = 1:23
    p_values(i) = ranksum(params_GhighF_ratio(:,i),params_GlowF_ratio(:,i));
end
for i=1:length(sim)
     low_F(i,1:length(sim{i}{1}(end,:))) = sim{i}{1}(end,:);
     high_F(i,1:length(sim{i}{2}(end,:))) = sim{i}{2}(end,:);            
end
A_high_F = [[1:1]' high_F (high_F(:,9)+high_F(:,8))./(high_F(:,8)+high_F(:,7))];
A_low_F = [[1:1]' low_F (low_F(:,9)+low_F(:,8))./(low_F(:,8)+low_F(:,7))];
B_high_F = [A_high_F A_low_F(:,11)-A_high_F(:,11) A_low_F(:,11)./A_high_F(:,11)];
B_low_F = [A_low_F A_low_F(:,11)-A_high_F(:,11) A_low_F(:,11)./A_high_F(:,11)];

save('median_low_hoksokratio2','sim_GhighF_ratio','params_GhighF_ratio','params_GlowF_ratio','median_lowhoksok','median_highhoksok','p_values','B_high_F','B_low_F','Gparams_ratio','sim')

%% 4. Boxplots parameter values giving low or high hok/sok ratio
load('median_low_hoksokratio2','params_GhighF_ratio','params_GlowF_ratio')
%Boxplots that visualize the distribution of the parameter values giving a
%low or high hok/sok ratio
figure(1);
subplot(2,4,1);
x = params_GlowF_ratio(:,1);
y = params_GhighF_ratio(:,1);
z = [x; y];
g = [zeros(length(x),1); ones(length(y),1)];
boxplot(z,g,'labels',{'high hoksok ratio','low hoksok ratio'});
set(gca,'yscale','log');
title('alpha T7 (a1)');
subplot(2,4,2);
x = params_GlowF_ratio(:,2);
y = params_GhighF_ratio(:,2);
z = [x; y];
g = [zeros(length(x),1); ones(length(y),1)];
boxplot(z,g,'labels',{'high hoksok ratio','low hoksok ratio'});
set(gca,'yscale','log');
title('alpha frm leakage (a2)');
subplot(2,4,3);
x = params_GlowF_ratio(:,3);
y = params_GhighF_ratio(:,3);
z = [x; y];
g = [zeros(length(x),1); ones(length(y),1)];
boxplot(z,g,'labels',{'high hoksok ratio','low hoksok ratio'});
set(gca,'yscale','log');
title('alpha frm max (a3)');
subplot(2,4,4);
x = params_GlowF_ratio(:,4);
y = params_GhighF_ratio(:,4);
z = [x; y];
g = [zeros(length(x),1); ones(length(y),1)];
boxplot(z,g,'labels',{'high hoksok ratio','low hoksok ratio'});
set(gca,'yscale','log');
title('lac promoter leakage (a4)');
subplot(2,4,5);
x = params_GlowF_ratio(:,5);
y = params_GhighF_ratio(:,5);
z = [x; y];
g = [zeros(length(x),1); ones(length(y),1)];
boxplot(z,g,'labels',{'high hoksok ratio','low hoksok ratio'});
set(gca,'yscale','log');
title('lac max transcription rate (a5)');
subplot(2,4,6);
x = params_GlowF_ratio(:,6);
y = params_GhighF_ratio(:,6);
z = [x; y];
g = [zeros(length(x),1); ones(length(y),1)];
boxplot(z,g,'labels',{'high hoksok ratio','low hoksok ratio'});
set(gca,'yscale','log');
title('degradation rate mFrmR (b1)');
subplot(2,4,7);
x = params_GlowF_ratio(:,7);
y = params_GhighF_ratio(:,7);
z = [x; y];
g = [zeros(length(x),1); ones(length(y),1)];
boxplot(z,g,'labels',{'high hoksok ratio','low hoksok ratio'});
set(gca,'yscale','log');
title('degradation rate FrmR (b2)');
subplot(2,4,8);
x = params_GlowF_ratio(:,8);
y = params_GhighF_ratio(:,8);
z = [x; y];
g = [zeros(length(x),1); ones(length(y),1)];
boxplot(z,g,'labels',{'high hoksok ratio','low hoksok ratio'});
set(gca,'yscale','log');
title('degradation rate mLacI (b3)');
sgtitle('parameter distribution giving a low (15-25) or high (5000-7000) ratio');

figure(2);
subplot(2,4,1);
x = params_GlowF_ratio(:,9);
y = params_GhighF_ratio(:,9);
z = [x; y];
g = [zeros(length(x),1); ones(length(y),1)];
boxplot(z,g,'labels',{'high hoksok ratio','low hoksok ratio'});
set(gca,'yscale','log');
title('degradation rate LacI (b4)');
subplot(2,4,2);
x = params_GlowF_ratio(:,10);
y = params_GhighF_ratio(:,10);
z = [x; y];
g = [zeros(length(x),1); ones(length(y),1)];
boxplot(z,g,'labels',{'high hoksok ratio','low hoksok ratio'});
set(gca,'yscale','log');
title('degradation rate sok (b5)');
subplot(2,4,3);
x = params_GlowF_ratio(:,11);
y = params_GhighF_ratio(:,11);
z = [x; y];
g = [zeros(length(x),1); ones(length(y),1)];
boxplot(z,g,'labels',{'high hoksok ratio','low hoksok ratio'});
set(gca,'yscale','log');
title('degradation rate hoksok (b6)');
subplot(2,4,4);
x = params_GlowF_ratio(:,12);
y = params_GhighF_ratio(:,12);
z = [x; y];
g = [zeros(length(x),1); ones(length(y),1)];
boxplot(z,g,'labels',{'high hoksok ratio','low hoksok ratio'});
set(gca,'yscale','log');
title('degradation rate hok (b7)');
subplot(2,4,5);
x = params_GlowF_ratio(:,13);
y = params_GhighF_ratio(:,13);
z = [x; y];
g = [zeros(length(x),1); ones(length(y),1)];
boxplot(z,g,'labels',{'high hoksok ratio','low hoksok ratio'});
set(gca,'yscale','log');
title('translation rate mfrmR (g1)');
subplot(2,4,6);
x = params_GlowF_ratio(:,14);
y = params_GhighF_ratio(:,14);
z = [x; y];
g = [zeros(length(x),1); ones(length(y),1)];
boxplot(z,g,'labels',{'high hoksok ratio','low hoksok ratio'});
set(gca,'yscale','log');
title('translation rate mlacI (g2)');
subplot(2,4,7);
x = params_GlowF_ratio(:,15);
y = params_GhighF_ratio(:,15);
z = [x; y];
g = [zeros(length(x),1); ones(length(y),1)];
boxplot(z,g,'labels',{'high hoksok ratio','low hoksok ratio'});
set(gca,'yscale','log');
title('unbinding formaldehyde-frmR complex (k1)');
subplot(2,4,8);
x = params_GlowF_ratio(:,16);
y = params_GhighF_ratio(:,16);
z = [x; y];
g = [zeros(length(x),1); ones(length(y),1)];
boxplot(z,g,'labels',{'high hoksok ratio','low hoksok ratio'});
set(gca,'yscale','log');
title('binding formaldehyde and frmR (k2)');
sgtitle('parameter distribution giving a low (15-25) or high (5000-7000) ratio');

figure(3);
subplot(2,3,1);
x = params_GlowF_ratio(:,17);
y = params_GhighF_ratio(:,17);
z = [x; y];
g = [zeros(length(x),1); ones(length(y),1)];
boxplot(z,g,'labels',{'high hoksok ratio','low hoksok ratio'});
set(gca,'yscale','log');
title('binding frmR in hill equation (k3)');
subplot(2,3,2);
x = params_GlowF_ratio(:,18);
y = params_GhighF_ratio(:,18);
z = [x; y];
g = [zeros(length(x),1); ones(length(y),1)];
boxplot(z,g,'labels',{'high hoksok ratio','low hoksok ratio'});
set(gca,'yscale','log');
title('binding hok and sok (k4)');
subplot(2,3,3);
x = params_GlowF_ratio(:,19);
y = params_GhighF_ratio(:,19);
z = [x; y];
g = [zeros(length(x),1); ones(length(y),1)];
boxplot(z,g,'labels',{'high hoksok ratio','low hoksok ratio'});
set(gca,'yscale','log');
title('unbinding hok and sok (k5)');
subplot(2,3,4);
x = params_GlowF_ratio(:,20);
y = params_GhighF_ratio(:,20);
z = [x; y];
g = [zeros(length(x),1); ones(length(y),1)];
boxplot(z,g,'labels',{'high hoksok ratio','low hoksok ratio'});
set(gca,'yscale','log');
title('binding LacI in hill equation (k6)');
subplot(2,3,5);
x = params_GlowF_ratio(:,21);
y = params_GhighF_ratio(:,21);
z = [x; y];
g = [zeros(length(x),1); ones(length(y),1)];
boxplot(z,g,'labels',{'high hoksok ratio','low hoksok ratio'});
set(gca,'yscale','log');
title('hill equation frmR (n1)');
subplot(2,3,6);
x = params_GlowF_ratio(:,22);
y = params_GhighF_ratio(:,22);
z = [x; y];
g = [zeros(length(x),1); ones(length(y),1)];
boxplot(z,g,'labels',{'high hoksok ratio','low hoksok ratio'});
set(gca,'yscale','log');
title('hill equation lacI (n2)');
sgtitle('parameter distribution giving a low (15-25) or high (5000-7000) ratio');

%% 5. scatter plots parameter value vs FA senstivity (all good p-sets)
load('median_low_hoksokratio2','params_GhighF_ratio','median_lowhoksok','median_highhoksok','p_values','B_high_F','B_low_F','Gparams_ratio','sim')

figure(4);
subplot(3,3,1);
scatter(Gparams_ratio(:,24),Gparams_ratio(:,2));
set(gca,'yscale','log','FontSize',13);
xlabel('FA sensitivity');
ylabel('a1');
title('transcription rate Pconst.');
subplot(3,3,2);
scatter(Gparams_ratio(:,24),Gparams_ratio(:,3));
set(gca,'yscale','log','FontSize',13);
xlabel('FA sensitivity');
ylabel('a2');
title('leakage transcription rate Pfrm');
subplot(3,3,3);
scatter(Gparams_ratio(:,24),Gparams_ratio(:,4));
set(gca,'yscale','log','FontSize',13);
xlabel('FA sensitivity');
ylabel('a3');
title('max. transcription rate Pfrm');
subplot(3,3,4);
scatter(Gparams_ratio(:,24),Gparams_ratio(:,5));
set(gca,'yscale','log','FontSize',13);
xlabel('FA sensitivity');
ylabel('a4');
title('leakage transcription rate Plac');
subplot(3,3,5);
scatter(Gparams_ratio(:,24),Gparams_ratio(:,6));
set(gca,'yscale','log','FontSize',13);
xlabel('FA sensitivity');
ylabel('a5');
title('max. transcription rate Plac');
subplot(3,3,6);
scatter(Gparams_ratio(:,24),Gparams_ratio(:,7));
set(gca,'yscale','log','FontSize',13);
xlabel('FA sensitivity');
ylabel('b1');
title('degradation rate mFrmR');
subplot(3,3,7);
scatter(Gparams_ratio(:,24),Gparams_ratio(:,8));
set(gca,'yscale','log','FontSize',13);
xlabel('FA sensitivity');
ylabel('b2');
title('degradation rate FrmR');
subplot(3,3,8);
scatter(Gparams_ratio(:,24),Gparams_ratio(:,9));
set(gca,'yscale','log','FontSize',13);
xlabel('FA sensitivity');
ylabel('b3');
title('degradation rate mLacI');
subplot(3,3,9);
scatter(Gparams_ratio(:,24),Gparams_ratio(:,10));
set(gca,'yscale','log','FontSize',13);
xlabel('FA sensitivity');
ylabel('b4');
title('degradation rate LacI');
sgtitle('formaldehyde sensitivity vs parameter value','Fontsize',22);

figure(5);
subplot(3,3,1);
scatter(Gparams_ratio(:,24),Gparams_ratio(:,11));
set(gca,'yscale','log','FontSize',13);
xlabel('FA sensitivity');
ylabel('b5');
title('degradation rate sok');
subplot(3,3,2);
scatter(Gparams_ratio(:,24),Gparams_ratio(:,12));
set(gca,'yscale','log','FontSize',13);
xlabel('FA sensitivity');
ylabel('b6');
title('degradation rate hok-sok');
subplot(3,3,3);
scatter(Gparams_ratio(:,24),Gparams_ratio(:,13));
set(gca,'yscale','log','FontSize',13);
xlabel('FA sensitivity');
ylabel('b7');
title('degradation rate hok');
subplot(3,3,4);
scatter(Gparams_ratio(:,24),Gparams_ratio(:,14));
set(gca,'yscale','log','FontSize',13);
xlabel('FA sensitivity');
ylabel('g1');
title('translation rate mFrmR');
subplot(3,3,5);
scatter(Gparams_ratio(:,24),Gparams_ratio(:,15));
set(gca,'yscale','log','FontSize',13);
xlabel('FA sensitivity');
ylabel('g2');
title('translation rate mLacI');
subplot(3,3,6);
scatter(Gparams_ratio(:,24),Gparams_ratio(:,16));
set(gca,'yscale','log','FontSize',13);
xlabel('FA sensitivity');
ylabel('k1');
title('unbinding rate FA-FrmR complex');
subplot(3,3,7);
scatter(Gparams_ratio(:,24),Gparams_ratio(:,17));
set(gca,'yscale','log','FontSize',13);
xlabel('FA sensitivity');
ylabel('k2');
title('binding rate FA and FrmR');
subplot(3,3,8);
scatter(Gparams_ratio(:,24),Gparams_ratio(:,18));
set(gca,'yscale','log','FontSize',13);
xlabel('FA sensitivity');
ylabel('k3');
title('binding FrmR in Hill function');
subplot(3,3,9);
scatter(Gparams_ratio(:,24),Gparams_ratio(:,19));
set(gca,'yscale','log','FontSize',13);
xlabel('FA sensitivity');
ylabel('k4');
title('binding rate hok and sok');
sgtitle('formaldehyde sensitivity vs parameter value','Fontsize',22);

figure(6);
subplot(3,3,1);
scatter(Gparams_ratio(:,24),Gparams_ratio(:,20));
set(gca,'yscale','log','FontSize',13);
xlabel('FA sensitivity');
ylabel('k5');
title('unbinding rate hok-sok');
subplot(3,3,2);
scatter(Gparams_ratio(:,24),Gparams_ratio(:,21));
set(gca,'yscale','log','FontSize',13);
xlabel('FA sensitivity');
ylabel('k6');
title('binding LacI in Hill function');
subplot(3,3,3);
scatter(Gparams_ratio(:,24),Gparams_ratio(:,22));
set(gca,'yscale','log','FontSize',13);
xlabel('FA sensitivity');
ylabel('n1');
title('Hill coefficient FrmR');
subplot(3,3,4);
scatter(Gparams_ratio(:,24),Gparams_ratio(:,23));
set(gca,'yscale','log','FontSize',13);
xlabel('FA sensitivity');
ylabel('n2');
title('Hill coefficient LacI');
sgtitle('formaldehyde sensitivity vs parameter value','Fontsize',22);