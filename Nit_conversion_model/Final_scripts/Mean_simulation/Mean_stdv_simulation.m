%% script to compute figure mean +/- st.dev simulation of 25,50 or 100 best parameter sets. 
% choose number of 'good_parameters'
I = 25;
% load parameter sets
load(sprintf('good_p_sets_gas_%i.mat',I))

%load data1 for YZN-001, data2 for XL-2
%columns data1: time, ammonia concentration, nitrate concentration, nitrite
%concentration, population volume, nitrogen gas
load('data1.mat')
%columns data2: time, ammonia concentration, nitrate concentration, nitrite
%concentration, population volume
load('data2.mat')

%load volume data (11) - YZN-001 and (22) - XL-2
%columns t_Vc_mu_1: time, population volume, growth rate
load('time_Volume_growth_matrix_11.mat','t_Vc_mu_1')
%columns t_Vc_mu_2: time, population volume, growth rate
load('time_Volume_growth_matrix_22.mat','t_Vc_mu_2')

%re-simulate the system to save the interpolated output as stacked
%matrices.
[SIM_tot] = mean_HNAD_sim(data1,data2,t_Vc_mu_1,t_Vc_mu_2,good_p_sets,I);

save(sprintf('Simulation_%i_Test_best_parameters',I),'SIM_tot')
