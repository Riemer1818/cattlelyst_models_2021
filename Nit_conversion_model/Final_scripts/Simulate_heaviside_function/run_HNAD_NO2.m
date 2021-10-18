%% Function to run heaviside model to match data2 for strain XL-2

%n = amount of parameters to be tested, set to 1 mostly for manipulations
%k = specify number of test

%data1 = data HNAD YZN per minute, columns: (time, NH3, NO3, NO2, Vc, N2)
%data2 = data Aerobic denitrification XL2 per minute, columns: (time, NH3,
%NO3, NO2, Vc)

function run_HNAD_NO2(k,n)
%load best parameter set (1 by 36)

load('p_sets_best.mat','p_sets_SA')
% convert to linear scale
p_sets = 10.^(p_sets_SA);

%load data1 and data2
load('data1.mat')  %data for YZN-001
load('data2.mat')  %data for XL-2

[sim, score_ij, score, output] = optimise_HNAD_sim_NO2(data1,data2,p_sets(1,:),n);
save(sprintf("Output_heaviside_%i.mat",k),'sim','score_ij','score','output');
end

