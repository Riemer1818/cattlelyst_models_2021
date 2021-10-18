%% mean_HNAD_sim is a function that computes stacked simulations in the third dimension
% inputs:
% data1: data for YZN-001, columns: Time, NH3, NO3, NO2, Vc, N2
% data2: data for XL-2, columns: Time, NH3, NO3, NO2, Vc
% t_Vc_mu_1 : time_volume_mu matrix for YZn-001 condition
% t_Vc_mu_2: time_volume_mu matrix for Xl-2 condition
% p_sets: parameter sets
% I: number of simulations needed

%output: SIM_tot: stacked interpolated simulation data for both conditions


function [SIM_tot] = mean_HNAD_sim(data1,data2,t_Vc_mu_1,t_Vc_mu_2,p_sets,I)
% time vector for simulations
time1 = 0:0.1:4320; %72 hours for YZN-001
time2 = 0:0.1:2160; %31 hours for XL-2

% Scores
% Set up empty vectors/matrices to store scoresdata1
score = zeros(I,1);
score_ij = zeros(I,2);

% Start testing parameter sets (from 1 to i)

for i= 1:I
%% Simulate different datasets
for j =1:2
    if j == 1 
%% condition 1, YZN-001 for HNAD and ammonia growth
        % initial conditions for YZN-001
        % columns: NH3, NH2OH, NO, NO3, NO2, N2O, N2, NH3_ex, NO3_ex,
        % NO2_ex, NO_ex, N2O_ex, N2_ex
    ics_YZN_001 = [0 0 0 0 0 0 0 5.294 0 0 0 0 0];    
    %Medium volume for YZN-001 condition
        Vm(1) = 0.125; %L 
    % simulate system
    [t,y] = ode15s(@HNAD_sim_Vol_optimize,time1,ics_YZN_001,[],Vm(1),t_Vc_mu_1,p_sets(i,:));
    
    %if statement to circumvent errors arising from warning: ode15s solver time-step below min value
    if t(end) < 4320
    score_ij(i,j) = 10000000;
    sim1 = zeros(9,14);
    SIM1(:,:,i)= [sim1];
    
    else
    % Interpolate output to compare same time-points
    sim1 = interp1(time1,[t,real(y)],data1(:,1));
    %store simulations in 3D 
    SIM1(:,:,i)=[sim1];
    % score sim1 versus data1
    score_ij(i,j) = sum((data1(:,2)-sim1(:,9)).^2)+ sum((data1(:,3)-sim1(:,10)).^2)+ sum((data1(:,4)-sim1(:,11)).^2)+ sum((data1(:,6)-sim1(:,14)).^2);
    end
%% Data2 = Data XL2, aerobic denitrification + growth on nitrate y(9)
    elseif j == 2
        % initial conditions for XL-2
        % columns: NH3, NH2OH, NO, NO3, NO2, N2O, N2, NH3_ex, NO3_ex,
        % NO2_ex, NO_ex, N2O_ex, N2_ex
    ics_XL2 = [0 0 0 0 0 0 0 0 1.5396 0.012 0 0 0]; 
       %Medium volume 
        Vm(2) = 0.250; %L
        %simulate system
    [t,y] = ode15s(@HNAD_sim_Vol_optimize,time2,ics_XL2,[],Vm(2),t_Vc_mu_2,p_sets(i,:));
    
    %if statement to circumvent errors arising from warning: ode15s solver time-step below min value
    if t(end) < 2160
    score_ij(i,j) = 10000000;
    sim2 = zeros(9,14);
    SIM2(:,:,i)= [sim2];
    else
        
    % Interpolate output to compare same time-points
    sim2 = interp1(time2,[t,real(y)],data2(:,1));
    % store sim
    SIM2(:,:,i)= [sim2];
    % score sim1 versus data1 + sim2 versus data2
    score_ij(i,j) = score_ij(i,j) + sum((data2(:,2)-sim2(:,9)).^2)+sum((data2(:,3)-sim2(:,10)).^2)+sum((data2(:,4)-sim2(:,11)).^2);
    end
    
      %% Merge output simulations for each condition
      SIM_tot = cat(1,SIM1,SIM2);
        
end
end 
end
end


