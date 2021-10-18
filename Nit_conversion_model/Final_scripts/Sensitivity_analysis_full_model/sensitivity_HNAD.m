%% Function sensitivity_HNAD calls the HNAD model and computes the score based on parameter sets
%inputs: parameter sets, best parameter set
%Output: normalized score



function [score] = sensitivity_HNAD(p_sets)

%load datafiles
% YZN-001 (columns: time, NH3, NO3, NO2, Vc, N2)
load('data1.mat')
% XL-2 (columns: time, NH3, NO3, NO2, Vc)
load('data2.mat')


%time vector for simulations
time1 = 0:0.1:4320; %72 hours for YZN-001 condition
time2 = 0:0.1:2160; %31 hours for XL-2 condition

% Scores
% Set up empty vectors/matrices to store scores
% only one parameter set is tested such that I = 1
I = 1;
score = zeros(I,1);
score_ij = zeros(I,2);

for i= I
%% Simulate different conditions
for j =1:2
    if j == 1 
%% condition 1, YZN-001 for HNAD and ammonia growth
        % initial conditions 
    % columns denote: NH3, NH2OH, NO, NO3, NO2, N2O, N2, NH3_ex, NO3_ex,
    % NO2_ex, NO_ex, N2O_ex, N2_ex,Vc
    ics_YZN_001 = [0 0 0 0 0 0 0 5.294 0 0 0 0 0 0.000019];    
    %Volume of the medium
        Vm(1) = 0.125; %L 
    
    % Simulate model, specify in model kg = 0.003 or kg = 0, J specifies
    % growth function
    [t,y] = ode15s(@HNAD_sim_Vol,time1,ics_YZN_001,[],Vm(1),p_sets(i,:),1);
    
    %if statement to circumvent errors arising from warning: ode15s solver time-step below min value
    if t(end) < 4320
    score_ij(i,j) = NaN;
    
    else
    % Store output simulations for each parameter set
    % columns denote:time, NH3, NH2OH, NO, NO3, NO2, N2O, N2, NH3_ex, NO3_ex,
    % NO2_ex, NO_ex, N2O_ex, N2_ex,Vc
      output{i}{j} = [t,real(y)];    
        
    % Interpolate output to compare same time-points
    sim1 = interp1(time1,[t,real(y)],data1(:,1));
    %store simulations
    sim{i}{j}=[sim1];
    % score sim1 versus data1 now volume is included as variable in score
    % function
    score_ij(i,j) = sum((data1(:,2)-sim1(:,9)).^2)+ sum((data1(:,3)-sim1(:,10)).^2)+ sum((data1(:,4)-sim1(:,11)).^2)+ sum((data1(:,6)-sim1(:,14)).^2)+100*sum((data1(:,5)-sim1(:,15)).^2);
    end
    elseif j == 2
%% Condition 2 XL2, aerobic denitrification + growth on nitrate y(9)
     % initial conditions
     % columns denote: NH3, NH2OH, NO, NO3, NO2, N2O, N2, NH3_ex, NO3_ex,
    % NO2_ex, NO_ex, N2O_ex, N2_ex,Vc
    ics_XL2 = [0 0 0 0 0 0 0 0 1.5396 0.012 0 0 0 0.000019]; 
    %Volume of the medium
        Vm(2) = 0.250; %L
       
        %simulate model specify whether kg = 0.003 or kg = 0, J specifies
        %growth function
    [t,y] = ode15s(@HNAD_sim_Vol,time2,ics_XL2,[],Vm(2),p_sets(i,:),2);
    
    %if statement to circumvent errors arising from warning: ode15s solver time-step below min value
    if t(end) < 2160
    score_ij(i,j) = NaN;
    
    else
    % Store output simulations for each parameter set
    % columns denote:time, NH3, NH2OH, NO, NO3, NO2, N2O, N2, NH3_ex, NO3_ex,
    % NO2_ex, NO_ex, N2O_ex, N2_ex,Vc
      output{i}{j} = [t,real(y)];    
        
    % Interpolate output to compare same time-points
    sim2 = interp1(time2,[t,real(y)],data2(:,1));
    % store sim
    sim{i}{j}= [sim2];
    % score sim1 versus data1 + sim2 versus data2
    score_ij(i,j) = score_ij(i,j) + sum((data2(:,2)-sim2(:,9)).^2)+sum((data2(:,3)-sim2(:,10)).^2)+sum((data2(:,4)-sim2(:,11)).^2)+ 100*sum((data2(:,5)-sim2(:,15)).^2);
    end

        
end
  
%     Sum scores for individual data sets. Needs to be normalized by the number of measured data
%     points used to calculate the score.
%     data1 is 9*6 matrix, measured data points are 9*5
    D1= 9*5;
%     data2 is 9*5 matrix, measured data points are 6*5
    D2 = 9*4;
    score(i,1) = (score_ij(i,1)+score_ij(i,2))/(D1+D2);
end  
end
end


