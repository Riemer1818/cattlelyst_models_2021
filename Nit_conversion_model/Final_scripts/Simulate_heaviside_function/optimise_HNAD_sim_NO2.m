%% optimise_HNAD_sim_NO2 is a function that calls the model to simulate the system
% inputs: 
% data1, for p. stutzeri YZN-001
% data2, for p. stutzeri Xl-2
% p_sets, parameter sets
% I denotes the number of simulations run before saving, will be set to 1

%outputs: 
% sim: interpolated output corresponding to time data1 and data2
% score_ij: score matrix, column 1 corresponds to condition 1, column two
% corresponds to condition 2
% score: normalized score
% Output: t,y matrix


function [sim, score_ij, score, output] = optimise_HNAD_sim_NO2(data1,data2,p_sets,I)

%time vector for simulations
time1 = 0:60:4320; %72 hours for data1
time2 = 0:60:2160; %31 hours for data2

% Scores
% Set up empty vectors/matrices to store scoresdata1
score = zeros(I,1);
score_ij = zeros(I,2);

% Start testing parameter sets (from 1 to i) will be 1

for i= 1:I
%% Simulate different conditions
for j =1:2
    if j == 1 
%% Condition 1, YZN_001, HNAD + growth in ammonia y(8)
        % initial conditions 
        % columns output: NH3, NH2OH, NO, NO3, NO2, N2O, N2, NH3_ex, NO3_ex,
        % NO2_ex, NO_ex, N2O_ex, N2_ex,Vc

    ics_YZN_001 = [0 0 0 0 0 0 0 5.294 0 0 0 0 0 0.000019]; 
    
    %Volume of the medium
        Vm(1) = 0.125; %L 
        
    % simulate system
    [t,y] = ode15s(@HNAD_sim_Vol_heaviside,time1,ics_YZN_001,[],Vm(1),p_sets(i,:),1);
    
    %if statement to circumvent errors arising from warning: ode15s solver time-step below min value
    if t(end) < 4320
    score_ij(i,j) = 10000000;
    
    else
    % Store output simulations for each parameter set
    % columns output: time, NH3, NH2OH, NO, NO3, NO2, N2O, N2, NH3_ex, NO3_ex,
    % NO2_ex, NO_ex, N2O_ex, N2_ex,Vc
      output{i}{j} = [t,real(y)];    
        
    % Interpolate output to compare same time-points
    sim1 = interp1(time1,[t,real(y)],data1(:,1));
    %store simulations
    sim{i}{j}=[sim1];
    % score sim1 versus data1
    score_ij(i,j) = sum((data1(:,2)-sim1(:,9)).^2)+ sum((data1(:,3)-sim1(:,10)).^2)+ sum((data1(:,4)-sim1(:,11)).^2)+ sum((data1(:,6)-sim1(:,14)).^2)+100*sum((data1(:,5)-sim1(:,15)).^2);
    
       end
    elseif j == 2
%% condition 2, XL-2 AD and nitrate growth
        % initial conditions
        % columns output:  NH3, NH2OH, NO, NO3, NO2, N2O, N2, NH3_ex, NO3_ex,
        % NO2_ex, NO_ex, N2O_ex, N2_ex,Vc

    ics_XL2 = [0 0 0 0 0 0 0 0 1.5396 0.012 0 0 0 0.000019]; 
    
        %Volume of the medium
        Vm(2) = 0.250; %L
 
        %simulate system
    [t,y] = ode15s(@HNAD_sim_Vol_heaviside,time2,ics_XL2,[],Vm(2),p_sets(i,:),2);
    
    %if statement to circumvent errors arising from warning: ode15s solver time-step below min value
    if t(end) < 2160
    score_ij(i,j) = 10000000;
    
    else
    % Store output simulations for each parameter set
    % columns output: time, NH3, NH2OH, NO, NO3, NO2, N2O, N2, NH3_ex, NO3_ex,
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
