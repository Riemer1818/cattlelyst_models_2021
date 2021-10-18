%% Function sensitivity_N2O is a function that computes the N2O score based on provided parameters
% inputs: parameter set
% outputs: score, but can be altered to score_ij(:,1) (condition 1), or score_ij(:,2)(condition 2) to
% compare conditions separately


function [score] = sensitivity_N2O(p_sets)

%time vector for simulations
time1 = 0:60:4320; %72 hours for data1
time2 = 0:60:2160; %31 hours for data2

% Scores
% Set up empty vectors/matrices to store scores only 1 per call
I = 1;
score = zeros(I,1);
score_ij = zeros(I,2);

for i= I
%% Simulate different conditions
for j =1:2
    if j == 1 
%% Condition 1, YZN-001 ammonia growth HNAD
        % initial conditions 
        % columns output: NH3, NH2OH, NO, NO3, NO2, N2O, N2, NH3_ex, NO3_ex,
        % NO2_ex, NO_ex, N2O_ex, N2_ex,Vc

    ics_YZN_001 = [0 0 0 0 0 0 0 5.294 0 0 0 0 0 0.000019];    
    %Medium volume
        Vm(1) = 0.125; %L 
    % simulate system
    [t,y] = ode15s(@HNAD_sim_Vol,time1,ics_YZN_001,[],Vm(1),p_sets(i,:),1);
    
    %if statement to circumvent errors arising from warning: ode15s solver time-step below min value
    if t(end) < 4320
    score_ij(i,j) = NaN;
    
    else
         
   % Store output simulations 
   % columns output: time, NH3, NH2OH, NO, NO3, NO2, N2O, N2, NH3_ex, NO3_ex,
   % NO2_ex, NO_ex, N2O_ex, N2_ex,Vc
    output{i}{j} = [t,real(y)]; 
    % compute score for YZN-001 condition (sum-N2O_ex)
    score_ij(i,j) = sum(output{i}{j}(:,13));
   
    end
%% Data2 = Data XL2, aerobic denitrification + growth on nitrate y(9)
    elseif j == 2
        % initial conditions
    % columns output: NH3, NH2OH, NO, NO3, NO2, N2O, N2, NH3_ex, NO3_ex,
    % NO2_ex, NO_ex, N2O_ex, N2_ex,Vc
    ics_XL2 = [0 0 0 0 0 0 0 0 1.5396 0.012 0 0 0 0.000019]; 
    %Medium volume
        Vm(2) = 0.250; %L
        %simulate system
    [t,y] = ode15s(@HNAD_sim_Vol,time2,ics_XL2,[],Vm(2),p_sets(i,:),2);
    
    %if statement to circumvent errors arising from warning: ode15s solver time-step below min value
    if t(end) < 2160
    score_ij(i,j) = NaN;
 
    else
        
   % Store output simulations for each parameter set
   % columns output: time, NH3, NH2OH, NO, NO3, NO2, N2O, N2, NH3_ex, NO3_ex,
   % NO2_ex, NO_ex, N2O_ex, N2_ex,Vc
      output{i}{j} = [t,real(y)];
  
    % compute score score XL-2
    score_ij(i,j) = score_ij(i,j) +  sum(output{i}{j}(:,13));
    end

        
end
  
%     Sum scores for individual data sets, normalized by dividing by 2
    score(i,1) = (score_ij(i,1)+score_ij(i,2))/2;
end
  
end
end


