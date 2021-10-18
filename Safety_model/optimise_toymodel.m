%% Goal: sample parameter space
% Function to sample parameter space and calculate scores for
% the toy model.

% data = data vector/matrix to compare simulations to
% N = number of parameter samples to simulate in the chunk
% q = number of the chunk
%% Optimise toy model function
function [score_ij,score,output] = optimise_toymodel(data,p_sets,N,q)

% Time vector for simulations
time = 0:0.01:1.1;

% I max level
I = (10^-3)*ones(N,1);%20*p_sets(:,5) OR 10*ones(N,1);
% t1on - time of I increase
ts1on = 0*ones(N,1);%60*p_sets(:,7) OR 30*ones(N,1);
% t1off - time of I decrease
ts1off = 1.1*ones(N,1);%80*p_sets(:,8) OR 40*ones(N,1);

% Set up empty vectors/matrices to store scores
score = zeros(N,1);
score_ij = zeros(N,length(data(1,1,:)));

% Start testing parameter sets (from 1 to N)
for i = 1:N
    
    % Set up initial conditions of model simulation
    ics=[0, 0, 0];
    
    % Simulate different datasets
    for j = 1:length(data(1,1,:))
        
        if j == 1 % 0% strength input
            [t,y] = ode15s(@hoksok_toymodel,time,ics,[],I(i)*0,ts1on(i),ts1off(i),p_sets(i,:));
            
        elseif j == 2 % 100% strength input
            [t,y] = ode15s(@hoksok_toymodel,time,ics,[],I(i),ts1on(i),ts1off(i),p_sets(i,:));
        
        end
        
        % Interpolate output so that you are comparing the same time-points.
        sim = interp1(t,[t,real(y)],data(:,1,j)/3600);
        
        % Store output simulations for each parameter set
         if t(end) < 1.1
           output{i}{j} = [0,0];
       else    
            %Store output for each parameter set
            output{i}{j} = [t,real(y)];
         end
        
        % Calculate score using "sum of squared differences"
        score_ij(i,j) = sum((data(end,2,j)-(sim(end,3)+sim(end,4))./(sim(end,2)+sim(end,3))).^2);
    end
    
    % Sum the scores for individual data sets so that each parameter set
    % has a single score. Normalised by the number of measured data
    % points used to calculate the score.
    score(i,1) = sum(score_ij(i,:))/length(data(:,2,1));
end

end