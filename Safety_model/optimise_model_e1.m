%% Goal: sample parameter space
% Function to sample parameter space for the model extension 1.

% N = number of parameter samples to simulate in the chunk
%% Optimise model extension 1 function
function [output] = optimise_model_e1(p_sets,N,i)

% Time vector for simulations
time = 0:0.01:48;

% max formaldehyde level
F = (18*(10^-6))*ones(N,1);
% t1on - time of high F
ts1on = 0*ones(N,1);
% t1off - time of low F
ts1off = 101*ones(N,1);

% Start testing parameter sets (from 1 to N)
for i = 1:N
    
    % Set up initial conditions of model simulation.
    ics=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    
    % There are two conditions, with low or high formaldehyde. 
    for j = 1:2
        
        if j == 1 % 50% strength input
            [t,y] = ode15s(@model_e1,time,ics,[],F(i)*0.5,ts1on(i),ts1off(i),p_sets(i,:));
            
        elseif j == 2 % 100% strength input
            [t,y] = ode15s(@model_e1,time,ics,[],F(i),ts1on(i),ts1off(i),p_sets(i,:));
            
        end
        
        if t(end) < 12
           output{i}{j} = [0,0];
       else    
            %Store output for each parameter set
            output{i}{j} = [t,real(y)];
            end
       
    end  
    
end