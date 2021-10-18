%% Goal: sample parameter space
% Function to sample parameter space for the model extension 1+7.

% N = number of parameter samples to simulate in the chunk
%% Optimise model extension 1+7 function
function [output] = optimise_model_e1_7(p_sets,N,i)

% Time vector for simulations (is concentrations are halved every x hours,
% end time of simulation needs to be x
time = 0:0.01:48; 

% max formaldehyde level
F = (18*(10^-6))*ones(N,1);
% t1on - time of high F
ts1on = 0*ones(N,1);
% t1off - time of low F
ts1off = 48*ones(N,1);

% Start testing parameter sets (from 1 to N)
for i = 1:N
    
    % Set up initial conditions of model simulation. Use this in case the
    % cell concentration is not halved and for the first run is the
    % concentrations are halved every x hours
    ics=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    
    % Set up initial conditions of model simulation when the concentrations
    % of compounds in the cell is halved every x hours. 
%     load('e1_7_conc.mat','c_highF','c_lowF'); 
%     ics_h = [(c_highF(1,1:13)).*0.5];
%     ics_l = [(c_lowF(1,1:13)).*0.5];
    
    
    % There are two conditions, with low or high formaldehyde. 
    for j = 1:2
        
        if j == 1 % 50% strength input
            [t,y] = ode15s(@model_e1_7,time,ics,[],F(i)*0.5,ts1on(i),ts1off(i),p_sets(i,:));
            %[t,y]=ode15s(@model_e1_7,time,ics_l,[],F(i)*0.5,ts1on(i),ts1off(i),p_sets(i,:));
            %first line is used if cell concentrations are not halved or
            %for the first run when the concentrations of the compounds in
            %the cell are halved. The second line is for the other cycles
            %when the cell concentration is halved every x hours.
            
        elseif j == 2 % 100% strength input
            [t,y] = ode15s(@model_e1_7,time,ics,[],F(i),ts1on(i),ts1off(i),p_sets(i,:));
            %[t,y] = ode15s(@model_e1_7,time,ics_h,[],F(i),ts1on(i),ts1off(i),p_sets(i,:));
            %first line is used if cell concentrations are not halved or
            %for the first run when the concentrations of the compounds in
            %the cell are halved. The second line is for the other cycles
            %when the cell concentration is halved every x hours.
        end
        
%         if t(end) < 12
%            output{i}{j} = [0,0];
%        else    
            %Store output for each parameter set
            output{i}{j} = [t,real(y)];
            end
       
    end  
    
end