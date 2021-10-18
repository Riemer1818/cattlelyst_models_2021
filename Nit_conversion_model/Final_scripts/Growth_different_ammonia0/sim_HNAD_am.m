%% function sim_HNAD_am is a function that simulates the HNAD model for different [ammonia]
% output: simulation output
% inputs: p_sets = best parameter set (36 colomns, 1 row)
%           I = number of simulations

function [output] = sim_HNAD_am(p_sets,I)


%time vector for simulations can be changed accordingly
time1 = 0:60:2880; 
time2 = 0:60:2880; 


for i= 1:I
%% Simulate different conditions
 % Ammonia concentrations to test, I corresponds to the length of this
 % vector. 
    Am_c = [0.25 0.5 1 2 4 5.01 10];
    
%% Specify function inputs; Initial conditions medium volume

% columns denote: NH3, NH2OH, NO, NO3-, NO2-, N2O, N2, NH3_ex, NO3_ex,
% NO2_ex, NO_ex, N2O_ex, N2_ex, Vc

 ics = [0 0 0 0 0 0 0 Am_c(i) 0 0 0 0 0 0.000019];    
       
%Volume of the medium
        Vm = 0.100; %L  

% simulate system to obtain output
    [t,y] = ode15s(@HNAD_sim_am,time1,ics,[],Vm,p_sets(1,:));

% Store output simulations for each condition

      output{i}{1} = [t,real(y)];                  
    end
     
end
      

