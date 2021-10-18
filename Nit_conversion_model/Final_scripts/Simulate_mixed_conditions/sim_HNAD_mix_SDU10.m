%% function sim_HNAD_mix_SDU10 is a function that simulates the HNAD model for different [ammonia]
% output: simulation output
% inputs: p_sets = best parameter set (36 colomns, 1 row)
% I = number of simulations, if more simulation rounds should be
% run, this can be done by specifying multiple initial conditions
% or providing a parameter matrix. 


function [output] = sim_HNAD_mix_SDU10(p_sets,I)

%time vectors for simulations based on the data for SDU10
% two conditions are simulated so specify two time-vectors
time1 = 0:60:1920; 
time2 = 0:60:1920; 


for i= 1:I
%% Simulate different conditions
for j =1:2
    %% simulate condition 1, growth on ammonia and nitrate
    if j == 1 
    
    % initial conditions 
    
    % columns denote: NH3, NH2OH, NO, NO3-, NO2-, N2O, N2, NH3_ex, NO3_ex,
    % NO2_ex, NO_ex, N2O_ex, N2_ex, Vc
    
    % specify initial ammonia, nitrate and population volume values based
    % on data
    ics_SDU10_1 = [0 0 0 0 0 0 0 5.871 1.6097 0 0 0 0 0.000019];    
    
    % specify volume of the medium
        Vm(1) = 0.100; %L 
   
    % simulate model to obtain output
    [t,y] = ode15s(@HNAD_sim_mix_SDU10,time1,ics_SDU10_1,[],Vm(1),p_sets(i,:));
    
    % Store output simulations for this condition
      output{i}{j} = [t,real(y)];       
  
    elseif j == 2
    %% simulate condition 2, growth on ammonia and nitrite
            % initial conditions 
    
    % columns denote: NH3, NH2OH, NO, NO3-, NO2-, N2O, N2, NH3_ex, NO3_ex,
    % NO2_ex, NO_ex, N2O_ex, N2_ex, Vc
    
    % specify initial ammonia, nitrate and population volume values based
    % on data
    ics_SDU10_2 = [0 0 0 0 0 0 0 5.8706 0 2.1804 0 0 0 0.000019]; 
    
    %Specify volume of the medium
        Vm(2) = 0.100; %L
         
    %simulate model to obtain output
    [t,y] = ode15s(@HNAD_sim_mix_SDU10,time2,ics_SDU10_2,[],Vm(2),p_sets(i,:));
    
    % Store output simulations for this condition
      output{i}{j} = [t,real(y)];    
           
    end
     
end
      
end
