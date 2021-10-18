%% function sim_HNAD_mix_ZN1 is a function that simulates the HNAD model for different [ammonia]
% output: simulation output
% inputs: p_sets = best parameter set (36 colomns, 1 row)
% I = number of simulations, if more simulation rounds should be
% run, this can be done by specifying multiple initial conditions
% or providing a parameter matrix. 


function [output] = sim_HNAD_mix_ZN1(p_sets,I)

%time vector for simulations based on the data for ZN1
time1 = 0:60:2880; 
time2 = 0:60:2880; 


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
        ics_ZN1_1 = [0 0 0 0 0 0 0 2.7176 0.8242 0 0 0 0 0.000019];    
    
        % specify volume of the medium
        Vm(1) = 0.100; %L 
        
    % simulate system to obtain output
    [t,y] = ode15s(@HNAD_sim_mix_ZN1,time1,ics_ZN1_1,[],Vm(1),p_sets(i,:));
    
    % Store output simulations for this condition
      output{i}{j} = [t,real(y)];       
    
    elseif j == 2
 %% simulate condition 2, growth on ammonia and nitrite
        % initial conditions
        
        % columns denote: NH3, NH2OH, NO, NO3-, NO2-, N2O, N2, NH3_ex, NO3_ex,
        % NO2_ex, NO_ex, N2O_ex, N2_ex, Vc
        
         % specify initial ammonia, nitrate and population volume values based
         % on data
     
    ics_ZN1_2 = [0 0 0 0 0 0 0 2.88 0 1.1000 0 0 0 0.000019]; 
    
    
        %Specify volume of the medium
        Vm(2) = 0.100; %L

           
        %simulate model to obtain output
    [t,y] = ode15s(@HNAD_sim_mix_ZN1,time2,ics_ZN1_2,[],Vm(2),p_sets(i,:));
    
     % Store output simulations for this condition
      output{i}{j} = [t,real(y)];    
           
    end
     
end
      
end
