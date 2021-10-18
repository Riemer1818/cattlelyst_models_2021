%% Model file for simulation N-dynamics for Pseudomonas putida ZN1

% Inputs:
% t = time
% y(1-13) = concentration of Nitrogen species
% y(14) = population volume
% Vm = volume medium
% p_sets = parameter set

% Outputs: 
% [t,y] = time-y matrix
%% Description of the species
    % Intracellular N-species
        % y(1)= [NH3] 
        % y(2)= [NH2OH] 
        % y(3)= [NO]
        % y(4)= [NO3-]
        % y(5)= [NO2-]
        % y(6)= [N2O]
        % y(7)= [N2]

    % extracellular N-species
        % y(8)= [NH3]ex
        % y(9)= [NO3-]ex
        % y(10)= [NO2-]ex
        % y(11)= [NO]ex
        % y(12)= [N2O]ex
        % y(13) = [N2]ex

    % population volume
        % y(14) = y(14) 

%% Parameters 
    %Max reaction rates
    % units: Vmax = mmol.min-1
        % Vmax1 = Vmax for R1 
        % Vmax2 = Vmax for R2 
        % Vmax3 = Vmax for R3
        % Vmax4a = Vmax for R4 forward
        % Vmax4b = Vmax for R4 backward
        % Vmax5a = Vmax for R5 forward
        % Vmax5b = Vmax for R5 backward
        % Vmax6 = Vmax for R6
        % Vmax7 = Vmax for R7 
        
     
    % Michaelis constants
    % units: Km = mmol
    
        % Km1 = Km for substrate NH3 for enzyme AMO
        % Km2 = Km for substrate NH2OH for enzyme HAO 
        % Km4a = Km for substrate NO2- for enzyme Nir
        % Km4b = Km for substrate NO for enzyme Nir
        % Km5a = Km for substrate NO3- for enzyme Nap
        % Km5b = Km for substrate NO2- for enzyme NxrAB
        % Km6 = Km for substrate NO for enzyme NOR
        % Km7 = Km for subsrate N2O for enzyme N2OR

    % Exchange reaction rates mmol.min-1
    
        % Tmax1a = uptake rate [NH3]ex
        % Tmax1b = excretion rate [NH3]
        % Tmax2a = uptake rate [NO3-]ex
        % Tmax2b = excretion rate [NO3-]
        % Tmax3a = uptake rate [NO2-]ex
        % Tmax3b = uptake rate [NO2-]
        
    % Sensitivites of transporters unit: mmol
        % KT1a = affinity of transporter for [NH3]ex
        % KT1b = affinity of transporter for [NH3}
        % KT2a = affinity of transporter for [NO3-]ex
        % KT2b = affinity of transporter for [NO3-]
        % KT3a = affinity of tranporter for [NO2-]ex
        % KT3b = affinity of transporter for [NO2-]
        
    % Gas diffusion constants unit: L.min-1    
        % kd1 = Diffusion constant NO
        % kd2 = Diffusion constant N2O
        % kd3 = diffusion constant N2

    %% Volume parameters
% mumax1 = maximum growth rate on [NH3]ex
% mumax2 = maximum growth rate on [NO3-]ex
% KNH3 = sensitivity constant for [NH3]ex
% KNO3- = sensitivity constant for [NO3-]ex



function [dydt] = HNAD_sim_mix_ZN1(t,y,Vm,p)
%HNAD_sim is a function describing the HNAD mechanism of Pseudomonas stutzeri.
% Input parameters (36)

    %  Intracellular reactions
Vmax1 = p(1);
Vmax2 = p(2);
Vmax3 = p(3);
Vmax4a = p(4);
Vmax4b = p(5);
Vmax5a = p(6);
Vmax5b = p(7);
Vmax6 = p(8);
Vmax7 = p(9);

Km1 = p(10);
Km2 = p(11);
Km4a = p(12);
Km4b = p(13);
Km5a = p(14);
Km5b = p(15);
Km6 = p(16);
Km7 = p(17);

    % Transport reactions
 % Transport rate
 
Tmax1a = 0.7*p(18);% specifically for ZN1
Tmax1b = p(19);
KT1a = p(20);
KT1b = p(21);

Tmax2a = p(22);
Tmax2b = p(23);
KT2a = p(24);
KT2b = p(25);

Tmax3a = p(26);
Tmax3b = p(27);
KT3a = p(28);
KT3b = p(29);

 % Gas release factors
kd1 = p(30);
kd2 = p(31);
kd3 = p(32);
 
    %Volume parameters 
    % ammonia
mumax1 = 1.15*p(33); 
KNH3 = p(34); 
    % nitrate
mumax2 = 0.85*p(35); 
KNO3 = p(36); 
    % nitrite
mumax3 = 0.8*p(35);
KNO2 = 3*p(36);

% from steady-state approximation growth factor
kg = 0.003;

% Model equations
% Matrix for y
dydt = zeros(size(y));


% Switch N-source according to preference
                        % Ammonia growth
                        if y(8) > 0
                        mu = mumax1*(y(8)/(KNH3+y(8)));
                        % Nitrate growth
                        elseif y(8)<= 0.2 && y(9) > 0.1
                        mu =  mumax2*(y(9)/(KNO3+y(9)));
                        % Nitrite growth
                        elseif y(8)<= 0.2 && y(10)> 0.1
                        mu =  mumax3*(y(10)/(KNO2+y(10)));
                        % No growth
                        else
                        mu = 0;
                        end 
                        
                        
% Switch N-uptake, only used for ZN1
                        if y(8) > 0.6
                            Tmax2a = 0;                             
                            Tmax3a = 0;                      
                        elseif y(8) < 0.6
                            Tmax2a = Tmax2a;                           
                            Tmax3a = Tmax3a;
                        end

  
% simulate volume
dydt(14) = mu*y(14);

%% Intracellular species
%dNH3/dt 
dydt(1) = ((Tmax1a*y(8)*Vm/KT1a-Tmax1b*y(1)*y(14)/KT1b)/(1+y(8)*Vm/KT1a+y(1)*y(14)/KT1b)-(Vmax1*y(1)*y(14)/(Km1+y(1)*y(14))))/y(14)-(mu+kg)*y(1);
%dNH2OH/dt
dydt(2) = (Vmax1*y(1)*y(14)/(Km1+y(1)*y(14))-Vmax2*y(2)*y(14)/(Km2+y(2)*y(14))-2*Vmax3*(y(2)^2)*(y(14)^2)/((Km2^2)+(y(2)^2)*(y(14)^2)))/y(14)-mu*y(2);
%dNO/dt
dydt(3) = (Vmax2*y(2)*y(14)/(Km2+y(2)*y(14))+(Vmax4a*y(5)*y(14)/Km4a-Vmax4b*y(3)*y(14)/Km4b)/(1+y(5)*y(14)/Km4a+y(3)*y(14)/Km4b)-2*Vmax6*(y(3)^2)*(y(14)^2)/((Km6^2)+(y(3)^2)*(y(14)^2))-kd1*(y(3)-y(11)))/y(14)-mu*y(3);
%dNO3-/dt
dydt(4)= ((Tmax2a*y(9)*Vm/KT2a-Tmax2b*y(4)*y(14)/KT2b)/(1+y(9)*Vm/KT2a+y(4)*y(14)/KT2b)-(Vmax5a*y(4)*y(14)/Km5a-Vmax5b*y(5)*y(14)/Km5b)/(1+y(4)*y(14)/Km5a+y(5)*y(14)/Km5b))/y(14)-mu*y(4);
%dNO2-/dt
dydt(5)= ((Tmax3a*y(10)*Vm/KT3a-Tmax3b*y(5)*y(14)/KT3b)/(1+y(10)*Vm/KT3a+y(5)*y(14)/KT3b)-(Vmax4a*y(5)*y(14)/Km4a-Vmax4b*y(3)*y(14)/Km4b)/(1+y(5)*y(14)/Km4a+y(3)*y(14)/Km4b)+(Vmax5a*y(4)*y(14)/Km5a-Vmax5b*y(5)*y(14)/Km5b)/(1+y(4)*y(14)/Km5a+y(5)*y(14)/Km5b))/y(14)-mu*y(5);
%N2O/dt
dydt(6)= (Vmax3*(y(2)^2)*(y(14)^2)/((Km2^2)+(y(2)^2)*(y(14)^2))+Vmax6*(y(3)^2)*(y(14)^2)/((Km6^2)+(y(3)^2)*(y(14)^2))-Vmax7*y(6)*y(14)/(Km7+y(6)*y(14))-kd2*(y(6)-y(12)))/y(14)-mu*y(6);
%N2/dt
dydt(7)= (Vmax7*y(6)*y(14)/(Km7+y(6)*y(14))-kd3*(y(7)-y(13)))/y(14)-mu*y(7);

%% extracellular species
%NH3ex/dt
dydt(8)= -((Tmax1a*y(8)*Vm/KT1a-Tmax1b*y(1)*y(14)/KT1b)/(1+y(8)*Vm/KT1a+y(1)*y(14)/KT1b))/Vm;
%NO3-ex/dt
dydt(9)= -((Tmax2a*y(9)*Vm/KT2a-Tmax2b*y(4)*y(14)/KT2b)/(1+y(9)*Vm/KT2a+y(4)*y(14)/KT2b))/Vm;
%NO2-ex/dt
dydt(10)= -((Tmax3a*y(10)*Vm/KT3a-Tmax3b*y(5)*y(14)/KT3b)/(1+y(10)*Vm/KT3a+y(5)*y(14)/KT3b))/Vm;
%NOex/dt
dydt(11)= kd1*(y(3)-y(11))/Vm;
%N2Oex/dt
dydt(12)= kd2*(y(6)-y(12))/Vm;
%N2ex/dt
dydt(13)= kd3*(y(7)-y(13))/Vm;

end
