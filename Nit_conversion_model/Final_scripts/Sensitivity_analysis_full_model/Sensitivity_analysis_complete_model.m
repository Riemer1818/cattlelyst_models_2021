%% Script to conduct sensitivity analysis (takes approx 3-5 min to run)
%y-axis relative sensitivity coefficients, numerically approached.
%x-axis , discrete, all parameters

%load best parameter set (log scale)
load('p_sets_best.mat','p_sets_SA')
% convert to linear scale to simulate model
p_sets = 10.^(p_sets_SA);

% Specify parameter names
A = ["Vmax1","Vmax2","Vmax3","Vmax4a","Vmax4b","Vmax5a","Vmax5b","Vmax6","Vmax7","Km1","Km2","Km4a","Km4b","Km5a","Km5b","Km6","Km7","Tmax1a","Tmax1b","KT1a","KT1b","Tmax2a","Tmax2b","KT2a","KT2b","Tmax3a","Tmax3b","KT3a","KT3b","kd1","kd2","kd3","\mumax1","KNH3","\mumax2","KNO3"];
%% figure to visualize the relative sensitivity coefficients

figure()
set(gcf,'DefaultAxesFontSize',14); 

%% Simulate system with best parameter set
k = [1:1:36];
% obtain baseline score for p_sets_best
 [score] = sensitivity_HNAD(p_sets);
   score_i = score;
%    score_plot = zeros(1,36);
%    plot(k,score_plot,'db','MarkerFaceColor','b') 
%         hold on

%% 1 percent increase and decrease of parameter values
        subplot(2,2,1)
                
for k = 1:36 % specifies number of parameters
    %Pertube one parameter
    for i = 1:2 % specifies increase (1), decrease (2)
        if i == 1
        p_pertube = 1.01*p_sets(1,k); % 1 % parameter change
        p_sets_pertubed = [p_sets(1,1:(k-1)),p_pertube,p_sets(1,(k+1):36)]; %parameter set with 1 perturbed parameter
        [score] = sensitivity_HNAD(p_sets_pertubed);% compute score with perturbed parameter set
        if isnan(score) %if simulation error occurs plot asterix
              plot(k,0,'*k');
            hold on
        else
    % Compute relative sensitivity (%) coefficient 
    score_r1 =(score-score_i)/abs(p_pertube-p_sets(1,k))*(p_sets(1,k)/score_i);
    % plot green diamonds for increased parameter values
    plot(k,score_r1,'dg','MarkerFaceColor','g')
        hold on
        end
        elseif i == 2 % start simulation decreased parameter set
        p_pertube = 0.99*p_sets(1,k);% 1% parameter change
        p_sets_pertubed = [p_sets(1,1:(k-1)),p_pertube,p_sets(1,(k+1):36)];% parameter set with 1 perturbed parameter
        [score] = sensitivity_HNAD(p_sets_pertubed);% compute score with pertrubed parameter set
        if isnan(score)% plot asterix if simulation error occurs
             plot(k,0,'*k');
            hold on
        else
    %Compute relative sensitivity (%) coefficient 
    score_r1 =(score-score_i)/abs(p_pertube-p_sets(1,k))*(p_sets(1,k)/score_i); 
    % plot red diamonds for decreased parameter values
    plot(k,score_r1,'dr','MarkerFaceColor','r')
    hold on
        end
        end
        hold on
        end
end
%% set figure specs
% axis titles
title('Sensitivity coefficients (1% \Deltap_0)')
xlabel('Parameters')
xticks(1:1:36)
xticklabels(A)
xtickangle(90)
ylabel('relative sensitivity (%)')
set(gca,'YScale','linear')
set(gca,'FontName','Book Antiqua');

     
%% 3 percent increase and decrease of parameter values    
    subplot(2,2,2)
           
for k = 1:36 % specifies number of parameters
    %Pertube one parameter
    for i = 1:2 %denotes increase (1) or decrease
        if i == 1
        p_pertube = 1.03*p_sets(1,k); % 3% parameter increase
        p_sets_pertubed = [p_sets(1,1:(k-1)),p_pertube,p_sets(1,(k+1):36)]; %parameter set with one perturbed parameter
        [score] = sensitivity_HNAD(p_sets_pertubed); %compute score with pertrubed parameter set
                if isnan(score) % if simulation error occurs plot asterix
            plot(k,0,'*k');
            hold on
                else
    % Compute relative sensitivity (%) coefficient 
    score_r3 =(score-score_i)/abs(p_pertube-p_sets(1,k))*(p_sets(1,k)/score_i); 
    % plot green diamonds for increased parameter values
    plot(k,score_r3,'dg','MarkerFaceColor','g')
        hold on
                end
        elseif i == 2 % start simulating decreased parameter values
        p_pertube = 0.97*p_sets(1,k);% 3% parameter change
        p_sets_pertubed = [p_sets(1,1:(k-1)),p_pertube,p_sets(1,(k+1):36)];% parameter set with one perturbed parameter
        [score] = sensitivity_HNAD(p_sets_pertubed);% compute score with perturbed parameter set
        if isnan(score)% plot asterix if simulation error occurs
             plot(k,0,'*k');
            hold on
        else
    %Compute relative sensitivity (%) coefficient 
    score_r3 = (score-score_i)/abs(p_pertube-p_sets(1,k))*(p_sets(1,k)/score_i);
    % plot red diamonds for decreased parameter sets
    plot(k,score_r3,'dr','MarkerFaceColor','r')
    hold on
        end
        end
        hold on
        end
end
%% set figure specs
% axis titles
title('Sensitivity coefficients (3% \Deltap_0)')
xlabel('Parameters')
xticks(1:1:36)
xticklabels(A)
xtickangle(90)
ylabel('relative sensitivity (%)')
set(gca,'YScale','linear')
set(gca,'FontName','Book Antiqua');

%% 5 percent increase and decrease of parameter values
    subplot(2,2,3)
            
for k = 1:36 % specifies number of parameters
    %Pertube one parameter
    for i = 1:2 % specifies increase (1) or decrease (2)
        if i == 1
        p_pertube = 1.05*p_sets(1,k); % 5% parameter change
        p_sets_pertubed = [p_sets(1,1:(k-1)),p_pertube,p_sets(1,(k+1):36)];% parameter set with 1 perturbed parmeter
        [score] = sensitivity_HNAD(p_sets_pertubed); %compute score with perturbed parameter set
        if isnan(score) % if simulation error occurs plot asterix
             plot(k,0,'*k');
            hold on
        else
    %compute relative sensitivity (%) coefficient 
    score_r5 =(score-score_i)/abs(p_pertube-p_sets(1,k))*(p_sets(1,k)/score_i); 
    %plot green diamonds for increased parameter values
    plot(k,score_r5,'dg','MarkerFaceColor','g')
        hold on
        end
        elseif i == 2 % start simulating negative parameter values
        p_pertube = 0.95*p_sets(1,k); % 5% parameter change
        p_sets_pertubed = [p_sets(1,1:(k-1)),p_pertube,p_sets(1,(k+1):36)]; % parameter set with 1 perturbed parameter
        [score] = sensitivity_HNAD(p_sets_pertubed);% compute score with perturbed parameter value
        if isnan(score)% if simulation error occurs plot asterix
              plot(k,0,'*k');
            hold on
        else
    %compute relative sensitivity (%) coefficient 
    score_r5 = (score-score_i)/abs(p_pertube-p_sets(1,k))*(p_sets(1,k)/score_i); 
    % plot red diamond for decreased parameter values
    plot(k,score_r5,'dr','MarkerFaceColor','r')
    hold on
        end
        end
        hold on
        end
end
%% set figure specs
% axis title
title('Sensitivity coefficients (5% \Deltap_0)')
xlabel('Parameters')
xticks(1:1:36)
xticklabels(A)
xtickangle(90)
ylabel('relative sensitivity (%)')
set(gca,'YScale','linear')
set(gca,'FontName','Book Antiqua');

 %% 10 percent increase and decrease of parameter value
    subplot(2,2,4)
           
for k = 1:36 % specifies number of parameters
    %Pertube one parameter
    for i = 1:2 % specifies increase (1), or decrease (2)
        if i == 1
        p_pertube = 1.1*p_sets(1,k);% 10% parameter change
        p_sets_pertubed = [p_sets(1,1:(k-1)),p_pertube,p_sets(1,(k+1):36)]; %parameter set with 1 perturbed parameter
        [score] = sensitivity_HNAD(p_sets_pertubed);% compute score with perturbed parameter set
        if isnan(score) % if simulation error occurs plot asterix
              plot(k,0,'*k');
            hold on
        else
    %compute relative sensitivity (%) coefficient 
    score_r10 =(score-score_i)/abs(p_pertube-p_sets(1,k))*(p_sets(1,k)/score_i);
    % plot green diamond for positive parameter changes
    plot(k,score_r10,'dg','MarkerFaceColor','g')
            hold on
        end
            
        elseif i == 2 % start simulating model with decreased parameter sets
        p_pertube = 0.9*p_sets(1,k);% 10% parameter change
        p_sets_pertubed = [p_sets(1,1:(k-1)),p_pertube,p_sets(1,(k+1):36)];% parameter set with 1 perturbed parameter
        [score] = sensitivity_HNAD(p_sets_pertubed);% compute score with perturbed parameter set
        if isnan(score)% plot asterix if simulation error occurs
              plot(k,0,'*k');
            hold on
        else
    %compute relative sensitivity (%) coefficient 
    score_r10 = (score-score_i)/abs(p_pertube-p_sets(1,k))*(p_sets(1,k)/score_i);
    % plot red diamond for parameter decrease
    plot(k,score_r10,'dr','MarkerFaceColor','r')
    hold on
        end
        end
        hold on
        end
end

%% set figure specs
% set axis titles
title('Sensitivity coefficients (10% \Deltap_0)')
xlabel('Parameters')
xticks(1:1:36)
xticklabels(A)
xtickangle(90)
ylabel('relative sensitivity (%)')
set(gca,'YScale','linear')
set(gca,'FontName','Book Antiqua');

 %% 5 percent increase and decrease parameter values separate figure
      figure()
set(gcf,'DefaultLineLineWidth',2);
set(gcf,'DefaultAxesFontSize',24);      
   
          
for k = 1:36 % specifies number of parameters
    %Pertube one parameter
    for i = 1:2 % denotes parameter increase (1) or decrease (2)
        if i == 1
        p_pertube = 1.05*p_sets(1,k);% 5% parameter change
        p_sets_pertubed = [p_sets(1,1:(k-1)),p_pertube,p_sets(1,(k+1):36)];% parameter set with 1 perturbed parameter
        [score] = sensitivity_HNAD(p_sets_pertubed);% compute score with perturbed parameter set
        if isnan(score)% plot asterix if simulation error occurs
              plot(k,0,'*k','MarkerSize',14);
            hold on
        else
    %compute relative sensitivity (%) coefficient 
    score_r5 =(score-score_i)/abs(p_pertube-p_sets(1,k))*(p_sets(1,k)/score_i);
    %plot green diamond for parameter increase
    plot(k,score_r5,'dg','MarkerFaceColor','g','MarkerSize',14)
            hold on
        end
            
        elseif i == 2 % start simulating model with parameter decrease
        p_pertube = 0.95*p_sets(1,k);% 5% parameter change
        p_sets_pertubed = [p_sets(1,1:(k-1)),p_pertube,p_sets(1,(k+1):36)];%parameter set with 1 perturbed parameter
        [score] = sensitivity_HNAD(p_sets_pertubed);% compute score with perturbed parameter set
        if isnan(score)% if simulation error occurs plot asterix
           plot(k,0,'*k','MarkerSize',14);
            hold on
        else
    %compute relative sensitivity (%) coefficient 
    score_r5 = (score-score_i)/abs(p_pertube-p_sets(1,k))*(p_sets(1,k)/score_i); 
    % plot red diamond for parameter decrease
    plot(k,score_r5,'dr','MarkerFaceColor','r','MarkerSize',14)
    hold on
        end
        end
        hold on
        end
end

%% set figure specs
% axis titles
title('Sensitivity coefficients (5% \Deltap_0)')
xlabel('Parameters')
xticks(1:1:36)
xticklabels(A)
xtickangle(90)
ylabel('Sensitivity coefficients')
set(gca,'YScale','linear')
set(gca,'FontName','Book Antiqua');

