%% Script for Sensitivity analysis to determine influential parameters to limit N2O production
% y-axis absolute scores
% x-axis all parameters 

%% in sensitivity N2O change output to score_ij(:,1) or score_ij(:,2) if separate conditions need to be simulated. Leave score as output in this script.

%load parameter set (log-scale)
load('p_sets_best.mat','p_sets_SA')
%convert to linear scale
p_sets = 10.^(p_sets_SA);

% A for 32 parameters excluding monod parameters
A = ["Vmax1","Vmax2","Vmax3","Vmax4a","Vmax4b","Vmax5a","Vmax5b","Vmax6","Vmax7","Km1","Km2","Km4a","Km4b","Km5a","Km5b","Km6","Km7","Tmax1a","Tmax1b","KT1a","KT1b","Tmax2a","Tmax2b","KT2a","KT2b","Tmax3a","Tmax3b","KT3a","KT3b","kd1","kd2","kd3"];

%% alter p_sets OAT by factor 10
figure()
k = [1:1:32];% k represents parameters
 [score] = sensitivity_N2O(p_sets); % compute score with best parameter sets
   score_i = score;%set that as initial score
   score_plot = ones(1,32)*score;% set baseline
      plot(k,score_plot,'db','MarkerFaceColor','b','MarkerSize',14); %plot blue diamonds on baseline
        hold on
  
for k = 1:32
    %Pertube one parameter
    p_pertube = 10*p_sets(1,k);
    %One perturbed parameter per parameter set
    p_sets_pertubed = [p_sets(1,1:(k-1)),p_pertube,p_sets(1,(k+1):32),p_sets(1,33:36)];
     
 %% compute score with altered parameter set
    [score] = sensitivity_N2O(p_sets_pertubed);
            if isnan(score)% plot asterix if simulation fails
            plot(k,0,'pk','MarkerFaceColor','k','MarkerSize',14); 
            hold on
        else
    %Compute absolute score 
    score_r1 = score; 
    % plot green diamonds for factor 10 increase of parameter value
    plot(k,score_r1,'dg','MarkerFaceColor','g','MarkerSize',14); 
    hold on
end
end
    hold on
for k = 1:32
    %Pertube one parameter, divide by 10
    p_pertube = 0.1*p_sets(1,k);
    %one perturbed parameter per parameter set
    p_sets_pertubed = [p_sets(1,1:(k-1)),p_pertube,p_sets(1,(k+1):32),p_sets(1,33:36)];
    
    
    %% calculate score with perturbed parameter set
    [score] = sensitivity_N2O(p_sets_pertubed);
           if isnan(score)
            plot(k,0,'pk','MarkerFaceColor','k','MarkerSize',14); 
            hold on
        else
    %Compute absolute score
    score_r2 = score; 
    plot(k,score_r2,'dr','MarkerFaceColor','r','MarkerSize',14); 
    hold on
end
end
%% figure specs
% axis specs
   title('Sensitivity N_2O production')
   xlabel('Parameters')
   xticks(1:1:32)
   xticklabels(A)
   xtickangle(90)
   ylabel('Total N_2O_e_x (mM)')
   set(gca,'YScale','linear')
% set lay out
set(gca,'Fontsize',24)
set(gca,'FontName','Book Antiqua')
   