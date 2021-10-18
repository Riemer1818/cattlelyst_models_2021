%tutorial: https://opencobra.github.io/cobratoolbox/latest/tutorials/tutorialSteadyCom.html
%% Initiate COBRA
%initCobraToolbox % if you want to update
initCobraToolbox(false) % if you dont need to update
%changeCobraSolver('glpk'); % should also work for most SteadyCom functions
changeCobraSolver('ibm_cplex', 'LP'); % recommanded solver note: The solver compatibility is not tested with MATLAB R2018b

%% Load the Paracoccus and Methylococcus models
% according to SteadyCom Tutorial

% Paracoccus with only modelseed namespaces
%modelP = readCbModel('Paracoccus_MinimalMedia_exchanges.mat');

% Paracoccus model with adjusted extracellular metabolite namespaces
modelP = readCbModel('Paracoccus_newNamespaces_cFBA.mat');

% Methylococcus model with old namespaces
%modelM = readCbModel('iMcBath_NH4.mat');

% Methylococcus model with new namespaces
modelM = readCbModel('iMcBath_newNamespaces.mat');

%% Polishing the old models (with old namespaces)

% Convert the compartment format from e.g., '_c' or '_c0[cytosol_0] to
% '[c]' (or c0) 

% for Methylococcus
modelM.mets = regexprep(modelM.mets, '_([^_]?)$', '\[$1\]');
modelM.mets = regexprep(modelM.mets, '_im$', '\[im\]');
modelM.rxns = regexprep(modelM.rxns, '_e$', '\[e\]');
modelM.rxns = regexprep(modelM.rxns, '_c$', '\[c\]');

% for Paracoccus
modelP.mets = regexprep(modelP.mets, '_([^e]+)$', '\[c0\]');
modelP.mets = regexprep(modelP.mets, '_([^_]+)$', '\[e\]');
modelP.rxns = regexprep(modelP.rxns, '_e0', '\[e\]');
modelP.rxns = regexprep(modelP.rxns, '_c0', '\[c0\]');

% Two exchange reactions in the Paracoccus model did not represent
% extracellular compounds, so to prevent them from being removed while
% reating a compartmentalized model, they were manually renamed into a
% demand reaction
% Changing bounds of renamed "demand" reactions in Paracoccus model
modelP = changeRxnBounds(modelP, 'DM_cpd15302[c0]', 0, 'l');
modelP = changeRxnBounds(modelP, 'DM_cpd02701[c0]', 0, 'l');

% The namespaces of extracellular metabolites present in both GSMMs were 
%manually altered into BiGG namespace in the P. denitrificans model.

% Saving the models with new namespaces
writeCbModel(modelP)
writeCbModel(modelM)

%% SteadyCom compartmentalized model creation
nameTagsModel = {'M'; 'P'};
modelJ = createMultipleSpeciesModel({modelM; modelP}, nameTagsModel);

%% Add the second objective to the compartmentalized model
modelJ = changeObjective(modelJ, {'MBIOMASS_Mcapsulatus'; 'Pbio1'} );
objectiveJ = checkObjective(modelJ);

%% Set medium for extracellular space 'u'
% Close all extracellular space 'u' exchange reactions
for i = 2406:2505
    modelJ =  changeRxnBounds(modelJ, modelJ.rxns{i}, 0, 'l');
end

% Reopen exchange reactions of compounds in a minimal medium
minimalMedia = {'EX_ca2[u]'; 'EX_cd2[u]'; 'EX_ch4[u]'; 'EX_cl[u]'; 'EX_co2[u]';
    'EX_co[u]'; 'EX_cobalt2[u]'; 'EX_cpd10516[u]'; 'EX_cpd11606[u]'; 'EX_cu2[u]';
    'EX_fe2[u]'; 'EX_h2o[u]'; 'EX_h[u]'; 'EX_k[u]'; 'EX_mg2[u]'; 'EX_mn2[u]';
    'EX_n2[u]'; 'EX_na1[u]'; 'EX_nh4[u]'; 'EX_ni2[u]'; 'EX_o2[u]'; 'EX_pi[u]';
    'EX_so4[u]'; 'EX_zn2[u]'};

for i = minimalMedia
    modelJ = changeRxnBounds(modelJ, i, -1000, 'l');
end 

%% Set organism-specific uptake rates to be the same as in the original model:
% Methylococcus
% Close all extracellular space '[u]tr' exchange reactions (lb)
for i = 864:917
    modelJ =  changeRxnBounds(modelJ, modelJ.rxns{i}, 0, 'l');
end

% Add specific methane uptake rate
modelJ =  changeRxnBounds(modelJ, 'MIEX_ch4[u]tr', -18.46, 'l');

% Reopen exchange reactions of organism (lb)
exM = {'MIEX_ca2[u]tr'; 'MIEX_cd2[u]tr'; 'MIEX_cl[u]tr'; 'MIEX_co2[u]tr'; 'MIEX_co[u]tr';
    'MIEX_cobalt2[u]tr'; 'MIEX_cu2[u]tr'; 'MIEX_fe2[u]tr'; 'MIEX_h2o[u]tr'; 'MIEX_h[u]tr';
    'MIEX_k[u]tr'; 'MIEX_mg2[u]tr'; 'MIEX_mn2[u]tr'; 'MIEX_n2[u]tr'; 'MIEX_na1[u]tr'; 'MIEX_nh4[u]tr';
    'MIEX_ni2[u]tr'; 'MIEX_o2[u]tr'; 'MIEX_pi[u]tr'; 'MIEX_so4[u]tr'; 'MIEX_zn2[u]tr'};

for i = exM
    modelJ = changeRxnBounds(modelJ, i, -1000, 'l');
end 

% Close exchange reactions of organism (ub)
modelJ = changeRxnBounds(modelJ, 'MIEX_ch4[u]tr', 0, 'u');
modelJ = changeRxnBounds(modelJ, 'MIEX_n2[u]tr', 0, 'u');

% Paracoccus >> original model on complete media (so almost as lb = -1000
% and ub = 1000)
% Close specific hydroxylamine uptake
modelJ =  changeRxnBounds(modelJ, 'PIEX_ham[u]tr', 0, 'l');
 
% close specific CO2 uptake (prevents growth on CO2 in unrealistic
% conditions) and restricting the total growth
modelJ =  changeRxnBounds(modelJ, 'PIEX_co2[u]tr', 0, 'l');
modelJ =  changeRxnBounds(modelJ, 'Pbio1', 0.2, 'u');

%% Auxotrophy scenario: run to get Methylococcus auxotrophic for L-lysine and Paracoccus to excrete it
% Knock out production L-tyrosine in Methylococcus
modelJ =  changeRxnBounds(modelJ, 'MTYRTA', 0, 'u');
% Open dual transport in Paracoccus
modelJ =  changeRxnBounds(modelJ, 'Prxn05301[c0]', -1000, 'l');
% Add uptake L-lysine from the environment
modelJ = addReaction(modelJ, 'Mrxn01', 'metaboliteList', {'cpd00069[u]', 'Mtyr__L[c]'}, 'stoichCoeffList', [-1; 1]);

%% Ammonia depletion: run to get growth rates and fluxes when ammonia is depleted
% Close environmental ammonia exchange in space '[u]'
modelJ =  changeRxnBounds(modelJ, 'EX_nh4[u]', 0, 'l');

%% SteadyCom retrieve information from the combined model
% Add field to joint model
modelJ.ctrs = modelM.ctrs;

% Retrieve names and ids for organism/community exchange
% reactions/metabolites
[modelJ.infoCom, modelJ.indCom] = getMultiSpeciesModelId(modelJ, nameTagsModel);
disp(modelJ.infoCom);

% Incorporate the names and indices for biomass reactions
rxnBiomassId = findRxnIDs(modelJ, {'MBIOMASS_Mcapsulatus'; 'Pbio1'});
modelJ.infoCom.spBm = {'MBIOMASS_Mcapsulatus'; 'Pbio1'};
modelJ.indCom.spBm = rxnBiomassId;

%% Running SteadyCom to calculate maximum growth and other fluxes
options = struct();
options.GRguess = 0.260; % initial guess for max. growth rate
options.GRtol = 1e-6; % tolerance for final growth rate
options.algorithm = 3; % use the default algorithm (simple guessing for bounds, followed by matlab fzero)
[sol, result] = SteadyCom(modelJ, options);
