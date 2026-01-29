%% Overview and start stuff
% Build strain specific models for each of the strains in Bloom2013 thanks
% to the expression data from Alberts et al. 2018, eLife. However that was
% only grown in YNB media. Therefore in the side need to check which genes
% will be differentially expressed between ethanol and glucose media. This
% data will be from Smith & Kruglyak 2008, Plos Biology. Enough talk
% 
clear
% addpath("H:/Rajeeva/Software/cobratoolbox/") % Adding COBRA to PATH
addpath("LocalGini")
addpath("SprintGapFiller")
initCobraToolbox(false)                             % Initializing
changeCobraSolver('gurobi', 'all')

FILTER_STRAINS = true; %If you want to filter some strains out
FILTER_STRAINS_PATH = "data/strains_0.5_ynb.csv";

%% Getting data

model = readCbModel('models/yeast8_ynb_final.mat');

expression_path = "data/expressionValues_imputed_median_antilog.csv";
gene_expression_matrix = readtable(expression_path);

% Filter out strains 
if FILTER_STRAINS
    strains = readtable(FILTER_STRAINS_PATH);
    strains = strains.Strain;

    common_strains = intersect(gene_expression_matrix.Properties.VariableNames, strains);
    columns = [{'Genes'}; common_strains];
    gene_expression_matrix = gene_expression_matrix(:, columns);
end


ge_data.value = table2array(gene_expression_matrix(:, 2:end));
ge_data.context = gene_expression_matrix.Properties.VariableNames(2:end);
ge_data.genes = gene_expression_matrix.Genes;


% Scaling bounds - To avoid floating point errors
model.lb = model.lb*10;
model.ub = model.ub*10;

optimizeCbModel(model)

%% Model params

MeM = 'GIMME';
contexts = ge_data.context;
lt =25;
ut=75;
ThS = 1; % impliying at gene level
core_reaction = [find(strcmp(model.rxns, 'r_2111'))]; % r_1761 is ethanol exchange
tol = 1e-8;
filename = 'models/strain_models_gimme/';
changeCobraSolverParams('LP','feasTol',1e-9);     % Higher values gave errors
cons_mod_rxn_id = sprintcc(model,tol);            % Blocked reactions

%% Building Models

[Models,RxnImp] = buildContextmodels(ge_data,model,MeM,contexts,ut,lt,ThS,core_reaction,filename,cons_mod_rxn_id,tol);





