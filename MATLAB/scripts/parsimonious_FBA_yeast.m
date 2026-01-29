% Parsimonious FBA to understand the important genes in the model. Obtain
% for high and low growing strains and compare the results.

clear
file_loc = 'models/strain_models_gimme/';
files = cellstr(ls(strcat(file_loc, '*.mat')));

albert_strains = split(files, '.'); % Strains with transcriptome data
albert_strains = albert_strains(1:end, 1);
results_array_genes = cell(size(albert_strains));
results_array_reactions = cell(size(albert_strains));

% load ynb_strains\concordant_strains_ynb.mat concordant_strains % Strains that show common phenotype in all the cases
concordant_strains = readtable("data/concordant_strains_ynb_gimme.csv").Strain;

% gene_class_strains = struct();
% rxn_class_strains = struct();

% environment = getEnvironment();

for i=1:numel(albert_strains)
    % restoreEnvironment(environment)
    strain_name = split(files{i}, '.');
    strain_name = strain_name{1}; % Extracting the name of the strain
    filename = strcat(file_loc, files(i));

    % if ~contains(concordant_strains, strain_name)
    %     continue
    % end
    disp(strain_name)
    
    % model = readCbModel(filename{1});
    % 
    % [GeneClasses, RxnClasses, ModelIrrevFM] = pFBA(model, 'tol', 1e-8);
    % 
    % results_array_genes{i} = GeneClasses;
    % results_array_reactions{i} = RxnClasses;

    try
        model = readCbModel(filename{1});

        [GeneClasses, RxnClasses, ModelIrrevFM] = pFBA(model, 'tol', 1e-8);

        results_array_genes{i} = GeneClasses;
        results_array_reactions{i} = RxnClasses;

    catch 
        warning('Problematic strain: %s', strain_name)
    end

    % break

end

gene_class_strains = cell2struct( results_array_genes, albert_strains, 1);
rxn_class_strains = cell2struct(results_array_reactions, albert_strains, 1);

save results/pFBA_genes_and_reactions_all_gimme.mat gene_class_strains rxn_class_strains