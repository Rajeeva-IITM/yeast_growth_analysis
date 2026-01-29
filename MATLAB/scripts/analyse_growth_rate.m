% %% Step 1
% Comparing all the models' growth rate and comparing it to data from
% Bloom2013
% clear

file_loc = 'models/strain_models_sprint/';
files = cellstr(ls(strcat(file_loc, '*.mat')));
pheno_data = readtable("data/BYxRM_PhenoData.txt");
strains = pheno_data.Strain; % Extracting Strain information

albert_strains = split(files, '.');  % Strains with transcriptome data
albert_strains = albert_strains(1:end, 1);

common_strains = intersect(albert_strains, strains); % Intersection of strains
common_indices_bloom = find(contains(strains, common_strains));
common_indices_albert = find(contains(albert_strains, common_strains));

growth_rate = zeros(size(common_strains)); % Will be used to store computational growth values

%% Loop

for i=1:numel(albert_strains)
    % i = common_indices_bloom(j);
    strain_name = split(files{i}, '.');
    strain_name = strain_name{1}; % Extracting the name of the strain
    filename = strcat(file_loc, files(i));
    % if ~contains(common_strains, strain_name)
    %     continue
    % end
     model = readCbModel(filename{1});
    
     % model = changeRxnBounds(model, 'r_1714', 0, 'b'); % Setting glucose to zero

    growth_rate(i) = optimizeCbModel(model).f;
end

%% Step 2

growth_rate = growth_rate(common_indices_albert);

growth_rate_normalized = (growth_rate - mean(growth_rate))/std(growth_rate);
histogram(growth_rate_normalized, 50)
title("Normalized Growth Rate from strain specific GSMM")
xlabel("Normalized Growth")
ylabel("Counts")

ynb = pheno_data.YNB(common_indices_bloom);
ynb = fillmissing(ynb, "constant", mean(ynb, "omitnan"));
ynb_normalized = (ynb - mean(ynb))/std(ynb);

figure
histogram(ynb_normalized, 50)
title("Normalized Experimental Growth Data from Bloom 2013 et al. - YNB")
xlabel("Normalized Growth Rate")
ylabel("Counts")

%% YPD
ypd = pheno_data.YPD(common_indices_bloom);
ypd_normalized = fillmissing(ypd, "constant", mean(ypd, "omitnan"));
ypd_normalized = (ypd_normalized - mean(ypd_normalized))/std(ypd_normalized);

figure
histogram(ypd_normalized, 50)
title("Normalized Experimental Growth Data from Bloom 2013 et al. - YPD")
xlabel("Normalized Growth")
ylabel("Counts")

%% Ethanol
eth = pheno_data.Ethanol(common_indices_bloom);
eth = fillmissing(eth, "constant", mean(eth, "omitnan"));
eth_normalized = (eth - mean(eth))/std(eth);

% figure
% histogram(eth_normalized, 50)
% title("Normalized Experimental Growth Data from Bloom 2013 et al. - Ethanol")
% xlabel("Normalized Growth")
% ylabel("Counts")

