% Going to make a YNB media

clear

model = readCbModel('models/yeast-GEM.mat');

%Vitamins
rxn_indices = [1141; 1041; 1238; 1384; 1390; 1083; 1431; 1440; 1466]; %Last one was 1446 before
rxns = model.rxns(rxn_indices);

model = changeRxnBounds(model, rxns, -1, 'l');
model = changeRxnBounds(model, rxns,  1000, 'u');

disp("Glucose")
optimizeCbModel(model)

model_2 = changeRxnBounds(model, "r_1714", 0, "l"); % No glucose
model_2 = changeRxnBounds(model_2, "r_1761", -1, "l");

disp("Ethanol")
optimizeCbModel(model_2)

save("models/yeast8_ynb_final.mat", "model", '-mat')


