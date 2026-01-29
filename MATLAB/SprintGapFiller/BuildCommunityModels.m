function [MicComModel,removedCoreRxns] = BuildCommunityModels(Folder_path, abbr, Umodel, lbs, ubs, ConsiderOtherTranRxn, TransferCore, bioCore, media, UmodelMap, tol)
%%INPUT
%       Folder_path: Matlab cell listing the paths to all the microbial
%                    model structures in .mat format. Lower (*lb) and upper
%                    (*ub) bounds can be provided for any reaction in the 
%                    microbial models other than default bounds. Reaction 
%                    ID (*rxns) and metabolite ID (*mets) should be in same
%                    format as in the corresponding Umodel. 
%                    Metabolite IDs(*mets) should include the compartment 
%                    info(Eg: glc_D[e], pyr[c])
%
%       abbr: Matlab cell listing model abbrevations. All rxns and mets 
%             will have this prefix. Must be same order as in Folder_path
%
%       Umodel: Matlab cell listing the COBRA model structure of Universal 
%               models. If only one model is provided, then it will 
%               considered as universal model for all the microbial models,
%               else the variable 'UmodelMap' has to be provided. 
%               All the exchange reactions ID should begin with 'EX_'
%               The following fields are required in the universal model:
%                   * S - `m x n` Stoichiometric matrix
%                   * b  - `m x 1` change in concentration with time
%                   * c  - `n x 1` Linear objective coefficients
%                   * lb - `n x 1` Lower bounds on net flux
%                   * ub - `n x 1` Upper bounds on net flux
%                   * mets - metabolite IDs
%                   * rxns - reaction IDs
%
%       lbs: vector of lower bounds to all microbial biomass reactions
%
%       ubs: vector of upper bounds to all microbial biomass reactions
%
%       ConsiderOtherTranRxn: Boolean value
%               1: Transport reactions that are not present in the given 
%                  microbial model but present in the universal model will
%                  also be considered for inclusion
%               0: Transport reactions that are present in the microbial 
%                  model will only be considered for inclusion
%
%       TransferCore: Boolen value
%               1: Transport reactions that are present in the microbial
%                  model will also be considered as core reactions
%               0: Transport reactions that are present in the microbial
%                  model will not be considered as core reactions
%
%       bioCore: Boolen value
%               1: Only individual biomass reactions are considered as core
%                  reactions
%               0: All the non-transport reactions in the individual models
%                  are considered as core reactions
%
%       media: matlab structure with fields
%              *exc_rxns: list of exchange reactions (should be the union 
%                         set of all the exchange reactions in the 
%                         universal models)
%              *lb: Lower bounds of corresponding reactions in exc_rxns
%              *ub: Upper bounds of corresponding reactions in exc_rxns
%
%       UmodelMap: Array of size equal to number of microbial models. The 
%                  values define the ids of the correspoding Universal 
%                  model that has to be used 
%
%       tol: minimum absolute flux value that every reaction in the 
%            community model should carry (default: 1e-4)

%%OUTPUT
%       ComModel: The consistent community model
%       
%       removedRxns: Reactions that were removed because of inconsistency 
%                    (or) inability to carry the minimum flux (tol) in the given
%                    community and media conditions

%%AUTHOR
%       Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras

if ~exist('tol', 'var') 
    tol = 1e-4;
end

if numel(Folder_path) ~= numel(abbr)
    error('Model names has to be of the same size as the number of models in Foler_path')
end

if numel(Umodel) > 1
    if ~exist('UmodelMap','var')
        error('The variable ''UmodelMap'' has to be provided')
    end
end

if exist('UmodelMap','var')
    if numel(UmodelMap) ~= numel(Folder_path)
        error('UmodelMap has to be of the same size as the number of models in Foler_path')
    end
else
    UmodelMap = ones(numel(Folder_path),1);
end

n_models = numel(Folder_path); % number of models
n_Umodels = numel(Umodel); % number of universal models
Umodel_new = Umodel;

exc_rxns={};exc_rxnFormulas={};
for j=1:n_Umodels
    Utemp = Umodel_new{j};
    Curr_exc=Utemp.rxns(startsWith(Utemp.rxns,'EX_'));
    exc_rxns=[exc_rxns;Curr_exc];
    form =printRxnFormula(Utemp,'rxnAbbrList',Curr_exc,'printFlag',false);
    exc_rxnFormulas = [exc_rxnFormulas;form];
    Utemp=removeRxns(Utemp,exc_rxns); % remove all exchange rxns
    Umodel_new{j} = Utemp;
end

[~,ia,~] = unique(exc_rxns);
exc_rxns = exc_rxns(ia);
exc_rxnFormulas = exc_rxnFormulas(ia);

S=[];lb=[];ub=[];c=[];b=[];rxns=[];mets=[];core=[];
for i=1:n_models
    load(Folder_path{i})
    CuUmodel = Umodel_new{UmodelMap(i)};
    if ConsiderOtherTranRxn
        UmodelTemp = CuUmodel;
    else
        % removing all the transport reactions that are not there in the
        % given microbial model
        UmodelTemp = removeRxns(CuUmodel,setdiff(CuUmodel.rxns(getTransRxns(CuUmodel)),...
            model.rxns(getTransRxns(model))));
    end
    
    if bioCore
        coreTemp = zeros(numel(UmodelTemp.rxns),1); % all reactions are non-core except biomass reactions
    elseif TransferCore
        coreTemp = ismember(UmodelTemp.rxns,model.rxns);
    else
        coreTemp = ismember(UmodelTemp.rxns,model.rxns(setdiff([1:numel(model.rxns)],getTransRxns(model))));
    end
    
    % Adding the biomass reaction
    BioForm = printRxnFormula(model,model.rxns(find(model.c)),0);
    UmodelTemp=addReaction(UmodelTemp,'bio','reactionFormula',BioForm{1},...
        'lowerBound',lbs(i),'upperBound',ubs(i));
    
    core = [core;coreTemp;1]; % one refers to the biomass reaction
    new_rxns = cellfun(@(x)rename_rxns(x,abbr{i}),UmodelTemp.rxns,'uni',false);
    rxns = [rxns;new_rxns];
    new_mets=cellfun(@(x)rename_mets(x,abbr{i}),UmodelTemp.mets,'uni',false);
    mets=[mets;new_mets];
    S = blkdiag(S,UmodelTemp.S);c=[c;UmodelTemp.c];b=[b;UmodelTemp.b];
    new_lb = UmodelTemp.lb;new_ub = UmodelTemp.ub;
    [loca,locb] = ismember(model.rxns,UmodelTemp.rxns);
    locb = locb(locb~=0);
    new_lb(locb)=model.lb(loca);new_ub(locb)=model.ub(loca);
    lb=[lb;new_lb];ub=[ub;new_ub];
end

% merging the extracellular metabolite rows and removing the extra
% metabolites
[Umets,~,ix] = unique(mets);
counts = accumarray(ix,1).';
counts = counts';
for j=1:numel(counts)
    if counts(j)>1
        ids = find(ismember(mets,Umets{j}));
        S(ids(1),:)=sum(S(ids,:),1);
        S(ids(2:end),:)=[];
        b(ids(2:end))=[];
        mets(ids(2:end))=[];
    end
end
    
ComModel=struct();
ComModel.S=S;ComModel.lb=lb;
ComModel.ub=ub;ComModel.c=c;
ComModel.b=b;ComModel.mets=mets;
ComModel.rxns=rxns;

for i=1:numel(media.exc_rxns)
    id = find(ismember(exc_rxns,media.exc_rxns{i}));
    if ~isempty(id)
        ComModel=addReaction(ComModel,exc_rxns{id},'reactionFormula',exc_rxnFormulas{id},...
            'lowerBound',media.lb(i),'upperBound',media.ub(i));
    end
end
core = [core;zeros(sum(ismember(exc_rxns,media.exc_rxns)),1)]; % exchange reactions are non-core

[MicComModel,removedCoreRxns] = SprintGapFiller(ComModel,core,tol);

end


function rxns = rename_rxns(a,ABR)
rxns = [ABR,'_',a];
end

function mets = rename_mets(a,ABR)
if ~strcmp(a(end-1),'e')&&~contains(a,'biomass')
    mets = [ABR,'_',a];
else
    mets = a;
end
end

function ids = getTransRxns(model)
% this part of code is adapted from findTransRxns.m in COBRA toolbox
[~,compSymbols]=arrayfun(@parseMetNames, model.mets);
ids=[];
for i = 1:numel(model.rxns)
    % get compartment symbols for each rxn
    compSymbolsTmp=compSymbols(model.S(:,i)~=0);
    % if there's more than 1 compartment involved, it's a transport rxn
    if length(unique(compSymbolsTmp))>1
        ids=[ids;i];
    end
end
end
