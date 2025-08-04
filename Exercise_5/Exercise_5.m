%% 📑 Exercise #5: NICEgame with biolog data
% In this exercise we use the NICEgame workflow to gapfill a model using biolog 
% growth data. Since the model used for this exercise is part of an unpublished 
% dataset, the model file will only be uploaded after publication.
%%   Load model, reactionn database, experimental data

% Load DB
load('keggModel.mat') % KEGG DB
% Load draft model
load('pscfSuccinateM9Vitconstrained.mat') % draft model to curate

% Load biolog data
biolog = readtable('BioLogSubset_pscf.xlsx');

% Check model structure, looks like a normal GEM non-TFA
open pscf
% Check DB structure, also looks like a normal GEM non-TFA
open keggmodel

%% FBA-based metabolite essentiality

% check growth of WT model
solWT = optimizeCbModel(pscf);

% Reaction ID to knock out
rxnToKO = 'EXC_BOTH_C00042_e';
% Find index of the reaction in the model
rxnIdx = find(strcmp(pscf.rxns, rxnToKO));
pscf.lb(rxnIdx) = 0;
%pscf.ub(rxnIdx) = 0;
solKO = optimizeCbModel(pscf);

% Create a loop where we close/open exchange for different carbon sources
% and check for growth, start with FBA
% Preallocate results
biolog.ModelGrowth = zeros(height(biolog), 1);
biolog.FBAGrowthRate = nan(height(biolog), 1);  % New column to store growth rate

for i = 1:height(biolog)
    model = pscf;  % Fresh model for each run
    metName = biolog{i, 1}{1};    % Metabolite name
    keggID  = biolog{i, 2}{1};    % KEGG ID like 'C00062'
    
    % Construct exchange reaction ID
    exchRxnID = ['EXC_BOTH_' keggID '_e'];

    % Find the reaction index
    rxnIdx = find(strcmp(model.rxns, exchRxnID), 1);

    if isempty(rxnIdx)
        fprintf('Reaction %s not found — skipping %s\n', exchRxnID, metName);
        continue;
    end

    % Enable uptake
    model.lb(rxnIdx) = -10;

    % Run FBA
    sol = optimizeCbModel(model);

    % Store results
    if ~isempty(sol.f) && sol.f > 1e-3
        biolog.ModelGrowth(i) = 1;
        biolog.FBAGrowthRate(i) = sol.f;
    else
        biolog.ModelGrowth(i) = 0;
        biolog.FBAGrowthRate(i) = sol.f;
    end

    fprintf('%-30s | Growth: %d | Rate: %.4f\n', ...
        metName, biolog.ModelGrowth(i), biolog.FBAGrowthRate(i));
end

%% Compare metabolite essentiality prediction against biolog data

% evaluate
% Create new column for label
biolog.PredictionLabel = strings(height(biolog), 1);  % initialize as string array
% Loop through each row and label
for i = 1:height(biolog)
    trueVal = biolog.pscf(i);          % Experimental value (0 or 1)
    predVal = biolog.ModelGrowth(i);   % Predicted value (0 or 1)

    if trueVal == 1 && predVal == 1
        biolog.PredictionLabel(i) = "TP";
    elseif trueVal == 0 && predVal == 0
        biolog.PredictionLabel(i) = "TN";
    elseif trueVal == 1 && predVal == 0
        biolog.PredictionLabel(i) = "FN";
    elseif trueVal == 0 && predVal == 1
        biolog.PredictionLabel(i) = "FP";
    else
        biolog.PredictionLabel(i) = "Unclear";  % Shouldn't happen if all values are 0 or 1
    end
end

% get number of FP/FN/TP/TN, calculate MCC
FP = height(biolog_original(strcmp(biolog_original.PredictionLabel,"FP"),:)); % 18 FP
FN = height(biolog_original(strcmp(biolog_original.PredictionLabel,"FN"),:)); % 2 FN
TP = height(biolog_original(strcmp(biolog_original.PredictionLabel,"TP"),:)); % 4 TP
TN = height(biolog_original(strcmp(biolog_original.PredictionLabel,"TN"),:)); % 4 TN
numerator = (TP * TN) - (FP * FN);
denominator = sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN));
MCC_original = numerator/denominator; % -0.1515, very bad original MCC!
%% Gapfill for false negative metabolite

% OK, lets gapfill for the metabolite D-Galacturonic-Acid, which the microbe 
% should be able to grow on but currently does not. 

% First we need to allow uptake of the FN carbon source
metID = 'C00333_e';
metIdx = find(strcmp(pscf.mets, metID));
exchangeRxns = pscf.rxns(rxnIdxs(contains(pscf.rxns(rxnIdxs), 'EX')));
exchRxnID = 'EXC_BOTH_C00333_e';
exchIdx = find(strcmp(pscf.rxns, exchRxnID));
pscf.lb(exchIdx) = -10;  

% % Check what are the biomass buidling blocks that cannot be produced
% % error when trying to run this code
% obj = find(pscf.c); % find objective function (biomass)
% bbbs_m = pscf.mets(find(pscf.S(:,find(pscf.c))<0)); % find metabolites consumed by growth
% r = length(pscf.rxns); % define number of rxns
% model_bbb = addDemandReaction(pscf,bbbs_m); % add demand reactions for bbbs in new model
% for i = 1:length(bbbs_m) % for each rxns
%     modeli = model_bbb; % define new model
% modeli.c = zeros(length(modeli.rxns),1); % clear objective function
% modeli.c(r+i) = 1; % define objective function to be one of newly added demand rxns
% sol = optimizeCbModel(modeli,'max'); % solve fba
% if sol.f < 10^-6 % assume no growth if value is very small
%     sol.f = 0;
% end
% sol_bbb(i,1) = sol.f; % check if bbb can be produced
% end
% problematic_BBBs = bbbs_m(sol_bbb == 0); %
% rxnID = strcat('DM_',problematic_BBBs);
% 
% %sourceModel=model;

% prepare for gapfilling
[GFmodel, conflict] = PrepareForGapFilling(pscf, {keggmodel},'c', 0,1,{},[],thermo_data);

% gapfilling
GFmodeli = GFmodel;

GFmodeli.var_lb(find(ismember(GFmodel.varNames,strcat('F_',pscf.rxns(find(pscf.c)))))) = 0.01;
GFmodeli.var_lb(find(ismember(GFmodel.varNames,strcat('NF_',pscf.rxns(find(pscf.c)))))) = 0.01;

%generate the alternative solutions
NumAlt = 10;
tagMin = 0; % solutions of min size
tagMinP1 = 0;  % solutions of min +1 size
tagSave = 1;
time = 600; % time limit for the solver in seconds
filename = './pscf_gf.mat';
indUSE = GFmodeli.indUSE;
mingrRate = 0.01; % force growth rate if you want to gap-fill for each bbb then mingrRate = 0
%rxnsToCheck = {'EXC_BOTH_C00042_e'}; % rxnsToCheck = GFmodeli.rxnRescued;
%[ActRxns, GFSumRes, DPsAll] = gapFilling(GFmodeli,indUSE,NumAlt,mingrRate,rxnsToCheck,tagMin,time,filename);
[resultStat,ActRxns,DPsAll] = gfbiomass(GFmodeli,GFmodeli.indUSE,NumAlt,time,tagMin,tagMinP1,tagSave,filename);

save(strcat('./GF_',organism_id,'.mat'),'ActRxns','DPsAll')
%% Evaluate gapfilling solution

%% Evaluate gapfilling solution #1

if ~exist('uptakeBound', 'var')
    uptakeBound = -10;
end

nAlts = size(ActRxns{1}, 1);
biolog_original = biolog
%biologResults = cell(nAlts, 1);

%for k = 1:nAlts
    fprintf('\nEvaluating Alternative #%d\n', k);

    % Apply alternative solution to model
    alt = ActRxns{1}{k,1};
    modelk = pscf;
    % remove rxnGeneMat from the model to prevent errors, should not affect
    % any results since no gene essentiality
    modelk = rmfield(modelk, 'rxnGeneMat');

    %add reaction 'R10846_c'
    [modelk_added] = addReaction(modelk,'R10846_c','C20889_c <=> C20896_c') % remove underscore
    % i think this is the command below that ran
    %[modelk_added] = addReaction(modelk,'R10846_c',{'C20889_c','C20896_c'},[-1 1],true,-10,10,0,'no_subsystem_info','')
    
    % use FBA for more appropriate comparison to baseline model prediction
    sol = optimizeCbModel(modelk_added) 

    
    %modelk.ub(ismember(modelk.rxns, alt))= 10;
    %modelk.lb(ismember(modelk.rxns, alt)) = -10;

    %modelk.var_ub(ismember(modelk.varNames, strcat('F_', alt))) = 50;
    %modelk.var_ub(ismember(modelk.varNames, strcat('R_', alt))) = 50;
    %modelk.var_ub(ismember(modelk.varNames, strcat('NF_', alt))) = 50;
    %modelk.var_lb(ismember(modelk.varNames, strcat('NF_', alt))) = -50;
    %solt = solveTFAmodelCplex(modelk)
%end


% First we need to close uptake of the FN carbon source
exchRxnID = 'EXC_BOTH_C00333_e';
exchIdx = find(strcmp(modelk_added.rxns, exchRxnID));
modelk_added.lb(exchIdx) = 0;
% update gapfilled reaction bounds from -10 to 10 to -50 to 50
% evaluate new essentiality with modelk_added

for i = 1:height(biolog)
    model = modelk_added;  % Fresh model for each run
    
    metName = biolog{i, 1}{1};    % Metabolite name
    keggID  = biolog{i, 2}{1};    % KEGG ID like 'C00062'
    
    % Construct exchange reaction ID
    exchRxnID = ['EXC_BOTH_' keggID '_e'];

    % Find the reaction index
    rxnIdx = find(strcmp(model.rxns, exchRxnID), 1);

    if isempty(rxnIdx)
        fprintf('Reaction %s not found — skipping %s\n', exchRxnID, metName);
        continue;
    end

    % Enable uptake
    model.lb(rxnIdx) = -10;

    % Run FBA
    sol = optimizeCbModel(model);

    % Store results
    if ~isempty(sol.f) && sol.f > 1e-3
        biolog.ModelGrowth(i) = 1;
        biolog.FBAGrowthRate(i) = sol.f;
    else
        biolog.ModelGrowth(i) = 0;
        biolog.FBAGrowthRate(i) = sol.f;
    end

    fprintf('%-30s | Growth: %d | Rate: %.4f\n', ...
        metName, biolog.ModelGrowth(i), biolog.FBAGrowthRate(i));
end

% evaluate
% Create new column for label
biolog.PredictionLabel = strings(height(biolog), 1);  % initialize as string array
% Loop through each row and label
for i = 1:height(biolog)
    trueVal = biolog.pscf(i);          % Experimental value (0 or 1)
    predVal = biolog.ModelGrowth(i);   % Predicted value (0 or 1)

    if trueVal == 1 && predVal == 1
        biolog.PredictionLabel(i) = "TP";
    elseif trueVal == 0 && predVal == 0
        biolog.PredictionLabel(i) = "TN";
    elseif trueVal == 1 && predVal == 0
        biolog.PredictionLabel(i) = "FN";
    elseif trueVal == 0 && predVal == 1
        biolog.PredictionLabel(i) = "FP";
    else
        biolog.PredictionLabel(i) = "Unclear";  % Shouldn't happen if all values are 0 or 1
    end
end

% get number of FP/FN/TP/TN, calculate MCC
FP = height(biolog(strcmp(biolog.PredictionLabel,"FP"),:)); % 18 FP
FN = height(biolog(strcmp(biolog.PredictionLabel,"FN"),:)); % 1 FN, fixed the gapfilling target
TP = height(biolog(strcmp(biolog.PredictionLabel,"TP"),:)); % 5 TP
TN = height(biolog(strcmp(biolog.PredictionLabel,"TN"),:)); % 4 TN
numerator = (TP * TN) - (FP * FN);
denominator = sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN));
MCC = numerator/denominator; % MCC very small 0.0162
% modest improvement in MCC, at least now the value is positive

%% Fix another false negative

% Next, try fixing FN tween-80 (Kegg compound ID: C11625, LCSB ID: 1467972114)
% save modelk_addded backup
%modelk_backup = modelk_added;
% First we need to allow uptake of the FN carbon source
metID = 'C11625_e';
%metIdx = find(strcmp(modelk_added.mets, metID));
%exchangeRxns = modelk_added.rxns(rxnIdxs(contains(modelk_added.rxns(rxnIdxs), 'EX')));
exchRxnID = 'EXC_BOTH_C11625_e';
exchIdx = find(strcmp(modelk_added.rxns, exchRxnID));
modelk_added.lb(exchIdx) = -10;  

% prepare for gapfilling
[GFmodel, conflict] = PrepareForGapFilling(modelk_added, {keggmodel},'c', 0,1,{},[],thermo_data);

% gapfilling
GFmodeli = GFmodel;

GFmodeli.var_lb(find(ismember(GFmodel.varNames,strcat('F_',modelk_added.rxns(find(modelk_added.c)))))) = 0.01;
GFmodeli.var_lb(find(ismember(GFmodel.varNames,strcat('NF_',modelk_added.rxns(find(modelk_added.c)))))) = 0.01;

%generate the alternative solutions
NumAlt = 10;
tagMin = 0; % solutions of min size
tagMinP1 = 0;  % solutions of min +1 size
tagSave = 1;
time = 600; % time limit for the solver in seconds
filename = './pscf_gf.mat';
indUSE = GFmodeli.indUSE;
mingrRate = 0.01; % force growth rate if you want to gap-fill for each bbb then mingrRate = 0
%rxnsToCheck = {'EXC_BOTH_C00042_e'}; % rxnsToCheck = GFmodeli.rxnRescued;
%[ActRxns, GFSumRes, DPsAll] = gapFilling(GFmodeli,indUSE,NumAlt,mingrRate,rxnsToCheck,tagMin,time,filename);
[resultStat,ActRxns2,DPsAll] = gfbiomass(GFmodeli,GFmodeli.indUSE,NumAlt,time,tagMin,tagMinP1,tagSave,filename);

%save(strcat('./GF_',organism_id,'.mat'),'ActRxns','DPsAll')

% looks like all of the alternative solutions involve degradation of other media
% compounds rather than directly degrading tween-80
% reactions for the direct degradation of this compound cannot be found on
% ATLASx
%% Fix a false positive

% as an exercise could try to further curate this model by fixing some FP
% instead. this involves blocking reactions 1 at a time and then evaluating
% which one(s) results in no growth under FP carbon source environment.
% Lets try with FP carbon source L-Arginine

% close uptake of tween-80 
exchRxnID = 'EXC_BOTH_C11625_e';
exchIdx = find(strcmp(modelk_added.rxns, exchRxnID));
modelk_added.lb(exchIdx) = 0;  

% open uptake of L-Arginine
exchRxnID = 'EXC_BOTH_C00062_e';
exchIdx = find(strcmp(modelk_added.rxns, exchRxnID));
modelk_added.lb(exchIdx) = -10;

% single reaction essentiality to suggest alternatives for KOing to fix FP
[grRatioFBA,grRateKOFBA,grRateWTFBA,hasEffectFBA,delRxnsFBA,fluxSolutionFBA] = singleRxnDeletion(modelk_added);
% Find indices of lethal knockouts
lethalIdx = find(grRatioFBA < 0.1);  %
% Get the corresponding reaction names
lethalRxns = delRxnsFBA(lethalIdx);
% Display
disp('Lethal reactions:');
printRxnFormula(modelk_added, lethalRxns, true, true,true);

% trivial solution is of course to block import of arginine, in a real
% situation perhaps some literature search and genomic analysis could
% reveal whether arginine transporters are present or not. to select the
% "best" of 171 alternative rxn KOs, would need to evaluate the metabolite
% essentiality for each KO and see which yields best MCC

results = struct();

for i = 1:length(lethalRxns)
    rxnID = lethalRxns{i};
    modelKO = modelk_added;

    % Knockout by setting bounds to zero
    rxnIdx = find(strcmp(modelKO.rxns, rxnID));
    if isempty(rxnIdx)
        warning('Reaction %s not found in model', rxnID);
        %continue;
    end
    modelKO.lb(rxnIdx) = 0;
    modelKO.ub(rxnIdx) = 0;

    TP = 0; TN = 0; FP = 0; FN = 0;

     for j = 1:height(biolog)
        metID = biolog.KEGGID{j};

        % Find exchange reaction
        exRxnIdx = find(contains(modelKO.rxns, metID) & startsWith(modelKO.rxns, 'EXC_'));
        if isempty(exRxnIdx)
            %continue; % skip if exchange rxn not found
        end

        % Temporarily allow uptake of this carbon source, block others
        modelTemp = modelKO;
        
        % Allow this carbon source
        modelTemp.lb(exRxnIdx) = -10;

        % Block arginine
        exchRxnID = 'EXC_BOTH_C00062_e';
        exchIdx = find(strcmp(modelTemp.rxns, exchRxnID));
        modelTemp.lb(exchIdx) = 0;

        % Solve the model
        sol = optimizeCbModel(modelTemp);

        grows = sol.f > 1e-6;
        predVal = double(grows);
        trueVal = biolog.pscf(j);
    
        % Update TP/TN/FP/FN
        if trueVal == 1 && predVal == 1
            TP = TP + 1;
        elseif trueVal == 0 && predVal == 0
            TN = TN + 1;
        elseif trueVal == 1 && predVal == 0
            FN = FN + 1;
        elseif trueVal == 0 && predVal == 1
            FP = FP + 1;
        end
     end
    
      % Compute MCC
    numerator = TP * TN - FP * FN;
    denominator = sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN));
    if denominator == 0
        MCC = NaN;
    else
        MCC = numerator / denominator;
    end

    % Store result
    results(i).rxn = rxnID;
    results(i).TP = TP;
    results(i).TN = TN;
    results(i).FP = FP;
    results(i).FN = FN;
    results(i).MCC = MCC;
end

% looks like most solutions result in very few FP but at the cost of more
% FN and less TP. manually verify KO results for R08162, indeed looks like
% it cannot grow on putrecine = TN

% 3 best solutions 'R00768', 'R02060', 'R05332' result in MCC 0.3685
best_solns = {'R00768', 'R02060', 'R05332'}
printRxnFormula(modelk_added, best_solns, true, true,true);

% in practice, continue and try fixing new false negatives, as exercise
% maybe end here. could also check glutamate is in the new FN, check if KO
% alternatives are part of arginine degradation

%%


%% Convert to TFA
% add description to avoid error 
pscf.description='pscf model from omar';

% fix metCharge to avoid error
if iscell(pscf.metCharge)
    newMetCharge = nan(length(pscf.metCharge), 1);  % Preallocate

    for i = 1:length(pscf.metCharge)
        val = pscf.metCharge{i};

        % Check if it's a number directly
        if isnumeric(val)
            newMetCharge(i) = val;

        % Try parsing if it's a string
        elseif ischar(val) || isstring(val)
            numVal = str2double(val);
            if ~isnan(numVal)
                newMetCharge(i) = numVal;
            end

        % Otherwise leave as NaN
        end
    end

    pscf.metCharge = newMetCharge;
end

% run functions to generate TFA-ready model
prepped_m = prepModelforTFA(pscf, ReactionDB, pscf.CompartmentData);
solFBA = optimizeCbModel(pscf);
min_obj = roundsd(0.5*solFBA.f, 2, 'floor');
tmp = convToTFA(prepped_m, ReactionDB, [], 'DGo', [], min_obj);
% Add net flux variables, which are equal to forwards flux - backwards flux
% NF_rxn = F_rxn - B_rxn
this_tmodel = addNetFluxVariables(tmp);
% Solve tFA
soltFA = solveTFAmodelCplex(this_tmodel);

% Make sure succinate is blocked
rxnID = 'EXC_BOTH_C00042_e';
varName = ['F_' rxnID];
vIdx = strcmp(this_tmodel.varNames, varName);
if any(vIdx)
    this_tmodel.var.lb(vIdx) = 0;  % Block uptake
    fprintf('Blocked uptake for %s\n', rxnID);
else
    warning('Could not find TFA flux variable for %s', rxnID);
end

% Create a loop where we close/open exchange for different carbon sources
% and check for growth, now with TFA
% Preallocate results
biolog.ModelGrowthTFA = zeros(height(biolog), 1);
flkBase = this_tmodel;
for i = 1:height(biolog)
    metName = biolog{i, 1}{1};    % Metabolite name (not used in ID)
    keggID  = biolog{i, 2}{1};    % KEGG ID like 'C00062'
    
    % Construct exchange reaction ID
    exchRxnID = ['EXC_BOTH_' keggID '_e'];

    % Create a fresh copy of the model for each test
    flk = flkBase;

    % Find the reaction index
    rxnIdx = find(strcmp(flk.rxns, exchRxnID), 1);

    if isempty(rxnIdx)
        fprintf('Reaction %s not found — skipping %s\n', exchRxnID, metName);
        continue;
    end

    % Enable uptake: set lower bound to -10 (or your model's default)
    flk.lb(rxnIdx) = -10;

    % Run FBA
    sol = solveTFAmodelCplex(flk);

    % Store whether the model grows
    if ~isempty(sol.val) && sol.val > 1e-6
        biolog.ModelGrowth(i) = 1;
    else
        biolog.ModelGrowth(i) = 0;
    end

    fprintf('%-30s | Growth: %d\n', metName, biolog.ModelGrowth(i));
end
%% Miscellaneous
% Code to check media composition across various models


% Check media composition
exRxnsIdx = find(startsWith(pscf.rxns, 'EXC_'));
fprintf('Exchange Reaction     | Uptake | Metabolite Name\n');
fprintf('--------------------------------------------------\n');

for i = 1:length(exRxnsIdx)
    idx = exRxnsIdx(i);
    rxnID = pscf.rxns{idx};
    lb = pscf.lb(idx);

    uptake = lb < 0;
    if ~uptake
        continue;  % Skip if uptake not allowed
    end

    % --- Extract metabolite ID from reaction ID ---
    parts = split(rxnID, '_');
    metID = join(parts(end-1:end), '_');  % e.g., 'C00243_e'
    metID = metID{1};

    % --- Find index in flk.mets ---
    metIdx = find(strcmp(pscf.mets, metID), 1);  % exact match

    % --- Get metabolite name from bigg field ---
    %if ~isempty(metIdx) && isfield(flk, 'bigg') && isfield(flk.bigg, 'metabolite')
        metName = pscf.bigg.metabolite{metIdx};
    %else
    %    metName = '[Unknown]';
    %end

    if iscell(metName)
        metNameStr = metName{1};
    else
        metNameStr = metName;  % already a char array or string
    end
    
    fprintf('%-22s |  %5s | %s\n', rxnID, 'YES', metNameStr);
end

%%

% Check media composition
exRxnsIdx = find(startsWith(modelk_added.rxns, 'EXC_'));
fprintf('Exchange Reaction     | Uptake | Metabolite Name\n');
fprintf('--------------------------------------------------\n');

for i = 1:length(exRxnsIdx)
    idx = exRxnsIdx(i);
    rxnID = modelk_added.rxns{idx};
    lb = modelk_added.lb(idx);

    uptake = lb < 0;
    if ~uptake
        continue;  % Skip if uptake not allowed
    end

    % --- Extract metabolite ID from reaction ID ---
    parts = split(rxnID, '_');
    metID = join(parts(end-1:end), '_');  % e.g., 'C00243_e'
    metID = metID{1};

    % --- Find index in flk.mets ---
    metIdx = find(strcmp(modelk_added.mets, metID), 1);  % exact match

    % --- Get metabolite name from bigg field ---
    %if ~isempty(metIdx) && isfield(flk, 'bigg') && isfield(flk.bigg, 'metabolite')
        metName = modelk_added.bigg.metabolite{metIdx};
    %else
    %    metName = '[Unknown]';
    %end

    if iscell(metName)
        metNameStr = metName{1};
    else
        metNameStr = metName;  % already a char array or string
    end
    
    fprintf('%-22s |  %5s | %s\n', rxnID, 'YES', metNameStr);
end

%%


% Check media composition modelTemp
exRxnsIdx = find(startsWith(modelTemp.rxns, 'EXC_'));
fprintf('Exchange Reaction     | Uptake | Metabolite Name\n');
fprintf('--------------------------------------------------\n');

for i = 1:length(exRxnsIdx)
    idx = exRxnsIdx(i);
    rxnID = modelTemp.rxns{idx};
    lb = modelTemp.lb(idx);

    uptake = lb < 0;
    if ~uptake
        continue;  % Skip if uptake not allowed
    end

    % --- Extract metabolite ID from reaction ID ---
    parts = split(rxnID, '_');
    metID = join(parts(end-1:end), '_');  % e.g., 'C00243_e'
    metID = metID{1};

    % --- Find index in flk.mets ---
    metIdx = find(strcmp(modelTemp.mets, metID), 1);  % exact match

    % --- Get metabolite name from bigg field ---
    %if ~isempty(metIdx) && isfield(flk, 'bigg') && isfield(flk.bigg, 'metabolite')
        metName = modelTemp.bigg.metabolite{metIdx};
    %else
    %    metName = '[Unknown]';
    %end

    if iscell(metName)
        metNameStr = metName{1};
    else
        metNameStr = metName;  % already a char array or string
    end
    
    fprintf('%-22s |  %5s | %s\n', rxnID, 'YES', metNameStr);
end


%%

% Check media composition model
exRxnsIdx = find(startsWith(model.rxns, 'EXC_'));
fprintf('Exchange Reaction     | Uptake | Metabolite Name\n');
fprintf('--------------------------------------------------\n');

for i = 1:length(exRxnsIdx)
    idx = exRxnsIdx(i);
    rxnID = model.rxns{idx};
    lb = model.lb(idx);

    uptake = lb < 0;
    if ~uptake
        continue;  % Skip if uptake not allowed
    end

    % --- Extract metabolite ID from reaction ID ---
    parts = split(rxnID, '_');
    metID = join(parts(end-1:end), '_');  % e.g., 'C00243_e'
    metID = metID{1};

    % --- Find index in flk.mets ---
    metIdx = find(strcmp(model.mets, metID), 1);  % exact match

    % --- Get metabolite name from bigg field ---
    %if ~isempty(metIdx) && isfield(flk, 'bigg') && isfield(flk.bigg, 'metabolite')
        metName = model.bigg.metabolite{metIdx};
    %else
    %    metName = '[Unknown]';
    %end

    if iscell(metName)
        metNameStr = metName{1};
    else
        metNameStr = metName;  % already a char array or string
    end
    
    fprintf('%-22s |  %5s | %s\n', rxnID, 'YES', metNameStr);
end

%% check media composition for gfmodel

exRxnsIdx = find(startsWith(GFmodel.rxns, 'EX_'));
fprintf('Exchange Reaction     | Uptake | Metabolite Name\n');
fprintf('--------------------------------------------------\n');

for i = 1:length(exRxnsIdx)
    idx = exRxnsIdx(i);
    rxnID = GFmodel.rxns{idx};
    lb = GFmodel.lb(idx);

    uptake = lb < 0;
    if ~uptake
        continue;  % Skip if uptake not allowed
    end

    % --- Extract metabolite ID from reaction ID ---
    parts = split(rxnID, '_');
    if length(parts) < 3
        metID = '[Unknown]';
    else
        metID = join(parts(end-1:end), '_');  % e.g., 'glc__D_e'
        metID = metID{1};
    end

    % --- Find index in mets ---
    metIdx = find(strcmp(GFmodel.mets, metID), 1);  % exact match

    if isempty(metIdx)
        metNameStr = '[Unknown]';
    else
        metName = GFmodel.mets{metIdx};
        if iscell(metName)
            metNameStr = metName{1};
        else
            metNameStr = metName;
        end
    end

    fprintf('%-22s |  %5s | %s\n', rxnID, 'YES', metNameStr);
end
%% 
%