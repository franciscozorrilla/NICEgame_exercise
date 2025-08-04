%% ⚖️ Exercise #4: Reproduction of result from NICEgame paper (AMAOTr)
% We previosuly identified AMAOTr as target for gapfilling. Now we knockout 
% this reaction in iML1515, merge with the ATLAS database, generate alternative 
% solutions, and evaluate them.
%% Load model and ATLAS subset reaction database

% iML1515_KEGG_final_thermo, mapped manually by Eva
load('../Exercise_3/iML1515_KEGG_final_thermo.mat') 
model = ttmodel;

% define organism_id for filenaming later on
organism_id = 'iML1515';

% load ATLAS subset
load('../Exercise_3/reducedATLAS_ecoli_yeast.mat') % ATLAS ecoli yeast metabolite subset
DBmodel = reducedATLAS;

% load also thermo database
tmp = load('../Exercise_2/thermo_data.mat');
ReactionDB = tmp.DB_AlbertyUpdate;
clear tmp

% load essentiality data
gene_table = readtable('../Exercise_3/gene_table.tsv', 'FileType', 'text', 'Delimiter', '\t');
essentiality_table = readtable('../Exercise_3/experimental_essentiality.tsv', 'FileType', 'text', 'Delimiter', '\t');
% remove ambiguous genes
essentiality_table = essentiality_table(strcmp(essentiality_table.Unclear, 'FALSE'), :);
% merge tables
essentiality_table = innerjoin(essentiality_table, gene_table, 'Keys', 'Gene');
% pull out essential & nonessential genes cell strucutre
essential = essentiality_table(strcmp(essentiality_table.Essential,"TRUE"),:);essential=essential.Genes;
nonEssential = essentiality_table(strcmp(essentiality_table.Essential,"FALSE"),:);nonEssential=nonEssential.Genes;

%% Knockout target gapfilling reaction AMAOTr

% check growth of WT model
solWT = optimizeCbModel(model);

rxnsToKO = {'AMAOTr'}; % to reproduce result from NICEgame paper fig 3a, also identified in exerrcise 3

% find reaction and block it
f = find(ismember(model.rxns,rxnsToKO));
model.lb(f) = 0;
model.ub(f) = 0;

% check growth of KO model
solKO = optimizeCbModel(model);

% looks like it can grow in TFA
solKOTFA = optimizeThermoModel(model);

% confirm by printing TFA-associated constraints for this reaction
rxnID = 'AMAOTr';
% Find indices of associated variables
F_idx = find(strcmp(model.varNames, ['F_', rxnID]));
R_idx = find(strcmp(model.varNames, ['R_', rxnID]));
NF_idx = find(strcmp(model.varNames, ['NF_', rxnID]));  % if net flux is present
% Print bounds
fprintf('Bounds for reaction %s:\n', rxnID);
if ~isempty(F_idx)
    fprintf('  F_ (forward):  [%g, %g]\n', model.var_lb(F_idx), model.var_ub(F_idx));
end
if ~isempty(R_idx)
    fprintf('  R_ (reverse):  [%g, %g]\n', model.var_lb(R_idx), model.var_ub(R_idx));
end
if ~isempty(NF_idx)
    fprintf('  NF_ (net flux):  [%g, %g]\n', model.var_lb(NF_idx), model.var_ub(NF_idx));
end

% block reaction TFA bounds as well, re-run loop to print bounds to confirm
model.var_ub(ismember(model.varNames, strcat('F_', rxnsToKO))) = 0;
model.var_ub(ismember(model.varNames, strcat('R_', rxnsToKO))) = 0;
model.var_ub(ismember(model.varNames, strcat('NF_', rxnsToKO))) = 0;
model.var_lb(ismember(model.varNames, strcat('NF_', rxnsToKO))) = 0;

% now the model cannot grow as the reaction has been KOd
solKOTFA = optimizeThermoModel(model);

%% 
% * Does the model grow after KO? Is the observed behavior expected?
%% Merge model with ATLAS subset

% run PrepareForGapFilling to merge model and run essentiality
% alternatively ru-use GFmodel from previous exercise
flagEss = 0; % run essentiality prediction on the original and merged model
flagTFA = 1; % model is already TFA-ready, takes longer but cannot reproduce gapfilling otherwise
[GFmodel, conflict] = PrepareForGapFilling(model, {DBmodel},'c', flagEss,flagTFA,{},[],ReactionDB);

% save output to avoid having to re-run merging step
save('GFmodel.mat','GFmodel');

% load in future
load('GFmodel.mat');

% check if output model can grow
optimizeCbModel(GFmodel)

% TFA
optimizeThermoModel(GFmodel)
optimizeThermoModel(GFmodel,true)

%% 
% * Does the output model grow? What explains the magnitude of the computed 
% growth rate?
%% Generate gapfilling alternatives

%% gapfilling
GFmodeli = GFmodel;

% force growth
GFmodeli.var_lb(find(ismember(GFmodel.varNames,strcat('F_',model.rxns(find(model.c)))))) = 0.01;
GFmodeli.var_lb(find(ismember(GFmodel.varNames,strcat('NF_',model.rxns(find(model.c)))))) = 0.01;

%generate the alternative solutions
NumAlt = 30;
tagMin = 0; % solutions of min size
tagMinP1 = 0;  % solutions of min +1 size
tagSave = 1;
time = 600; % time limit for the solver in seconds
filename = './gf.mat';
indUSE = GFmodeli.indUSE;
mingrRate = 0.01; % force growth rate if you want to gap-fill for each bbb then mingrRate = 0
rxnsToCheck = {'AMAOTr'};
GFmodeli.rxnRescued = rxnsToCheck;
[ActRxns, GFSumRes, DPsAll] = gapFilling(GFmodeli,indUSE,NumAlt,mingrRate,rxnsToCheck,tagMin,time,filename);
save(strcat('./GF_',organism_id,'.mat'),'ActRxns','DPsAll')

%% Rank alternative gapfilling solutions

% save model backup before KOing (original Esseniality() run with equal output for all alternatives)
%model_backup = model;
%model = model_backup;
model = GFmodel;

% check growth of KO model
solWT = optimizeCbModel(model);

% check original bounds
rxnsToKO = {'AMAOTr'}; % to reproduce result from NICEgame paper fig 3a
f = find(ismember(model.rxns,rxnsToKO));
model.lb(f)
model.ub(f)
model.var_ub(ismember(model.varNames, strcat('F_', rxnsToKO)))
model.var_ub(ismember(model.varNames, strcat('R_', rxnsToKO)))
model.var_ub(ismember(model.varNames, strcat('NF_', rxnsToKO)))
model.var_lb(ismember(model.varNames, strcat('NF_', rxnsToKO)))

% block/unblock AMAOTr rxn. in principle we are trying to see if the model
% is now better at essentiality by fixing the false negative rxn AMAOTr;
% therefore the reaction should be opened for the evaluation stage
% FBA bounds
model.lb(f) = -50;
model.ub(f) = 50;
% TFA bounds
model.var_ub(ismember(model.varNames, strcat('F_', rxnsToKO))) = 100;
model.var_ub(ismember(model.varNames, strcat('R_', rxnsToKO))) = 100;
model.var_ub(ismember(model.varNames, strcat('NF_', rxnsToKO))) = 100;
model.var_lb(ismember(model.varNames, strcat('NF_', rxnsToKO))) = -100;

% check growth of KO model
solKO = optimizeCbModel(model);

%% block all the reactions from the database before running essentiality
model.genes = ttmodel.genes;
model.rxnGeneMat = ttmodel.rxnGeneMat;
model = addNetFluxVariables(model);
indNF = getAllVar(model,{'NF'});
fromDB=(find(ismember(model.rxnFrom,{'DBmodel_1'})));
modelk.var_lb(ismember(model.varNames, strcat( 'NF_',model.rxns(fromDB))))=0;
model.var_ub(ismember(model.varNames, strcat( 'NF_',model.rxns(fromDB))))=0;
model.var_ub(ismember(model.varNames, strcat( 'F_',model.rxns(fromDB))))=0;
model.var_ub(ismember(model.varNames, strcat( 'R_',model.rxns(fromDB))))=0;
model.ub(fromDB)=0;
model.lb(fromDB)=0;

% run essentiality evaluation, i cant actually get this code to run, see
% error message commented below
[rEssentiality]=Essentiality(model,rxnsToCheck,ActRxns,nonEssential,essential)

% No field cplex.Solution.x
% Warning: The solver does not return a solution! 
% > In solveTFBAmodelCplex>x_solveCplex (line 93)
% In solveTFBAmodelCplex (line 50)
% In optimizeThermoModel (line 27)
% In thermoSingleGeneDeletion (line 149)
% In Essentiality (line 37) 
% Index exceeds the number of array elements (0).
% 
% Error in thermoSingleGeneDeletion (line 150)
%     Jzrxns = model.rxns(solWTtfa.x(indNF)<1E-8);
% 
% Error in Essentiality (line 37)
%         [grRatio_genetfa,grRateKO_genetfa] = thermoSingleGeneDeletion(modelk, 'TFA', modelk.genes, 0, 0, 0, essThr);


% in S2 Table solutions for AMAOTr do all have very similar values
% Test biomass and BBB yield FBA
[resultFBA]=TestFBA_growth(model,ActRxns{1});
% Test biomass and BBB yield TFBA
[resultTFBA]=TestTFBA_growth(model,ActRxns{1});
% Length
[resultLength]=TestLength_growth(ActRxns{1});   

% At this point we have generated alternative gapfilling solutions for
% target reaction AMAOTr (FN in essentiality + rescued in merged model).
% Most alternatives have about the same essentiality score but I see
% differences in FBA solutions. For some reason TFA solution values are
% massive and there is very little difference between alternatives,
% probably something not properly constrained in the model.