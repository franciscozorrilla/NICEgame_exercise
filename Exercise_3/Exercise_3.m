%% 💻 Exercise #3: Essentiality prediction & evaluation
% The goal of this exercise is to generate gene and reaction essentiality predictions, 
% using FBA and TFA-based approaches.
% 
% Finally, we also predict reaction essentiality on the merged _E. coli_ and 
% ATLAS subset model to identify rescued reactions.
%% Load model and ATLAS subset reaction database

% iML1515_KEGG_final_thermo, mapped manually by Eva
load('iML1515_KEGG_final_thermo.mat') 
model = ttmodel;

% define organism_id for filenaming later on
organism_id = 'iML1515';

% load ATLAS subset
load('reducedATLAS_ecoli_yeast.mat') % ATLAS ecoli yeast metabolite subset
DBmodel = reducedATLAS;

% load also thermo database
tmp = load('../Exercise_2/thermo_data.mat');
ReactionDB = tmp.DB_AlbertyUpdate;
clear tmp
%% 
% * Is the loaded model TFA ready?
% * How large is the reaction database?
%% FBA-based single reaction deletion

% make sure solver is set
changeCobraSolver('cplex_direct');

% Run in single command
[grRatioFBA,grRateKOFBA,grRateWTFBA,hasEffectFBA,delRxnsFBA,fluxSolutionFBA] = singleRxnDeletion(model);
% Lets count the number of KOs where growth is less than 10% of original
Rxns = model.rxns;
TFBA_rxns = table(Rxns, grRatioFBA, 'VariableNames', {'Genes', 'grRatio'});
% 355 lethal rxn KOs with FBA
sum(TFBA_rxns.grRatio < 0.1)
%% 
% * How many lethal reactions are predicted by FBA?
%% TFA-based single reaction deletion

% Single command rxn TFA essentiality
[grRatioR, grRateKOR, grRateWTR, hasEffectR, delRxnsR, fluxSolutionR, impactTasksR] = thermoSingleRxnDeletion(model);
Rxns = model.rxns;
TTFA_rxns = table(Rxns, grRatioR, 'VariableNames', {'Reaction', 'grRatio_original'});
% 324 lethal rxn KOs with TFA
sum(TTFA_rxns.grRatio_original < 0.1)
%% 
% * How many lethal reactions are predicted by TFA?
%% FBA-based single gene deletion

% Run single command gene essentiality - FBA based
[grRatioFBA,grRateKOFBA,grRateWTFBA,hasEffectFBA,delRxnsFBA,fluxSolutionFBA] = singleGeneDeletion(model);
Genes = model.genes;
TFBA_genes = table(Genes, grRatioFBA, 'VariableNames', {'Genes', 'grRatio'});
% 258 lethal gene KOs with FBA
sum(TFBA_genes.grRatio < 0.1)
%% 
% * How many lethal genes are predicted by FBA?
%% TFA-based single gene deletion

% Run single command gene essentiality - TFA based
[grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution, impactTasks] = thermoSingleGeneDeletion(model);
TTFA_genes = table(Genes, grRatio, 'VariableNames', {'Genes', 'grRatio'});
% 252 lethal gene KOs with TFA
sum(TTFA_genes.grRatio < 0.1)
%% 
% * How many lethal genes are predicted by TFA?
%% Merge model with ATLAS subset

% run PrepareForGapFilling to merge model and run essentiality
flagEss = 1; % run essentiality prediction on the original and merged model
flagTFA = 1; % ideally this flag should be set to 1 so that TFA-based essentiality is computed, however this takes a few hours to run
[GFmodel, conflict] = PrepareForGapFilling(model, {DBmodel},'c', flagEss,flagTFA,{},[],ReactionDB);
% save to avoid having to re-run
save('GFmodel.mat','GFmodel')
%% 
% * What are the dimension of the merged model S matrix?
%% Merged model TFA-based single reaction deletion

% Run single command rxn essentiality - TFA based
[grRatioM, grRateKOM, grRateWTM, hasEffectM, delRxnsM, fluxSolutionM, impactTasksM] = thermoSingleRxnDeletion(GFmodel);

% when trying to reproduce this result I now get the following error:
% No field cplex.Solution.x
% Warning: The solver does not return a solution! 
% > In solveTFBAmodelCplex>x_solveCplex (line 93)
% In solveTFBAmodelCplex (line 50)
% In optimizeThermoModel (line 27)
% In thermoSingleRxnDeletion (line 81) 
% Index exceeds the number of array elements (0).
% 
% Error in thermoSingleRxnDeletion (line 82)
% Jzrxns = model.rxns(solWTtfa.x(indNF)<1E-8);

Rxns_merged = GFmodel.rxns;
TTFA_rxnsM = table(Rxns_merged, grRatioM, 'VariableNames', {'Reaction', 'grRatio'});
%  lethal reaction KOs with TFA
sum(TTFA_rxnsM.grRatio < 0.1) % 119
%% 
% * How many lethal reactions are predicted by TFA in the merged model?
%% Identify rescued reactions as gapfilling targets

% TFA-based
% Compute rescued reactions by combining tables and filtering
merged_essentiality = innerjoin(TTFA_rxnsM,TTFA_rxns, 'Keys','Reaction');
rescued_reactions = merged_essentiality(merged_essentiality.grRatio_original<0.1 & merged_essentiality.grRatio>0.1,:) % 205 rescued rxns

% in any case, in the previous code chunk we ran PrepareForGapFilling() with TFA and essentiality flags
% so we can look at GFmodel.rxnRescued structure and see that AMAOTr is present

% if PrepareForGapFilling() was run with essentiality flag
% set to 1, the model will have structure rxnRescued = rxnEssentialWT - rxnEssentialMerged
% with the thermo flag set to 0 these predictions are FBA based, with
% thermo flag set to 1 these predictions are TFA based
height(GFmodel.rxnRescued); % 223 rescued reactions
sum(strcmp(GFmodel.rxnRescued,'AMAOTr')); % AMAOTr identified as gapfilling target

%% 
% 
%% Compare essentiality predictions against experimental data
% In the NICEgame paper, they use essentiality data from <https://journals.asm.org/doi/10.1128/mbio.02096-17 
% this publication>.

% load table to translate gene IDs to gene names from essentiality data
gene_table = readtable('gene_table.tsv', 'FileType', 'text', 'Delimiter', '\t');

% load table with experiemntal essentiality (table S1,
% https://doi.org/10.1128/mbio.02096-17)
essentiality_table = readtable('experimental_essentiality.tsv', 'FileType', 'text', 'Delimiter', '\t');
% sanity check
height(essentiality_table(strcmp(essentiality_table.Essential,'TRUE'),:)); % 358 essential genes
height(essentiality_table(strcmp(essentiality_table.Essential,'FALSE'),:)); % 3955 non-essential

% filter out genes with ambiguous essentiality
essentiality_table = essentiality_table(strcmp(essentiality_table.Unclear, 'FALSE'), :);
% sanity check
height(essentiality_table(strcmp(essentiality_table.Essential,'TRUE'),:)); % 358 essential genes
height(essentiality_table(strcmp(essentiality_table.Essential,'FALSE'),:)); % 3793 non-essential

% merge tables
essentiality_table = innerjoin(essentiality_table, gene_table, 'Keys', 'Gene');
% sanity check
height(essentiality_table(strcmp(essentiality_table.Essential,'TRUE'),:)); % 146 essential genes
height(essentiality_table(strcmp(essentiality_table.Essential,'FALSE'),:)); % 1278 non-essential

% pull out essential & nonessential genes cell strucutre
essential = essentiality_table(strcmp(essentiality_table.Essential,"TRUE"),:);essential=essential.Genes;
nonEssential = essentiality_table(strcmp(essentiality_table.Essential,"FALSE"),:);nonEssential=nonEssential.Genes;
%% Identify false negatives as targets for gapfilling
% These are genes that are predicted to be lethal when KO'd, yet experiments 
% show otherwise. This means that there must be missing pathways that can compensate 
% for gene KO. To gap-fill for these, we KO the reaction and then generate alternatives 
% that can restore growth.

% merge TFA-based gene essentiality predictions with essentiality experiment data
essentiality_eval = innerjoin(TTFA_genes, essentiality_table, 'Keys', 'Genes');

% calculate FN, TN, FP, TP
FN = essentiality_eval(essentiality_eval.grRatio < 0.1 & strcmp(essentiality_eval.Essential, 'FALSE'), :); % 136 FN
TN = essentiality_eval(essentiality_eval.grRatio < 0.1 & strcmp(essentiality_eval.Essential, 'TRUE'), :); % 103 TN
FP = essentiality_eval(essentiality_eval.grRatio > 0.1 & strcmp(essentiality_eval.Essential, 'TRUE'), :); % 43 FP
TP = essentiality_eval(essentiality_eval.grRatio > 0.1 & strcmp(essentiality_eval.Essential, 'FALSE'), :); % 1142 TP

% I can see that the gene bioA (b0774) is found in the FN structure, meaning that
% it is experimentally non-essential but computationally essential.
% Furthermore, we find associated reaction AMAOTr in the set of rescued
% reactions, meaning that it is computationally essential in the original
% model but no longer in the merged model. Lets generate some alternative
% solutions for gapfilling next. This is a gapfilling target in figure 3 of
% the NICEgame paper