%% Exercise #2
% Note: This tutorial is directly borrowed from the NICEgame repo, with some 
% minor modifications (see <https://github.com/EPFL-LCSB/NICEgame/blob/master/matlab/scripts/tutorial.m 
% script>).
% 
% The goal is to learn the basics of running NICEgame, including key functions 
% used, as well as required input files.
% 
% 
%% Loading the input files
%% 
% * BiGG reaction database (King et al., Nucleic Acids Research, 2016)
% * _E. coli_ GEM iML1515 (Monk et al.,Nat Biotechnol, 2017)
% * Thermodynamic properties DB (Alberty et al., John Wiley & Sons, 2006)

% load BiGG reaction database
load('BiGG_DB.mat')
% create a copy of BiGG db
DBmodel = model;

% load GEM
load('iML1515.mat')
% create a copy of the original model
sourceModel = model;
% define organism id
organism_id = 'iML1515';

% load thermo database
load('thermo_data.mat')
thermo_data = DB_AlbertyUpdate; 
% thermo_data = {} % if no thermo
%% 
% As an exercise, inspect the strucutre of the GEM, BiGG database, and thermo 
% database (e.g. run |open model| in the command window). 
%% 
% * What are the dimensions of the S matrix in the GEM and BiGG database? What 
% do these dimensions represent? 
% * What kind of data is contained in the thermo database?
%% 
% 
% Knockout reaction for gapfilling
%% 
% * Simulate wild-type model growth using the |optimizeCbModel()| function
% * Knockout reaction by closing bounds
% * Simulate knockout model growth

% make sure CPLEX solver is set before optimizing
changeCobraSolver('cplex_direct','LP');

% check growth of wild-type model, i.e. before knocking out any reaction
solWT = optimizeCbModel(model);

% define rxn ID to KO
rxnsToKO = {'5DOAN'};

% find index for KO rxn
f = find(ismember(model.rxns,rxnsToKO));

% set upper and lower bounds to 0, effectively blocking flux through rxn
model.lb(f) = 0;
model.ub(f) = 0;

% check growth of KO model
solKO = optimizeCbModel(model);
%% 
% As an exercise, inspect the |optimizeCbModel()| function and output (e.g. 
% |open optimizeCbModel.m|)
%% 
% * What is the function doing with the model? What is |S|, |b|, |c|, |ub|, 
% and |lb|?
% * What do the different fields in the output structure represent? What is 
% |x|, |f|, and |stat|?
% * What effect did knocking out a reaction have in the model?
%% 
% 
% Check biomass building blocks that cannot be produced
%% 
% * Find biomass reaction index through the objective function
% * Find biomass building blocks (bbb) and add a demand reaction for each
% * Check which bbb's can be produced and which cannot

% find biomass reaction index
obj = find(model.c);

% find biomass building blocks, i.e. substrates in the biomass reaction
bbbs_m = model.mets(find(model.S(:,find(model.c))<0));

% add demand reactions
model_bbb = addDemandReaction(model,bbbs_m);

% get size of model reactions
r = length(model.rxns);

% initialize loop through each biomass building block
for i = 1:length(bbbs_m)
    % create new copy of model upon starting loop
    modeli = model_bbb;
    % clear objective function
    modeli.c = zeros(length(modeli.rxns),1);
    % set objective function to one of newly added demand reactions
    modeli.c(r+i) = 1;
    % solve optimization problem to check if model can produce each bbb
    sol = optimizeCbModel(modeli,'max');
    % set production to 0 if value is very small
    if sol.f < 10^-4 % units of mmol/gDCW*hr
        sol.f = 0;
    end
    % save how much flux of each bbb can be produced
    sol_bbb(i,1) = sol.f;
end

% identify problematic bbb's as those that cannot be produced after
% knocking out rxn 5DOAN
problematic_BBBs = bbbs_m(sol_bbb == 0); %

% store reaction IDs
rxnID = strcat('DM_',problematic_BBBs);
%% 
% Note: we can gapfill in a more targeted way by filling in gaps for each problematic 
% biomass building block, or we can do a more general gapfilling by directly targeting 
% biomass.
%% 
% * What does the |addDemandReaction()| function do? What is a demand reaction?
% * How many biomass building blocks does the model have?
% * Which biomass building blocks cannot be produced?
% 
% Prepare for gapfilling
%% 
% * Set parameters for the |PrepareForGapFilling()| function
% * Merge model and reaction database

% perform an essentiality analysis for the WT and the merged model, skip
% for now, we will carry out essentiality predictions in exercise 3
tagEssentiality = 0;

% perform conversion of model for TFA
flagTFA = 1;

% merge model with reaction db, use also the thermo db to create TFA ready
% model.
[GFmodel, conflict] = PrepareForGapFilling(sourceModel, {DBmodel},'', tagEssentiality,flagTFA,{},[],thermo_data);
%% 
% Inspect the |PrepareForGapFilling()| function and output |GFmodel|.
%% 
% * What are the inputs and what does the function do?
% * What are the dimension of the S matrix in the merged model?
% 
% Check media composition
%% 
% * Print exchange reactions with allowed uptake and secretion
% * Optionally modify media composition

% In this model with BiGG reaction IDs we can identify exchanges with the
% prefix 'EX_'. The following code will print all exchange reactions with
% an open uptake or secretion

exRxnsIdx = find(startsWith(GFmodel.rxns, 'EX_'));
fprintf('Exchange Reaction     | Uptake | Secretion |    LB   |    UB   | Metabolite Name\n');
fprintf('----------------------------------------------------------------------------------------\n');

uptakeCount = 0;
secretionCount = 0;

for i = 1:length(exRxnsIdx)
    idx = exRxnsIdx(i);
    rxnID = GFmodel.rxns{idx};
    lb = GFmodel.lb(idx);
    ub = GFmodel.ub(idx);

    uptake = lb < 0;
    secretion = ub > 0;

    if ~uptake && ~secretion
        continue;
    end

    if uptake
        uptakeCount = uptakeCount + 1;
        uptakeStr = 'YES';
    else
        uptakeStr = 'NO';
    end

    if secretion
        secretionCount = secretionCount + 1;
        secretionStr = 'YES';
    else
        secretionStr = 'NO';
    end

    if startsWith(rxnID, 'EX_')
        metID = extractAfter(rxnID, 'EX_');
    else
        metID = '[Unknown]';
    end

    % Use metabolite names instead of IDs
    metIdx = find(strcmp(GFmodel.mets, metID), 1);

    if isempty(metIdx)
        metNameStr = '[Unknown]';
    else
        metName = GFmodel.metNames{metIdx};
        if iscell(metName)
            metNameStr = metName{1};
        else
            metNameStr = metName;
        end
    end

    fprintf('%-22s |  %5s  |   %5s   | %7.2f | %7.2f | %s\n', ...
        rxnID, uptakeStr, secretionStr, lb, ub, metNameStr);
end

fprintf('----------------------------------------------------------------------------------------\n');
fprintf('Total exchange reactions open for uptake:   %d\n', uptakeCount);
fprintf('Total exchange reactions open for secretion: %d\n', secretionCount);

% it is not required for this exercise, but media composition can be
% changed by modifying list of metabolites

% first find all exchange reactions
EX = find(contains(GFmodel.rxns,'EX_'));

% close the uptake of all metabolites
GFmodel.lb(EX) = 0;

% verify that the model should not be able to grow without uptaking
optimizeCbModel(GFmodel)

% define list of metabolites desired to be added for uptake
media = {'EX_pi_e','EX_h_e','EX_fe3_e','EX_mn2_e','EX_co2_e','EX_fe2_e', ...
    'EX_glc__D_e','EX_zn2_e','EX_cbl1_e','EX_mg2_e','EX_ca2_e','EX_ni2_e', ...
    'EX_cu2_e','EX_cobalt2_e','EX_sel_e','EX_h2o_e','EX_nh4_e','EX_mobd_e', ...
    'EX_so4_e','EX_k_e','EX_na1_e','EX_o2_e','EX_cl_e','EX_tungs_e','EX_slnt_e'};

% get reaction indexes for metabolites
f = find(ismember(GFmodel.rxns,media));

% re-open closed reactions
GFmodel.lb(f) = -50;

% verify that the model can now grow again
optimizeCbModel(GFmodel)
%% 
% It is important to check the bounds on exchange reactions when running simulations 
% and interpreting results.
% 
% 
% Generate gapfilling alternatives
%% 
% * Set parameters and run |gapFilling()|

% gapfill for biomass
NumAlt = 10; % number of alternatives to generate
tagMin = 1; % solutions of min size
tagMinP1 = 1;  % solutions of min +1 size
tagSave = 1;
time = 600; % time limit for the solver in seconds
filename = './gf.mat';
indUSE = GFmodeli.indUSE;
mingrRate = 0.01; % force growth rate if you want to gapfill for biomass
rxnsToCheck = {};
[ActRxns, GFSumRes, DPsAll] = gapFilling(GFmodeli,indUSE,NumAlt,mingrRate,rxnsToCheck,tagMin,time,filename);
save(strcat('./GF_',organism_id,'.mat'),'ActRxns','DPsAll')

%% 
% Inspect the |gapfilling()| function's code and inspect the input and output 
% structures.
%% 
% * Inspect the generated gapfilling solutions stored in |ActRxns|