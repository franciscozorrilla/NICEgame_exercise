%% Exercise 1
% Note: This tutorial is directly borrowed from the matTFA repo, with some minor 
% modifications (see <https://github.com/EPFL-LCSB/matTFA/blob/master/tutorials/tutorial.m 
% script>).
% 
% The goal is to learn the basics of running thermodynamic-based flux analysis 
% (TFA), including key functions used, as well as required input files & data.
% 
% 
%% Loading the input files
%% 
% * Small _E. coli_ model
% * Thermodynamic properties DB

% make sure solver is set
changeCobraSolver('cplex_direct')

%% load the model
load('smallEcoli.mat');
mymodel = smallEcoli;

% Limit the bounds of the fluxes that are higher than 100 or lower than
% -100 mmol/(gDW * h)
if any(mymodel.lb<-100) || any(mymodel.ub>100)
    mymodel.lb(mymodel.lb<-100) = -100;
    mymodel.ub(mymodel.ub>+100) = +100;
end

%% Load the thermodynamics database
tmp = load('../Exercise_2/thermo_data.mat');
ReactionDB = tmp.DB_AlbertyUpdate;
clear tmp

%% 
% 
%% FBA, FVA, and blocked reactions
%% 
% * Run FBA simulation with the |optimizeCbModel()| function
% * Run FVA with the |runMinaMax()| function
% * Remove blocked reactions with the |removeRxns()| function

% Test simple fba
solFBA = optimizeCbModel(mymodel);

% We can set a lower bound for growth (e.g. 50% of maximal growth)
min_obj = roundsd(0.5*solFBA.f, 2, 'floor');
mymodel.lb(mymodel.c==1) = min_obj;

%% Perform FVA 
fva = runMinMax(mymodel);

% Are there any blocked reactions?
% solver tolerance is 1e-9
SolTol = 1e-9;
id_Blocked_in_FBA = find( (fva(:,1)>-SolTol & fva(:,1)<SolTol) & ...
                          (fva(:,2)>-SolTol & fva(:,2)<SolTol) );

% If there exist block reactions
while ~isempty(id_Blocked_in_FBA)
    % remove them
    mymodel = removeRxns(mymodel, mymodel.rxns(id_Blocked_in_FBA));
    fva = runMinMax(mymodel);
    id_Blocked_in_FBA = find( (fva(:,1)>-SolTol & fva(:,1)<SolTol) & ...
                              (fva(:,2)>-SolTol & fva(:,2)<SolTol) );
end

%% 
%% TFA & TFA with concentration data
%% 
% * |prepModelforTFA()|: inputs model + thermo_data, calculates gibbs energies 
% of reactions, outputs prepared model
% * |convToTFA()|: inputs prepared model, adds new thermodynamic constraints 
% based on gibbs energies, outputs converted model
% * |addNetFluxVariables()|: inputs converted model, adds net flux variables 
% and constraints to model, outputs TFA-ready model
% * |solveTFAmodelCplex()|: inputs TFA-ready model, outputs TFA solution

% prepare model for TFA
prepped_m = prepModelforTFA(mymodel, ReactionDB, mymodel.CompartmentData);

%% Convert to TFA
tmp = convToTFA(prepped_m, ReactionDB, [], 'DGo', [], min_obj);

% Add net flux variables, which are equal to forwards flux - backwards flux
% NF_rxn = F_rxn - B_rxn
this_tmodel = addNetFluxVariables(tmp);

%% Solve tFA
soltFA = solveTFAmodelCplex(this_tmodel);

%% Perform TVA
% Get the variables representing the net fluxes
NF_ix = getAllVar(this_tmodel,{'NF'});
tva = runTMinMax(this_tmodel, this_tmodel.varNames(NF_ix));

%% We add some generic data for cofactor concentrations
metNames = {'adp_c', 'amp_c', 'atp_c'};
C_lb = [1e-06, 2e-04, 1e-03]';
C_ub = [7e-04, 3e-04, 5e-02]';
LC_varNames = {'LC_adp_c', 'LC_amp_c', 'LC_atp_c'};
% find the indices of these variables in the variable names of the tfa
id_LC_varNames = find_cell(LC_varNames, this_tmodel.varNames);
% Set to the model these log-concentration values
this_tmodel.var_lb(id_LC_varNames) = log(C_lb);
this_tmodel.var_ub(id_LC_varNames) = log(C_ub);
% Run another tva with the data
tva_wData = runTMinMax(this_tmodel, this_tmodel.varNames(NF_ix));

%% 
% 
%% Visualize bidirectional reactions that become unidirectional with thermodynamic constraints
% 


%% Plot the differences
% Inline function definitions to get:
% - values of flux ranges from a two column vector of lower and upper bounds
f = @(x) abs(x(:,2) - x(:,1));
% - logical vector with indices for bidirectional reactions (flux ranges crossing zero)
n = @(x) x(:,1)<-1e-9 & x(:,2)>1e-9;
% - scoring metric of relative changes in ranges, for ranges greater than 0.01
s = @(x,y) ( abs(f(x)-f(y))./f(x) ) .* (f(x) > 10);

% Calculate this scoring metric for the differences between
% - tfa and tva
loss = s(fva, tva);
% - tva witout and with concentration data
loss_tva_wData = s(tva, tva_wData);

% find bidirectional reactions based on
% - (1) fva
is_bd_fva = (n(fva));
% - (2) tva without concentration data
is_bd_tva = (n(tva));
% - (3) tva with concentration data
is_bd_tvawData = (n(tva_wData));

% Reactions that become unidirectional from fva (1) To tva without
% concentration data (2)
criterion1 = setdiff(find(is_bd_fva), find(is_bd_tva));
% Reactions that become unidirectional tva without (2) to tva with
% concentration data (3)
criterion2 = setdiff(find(is_bd_tva), find(is_bd_tvawData));

% We plot the ranges of these reactions
id_rxns_to_plot = [criterion1; criterion2];
xlabels = this_tmodel.rxns(id_rxns_to_plot);
figure
hold on
p1 = p_patch_va(fva(id_rxns_to_plot,:),0,[34 29 35]./255, 0.2);
p2 = p_patch_va(tva(id_rxns_to_plot,:),0.20,[61 90 128]./255, 0.2);
p3 = p_patch_va(tva_wData(id_rxns_to_plot,:),0.4,[152 193 217]./255, 0.2);
plot([0 size(id_rxns_to_plot,1)+1],[0 0],'k--','linewidth',2)
set(gca,'XTick',1:size(id_rxns_to_plot,1),'XTickLabel',xlabels,'XTickLabelRotation',-45,'TickLabelInterpreter','none')
legend([p1 p2 p3] , 'FBA', 'TFA & default conc. ranges','TFA & conc. data','Location','southwest')
ylabel ('flux [mmol.gDw^{-1}]')
title('Bi-directional reactions become unidirectional upon imposing thermodynamic constraints and data')%,'fontsize',20)

%% 
% 
%% Reactions with flux ranges affected by thermondynamic constraints
% 

% To illustrate the impact of thermodynamics on the network we select to
% plot also the ranges of reactions that are not necessarily bidirectional
% but are affected by the addition of thermodynamics
% (1) the bi-directionals scoring higher than 0.2 (>20% difference) in the
% relative difference metric between tfa and tva
criterion3 = find(loss>0.2);
% (2) the reactions scoring higher tnan 0.2 (>20% difference) in the
% relative difference metric between tva witout and tva with concentration
% data
criterion4 = find(loss_tva_wData>0.2);

id_rxns_to_plot = union(criterion3, criterion4);
xlabels = this_tmodel.rxns(id_rxns_to_plot);
figure
hold on
p1 = p_patch_va(fva(id_rxns_to_plot,:),0,[34 29 35]./255, 0.2);
p2 = p_patch_va(tva(id_rxns_to_plot,:),0.20,[61 90 128]./255, 0.2);
p3 = p_patch_va(tva_wData(id_rxns_to_plot,:),0.4,[152 193 217]./255, 0.2);
plot([0 size(id_rxns_to_plot,1)+1],[0 0],'k--','linewidth',2)
set(gca,'XTick',1:size(id_rxns_to_plot,1),'XTickLabel',xlabels,'XTickLabelRotation',-45,'TickLabelInterpreter','none')
legend([p1 p2 p3] , 'FBA', 'TFA & default conc. ranges','TFA & conc. data','Location','southwest')
ylabel ('flux [mmol.gDw^{-1}]')
title('Impact of thermodynamic constraints and data on network flux ranges','fontsize',20)
%% 
%