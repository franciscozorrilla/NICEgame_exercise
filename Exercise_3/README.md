### 💻 Exercise 3: Essentiality prediction & evaluation
- Goal
  - Run single reaction essentiality (FBA & TFA-based) and single gene essentiality (FBA & TFA-based). Also run TFA-based single reaction essentiality on merged model to identify rescued reactions. These are used to identify targets for gapfilling.
  - Obtain experimental data and compare against predictions. False negatives are experimentally dispensible genes predicted to be lethal by gene knockout; these are targets for gapfilling, as this indicates that there must be some other reaction(s) compensating for the knockout. Also compare merged model essentiality to identify rescued reactions (e.g. AMAOTr)
- Files
  - `iML1515_KEGG_final_thermo.mat`: TFA-ready *E. coli* model with mapped metabolite IDs to ATLAS
  - `experimental_essentiality.tsv`: gene knockout experiment results from [source publication](https://journals.asm.org/doi/10.1128/mbio.02096-17).
  - `gene_table.tsv`: lookup table to convert from publication gene IDs (e.g. glyA) to model gene IDs (e.g. b2551)
- Functions
  - `singleRxnDeletion()`: inputs model, runs FBA-based single reaction knockouts, returns lethal reaction knockouts
  - `thermoSingleRxnDeletion()`: inputs TFA model, runs TFA-based single reaction knockouts, returns lethal reaction knockouts
  - `singleGeneDeletion()`: inputs model, runs FBA-based single gene knockout, returns lethal gene knckouts
  - `thermoSingleGeneDeletion()`: inputs model, runs TFA-based single gene knockout, returns lethal gene knockouts
  - `strcmp()`: used to compare values from predictions and experimental data
  - `readtable()`: used to load data files into MATLAB

### 🐛 Bug report
- Issue
  - Merged model triggers error when running `thermoSingleRxnDeletion()`, possibly related to overly constrained thermo parameters?
   ```MATLAB
   % Run single command rxn essentiality - TFA based
  [grRatioM, grRateKOM, grRateWTM, hasEffectM, delRxnsM, fluxSolutionM, impactTasksM] = thermoSingleRxnDeletion(GFmodel);
  
  No field cplex.Solution.x
  Warning: The solver does not return a solution! 
  > In solveTFBAmodelCplex>x_solveCplex (line 93)
  In solveTFBAmodelCplex (line 50)
  In optimizeThermoModel (line 27)
  In thermoSingleRxnDeletion (line 81) 
  Index exceeds the number of array elements (0).
  
  Error in thermoSingleRxnDeletion (line 82)
  Jzrxns = model.rxns(solWTtfa.x(indNF)<1E-8);
   ```
- Workaround
  - Turn on essentiality and thermo flags when running `PrepareForGapFilling`, this runs TFA-based reaction essentiality in original and merged models, based on which it identifies rescued reactions.

[Next exercise](https://github.com/franciscozorrilla/NICEgame_exercise/tree/master/Exercise_4)
