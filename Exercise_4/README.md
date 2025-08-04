### 📑 Exercise 4: Reproduction of result from NICEgame paper (AMAOTr)
- Goal
  - We previosuly identified AMAOTr as target for gapfilling. Now we knockout this reaction in iML1515, merge with the ATLAS database, generate alternative solutions, and evaluate them (e.g. confusion matrix components, growth rate)
  - We obtain slightly different values (difference likely due to gene ID matching, publication uses total of 1470 while we use 1424), but overall very similar MCC scores across alternatives (~0.486) and compared to publication results (~0.487).
- Files
  - `reducedATLAS_ecoli_yeast.mat`: ATLAS subset of reactions involving metabolites in *E. coli* and *S. cerevisiae*
- Functions
  - `PrepareForGapFilling()`: inputs GEM and reaction database, outputs merged model
  - `gapFilling()`: inputs merged model, outputs alternative solutions for gapfilling 
  - `Essentiality()`:  inputs alternative gapfilling solutions, model, and essentiality experimental data, outputs confusion matrix components (FP, FN, TP, TN) and Matthews correlation coefficient (MCC) for each alternative
  - `TestFBA_growth()`: inputs merged model and alternative gapfilling solutions, outputs FBA-based growth rates for each alternative
  - `TestTFBA_growth()`: inputs merged model and alternative gapfilling solutions, outputs TFA-based growth rates for each alternative

[Next exercise](https://github.com/franciscozorrilla/NICEgame_exercise/tree/master/Exercise_5)

### 🐛 Bug report
- Issue
  - Problem when trying to evaluate essentiality predictions for alternative gapfilling solutions, similar to bug in notebook 3 looks like the error is possibly related to thermo constraints?
   ```MATLAB
    > [rEssentiality]=Essentiality(model,rxnsToCheck,ActRxns,nonEssential,essential)
    
     No field cplex.Solution.x
     Warning: The solver does not return a solution! 
     > In solveTFBAmodelCplex>x_solveCplex (line 93)
     In solveTFBAmodelCplex (line 50)
     In optimizeThermoModel (line 27)
     In thermoSingleGeneDeletion (line 149)
     In Essentiality (line 37) 
     Index exceeds the number of array elements (0).
     
     Error in thermoSingleGeneDeletion (line 150)
         Jzrxns = model.rxns(solWTtfa.x(indNF)<1E-8);
     
     Error in Essentiality (line 37)
             [grRatio_genetfa,grRateKO_genetfa] = thermoSingleGeneDeletion(modelk, 'TFA', modelk.genes, 0, 0, 0, essThr);
   ```
