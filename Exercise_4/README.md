### 📑 Exercise 4: Reproduction of result from NICEgame paper (AMAOTr)
- Goal
  - We previosuly identified AMAOTr as target for gapfilling. Now we knockout this reaction in iML1515, merge with the ATLAS database, generate alternative solutions, and evaluate them (e.g. confusion matrix components, growth rate)
  - We obtain slightly different values (difference likely due to gene ID matching, publication uses 1470 while we use 1424), but overall very similar MCC scores across alternatives (~0.486) and compared to publication results (~0.487).
- Files
  - `reducedATLAS_ecoli_yeast.mat`: ATLAS subset of reactions involving metabolites in *E. coli* and *S. cerevisiae*
- Functions
  - `PrepareForGapFilling()`: inputs GEM and reaction database, outputs merged model
  - `gapFilling()`: inputs merged model, outputs alternative solutions for gapfilling 
  - `Essentiality()`:  inputs alternative gapfilling solutions, model, and essentiality experimental data, outputs confusion matrix components (FP, FN, TP, TN) and Matthews correlation coefficient (MCC) for each alternative
  - `TestFBA_growth()`: inputs merged model and alternative gapfilling solutions, outputs FBA-based growth rates for each alternative
  - `TestTFBA_growth()`: inputs merged model and alternative gapfilling solutions, outputs TFA-based growth rates for each alternative

[Next exercise](https://github.com/franciscozorrilla/NICEgame_exercise/tree/master/Exercise_5)
