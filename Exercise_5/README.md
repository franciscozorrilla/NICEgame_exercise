### 🧫 Exercise 5: NICEgame with biolog data
- Goal
  - Obtain model with associated biolog growth data on different carbon sources, evaluate metabolite essentiality by changing the carbon source in base media and comparing to experimental results.
  - Fixing false negatives: FN compounds are experimentally viable carbon sources for growth, yet the model is unable to produce biomass, indicating they are targets for gapfilling. Choose false negative (e.g. D-Galacturonic-Acid), replace as carbon source, and suggest gapfilling alternatives. For one alternative, add the reaction to the model and re-evaluate confusion matrix components showing improvement in MCC score.
  - Fixing false positives: FP compounds are not experimentally viable carbon sources for growth, yet the model still produces biomass, indicating that some reaction(s) should be turned off. Fix a false positive target (e.g. Arginine); the model can grow using Arginine as a carbon source but experimental data shows otherwise. In this case we have to identify candidate reactions to knockout that resolve false positive growth without introducing new errors. Find lethal reaction knockouts with Arginine as a carbon source, evaluate confusion matrix components, and identify knockouts that maximise MCC.
- Files
  - `keggModel.mat`: KEGG reaction database 
- Functions
  - `PrepareForGapFilling()`: inputs GEM and reaction database, outputs merged model
  - `gfbiomass()`: inputs merged model, outputs alternative solutions for gapfilling
  - `addReaction()`: inputs model and desired reaction to be added, outputs expanded model

[Next exercise](https://github.com/franciscozorrilla/NICEgame_exercise/tree/master/Exercise_6)

### ⚠️ Note
- the model file and experimental data used in this exercise will be released once the associated paper is submitted for publication
