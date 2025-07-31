### 🌡️ Exercise 1: matTFA tutorial
- Goal
  - Learn the basics of running TFA: prepare GEM for TFA, compare TFA and FBA results to show that TFA is more restrictive
- Files
  - `thermo_data.mat`: database containing chemical information for ~18k compounds (e.g. mass, charge, pKa, deltaG_formation)
- Functions
  - `prepModelforTFA()`: inputs model + thermo_data, calculates gibbs energies of reactions, outputs prepared model
  - `convToTFA()`: inputs prepared model, adds new thermodynamic constraints based on gibbs energies, outputs converted model
  - `addNetFluxVariables()`: inputs converted model, adds net flux variables and constraints to model, outputs TFA-ready model
  - `solveTFAmodelCplex()`: inputs TFA-ready model, outputs TFA solution
  - `optimizeCbModel()`: inputs FBA model, outputs FBA solution

[Next exercise](https://github.com/franciscozorrilla/NICEgame_exercise/tree/master/Exercise_2)
