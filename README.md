# 🏓 NICEgame exercise 🏓

This repo contains code and files used during a 2 week visit to the LCSB lab. The goal was to get hands on experience using the NICEgame workflow to curate metabolic models based on experimental data. Through these exercises we also learn about thermodynamic-based flux analysis (TFA), BridgIT, and ATLASx.

### 🌡️ Exercise 1: matTFA tutorial
- Goal
  - Learn the basics of TFA: prepare GEM for TFA, compare TFA and FBA results to show that TFA is more restrictive
- Files
    - thermo_data.mat: database containing chemical information for ~18k compounds (e.g. mass, charge, pKa, deltaG_formation)
- Functions
    - prepModelforTFA(): inputs model + thermo_data, calculates gibbs energies of reactions, outputs prepared model
    - convToTFA(): inputs prepared model, adds new thermodynamic constraints based on gibbs energies, outputs converted model
    - addNetFluxVariables(): inputs converted model, adds net flux variables and constraints to model, outputs TFA-ready model
    - solveTFAmodelCplex(): inputs TFA-ready model, outputs TFA solution

### 🐤 Exercise 2: NICEgame tutorial
- Goal
  - Learn the basics of NICEgame: KO reaction in iML1515 model, merge model with BiGG database, generate alternative gapfilling solutions
- Files
  - iML1515.mat: *E. coli* model
  - BiGG_DB.mat: BiGG database of reactions
- Functions
  - PrepareForGapFilling(): inputs GEM and reaction database, outputs merged model
  - gapFilling(): inputs merged model, outputs alternative solutions for gapfilling

### 💻 Exercise 3: Essentiality prediction
- Goal
  - Single reaction essentiality (FBA & TFA), single gene essentiality (FBA & TFA). These are used to identify targets for gapfilling
- Files
- Functions 

### ⚖️ Exercise 4: Essentiality evaluation
- Goal
  - Obtain experimental data and compare against predictions (FN = targets for gapfilling, e.g. AMAOTr). Also compare merged model essentiality to identify rescued reactions (e.g. AMAOTr)
- Files
- Functions
  
### 📑 Exercise 5: Results reproduction in NICEgame
- Goal
  - Identified AMAOTr as target for gapfilling, KO in iML1515, merge with ATLAS database, generate alternative solutions and evaluate them (e.g. confusion matrix, growth rate)
- Files
- Functions
  
### 🧫 Exercise 6: NICEgame with BIOLOG data
- Goal
  - Obtained model + processed biolog growth data on different carbon sources, evaluate metabolite essentiality, choose false negative target for gapfilling (D-Galacturonic-Acid), replace as carbon source, suggest gapfilling alternatives. For one alternative, add the reaction to the model and re-evaluate TP/FP/TN/FN, improve MCC score. Repeat for FP gapfilling target (Arginine).
- Files
- Functions
  
### 🌉 Exercise 7: Submitting a reaction query to the BridgIT server
- Goal
  - Identify E.C. number for KEGG reaction R00919
- Files
- Functions
  
### 🗺️ Exercise 8: Submitting a metabolite query to the ATLASx server
- Goal
  - Identify reactions to degrade tween-80
- Files
- Functions 

## Requirements

- MATLAB: v2020b recommended for apple silicon
- CPLEX: v12.9 (full version free for academic use)

## Setup

After you have installed MATLAB and CPLEX, make sure that you add CPLEX to the MATLAB search path (accessed via the `Set path` button). Next, follow the installation instruction from the respective repos to install the following tools:

- [NICEgame](https://github.com/EPFL-LCSB/NICEgame)
- [matTFA](https://github.com/EPFL-LCSB/matTFA)

After cloning and pulling the repos, make sure to also add them in the MATLAB search path.

## License
The software in this repository is put under an APACHE licensing scheme

## Reference
If you use NICEgame (or any code in this repo) for your research, please consider citing the original publciation:

> E. Vayena, A. Chiappino-Pepe, H. MohammadiPeyhani, Y. Francioli, N. Hadadi, M. Ataman, J. Hafner, S. Pavlou, & V. Hatzimanikatis, A workflow for annotating the knowledge gaps in metabolic reconstructions using known and hypothetical reactions, Proc. Natl. Acad. Sci. U.S.A. 119 (46) e2211197119, https://doi.org/10.1073/pnas.2211197119 (2022).
