# 🏓 NICEgame exercise 🏓

This repo contains code and files used during a 2 week visit to the LCSB lab. The goal was to get hands on experience using the [NICEgame](https://www.pnas.org/doi/10.1073/pnas.2211197119) workflow to curate metabolic models based on experimental data. Through these exercises we also learn about thermodynamic-based flux analysis (TFA), BridgIT, and ATLASx.

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

### 🐤 Exercise 2: NICEgame tutorial
- Goal
  - Learn the basics of NICEgame: KO reaction in iML1515 model, merge model with BiGG database, generate alternative gapfilling solutions
- Files
  - `iML1515.mat`: *E. coli* model
  - `BiGG_DB.mat`: BiGG database of reactions
- Functions
  - `PrepareForGapFilling()`: inputs GEM and reaction database, outputs merged model
  - `gapFilling()`: inputs merged model, outputs alternative solutions for gapfilling

### 💻 Exercise 3: Essentiality prediction
- Goal
  - Run single reaction essentiality (FBA & TFA-based) and single gene essentiality (FBA & TFA-based). Also run TFA-based single reaction essentiality on merged model to identify rescued reactions. These are used to identify targets for gapfilling.
- Files
  - `iML1515_KEGG_final_thermo.mat`: TFA-ready *E. coli* model with mapped metabolite IDs to ATLAS
- Functions
  - `singleRxnDeletion()`: inputs model, runs FBA-based single reaction knockouts, returns lethal reaction knockouts
  - `thermoSingleRxnDeletion()`: inputs TFA model, runs TFA-based single reaction knockouts, returns lethal reaction knockouts
  - `singleGeneDeletion()`: inputs model, runs FBA-based single gene knockout, returns lethal gene knckouts
  - `thermoSingleGeneDeletion()`: inputs model, runs TFA-based single gene knockout, returns lethal gene knockouts

### ⚖️ Exercise 4: Essentiality evaluation
- Goal
  - Obtain experimental data and compare against predictions (false negatives = targets for gapfilling, e.g. AMAOTr). Also compare merged model essentiality to identify rescued reactions (e.g. AMAOTr)
- Files
  - `experimental_essentiality.tsv`: gene knockout experiment results from [source publication](https://journals.asm.org/doi/10.1128/mbio.02096-17).
  - `gene_table.tsv`: lookup table to convert from publication gene IDs (e.g. glyA) to model gene IDs (e.g. b2551)
- Functions
  - `strcmp()`: used to compare values from predictions and experimental data
  - `readtable()`: used to load data files into MATLAB
  
### 📑 Exercise 5: Results reproduction in NICEgame (AMAOTr)
- Goal
  - We previosuly identified AMAOTr as target for gapfilling. Now we knockout this reaction in iML1515, merge with the ATLAS database, generate alternative solutions, and evaluate them (e.g. confusion matrix, growth rate)
- Files
  - `reducedATLAS_ecoli_yeast.mat`: ATLAS subset of reactions involving metabolites in *E. coli* and *S. cerevisiae*
- Functions
  
### 🧫 Exercise 6: NICEgame with BIOLOG data
- Goal
  - Obtained model + processed biolog growth data on different carbon sources, evaluate metabolite essentiality, choose false negative target for gapfilling (D-Galacturonic-Acid), replace as carbon source, suggest gapfilling alternatives. For one alternative, add the reaction to the model and re-evaluate TP/FP/TN/FN, improve MCC score. Repeat for FP gapfilling target (Arginine).
- Files
  - `keggModel.mat`: KEGG reaction database 
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
