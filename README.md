# 🏓 NICEgame exercise 🏓

This repo contains code and files used during a 2 week visit to the [EPFL LCSB](https://www.epfl.ch/labs/lcsb/). The goal was to get hands on experience using the [NICEgame](https://www.pnas.org/doi/10.1073/pnas.2211197119) workflow to curate metabolic models based on experimental data. Through these exercises we also learn about thermodynamic-based flux analysis (TFA), BridgIT, and ATLASx.

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

### 🐤 Exercise 2: NICEgame tutorial
- Goal
  - Learn the basics of NICEgame: knockout reaction in iML1515 model (e.g. 5DOAN), merge model with BiGG database, generate alternative gapfilling solutions
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
  - Obtain experimental data and compare against predictions. False negatives are experimentally dispensible genes predicted to be lethal by gene knockout; these are targets for gapfilling, as this indicates that there must be some other reaction(s) compensating for the knockout. Also compare merged model essentiality to identify rescued reactions (e.g. AMAOTr)
- Files
  - `experimental_essentiality.tsv`: gene knockout experiment results from [source publication](https://journals.asm.org/doi/10.1128/mbio.02096-17).
  - `gene_table.tsv`: lookup table to convert from publication gene IDs (e.g. glyA) to model gene IDs (e.g. b2551)
- Functions
  - `strcmp()`: used to compare values from predictions and experimental data
  - `readtable()`: used to load data files into MATLAB
  
### 📑 Exercise 5: Reproduction of result from NICEgame paper (AMAOTr)
- Goal
  - We previosuly identified AMAOTr as target for gapfilling. Now we knockout this reaction in iML1515, merge with the ATLAS database, generate alternative solutions, and evaluate them (e.g. confusion matrix components, growth rate)
- Files
  - `reducedATLAS_ecoli_yeast.mat`: ATLAS subset of reactions involving metabolites in *E. coli* and *S. cerevisiae*
- Functions
  - `PrepareForGapFilling()`: inputs GEM and reaction database, outputs merged model
  - `gapFilling()`: inputs merged model, outputs alternative solutions for gapfilling 
  - `Essentiality()`:  inputs alternative gapfilling solutions, model, and essentiality experimental data, outputs confusion matrix components (FP, FN, TP, TN) and Matthews correlation coefficient (MCC) for each alternative
  - `TestFBA_growth()`: inputs merged model and alternative gapfilling solutions, outputs FBA-based growth rates for each alternative
  - `TestTFBA_growth()`: inputs merged model and alternative gapfilling solutions, outputs TFA-based growth rates for each alternative
  
### 🧫 Exercise 6: NICEgame with BIOLOG data
- Goal
  - Obtain model with associated biolog growth data on different carbon sources, evaluate metabolite essentiality by changing the carbon source in base media and comparing to experimental results.
  - False negatives compounds are experimentally viable carbon sources for growth, yet the model is unable to produce biomass, indicating they are targets for gapfilling. Choose false negative (e.g. D-Galacturonic-Acid), replace as carbon source, and suggest gapfilling alternatives. For one alternative, add the reaction to the model and re-evaluate confusion matrix components showing improvement in MCC score.
  - False positive compounds are not experimentally viable carbon sources for growth, yet the model still produces biomass, indicating that some reaction(s) should be turned off. Fix a false positive target (e.g. Arginine); the model can grow using Arginine as a carbon source but experimental data shows otherwise. In this case we have to identify candidate reactions to knockout that resolve false positive growth without introducing new errors. Find lethal reaction knockouts with Arginine as a carbon source, evaluate confusion matrix components, and identify knockouts that maximise MCC.
- Files
  - `keggModel.mat`: KEGG reaction database 
- Functions
  - `PrepareForGapFilling()`: inputs GEM and reaction database, outputs merged model
  - `gfbiomass()`: inputs merged model, outputs alternative solutions for gapfilling
  - `addReaction()`: inputs model and desired reaction to be added, outputs expanded model
  
### 🌉 Exercise 7: Submitting a reaction query to the BridgIT server
- Goal
  - Assign E.C. number for a given reaction, e.g. [KEGG reaction R00919](https://www.genome.jp/dbget-bin/www_bget?rn:R00919)
- Files
  - `R00919.zip`: zipfile for submission to BridgIT server, contains the following
     - `molfiles/`: folder with a .mol file for each metabolite involved in the query reaction, downloaded from KEGG (e.g. [Propanoyl-CoA](https://www.genome.jp/entry/C00100))
         - `C00005.mol`: NADPH
         - `C00006.mol`: NADP+
         - `C00080.mol`: H+
         - `C00100.mol`: Propanoyl-CoA
         - `C00894.mol`: Propenoyl-CoA
     - `test_systemfile.txt`: textfile specifying reaction string
        ```
        COMPOUNDS
        ENTRY
        reactionsS
        ENTRY;KEGG;EQUATION;OPERATORS
        1;;C00100+C00006<=>C00894+C00005+C00080;
        ```
- Server
   - https://lcsb-databases.epfl.ch/Bridgit/addProject
  
### 🗺️ Exercise 8: Submitting a metabolite query to the ATLASx server
- Goal
  - One false negative target from Exercise 6, Tween-80, cannot be fixed since there are no reactions that catalyze its degradation, submit this metabolite to ATLASx to suggest degradation routes in Exercise 8.
- Input
   - Tween-80 LCSB compound ID: 1467972114 
- Server
   - https://lcsb-databases.epfl.ch/Atlas2/Analysis
- Status
   - submitted to KEGG: blank result
   - submitted to bioATLAS: https://lcsb-databases.epfl.ch/Atlas2/GetResults/7828591
   - submitted to chemATLAS: https://lcsb-databases.epfl.ch/Atlas2/GetResults/4430457

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
