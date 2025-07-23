# 🏓 NICEgame exercise 🏓

- Exercise 1: matTFA tutorial
    - Basics of converting GEM for TFA, compared TFA and FBA results to show that TFA is more restrictive
- Exercise 2: NICEgame tutorial
    - KO reaction in iML1515, merge model with BiGG database, generate alternative gapfilling solutions
- Exercise 3: Essentiality prediction
    - Single reaction essentiality (FBA & TFA), single gene essentiality (FBA & TFA)
- Exercise 4: Essentiality evaluation
    - Obtain experimental data and compare against predictions (FN = targets for gapfilling, e.g. AMAOTr). Also compare merged model essentiality to identify rescued reactions (e.g. AMAOTr)
- Exercise 5: Results reproduction in NICEgame
     - Identified AMAOTr as target for gapfilling, KO in iML1515, merge with ATLAS database, generate alternative solutions and evaluate them (e.g. confusion matrix, growth rate)
- Exercise 6: NICEgame with BIOLOG data
     - Obtained model + processed biolog growth data on different carbon sources, evaluate metabolite essentiality, choose false negative target for gapfilling (D-Galacturonic-Acid), replace as carbon source, suggest gapfilling alternatives. For one alternative, add the reaction to the model and re-evaluate TP/FP/TN/FN, improve MCC score. Repeat for FP gapfilling target (Arginine).
- Exercise 7: Submitting a reaction query to the BridgIT server
     - Identify E.C. number for KEGG reaction R00919 
- Exercise 8: Submitting a metabolite query to the ATLASx server
     - Identify reactions to degrade tween-80

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
