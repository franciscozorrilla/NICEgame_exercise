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
