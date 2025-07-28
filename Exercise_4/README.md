### ⚖️ Exercise 4: Essentiality evaluation
- Goal
  - Obtain experimental data and compare against predictions. False negatives are experimentally dispensible genes predicted to be lethal by gene knockout; these are targets for gapfilling, as this indicates that there must be some other reaction(s) compensating for the knockout. Also compare merged model essentiality to identify rescued reactions (e.g. AMAOTr)

- Files
  - `experimental_essentiality.tsv`: gene knockout experiment results from [source publication](https://journals.asm.org/doi/10.1128/mbio.02096-17).
  - `gene_table.tsv`: lookup table to convert from publication gene IDs (e.g. glyA) to model gene IDs (e.g. b2551)
- Functions
  - `strcmp()`: used to compare values from predictions and experimental data
  - `readtable()`: used to load data files into MATLAB
