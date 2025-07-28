### 🐤 Exercise 2: NICEgame tutorial
- Goal
  - Learn the basics of NICEgame: knockout reaction in iML1515 model (e.g. 5DOAN), merge model with BiGG database, generate alternative gapfilling solutions
- Files
  - `iML1515.mat`: *E. coli* model
  - `BiGG_DB.mat`: BiGG database of reactions
- Functions
  - `PrepareForGapFilling()`: inputs GEM and reaction database, outputs merged model
  - `gapFilling()`: inputs merged model, outputs alternative solutions for gapfilling
