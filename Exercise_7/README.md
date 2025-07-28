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
