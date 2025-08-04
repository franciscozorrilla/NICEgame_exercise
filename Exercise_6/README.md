### 🌉 Exercise 6: Submitting a reaction query to the BridgIT server
- Goal
  - Submit a job to assign an EC number for a given reaction, e.g. [KEGG reaction R00919 associated with EC 1.3.1.84](https://www.genome.jp/dbget-bin/www_bget?rn:R00919)

1. Navigate to the [LCSB website](https://lcsb-databases.epfl.ch/Home) and request an account (free for academia) in order to access servers and databases
2. Once you have been granted access, click on the `BridgIT` tab and then click `Analyze`
3. Submit the `R00919.zip` file stored in this repo, or create your own input file for a reaction of interest based on the input file descriptions below
- Input files
  - `R00919.zip`: zipfile for submission to BridgIT server, contains the following
     - `molfiles/`: folder with a `.mol` file for each metabolite involved in the query reaction, downloaded from KEGG (e.g. [Propanoyl-CoA](https://www.genome.jp/entry/C00100)). MOL files (`.mol`) are a text-based chemical file format that describes the atoms, bonds, and connectivity of a given metabolite.
         - `C00005.mol`: NADPH
         - `C00006.mol`: NADP+
         - `C00080.mol`: H+
         - `C00100.mol`: Propanoyl-CoA
         - `C00894.mol`: Propenoyl-CoA
     - `test_systemfile.txt`: textfile specifying reaction string ⚠️ no spaces between metabolite IDs
        ```
        COMPOUNDS
        ENTRY
        reactionsS
        ENTRY;KEGG;EQUATION;OPERATORS
        1;;C00100+C00006<=>C00894+C00005+C00080;
        ```
4. Wait a few moments until job runs and output file is produced and begins downloading
- Output file
  - `ResultsReactionsBridgIt_6731978/Tanimoto_Atom_reactions_1.txt`: downloaded output from BridgIT server with the following columns
     - `reactionsA/ECA`: KEGG reaction/EC number hit
     - `reactionsB/ECB`: query reaction/EC number as determined from input files
     - `Tanimoto_FBI_Scores`: overall reaction similarity score between query and hit
     - `TL0` - `TL7`: similarity score at increasing bond distances from the active site (see [methods](https://www.pnas.org/doi/10.1073/pnas.1818877116#sec-2) for more details)
5. Inspect the output file
- We can see that there are two top hits with perfect reaction similarity to the query, corresponding to R10161/1.3.1.95 and R00919/1.3.1.84
    ```bash
    $ head ~/path/to/repo/NICEgame_exercise/Exercise_6/ResultsReactionsBridgIt_6731978/Tanimoto_Atom_reactions_1.txt
    
    reactionsA/ECA	reactionsB/ECB	Tanimoto_FBI_Scores	 TL0	TL1	TL2	TL3	TL4	TL5	TL6	TL7
    R10161(r)/1.3.1A1(rev)/1.3.1.95;	1/1.3.1A1(rev);	1	1	1	1	1	1	1	1	1
    R00919(r)/1.3.1A2(rev)/1.3.1.84;	1/1.3.1A1(rev);	1	1	1	1	1	1	1	1	1
    R09738(r)/1.3.1A2(rev)/1.3.1.86;	1/1.3.1A1(rev);	0.968992	1	0.909091	0.941176	0.96	0.942857	0.978261	0.981132	0.984127
    R10828(r)/1.3.1A2(rev)/1.3.1.93;	1/1.3.1A1(rev);	0.968992	1	0.909091	0.941176	0.96	0.942857	0.978261	0.981132	0.984127
    R01171(r)/1.3.1A1(rev)/1.3.1.44;	1/1.3.1A1(rev);	0.968992	1	0.909091	0.941176	0.96	0.942857	0.978261	0.981132	0.984127
    R06985/1.3.1A2(rev)/1.3.1.38,1.3.1.8;	1/1.3.1A1(rev);	0.922509	1	0.909091	0.888889	0.888889	0.891892	0.918367	0.928571	0.953846
    R03776(r)/1.3.1A2(rev)/1.3.1.38,1.3.1.8;	1/1.3.1A1(rev);	0.896057	1	0.909091	0.888889	0.888889	0.868421	0.882353	0.896552	0.911765
    R04753(r)/1.3.1A2(rev)/1.3.1.38,1.3.1.8;	1/1.3.1A1(rev);	0.886525	1	0.909091	0.888889	0.888889	0.868421	0.882353	0.881356	0.885714
    R07761/1.3.1A2(rev)/1.3.1.38,1.3.1.93;	1/1.3.1A1(rev);	0.886525	1	0.909091	0.888889	0.888889	0.868421	0.882353	0.881356	0.885714
    ```
- As expected, BridgIT is able to correctly assign the EC number 1.3.1.84 to the KEGG reaction [R00919](https://www.genome.jp/dbget-bin/www_bget?rn:R00919) (this was the input). However, we also see that there is an equally confident hit called for the EC number 1.3.1.95 to the KEGG reaction R10161. Visually inspect and compare the [two reactions on KEGG](https://www.genome.jp/dbget-bin/www_bget?rn:R10161+R00919). How do the reactions differ? Does it make sense to you that these two reactions have such a high and identical score?
   
[Next exercise](https://github.com/franciscozorrilla/NICEgame_exercise/tree/master/Exercise_7)
