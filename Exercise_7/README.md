### 🗺️ Exercise 7: Submitting a metabolite query to the ATLASx server
- Goal
  - One false negative target from Exercise 5, Tween-80, cannot be fixed since there are no reactions that catalyze its degradation, submit this metabolite to ATLASx to suggest degradation routes in Exercise 8.
1. Navigate to the [LCSB website](https://lcsb-databases.epfl.ch/Home) and request an account (free for academia) in order to access servers and databases
2. Once you have been granted access, click on the ATLASx tab and then click `Search`
3. Find the LCSB compound ID for your metabolite of interest ⚠️ Note that you must manually refresh the website in between searches of metabolites
4. Next, click on the `Analysis` tab and input the LCSB metabolite ID
- Input
   - LCSB compound ID (Tween-80): 1467972114, or other metabolite of interest
   - Number of generations: 1, this specifies the number of rounds that BNICE rules should be applied query input and its predicted degradation products; use a low value to avoid long runtimes
   - Search scope: different databases to query depending on compound and research question
- Status
   - Search scope KEGG: no hits (empty output)
   - Search scope bioATLAS: https://lcsb-databases.epfl.ch/Atlas2/GetResults/7828591
   - Search scope chemATLAS: https://lcsb-databases.epfl.ch/Atlas2/GetResults/4430457
