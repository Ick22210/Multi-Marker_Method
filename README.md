# Multi-Marker_Method
Code for simulations from our multi-marker method paper (JST)

You will have to change the directories and locations of the data and scripts, but otherwise everything should run fine.

Directories are organized as follows:

1) Data Simulation Code: Contains all code used to simulate genetic data from the genes NAT2, CHI3L2, and ASAH1. Data is subset from 1000 Genomes reference using the program Hapgen2. Code requires plink2 and hapgen2 software (located in the folder). The Subset1KGenoSNPs.sh code has already been run to produce the .haps and .legend files for each of the 3 genes. If you want to use other gene regions, you will have to subset for them with this code. 1000 Genomes reference data can be found online, or I have it stored if you need other chromosomes. It is too large to upload to github right now but I think I downloaded them from here: https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html.

2) JST_sPVE_Code: Contains the shell scripts and R code to run TIE and power simulations for JST, where s = PVE. All of this code is dependent on the code in the Sourced_Functions folder. 

3) SKAT_Code: Contains the shell scripts and R code to run TIE and power simulations for SKAT. The files containing "DoF" give a calculation of estimated degrees of freedom for SKAT, since that is not provided by the original package. This code is also dependent on the code in the Sourced_Functions folder.

4) Sourced_Functions: The main functions and full pipelines that are run using the shell scripts and code in the other folders. Be sure to check for sourced R libraries. For my simulations, code was either run on R/3.5.3 or R/4.0 since that is what was present on our computing cluster. Be sure to install the packages you need for the correct version of R. You can also install any R package needed to your local directory on a cluster computer and then source it from there. 

5) TIE_And_Power_Code_GLMs: Code for the comparison with generalized linear model. Similarly to the other directories, it contains the R files and shell scripts for running TIE and Power simulations. 