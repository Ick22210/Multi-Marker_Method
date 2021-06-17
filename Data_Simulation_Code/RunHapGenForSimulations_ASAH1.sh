#!/bin/bash
#will have to change the dl value if using for other genes, but figure that out later

#For ASAH1
/home/vlynn/Paper_II_Sims/HapGen_Files/hapgen2 \
-h /home/vlynn/Paper_II_Sims/HapGen_Files/1000GP_Phase3/${GENE}_SNPs.haps \
-l /home/vlynn/Paper_II_Sims/HapGen_Files/1000GP_Phase3/${GENE}_SNPs.legend \
-m /home/vlynn/Paper_II_Sims/HapGen_Files/1000GP_Phase3/genetic_map_chr${CHR}_combined_b37.txt \
-n ${SAMPLES} 0 \
-o /home/vlynn/Paper_II_Sims/HapGen_Files/${GENE}_Results_${SS}Pairs/${GENE}_${SS}Pairs_SimNum${j}.out \
-t /home/vlynn/Paper_II_Sims/HapGen_Files/1000GP_Phase3/${GENE}_Subset.tags \
-dl 17913940 1 1.5 2.25 