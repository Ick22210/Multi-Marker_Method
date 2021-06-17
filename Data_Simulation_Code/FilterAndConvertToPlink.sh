#!/bin/bash

/home/vlynn/Paper_II_Sims/HapGen_Files/plink2 \
--gen /home/vlynn/Paper_II_Sims/HapGen_Files/${GENE}_Results_${SS}Pairs/${GENE}_${SS}Pairs_SimNum${j}.out.controls.gen ref-first \
--oxford-single-chr ${CHR} \
--sample /home/vlynn/Paper_II_Sims/HapGen_Files/${GENE}_Results_${SS}Pairs/${GENE}_${SS}Pairs_SimNum${j}.out.controls.sample \
--maf 0.001 \
--geno 0.25 \
--hwe 0.001 \
--make-bed \
--out /home/vlynn/Paper_II_Sims/HapGen_Files/${GENE}_Results_${SS}Pairs/${GENE}_${SS}Pairs_SimNum${j}
