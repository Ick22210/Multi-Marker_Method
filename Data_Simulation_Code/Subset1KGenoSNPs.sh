#!/bin/bash
# Subsets Genes based on bp position from hap/legend files
# Vickie Arthur - June 4, 2020

#NAT2
/home/vlynn/Paper_II_Sims/HapGen_Files/plink2 --haps /home/vlynn/Paper_II_Sims/HapGen_Files/1000GP_Phase3/1000GP_Phase3_chr8.hap.gz \
--legend /home/vlynn/Paper_II_Sims/HapGen_Files/1000GP_Phase3/1000GP_Phase3_chr8.legend.gz 8 \
--chr 8 --from-bp 18386585 --to-bp 18401219 --export hapslegend -out NAT2_SNPs

# #CHI3L2
# /home/vlynn/Paper_II_Sims/HapGen_Files/plink2 --haps /home/vlynn/Paper_II_Sims/HapGen_Files/1000GP_Phase3/1000GP_Phase3_chr1.hap.gz \
# --legend /home/vlynn/Paper_II_Sims/HapGen_Files/1000GP_Phase3/1000GP_Phase3_chr1.legend.gz 1 \
# --chr 1 --from-bp 111743393 --to-bp 111786062 --export hapslegend -out CHI3L2_SNPs

# #ASAH1
# /home/vlynn/Paper_II_Sims/HapGen_Files/plink2 --haps /home/vlynn/Paper_II_Sims/HapGen_Files/1000GP_Phase3/1000GP_Phase3_chr8.hap.gz \
# --legend /home/vlynn/Paper_II_Sims/HapGen_Files/1000GP_Phase3/1000GP_Phase3_chr8.legend.gz 8 \
# --chr 8 --from-bp 17913925 --to-bp 17942507 --export hapslegend -out ASAH1_SNPs
