#!/bin/bash

CHR=${1?Error: no chr given}
SS=${2?Error: no sample size given}
NSIMS=${3?Error: no number of simulations given}
GENE=${4?Error: no gene name given}
SAMPLES=${5?Error:no number of samples given}

export CHR
export SS
export NSIMS
export GENE
export SAMPLES

#create directory to hold files that should not be erased
mkdir /home/vlynn/Paper_II_Sims/HapGen_Files/${GENE}_Results_${SS}Pairs

for ((j=1;j<=NSIMS;j++));
do		
	export j 
	
	#use hapgen to create genotypes
	if [ ${GENE} == "NAT2" ]; then
		sh /home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/RunHapGenForSimulations_NAT2.sh
	elif [ ${GENE} == "CHI3L2" ]; then
		sh /home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/RunHapGenForSimulations_CHI3L2.sh
	elif [ ${GENE} == "ASAH1" ]; then
		sh /home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/RunHapGenForSimulations_ASAH1.sh
	elif [ ${GENE} == "ABLIM3" ]; then
		sh /home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/RunHapGenForSimulations_ABLIM3.sh
	else 
		sh /home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/RunHapGenForSimulations_ADAMTS12.sh
	fi
	
	#remove cases files
	rm /home/vlynn/Paper_II_Sims/HapGen_Files/${GENE}_Results_${SS}Pairs/*.cases*
	
	#convert to plink file format
	sh /home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/FilterAndConvertToPlink.sh
	
	#remove files that aren't bed/bim/fam or plink logs, legend or summary files
	rm /home/vlynn/Paper_II_Sims/HapGen_Files/${GENE}_Results_${SS}Pairs_RareMAFs/*.controls.*
done