#!/bin/bash
if [ -f /etc/profile.d/modules.sh ]; then
   source /etc/profile.d/modules.sh 
fi

CHR=${1?Error: no chromosome given}
SS=${2?Error: no sample size given}
SEED=${3?Error: no random seed given}
GENE=${4?Error: no gene name given}
YPREV=${5?Error: no outcome prevalence given}
ORSize=${6?Error: no OR effect size given}
percentageAssoc=${7?Error: no % of associated SNPs given}
lowLDTF=${8?Error: no true or false given for low LD}
score=${9?Error: no true score given}
startVal=${10?Error: no start value given}
stopVal=${11?Error: no stop value given}

export CHR
export SS
export SEED
export GENE
export YPREV
export ORSize
export percentageAssoc
export lowLDTF
export score
export startVal
export stopVal

path="/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts"
cd $path

module load R/4.0
Rscript MainPipeline_Power_GLM_RSNPs_Cat.R $CHR $SS $SEED $GENE $YPREV $ORSize $percentageAssoc $lowLDTF $score $startVal $stopVal --nosave