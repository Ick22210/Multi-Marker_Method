#!/bin/bash
if [ -f /etc/profile.d/modules.sh ]; then
   source /etc/profile.d/modules.sh 
fi

CHR=${1?Error: no chromosome given}
SS=${2?Error: no sample size given}
SEED=${3?Error: no random seed given}
GENE=${4?Error: no gene name given}
YPREV=${5?Error: no outcome prevalence given}
S=${6?Error: no percent var explained given}

export CHR
export SS
export SEED
export GENE
export YPREV
export S

path="/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts"
cd $path

module load R/3.5.3
Rscript MainPipeline_TIE_Ghat.R $CHR $SS $SEED $GENE $YPREV $S --nosave