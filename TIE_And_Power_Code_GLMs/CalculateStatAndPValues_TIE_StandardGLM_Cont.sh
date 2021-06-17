#!/bin/bash
if [ -f /etc/profile.d/modules.sh ]; then
   source /etc/profile.d/modules.sh 
fi

CHR=${1?Error: no chromosome given}
SS=${2?Error: no sample size given}
SEED=${3?Error: no random seed given}
GENE=${4?Error: no gene name given}
SCORE=${5?Error: no score given}
START=${6?Error: no start value given}
STOP=${7?Error: no stop value given}

export CHR
export SS
export SEED
export GENE
export SCORE
export START
export STOP

path="/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts"
cd $path

module load R/3.5.3
Rscript MainPipeline_TIE_StandardGLM_Cont.R $CHR $SS $SEED $GENE $SCORE $START $STOP --nosave