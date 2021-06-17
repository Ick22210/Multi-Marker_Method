############################################################
### Type I Error Simulation code
############################################################
#edit 5/21/20 added in ability to standardize scores and weight scores

#to obtain the arguments from the bash files (chr, ss)
args <- commandArgs()

chr <- args[6]
numPairs <- args[7]
seed <- args[8]
gene <- args[9]
YPrev <- args[10]
s <- args[11]

chr = as.integer(chr)
seed = as.integer(seed)
numPairs = as.integer(numPairs)
YPrev = as.integer(YPrev)/100
s = as.integer(s)

#set random seed for reproducibility
set.seed(seed)
#source needed functions
source('/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/MainPipeLineFunctionsToSource_v3.R')

#define score weights if wanted
scoreWeights = c()
############################################################
#for unweighted and unstandardized scores
RunTIEPipelineLocalGHat(chr = chr, gene = gene, numPairs = numPairs, YPrev = YPrev, s = s, standardizeScores = FALSE, weightedScores = FALSE)
#for standardized scores
#RunTIEPipelineLocalGHat(chr = chr, gene = gene, numPairs = numPairs, YPrev = YPrev, s = s, standardizeScores = TRUE, weightedScores = FALSE)
#for weighted scores
#RunTIEPipelineLocalGHat(chr = chr, gene = gene, numPairs = numPairs, YPrev = YPrev, s = s, standardizeScores = FALSE, weightedScores = FALSE, scoreWeights = scoreWeights)
