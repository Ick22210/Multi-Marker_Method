############################################################
### Type I Error Simulation code
############################################################
#to obtain the arguments from the bash files (chr, ss)
args <- commandArgs()

chr <- args[6]
numPairs <- args[7]
seed <- args[8]
gene <- args[9]
score <- args[10]
startVal <- args[11]
stopVal <- args[12]

chr = as.integer(chr)
seed = as.integer(seed)
numPairs = as.integer(numPairs)
startVal = as.integer(startVal)
stopVal = as.integer(stopVal)

#set random seed for reproducibility
set.seed(seed)
#source needed functions
source('/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/MainPipeLineFunctionsToSource_LinHypTest.R')

#define score weights if wanted
scoreWeights = c()
############################################################
#for unweighted and unstandardized scores
RunTIEPipelineGLMCont(chr = chr, gene = gene, numPairs = numPairs, standardizeScores = FALSE, weightedScores = FALSE, score = score, start = startVal, stop = stopVal)
