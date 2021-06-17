############################################################
### Power Simulation code
############################################################
#edit 5/21/20 added in ability to standardize scores and weight scores

#to obtain the arguments from the bash files (chr, ss)
args <- commandArgs()

chr <- args[6]
numPairs <- args[7]
seed <- args[8]
gene <- args[9]
YPrev <- args[10]
ORSize <- args[11]
percentageAssoc <- args[12]
lowLDTF <- args[13]
score <- args[14]
startVal <- args[15]
stopVal <- args[16]

chr = as.integer(chr)
seed = as.integer(seed)
numPairs = as.integer(numPairs)
YPrev = as.integer(YPrev)/100
percentageAssoc = as.integer(percentageAssoc)
startVal = as.integer(startVal)
stopVal = as.integer(stopVal)

#set random seed for reproducibility
set.seed(seed)
#source needed functions
source('/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/MainPipeLineFunctionsToSource_LinHypTest.R')

#make sure low LD is T or FALSE
if(lowLDTF == "TRUE"){
	LowLD = TRUE
} else {
	LowLD = FALSE
}

#define score weights if wanted
scoreWeights = c()

#define gammas
null_Gamma = c(0)

gammaEff = c(0)
############################################################
#for original, unstandardized and unweighted tests
#Low LD SNPs associated
RunPowerPipelineGLMCat_RSNP(chr = chr, gene = gene, numPairs = numPairs, YPrev = YPrev, Gamma = gammaEff, TrueScore = score, ORSize = ORSize, standardizeScores = FALSE, weightedScores = FALSE, start = startVal, stop = stopVal, percentageAssoc = percentageAssoc, LowLD = LowLD)
