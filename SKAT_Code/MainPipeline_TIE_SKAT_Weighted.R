############################################################
### Type I Error Simulation code
############################################################
#to obtain the arguments from the bash files (chr, ss)
args <- commandArgs()

chr <- args[6]
numPairs <- args[7]
seed <- args[8]
gene <- args[9]
YPrev <- args[10]

chr = as.integer(chr)
seed = as.integer(seed)
numPairs = as.integer(numPairs)
YPrev = as.integer(YPrev)/100

#set kernel 
kern = "linear.weighted"

#set random seed for reproducibility
set.seed(seed)
#source needed functions
source('/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/MainPipeLineFunctionsToSource_ScoreOrRGenoWeighting.R')

#define score weights if wanted
NAT2_MAF = c(0.085, 0.062, 0.077, 0.072, 0.071,
             0.095, 0.233, 0.249, 0.077, 0.062,
             0.067, 0.11, 0.188, 0.187, 0.19,
             0.419, 0.188, 0.238, 0.236)

CHI3L2_MAF = c(0.204, 0.208, 0.137, 0.311, 0.334,
               0.334, 0.347, 0.457, 0.139, 0.288,
               0.449, 0.061, 0.055, 0.07, 0.055,
               0.27, 0.056, 0.055, 0.149, 0.147,
               0.121, 0.121, 0.06, 0.132, 0.28,
               0.281, 0.254, 0.172, 0.056, 0.256,
               0.309, 0.256, 0.33, 0.066, 0.405,
               0.411, 0.072, 0.352, 0.069, 0.069,
               0.073, 0.353, 0.158, 0.24)

ASAH1_MAF = c(0.08, 0.106, 0.106, 0.058, 0.106,
              0.082, 0.087, 0.418, 0.364, 0.393,
              0.078, 0.076, 0.077, 0.079, 0.123,
              0.111, 0.113, 0.081, 0.113, 0.113,
              0.096, 0.096, 0.11, 0.11, 0.106,
              0.11, 0.11, 0.11, 0.106, 0.11,
              0.478, 0.105, 0.108, 0.108, 0.071,
              0.108, 0.108, 0.108, 0.101, 0.454,
              0.454, 0.456 ,0.448, 0.446, 0.443,
              0.453, 0.443, 0.13, 0.439, 0.446, 
              0.446, 0.444, 0.444, 0.45, 0.311,
              0.448, 0.07, 0.454, 0.453, 0.069,
              0.446, 0.468, 0.096, 0.115, 0.076,
              0.06, 0.069, 0.069)

NAT2_MAF_Weights = 1/sqrt(NAT2_MAF)
CHI3L2_MAF_Weights = 1/sqrt(CHI3L2_MAF)
ASAH1_MAF_Weights = 1/sqrt(ASAH1_MAF)

NAT2_MAF_Weights_SKAT = sqrt(NAT2_MAF_Weights)
CHI3L2_MAF_Weights_SKAT = sqrt(CHI3L2_MAF_Weights)
ASAH1_MAF_Weights_SKAT = sqrt(ASAH1_MAF_Weights)

#define score weights if wanted
if(gene == "NAT2"){
	scoreWeights = NAT2_MAF_Weights_SKAT
} else if(gene == "CHI3L2"){
	scoreWeights = CHI3L2_MAF_Weights_SKAT
} else {
	scoreWeights = ASAH1_MAF_Weights_SKAT
}
############################################################
#for weighted scores
#RunTIEPipelineSKAT(chr = chr, gene = gene, numPairs = numPairs, YPrev = YPrev, kernel = kern, kernelWeights = scoreWeights, standardizeScores = FALSE, weightedScores = TRUE, weightedME = FALSE)
RunTIEPipelineSKAT(chr = chr, gene = gene, numPairs = numPairs, YPrev = YPrev, kernel = kern, kernelWeights = scoreWeights, standardizeScores = FALSE, weightedScores = FALSE, weightedME = TRUE)
#RunTIEPipelineSKAT(chr = chr, gene = gene, numPairs = numPairs, YPrev = YPrev, kernel = kern, kernelWeights = scoreWeights, standardizeScores = FALSE, weightedScores = TRUE, weightedME = TRUE)
