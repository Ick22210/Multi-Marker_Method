######################################################################
## Code Examples - Joint Testing 
######################################################################

#define path to where data is located
#also where the file with source functions is located
path = paste0("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Project_II/Real_Data_Analysis/Real_Data_Analysis_Code/Code_To_Send")

#source pipeline functions
source(paste0(path,"/RealDataPipeline_v2.R"))

#Binary Phenotypes
datafile = "Example_Datafile.csv"
matchedPairs = "Example_PairIds.txt"
covPhenoFile = "Example_Covariates_BinPhenos.txt"
s = 1
score = "IBS"
kernel = "linear" # "IBS"
outFileName1 = "Example_Outfile_BinPhenos.csv"
outFileName1_SKAT = "Example_Outfile_BinPhenos_SKAT.csv"

# #unweighted
# All 4 Scores at once
realDataPipeline_JST(datafile = datafile, matchedPairs = matchedPairs, covPhenoFile = covPhenoFile, dataPath = dataPath,
                     standardizeScores = FALSE, weightedScores = FALSE, phenosBinary = TRUE, s = s, score = "All", outFile = outFileName1)
# testing the 4 scores separately
realDataPipeline_JST(datafile = datafile, matchedPairs = matchedPairs, covPhenoFile = covPhenoFile, dataPath = dataPath,
                     standardizeScores = FALSE, weightedScores = FALSE, phenosBinary = TRUE, s = s, score = "IBS", outFile = "Example_Outfile_BinPhenos_IBS.csv")
realDataPipeline_JST(datafile = datafile, matchedPairs = matchedPairs, covPhenoFile = covPhenoFile, dataPath = dataPath,
                     standardizeScores = FALSE, weightedScores = FALSE, phenosBinary = TRUE, s = s, score = "Incomp", outFile = "Example_Outfile_BinPhenos_Incomp.csv")
realDataPipeline_JST(datafile = datafile, matchedPairs = matchedPairs, covPhenoFile = covPhenoFile, dataPath = dataPath,
                     standardizeScores = FALSE, weightedScores = FALSE, phenosBinary = TRUE, s = s, score = "AMS", outFile = "Example_Outfile_BinPhenos_AMS.csv")
realDataPipeline_JST(datafile = datafile, matchedPairs = matchedPairs, covPhenoFile = covPhenoFile, dataPath = dataPath,
                     standardizeScores = FALSE, weightedScores = FALSE, phenosBinary = TRUE, s = s, score = "Bin MM", outFile = "Example_Outfile_BinPhenos_BinMM.csv")

# testing the SKAT pipeline
realDataPipeline_SKAT(datafile = datafile, matchedPairs = matchedPairs, covPhenoFile = covPhenoFile, dataPath = dataPath,
                      standardizeScores = FALSE, weightedScores = FALSE, phenosBinary = TRUE, kernel = kernel, score = score, outFile = outFileName1_SKAT)

# IBS kernel
realDataPipeline_SKAT(datafile = datafile, matchedPairs = matchedPairs, covPhenoFile = covPhenoFile, dataPath = dataPath,
                      standardizeScores = FALSE, weightedScores = FALSE, phenosBinary = TRUE, kernel = "IBS", score = score, outFile = "Example_Outfile_BinPhenos_SKAT_IBS.csv")

#Continuous Phenotypes
covPhenoFileCont = "Example_Covariates_ContPhenos.txt"
outFileName_cont = "Example_Outfile_ContPhenos.csv"
outFileName_cont_SKAT = "Example_Outfile_ContPhenos_SKAT.csv"

# All 4 Scores at once
realDataPipeline_JST(datafile = datafile, matchedPairs = matchedPairs, covPhenoFile = covPhenoFileCont, dataPath = dataPath,
                     standardizeScores = FALSE, weightedScores = FALSE, phenosBinary = FALSE, s = s, score = "All", outFile = outFileName_cont)
# SKAT pipeline
realDataPipeline_SKAT(datafile = datafile, matchedPairs = matchedPairs, covPhenoFile = covPhenoFileCont, dataPath = dataPath,
                      standardizeScores = FALSE, weightedScores = FALSE, phenosBinary = FALSE, kernel = kernel, score = "All", outFile = outFileName_cont_SKAT)

######################################################################
## Code Examples - Testing Gene-Score Only (Gamma = 0)
######################################################################
#Binary Phenotypes
datafile = "Example_Datafile.csv"
matchedPairs = "Example_PairIds.txt"
covPhenoFile = "Example_Covariates_BinPhenos.txt"
dataPath = path

outfileName1_gamma = "Example_Outfile_BinPhenos_GammaOnly.csv"
outfileName2_gamma = "Example_Outfile_BinPhenos_GammaOnly_SingleScore.csv"
outfileName3_gamma = "Example_Outfile_BinPhenos_GammaOnly_NoIntercept.csv"
outfileName4_gamma = "Example_Outfile_BinPhenos_GammaOnly_SpecifyUnpen.csv"

#defaults
realDataPipeline_GeneScoreTest(datafile = datafile, matchedPairs = matchedPairs, covPhenoFile = covPhenoFile, dataPath = dataPath,
                               standardizeScores = FALSE, weightedScores = FALSE, phenosBinary = TRUE, score = "All", fitIntercept = TRUE, 
                               unpen = c(), outFile = outfileName1_gamma)

#single score
realDataPipeline_GeneScoreTest(datafile = datafile, matchedPairs = matchedPairs, covPhenoFile = covPhenoFile, dataPath = dataPath,
                               standardizeScores = FALSE, weightedScores = FALSE, phenosBinary = TRUE, score = "AMS", fitIntercept = TRUE, 
                               unpen = c(), outFile = outfileName2_gamma)
#no intercept
realDataPipeline_GeneScoreTest(datafile = datafile, matchedPairs = matchedPairs, covPhenoFile = covPhenoFile, dataPath = dataPath,
                               standardizeScores = FALSE, weightedScores = FALSE, phenosBinary = TRUE, score = "Incomp", fitIntercept = FALSE, 
                               unpen = c(), outFile = outfileName3_gamma)
#specify which covariates remain unpenalized (age, gender and intercept)
realDataPipeline_GeneScoreTest(datafile = datafile, matchedPairs = matchedPairs, covPhenoFile = covPhenoFile, dataPath = dataPath,
                               standardizeScores = FALSE, weightedScores = FALSE, phenosBinary = TRUE, score = "IBS", fitIntercept = FALSE, 
                               unpen = c(8), outFile = outfileName4_gamma)

#Continuous Phenotypes
covPhenoFileCont = "Example_Covariates_ContPhenos.txt"
outFileName_gamma_cont = "Example_Outfile_ContPhenos_GammaOnly.csv"
outfileName2_gamma_cont = "Example_Outfile_ContPhenos_GammaOnly_SingleScore.csv"
outfileName3_gamma_cont = "Example_Outfile_ContPhenos_GammaOnly_NoIntercept.csv"
outfileName4_gamma_cont = "Example_Outfile_ContPhenos_GammaOnly_SpecifyUnpen.csv"

#defaults
realDataPipeline_GeneScoreTest(datafile = datafile, matchedPairs = matchedPairs, covPhenoFile = covPhenoFileCont, dataPath = dataPath,
                               standardizeScores = FALSE, weightedScores = FALSE, phenosBinary = FALSE, score = "All", fitIntercept = TRUE, 
                               unpen = c(), outFile = outFileName_gamma_cont)

#single score
realDataPipeline_GeneScoreTest(datafile = datafile, matchedPairs = matchedPairs, covPhenoFile = covPhenoFileCont, dataPath = dataPath,
                               standardizeScores = FALSE, weightedScores = FALSE, phenosBinary = FALSE, score = "AMS", fitIntercept = TRUE, 
                               unpen = c(), outFile = outfileName2_gamma_cont)
#no intercept
realDataPipeline_GeneScoreTest(datafile = datafile, matchedPairs = matchedPairs, covPhenoFile = covPhenoFileCont, dataPath = dataPath,
                               standardizeScores = FALSE, weightedScores = FALSE, phenosBinary = FALSE, score = "Incomp", fitIntercept = FALSE, 
                               unpen = c(), outFile = outfileName3_gamma_cont)
#specify which covariates remain unpenalized (age, gender and intercept)
realDataPipeline_GeneScoreTest(datafile = datafile, matchedPairs = matchedPairs, covPhenoFile = covPhenoFileCont, dataPath = dataPath,
                               standardizeScores = FALSE, weightedScores = FALSE, phenosBinary = FALSE, score = "IBS", fitIntercept = FALSE, 
                               unpen = c(8), outFile = outfileName4_gamma_cont)
