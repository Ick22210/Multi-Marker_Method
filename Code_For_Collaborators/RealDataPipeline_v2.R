##########################################################
## Real Data Analysis Pipeline
## 
## Written by Victoria Arthur
## Last Edited: 05/26/2020
## Dependent on: RealDataSourceFunctions.R
##########################################################
#Input data formats:
#
## datafile:
# Data should be in matrix-type form, any file format is fine (text file, csv, etc) 
# with number of rows = 2*number of D/R pairs (ie all Donors and Recips from the pairs)
# and num of columns should correspond to SNPs, 
# columns should have values between 0-2 for additive model
# Should also have rownames that are the Ids of the Donor and Recips for matching
#
## matchedPairs:
# A file of Ids, first column is Recipient ID, second column is matched Donor ID
# This file should be in ".txt" format
#
## covPhenoFile:
# Additional covariates should be in ".txt" format as well
# first column: Recipient ID (rownames)
# second column: Phenotype of interest (Y), either binary or continous
# remaining columns: all other covariates (must be numeric - so gender should be coded as 0/1 for male/female for example)

#function to run joint score test (JST) analyses on real data
realDataPipeline_JST = function(datafile = "genotypeData.csv", matchedPairs = "pairedIDs.txt", covPhenoFile = "covariatesAndPhenotypes.txt", dataPath = path,
                            standardizeScores = FALSE, weightedScores = FALSE, scoreWeights, phenosBinary = TRUE, s, score = "All", outFile = "out.csv"){
  #Inputs: 
  ## datafile: the name of the file containing the genotype data
  ## matchedPairs: a text file containing the information on which donor id is matched to which recipient id
  ## covPhenoFile: a text file containing additional covariate information, the first column of the covariate file should be the Recipient ids,
  ##                the second column should be the phenotype of interest, additional columns are other covariates (all numeric)
  ## dataPath: the path to where the data is located
  ## standardizeScores: TRUE or FALSE for whether scores should be standardized
  ## weightedScores: TRUE or FALSE for whether scores should be weighted
  ## scoreWeights: Weights for each SNP if weightedScores = TRUE
  ## phenosBinary: TRUE or FALSE for whether the phenotype of interest is binary, if FALSE phenotype must be continuous
  ## s: the number of principal components you want to keep, should be an integer of 1 or greater.
  ## score: the gene-based score used in the model, either 'IBS', 'Incomp', 'AMS', 'Bin MM', or 'All' to run all 4 at the same time
  ## outFile: name of file to write output to
  
  #Output: Output file containing the score values and p-values for the R Geno/IBS Score, R Geno/Incomp Score, R Geno/AMS, and R Geno/Bin MM Score
  
  #set working directory to where data is
  setwd(path)
  
  #source the needed functions
  source(paste0(path,"/RealDataSourceFunctions.R"))
  library(stringr)
  
  if(score == "All"){
    #define matrix to hold all Stats and Pvalues
    statsAndPVals = matrix(NA, nrow = 4, ncol = 2)
    
    #Read in the data
    #determine what file type
    fileName = read.table(text = datafile, sep = ".", as.is = TRUE)$V2
    if(fileName == "csv"){ #if it's a csv file, read in with read.csv
      geneticData = read.csv(datafile, header = TRUE, row.names = 1)
    } else { #otherwise should be good with a read.table
      geneticData = read.table(datafile, header = TRUE, row.names = 1)
    }
    
    geneticData = na.omit(geneticData)
    
    #if only 1 SNP, need to add in a column so it doesn't screw everything up
    if(dim(geneticData)[2] == 1){
      geneticData$Test = 0
    }
    ###############################################
    #Determine D/R pairs and match
    ###############################################  
    #read in the file with the donor/recipient pair information
    pairings = read.table(matchedPairs)
    #make sure column names are specified
    colnames(pairings) = c("R_Id", "D_Id")
    #assign each pair a number
    pairings$PairNumber = seq(1:length(pairings$R_Id))
    #split the original data into donor and recipient subsets
    donorData = geneticData[rownames(geneticData) %in% pairings$D_Id,]
    recipData = geneticData[rownames(geneticData) %in% pairings$R_Id,]
    #add an id column based on row names for merging
    donorData$D_Id = rownames(donorData)
    recipData$R_Id = rownames(recipData)
    #now merge the genotype data with the ids so that the genotype data has the pair numbers
    donorGenotypes_numbered = merge(donorData, pairings, by.x = "D_Id", by.y = "D_Id")
    recipGenotypes_numbered = merge(recipData, pairings, by.x = "R_Id", by.y = "R_Id")
    #make everything a data frame 
    donorGenotypes_numbered.df = as.data.frame(donorGenotypes_numbered)
    recipGenotypes_numbered.df = as.data.frame(recipGenotypes_numbered)
    #subset to only include Pairs with both D and R info
    donorGenotypes_numbered.df.subset = donorGenotypes_numbered.df[(donorGenotypes_numbered.df$PairNumber %in% recipGenotypes_numbered.df$PairNumber),]
    recipGenotypes_numbered.df.subset = recipGenotypes_numbered.df[(recipGenotypes_numbered.df$PairNumber %in% donorGenotypes_numbered.df$PairNumber),]
    #order by Pair number
    donorGenotypes_numbered.df.ordered = donorGenotypes_numbered.df.subset[order(donorGenotypes_numbered.df.subset$PairNumber),]
    recipGenotypes_numbered.df.ordered = recipGenotypes_numbered.df.subset[order(recipGenotypes_numbered.df.subset$PairNumber),]
    #make rownames into R_Id or D_Id
    rownames(donorGenotypes_numbered.df.ordered) = donorGenotypes_numbered.df.ordered$D_Id
    rownames(recipGenotypes_numbered.df.ordered) = recipGenotypes_numbered.df.ordered$R_Id
    #remove non-SNP columns
    DGenosData = donorGenotypes_numbered.df.ordered[,2:(ncol(donorGenotypes_numbered.df.ordered)-2)]
    RGenosData = recipGenotypes_numbered.df.ordered[,2:(ncol(recipGenotypes_numbered.df.ordered)-2)]
    
    ###############################################
    # Read in the covariate/phenotype file
    ############################################### 
    covariates = read.table(covPhenoFile, header = TRUE)
    #remove nas
    covariates = na.omit(covariates)
    #pull only those Rs from the R geno matrix
    covariates = covariates[covariates[,1] %in% rownames(RGenosData),]
    #set row names to be Ids from 1st row
    rownames(covariates) = covariates[,1]
    #remove the ids column
    covariates = covariates[,-1]
    
    if(dim(RGenosData)[1] > dim(covariates)[1]){
      RGenosData = RGenosData[rownames(RGenosData) %in% rownames(covariates),]
      DGenosData = DGenosData[which(rownames(RGenosData) %in% rownames(covariates)),]
    } else if (dim(RGenosData)[1] < dim(covariates)[1]){
      covariates = covariates[rownames(covariates) %in% rownames(RGenosData),]
    }
    
    #make matrices
    #unlist and make into matrices
    RGenosMat = matrix(unlist(RGenosData), ncol = ncol(RGenosData), byrow = F)
    DGenosMat = matrix(unlist(DGenosData), ncol = ncol(DGenosData), byrow = F)
    #get rid of the Test column if it exists
    if(("Test" %in% colnames(RGenosData)) == TRUE){
      RGenosMat = matrix(RGenosMat[,1], ncol = 1)
      DGenosMat = matrix(DGenosMat[,1], ncol = 1)
    }
    #keep number of pairs for later use
    numPairs = nrow(RGenosMat)
    numSNPs = ncol(RGenosMat)
    
    #output number of pairs
    strings = str_split(datafile, ".csv", simplify = TRUE)
    # filenameNumPairs = paste0(strings[1],"_numPairs.txt")
    # write.table(numPairs, filenameNumPairs)
    
    #write out number of SNPs
    filenameNumSNPs = paste0(strings[1],"_numSNPs.txt")
    write.table(numSNPs, filenameNumSNPs)
    
    #make a matrix
    #pull the phenotypes from the covariates file
    phenos = as.matrix(covariates[,1],nrow=numPairs, ncol = 1)
    #remove the phenotypes column from the covariates
    CovData.list = covariates[,-1]
    CovData = matrix(unlist(CovData.list), ncol = ncol(CovData.list), byrow = F)
    
    ###############################################
    # Calculate Gene-based Scores
    ###############################################   
    # First: Calculate single snp scores
    IBS.snp = calcIBSMismatch(RGenosMat, DGenosMat)
    Incomp.snp = calcIncompatibilityScore(RGenosMat, DGenosMat)
    AMS.snp = calcAMS(RGenosMat, DGenosMat)
    BinMM.snp = calcBinaryMM(RGenosMat, DGenosMat)
    
    # Then: Calculate gene based scores
    ##Check to see if weights are used for scores
    if(weightedScores == FALSE){
      #check to see if scores should be standardized (can't be weighted and standardized)
      if(standardizeScores == FALSE){
        IBS.gene = calcGeneScore(SingleSNPKernel = IBS.snp, standardize = FALSE, useWeights = FALSE)
        Incomp.gene = calcGeneScore(SingleSNPKernel = Incomp.snp, standardize = FALSE, useWeights = FALSE)
        AMS.gene = calcGeneScore(SingleSNPKernel = AMS.snp, standardize = FALSE, useWeights = FALSE)
        BinMM.gene = calcGeneScore(SingleSNPKernel = BinMM.snp, standardize = FALSE, useWeights = FALSE)
      } else {
        IBS.gene = calcGeneScore(SingleSNPKernel = IBS.snp, standardize = TRUE, useWeights = FALSE)
        Incomp.gene = calcGeneScore(SingleSNPKernel = Incomp.snp, standardize = TRUE, useWeights = FALSE)
        AMS.gene = calcGeneScore(SingleSNPKernel = AMS.snp, standardize = TRUE, useWeights = FALSE)
        BinMM.gene = calcGeneScore(SingleSNPKernel = BinMM.snp, standardize = TRUE, useWeights = FALSE)
      }
    } else {
      IBS.gene = calcGeneScore(SingleSNPKernel = IBS.snp, standardize = FALSE, useWeights = TRUE, scoreWeights)
      Incomp.gene = calcGeneScore(SingleSNPKernel = Incomp.snp, standardize = FALSE, useWeights = TRUE, scoreWeights)
      AMS.gene = calcGeneScore(SingleSNPKernel = AMS.snp, standardize = FALSE, useWeights = TRUE, scoreWeights)
      BinMM.gene = calcGeneScore(SingleSNPKernel = BinMM.snp, standardize = FALSE, useWeights = TRUE, scoreWeights)
    }
    
    ###############################################
    # Calculate UR and US scores, and Q values
    # Need separate US for each score type
    ############################################### 
    #For binary phenotypes
    if(phenosBinary == TRUE){
      #Recipient genotype UR values:
      UR = CalcUScore(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = TRUE, RGenoData = RGenosMat, Phenos = phenos, BinPhenos = TRUE)
      #US for IBS Score
      US_IBS = CalcUScore(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = IBS.gene, Phenos = phenos, BinPhenos = TRUE)
      #US for Incomp Score
      US_Incomp = CalcUScore(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = Incomp.gene, Phenos = phenos, BinPhenos = TRUE)
      #US for AMS Score
      US_AMS = CalcUScore(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = AMS.gene, Phenos = phenos, BinPhenos = TRUE)
      #US for Bin MM Score
      US_BinMM = CalcUScore(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = BinMM.gene, Phenos = phenos, BinPhenos = TRUE)
      
      #Need to calculate Q values as well
      #R geno QR values
      QR = CalcQValues(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = TRUE, RGenoData = RGenosMat, Phenos = phenos, BinPhenos = TRUE)
      #QS for IBS Score
      QS_IBS = CalcQValues(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = IBS.gene, Phenos = phenos, BinPhenos = TRUE)
      #QS for Incomp Score
      QS_Incomp = CalcQValues(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = Incomp.gene, Phenos = phenos, BinPhenos = TRUE)
      #QS for AMS score
      QS_AMS = CalcQValues(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = AMS.gene, Phenos = phenos, BinPhenos = TRUE)
      #QS for Bin MM Score
      QS_BinMM = CalcQValues(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = BinMM.gene, Phenos = phenos, BinPhenos = TRUE)
    } 
    #for continuour phenotypes
    else {
      #Recipient genotype UR values:
      UR = CalcUScore(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = TRUE, RGenoData = RGenosMat, Phenos = phenos, BinPhenos = FALSE)
      #US for IBS Score
      US_IBS = CalcUScore(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = IBS.gene, Phenos = phenos, BinPhenos = FALSE)
      #US for Incomp Score
      US_Incomp = CalcUScore(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = Incomp.gene, Phenos = phenos, BinPhenos = FALSE)
      #US for AMS Score
      US_AMS = CalcUScore(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = AMS.gene, Phenos = phenos, BinPhenos = FALSE)
      #US for Bin MM Score
      US_BinMM = CalcUScore(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = BinMM.gene, Phenos = phenos, BinPhenos = FALSE)
      
      #Need to calculate Q values as well
      #R geno QR values
      QR = CalcQValues(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = TRUE, RGenoData = RGenosMat, Phenos = phenos, BinPhenos = FALSE)
      #QS for IBS Score
      QS_IBS = CalcQValues(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = IBS.gene, Phenos = phenos, BinPhenos = FALSE)
      #QS for Incomp Score
      QS_Incomp = CalcQValues(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = Incomp.gene, Phenos = phenos, BinPhenos = FALSE)
      #QS for AMS score
      QS_AMS = CalcQValues(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = AMS.gene, Phenos = phenos, BinPhenos = FALSE)
      #QS for Bin MM Score
      QS_BinMM = CalcQValues(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = BinMM.gene, Phenos = phenos, BinPhenos = FALSE)
    }
    
    ###############################################
    # Create combined Qs (QR, QS)
    ############################################### 
    #R geno and IBS score
    Q_IBS = cbind(QR, QS_IBS)
    #R geno and Incomp Score
    Q_Incomp = cbind(QR, QS_Incomp)
    #R geno and AMS score
    Q_AMS = cbind(QR, QS_AMS)
    #R geno and Bin MM score
    Q_BinMM = cbind(QR, QS_BinMM)
    
    ###############################################
    # Calculate full variance
    # Each set of Qs will have a variance calculation
    ############################################### 
    #R genos + IBS Score
    OrgVar_Q_IBS = CalcVariance(SampleSize = numPairs, QValues = Q_IBS)
    #R genos + Incomp Score
    OrgVar_Q_Incomp = CalcVariance(SampleSize = numPairs, QValues = Q_Incomp)
    #R genos + AMS Score
    OrgVar_Q_AMS = CalcVariance(SampleSize = numPairs, QValues = Q_AMS)
    #R genos + Bin MM Score
    OrgVar_Q_BinMM = CalcVariance(SampleSize = numPairs, QValues = Q_BinMM)
    
    ###############################################
    # Calculate final statistics and p-values
    ############################################### 
    #R geno and IBS score
    Stat_IBS = CalcStatisticPVal(SampleSize = numPairs, Variance = OrgVar_Q_IBS, UscoresR = UR, UscoreS = US_IBS, s = s)
    #R geno and Incomp score
    Stat_Incomp = CalcStatisticPVal(SampleSize = numPairs, Variance = OrgVar_Q_Incomp, UscoresR = UR, UscoreS = US_Incomp, s = s)
    #R geno and AMS Score
    Stat_AMS = CalcStatisticPVal(SampleSize = numPairs, Variance = OrgVar_Q_AMS, UscoresR = UR, UscoreS = US_AMS, s = s)
    #R geno and Bin MM score
    Stat_BinMM = CalcStatisticPVal(SampleSize = numPairs, Variance = OrgVar_Q_BinMM, UscoresR = UR, UscoreS = US_BinMM, s = s)
    
    #fill columns in order
    ##IBS, Incomp, AMS, Bin MM , 
    statsAndPVals[1,] = Stat_IBS
    statsAndPVals[2,] = Stat_Incomp
    statsAndPVals[3,] = Stat_AMS
    statsAndPVals[4,] = Stat_BinMM
    
    #labeling
    colnames(statsAndPVals) = c("Multi-marker Score", "p-value")
    rownames(statsAndPVals) = c("IBS Score", "Incomp Score", "AMS Score", "Bin MM Score")
    
    write.csv(statsAndPVals, file = paste0(outFile))
  } else {
    #define matrix to hold all Stats and Pvalues
    statsAndPVals = matrix(NA, nrow = 1, ncol = 2)
    
    #Read in the data
    #determine what file type
    fileName = read.table(text = datafile, sep = ".", as.is = TRUE)$V2
    if(fileName == "csv"){ #if it's a csv file, read in with read.csv
      geneticData = read.csv(datafile, header = TRUE, row.names = 1)
    } else { #otherwise should be good with a read.table
      geneticData = read.table(datafile, header = TRUE, row.names = 1)
    }
    
    geneticData = na.omit(geneticData)
    
    #if only 1 SNP, need to add in a column so it doesn't screw everything up
    if(dim(geneticData)[2] == 1){
      geneticData$Test = 0
    }
    ###############################################
    #Determine D/R pairs and match
    ###############################################  
    #read in the file with the donor/recipient pair information
    pairings = read.table(matchedPairs)
    #make sure column names are specified
    colnames(pairings) = c("R_Id", "D_Id")
    #assign each pair a number
    pairings$PairNumber = seq(1:length(pairings$R_Id))
    #split the original data into donor and recipient subsets
    donorData = geneticData[rownames(geneticData) %in% pairings$D_Id,]
    recipData = geneticData[rownames(geneticData) %in% pairings$R_Id,]
    #add an id column based on row names for merging
    donorData$D_Id = rownames(donorData)
    recipData$R_Id = rownames(recipData)
    #now merge the genotype data with the ids so that the genotype data has the pair numbers
    donorGenotypes_numbered = merge(donorData, pairings, by.x = "D_Id", by.y = "D_Id")
    recipGenotypes_numbered = merge(recipData, pairings, by.x = "R_Id", by.y = "R_Id")
    #make everything a data frame 
    donorGenotypes_numbered.df = as.data.frame(donorGenotypes_numbered)
    recipGenotypes_numbered.df = as.data.frame(recipGenotypes_numbered)
    #subset to only include Pairs with both D and R info
    donorGenotypes_numbered.df.subset = donorGenotypes_numbered.df[(donorGenotypes_numbered.df$PairNumber %in% recipGenotypes_numbered.df$PairNumber),]
    recipGenotypes_numbered.df.subset = recipGenotypes_numbered.df[(recipGenotypes_numbered.df$PairNumber %in% donorGenotypes_numbered.df$PairNumber),]
    #order by Pair number
    donorGenotypes_numbered.df.ordered = donorGenotypes_numbered.df.subset[order(donorGenotypes_numbered.df.subset$PairNumber),]
    recipGenotypes_numbered.df.ordered = recipGenotypes_numbered.df.subset[order(recipGenotypes_numbered.df.subset$PairNumber),]
    #make rownames into R_Id or D_Id
    rownames(donorGenotypes_numbered.df.ordered) = donorGenotypes_numbered.df.ordered$D_Id
    rownames(recipGenotypes_numbered.df.ordered) = recipGenotypes_numbered.df.ordered$R_Id
    #remove non-SNP columns
    DGenosData = donorGenotypes_numbered.df.ordered[,2:(ncol(donorGenotypes_numbered.df.ordered)-2)]
    RGenosData = recipGenotypes_numbered.df.ordered[,2:(ncol(recipGenotypes_numbered.df.ordered)-2)]
    
    ###############################################
    # Read in the covariate/phenotype file
    ############################################### 
    covariates = read.table(covPhenoFile, header = TRUE)
    #remove nas
    covariates = na.omit(covariates)
    #pull only those Rs from the R geno matrix
    covariates = covariates[covariates[,1] %in% rownames(RGenosData),]
    #set row names to be Ids from 1st row
    rownames(covariates) = covariates[,1]
    #remove the ids column
    covariates = covariates[,-1]
    
    if(dim(RGenosData)[1] > dim(covariates)[1]){
      RGenosData = RGenosData[rownames(RGenosData) %in% rownames(covariates),]
      DGenosData = DGenosData[which(rownames(RGenosData) %in% rownames(covariates)),]
    } else if (dim(RGenosData)[1] < dim(covariates)[1]){
      covariates = covariates[rownames(covariates) %in% rownames(RGenosData),]
    }
    
    #make matrices
    #unlist and make into matrices
    RGenosMat = matrix(unlist(RGenosData), ncol = ncol(RGenosData), byrow = F)
    DGenosMat = matrix(unlist(DGenosData), ncol = ncol(DGenosData), byrow = F)
    #get rid of the Test column if it exists
    if(("Test" %in% colnames(RGenosData)) == TRUE){
      RGenosMat = matrix(RGenosMat[,1], ncol = 1)
      DGenosMat = matrix(DGenosMat[,1], ncol = 1)
    }
    #keep number of pairs for later use
    numPairs = nrow(RGenosMat)
    numSNPs = ncol(RGenosMat)
    
    #output number of pairs
    strings = str_split(datafile, ".csv", simplify = TRUE)
    # filenameNumPairs = paste0(strings[1],"_numPairs.txt")
    # write.table(numPairs, filenameNumPairs)
    
    #write out number of SNPs
    filenameNumSNPs = paste0(strings[1],"_numSNPs.txt")
    write.table(numSNPs, filenameNumSNPs)
    
    #make a matrix
    #pull the phenotypes from the covariates file
    phenos = as.matrix(covariates[,1],nrow=numPairs, ncol = 1)
    #remove the phenotypes column from the covariates
    CovData.list = covariates[,-1]
    CovData = matrix(unlist(CovData.list), ncol = ncol(CovData.list), byrow = F)
    
    ###############################################
    # Calculate Gene-based Scores
    ###############################################   
    # First: Calculate single snp scores
    # Need to calculate score based on score variable that was input
    if(score == "IBS"){
      Score.snp = calcIBSMismatch(RGenosMat, DGenosMat)
    } else if(score == "Incomp"){
      Score.snp = calcIncompatibilityScore(RGenosMat, DGenosMat)
    } else if(score == "AMS"){
      Score.snp = calcAMS(RGenosMat, DGenosMat)
    } else if(score == "Bin MM"){
      Score.snp = calcBinaryMM(RGenosMat, DGenosMat)
    } else {
      print("Please enter a valid score type from the following: 'IBS', 'Incomp', 'AMS', 'Bin MM'.")
    }
    
    # Then: Calculate gene based scores
    ##Check to see if weights are used for scores
    if(weightedScores == FALSE){
      #check to see if scores should be standardized (can't be weighted and standardized)
      if(standardizeScores == FALSE){
        Score.gene = calcGeneScore(SingleSNPKernel = Score.snp, standardize = FALSE, useWeights = FALSE)
      } else {
        Score.gene = calcGeneScore(SingleSNPKernel = Score.snp, standardize = TRUE, useWeights = FALSE)
      } 
    } else {
      Score.gene = calcGeneScore(SingleSNPKernel = Score.snp, standardize = FALSE, useWeights = TRUE, scoreWeights)
    }
    
    ###############################################
    # Calculate UR and US scores, and Q values
    # Need separate US for each score type
    ############################################### 
    #For binary phenotypes
    if(phenosBinary == TRUE){
      #Recipient genotype UR values:
      UR = CalcUScore(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = TRUE, RGenoData = RGenosMat, Phenos = phenos, BinPhenos = TRUE)
      #US for Score
      US_Score = CalcUScore(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = Score.gene, Phenos = phenos, BinPhenos = TRUE)
      
      #Need to calculate Q values as well
      #R geno QR values
      QR = CalcQValues(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = TRUE, RGenoData = RGenosMat, Phenos = phenos, BinPhenos = TRUE)
      #QS for Score
      QS_Score = CalcQValues(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = Score.gene, Phenos = phenos, BinPhenos = TRUE)
    } 
    #for continuour phenotypes
    else {
      #Recipient genotype UR values:
      UR = CalcUScore(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = TRUE, RGenoData = RGenosMat, Phenos = phenos, BinPhenos = FALSE)
      #US for IBS Score
      US_Score = CalcUScore(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = Score.gene, Phenos = phenos, BinPhenos = FALSE)
      
      #Need to calculate Q values as well
      #R geno QR values
      QR = CalcQValues(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = TRUE, RGenoData = RGenosMat, Phenos = phenos, BinPhenos = FALSE)
      #QS for IBS Score
      QS_Score = CalcQValues(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = Score.gene, Phenos = phenos, BinPhenos = FALSE)
    }
    
    ###############################################
    # Create combined Qs (QR, QS)
    ############################################### 
    #R geno and score
    Q_Score = cbind(QR, QS_Score)
    
    ###############################################
    # Calculate full variance
    # Each set of Qs will have a variance calculation
    ############################################### 
    #R genos + IBS Score
    OrgVar_Q_Score = CalcVariance(SampleSize = numPairs, QValues = Q_Score)
    
    ###############################################
    # Calculate final statistics and p-values
    ############################################### 
    #R geno and score
    Stat_Score = CalcStatisticPVal(SampleSize = numPairs, Variance = OrgVar_Q_Score, UscoresR = UR, UscoreS = US_Score, s = s)
    
    #fill columns in order
    statsAndPVals[1,] = Stat_Score
    
    #labeling
    colnames(statsAndPVals) = c("Multi-marker Score", "p-value")
    rownames(statsAndPVals) = c(paste0(score," Score"))
    
    write.csv(statsAndPVals, file = paste0(outFile)) 
  }
}

#function to run SKAT analyses on real data
realDataPipeline_SKAT = function(datafile = "genotypeData.csv", matchedPairs = "pairedIDs.txt", covPhenoFile = "covariatesAndPhenotypes.txt", dataPath = path,
                                 standardizeScores = FALSE, weightedScores = FALSE, scoreWeights, phenosBinary = TRUE, kernel, score = "All", outFile = "out.csv"){
  #Inputs: 
  ## datafile: the name of the file containing the genotype data
  ## matchedPairs: a text file containing the information on which donor id is matched to which recipient id
  ## covPhenoFile: a text file containing additional covariate information, the first column of the covariate file should be the Recipient ids,
  ##                the second column should be the phenotype of interest, additional columns are other covariates (all numeric)
  ## dataPath: the path to where the data is located
  ## standardizeScores: TRUE or FALSE for whether scores should be standardized
  ## weightedScores: TRUE or FALSE for whether scores should be weighted
  ## scoreWeights: Weights for each SNP if weightedScores = TRUE
  ## phenosBinary: TRUE or FALSE for whether the phenotype of interest is binary, if FALSE phenotype must be continuous
  ## kernel: the SKAT kernel to use, for simulations I used unweighted linear or unweighted IBS
  ## score: the gene-based score used in the model, either 'IBS', 'Incomp', 'AMS', 'Bin MM' or 'All' if you want to do all 4
  ## outFile: name of file to write output to
  
  #Output: Output file containing the score values and p-values for the R Geno/IBS Score, R Geno/Incomp Score, R Geno/AMS, and R Geno/Bin MM Score
  
  #set working directory to where data is
  setwd(dataPath)
  
  #source the needed functions
  source(paste0(dataPath,"/RealDataSourceFunctions.R"))
  library(stringr)
  library(SKAT)
  
  if(score == "All"){
    #define matrix to hold all Stats and Pvalues
    statsAndPVals = matrix(NA, nrow = 4, ncol = 2)
    
    #Read in the data
    #determine what file type
    fileName = read.table(text = datafile, sep = ".", as.is = TRUE)$V2
    if(fileName == "csv"){ #if it's a csv file, read in with read.csv
      geneticData = read.csv(datafile, header = TRUE, row.names = 1)
    } else { #otherwise should be good with a read.table
      geneticData = read.table(datafile, header = TRUE, row.names = 1)
    }
    
    geneticData = na.omit(geneticData)
    
    #if only 1 SNP, need to add in a column so it doesn't screw everything up
    if(dim(geneticData)[2] == 1){
      geneticData$Test = 0
    }
    ###############################################
    #Determine D/R pairs and match
    ###############################################  
    #read in the file with the donor/recipient pair information
    pairings = read.table(matchedPairs)
    #make sure column names are specified
    colnames(pairings) = c("R_Id", "D_Id")
    #assign each pair a number
    pairings$PairNumber = seq(1:length(pairings$R_Id))
    #split the original data into donor and recipient subsets
    donorData = geneticData[rownames(geneticData) %in% pairings$D_Id,]
    recipData = geneticData[rownames(geneticData) %in% pairings$R_Id,]
    #add an id column based on row names for merging
    donorData$D_Id = rownames(donorData)
    recipData$R_Id = rownames(recipData)
    #now merge the genotype data with the ids so that the genotype data has the pair numbers
    donorGenotypes_numbered = merge(donorData, pairings, by.x = "D_Id", by.y = "D_Id")
    recipGenotypes_numbered = merge(recipData, pairings, by.x = "R_Id", by.y = "R_Id")
    #make everything a data frame 
    donorGenotypes_numbered.df = as.data.frame(donorGenotypes_numbered)
    recipGenotypes_numbered.df = as.data.frame(recipGenotypes_numbered)
    #subset to only include Pairs with both D and R info
    donorGenotypes_numbered.df.subset = donorGenotypes_numbered.df[(donorGenotypes_numbered.df$PairNumber %in% recipGenotypes_numbered.df$PairNumber),]
    recipGenotypes_numbered.df.subset = recipGenotypes_numbered.df[(recipGenotypes_numbered.df$PairNumber %in% donorGenotypes_numbered.df$PairNumber),]
    #order by Pair number
    donorGenotypes_numbered.df.ordered = donorGenotypes_numbered.df.subset[order(donorGenotypes_numbered.df.subset$PairNumber),]
    recipGenotypes_numbered.df.ordered = recipGenotypes_numbered.df.subset[order(recipGenotypes_numbered.df.subset$PairNumber),]
    #make rownames into R_Id or D_Id
    rownames(donorGenotypes_numbered.df.ordered) = donorGenotypes_numbered.df.ordered$D_Id
    rownames(recipGenotypes_numbered.df.ordered) = recipGenotypes_numbered.df.ordered$R_Id
    #remove non-SNP columns
    DGenosData = donorGenotypes_numbered.df.ordered[,2:(ncol(donorGenotypes_numbered.df.ordered)-2)]
    RGenosData = recipGenotypes_numbered.df.ordered[,2:(ncol(recipGenotypes_numbered.df.ordered)-2)]
    
    ###############################################
    # Read in the covariate/phenotype file
    ############################################### 
    covariates = read.table(covPhenoFile, header = TRUE)
    #remove nas
    covariates = na.omit(covariates)
    #pull only those Rs from the R geno matrix
    covariates = covariates[covariates[,1] %in% rownames(RGenosData),]
    #set row names to be Ids from 1st row
    rownames(covariates) = covariates[,1]
    #remove the ids column
    covariates = covariates[,-1]
    
    if(dim(RGenosData)[1] > dim(covariates)[1]){
      RGenosData = RGenosData[rownames(RGenosData) %in% rownames(covariates),]
      DGenosData = DGenosData[which(rownames(RGenosData) %in% rownames(covariates)),]
    } else if (dim(RGenosData)[1] < dim(covariates)[1]){
      covariates = covariates[rownames(covariates) %in% rownames(RGenosData),]
    }
    
    #make matrices
    #unlist and make into matrices
    RGenosMat = matrix(unlist(RGenosData), ncol = ncol(RGenosData), byrow = F)
    DGenosMat = matrix(unlist(DGenosData), ncol = ncol(DGenosData), byrow = F)
    #get rid of the Test column if it exists
    if(("Test" %in% colnames(RGenosData)) == TRUE){
      RGenosMat = matrix(RGenosMat[,1], ncol = 1)
      DGenosMat = matrix(DGenosMat[,1], ncol = 1)
    }
    #keep number of pairs for later use
    numPairs = nrow(RGenosMat)
    numSNPs = ncol(RGenosMat)
    
    #make a matrix
    #pull the phenotypes from the covariates file
    phenos = as.matrix(covariates[,1],nrow=numPairs, ncol = 1)
    #remove the phenotypes column from the covariates
    CovData.list = covariates[,-1]
    CovData = matrix(unlist(CovData.list), ncol = ncol(CovData.list), byrow = F)
    
    ###############################################
    # Calculate Gene-based Scores
    ###############################################   
    # First: Calculate single snp scores
    # Need to calculate score based on score variable that was input
    IBS.snp = calcIBSMismatch(RGenosMat, DGenosMat)
    Incomp.snp = calcIncompatibilityScore(RGenosMat, DGenosMat)
    AMS.snp = calcAMS(RGenosMat, DGenosMat)
    BinMM.snp = calcBinaryMM(RGenosMat, DGenosMat)
    
    # Then: Calculate gene based scores
    #check to see if weights are used for scores
    if(weightedScores == FALSE){
      #check to see if scores should be standardized (can't be weighted and standardized)
      if(standardizeScores == FALSE){
        IBS.gene = calcGeneScore(SingleSNPKernel = IBS.snp, standardize = FALSE, useWeights = FALSE)
        Incomp.gene = calcGeneScore(SingleSNPKernel = Incomp.snp, standardize = FALSE, useWeights = FALSE)
        AMS.gene = calcGeneScore(SingleSNPKernel = AMS.snp, standardize = FALSE, useWeights = FALSE)
        BinMM.gene = calcGeneScore(SingleSNPKernel = BinMM.snp, standardize = FALSE, useWeights = FALSE)
      } else {
        IBS.gene = calcGeneScore(SingleSNPKernel = IBS.snp, standardize = TRUE, useWeights = FALSE)
        Incomp.gene = calcGeneScore(SingleSNPKernel = Incomp.snp, standardize = TRUE, useWeights = FALSE)
        AMS.gene = calcGeneScore(SingleSNPKernel = AMS.snp, standardize = TRUE, useWeights = FALSE)
        BinMM.gene = calcGeneScore(SingleSNPKernel = BinMM.snp, standardize = TRUE, useWeights = FALSE)
      }
    } else {
      IBS.gene = calcGeneScore(SingleSNPKernel = IBS.snp, standardize = FALSE, useWeights = TRUE, scoreWeights)
      Incomp.gene = calcGeneScore(SingleSNPKernel = Incomp.snp, standardize = FALSE, useWeights = TRUE, scoreWeights)
      AMS.gene = calcGeneScore(SingleSNPKernel = AMS.snp, standardize = FALSE, useWeights = TRUE, scoreWeights)
      BinMM.gene = calcGeneScore(SingleSNPKernel = BinMM.snp, standardize = FALSE, useWeights = TRUE, scoreWeights)
    }
    
    #Combine R geno and Scores into 4 separate datasets, size: N x (m+1)
    RGeno.IBS.gene = cbind(RGenosMat, IBS.gene)
    RGeno.Incomp.gene = cbind(RGenosMat, Incomp.gene)
    RGeno.AMS.gene = cbind(RGenosMat, AMS.gene)
    RGeno.BinMM.gene = cbind(RGenosMat, BinMM.gene)
    
    ## Generate SKAT Null Models
    # formulas will be Y ~ covariates for continuous and dichotomous Y
    if(phenosBinary == TRUE){
      obj=SKAT_Null_Model(phenos~CovData, out_type="D", Adjustment = FALSE)
    } else {
      obj=SKAT_Null_Model(phenos~CovData, out_type="C", Adjustment = FALSE)
    }
    
    #perform SKAT analysis
    Stat_IBS = SKAT(RGeno.IBS.gene, obj, kernel = kernel, is_check_genotype = FALSE)
    Stat_Incomp = SKAT(RGeno.Incomp.gene, obj, kernel = kernel, is_check_genotype = FALSE)
    Stat_AMS = SKAT(RGeno.AMS.gene, obj, kernel = kernel, is_check_genotype = FALSE)
    Stat_BinMM = SKAT(RGeno.BinMM.gene, obj, kernel = kernel, is_check_genotype = FALSE)
    
    #pull stat and p value
    statsAndPVals[1,] = c(Stat_IBS$Q, Stat_IBS$p.value)
    statsAndPVals[2,] = c(Stat_Incomp$Q, Stat_Incomp$p.value)
    statsAndPVals[3,] = c(Stat_AMS$Q, Stat_AMS$p.value)
    statsAndPVals[4,] = c(Stat_BinMM$Q, Stat_BinMM$p.value)
    
    #labeling
    colnames(statsAndPVals) = c("SKAT Q Stat", "p-value")
    rownames(statsAndPVals) = c("IBS Score", "Incomp Score", "AMS Score", "Binary Mismatch Score")
    
    write.csv(statsAndPVals, file = paste0(outFile))
  } else {
    #define matrix to hold all Stats and Pvalues
    statsAndPVals = matrix(NA, nrow = 1, ncol = 2)
    
    #Read in the data
    #determine what file type
    fileName = read.table(text = datafile, sep = ".", as.is = TRUE)$V2
    if(fileName == "csv"){ #if it's a csv file, read in with read.csv
      geneticData = read.csv(datafile, header = TRUE, row.names = 1)
    } else { #otherwise should be good with a read.table
      geneticData = read.table(datafile, header = TRUE, row.names = 1)
    }
    
    geneticData = na.omit(geneticData)
    
    #if only 1 SNP, need to add in a column so it doesn't screw everything up
    if(dim(geneticData)[2] == 1){
      geneticData$Test = 0
    }
    ###############################################
    #Determine D/R pairs and match
    ###############################################  
    #read in the file with the donor/recipient pair information
    pairings = read.table(matchedPairs)
    #make sure column names are specified
    colnames(pairings) = c("R_Id", "D_Id")
    #assign each pair a number
    pairings$PairNumber = seq(1:length(pairings$R_Id))
    #split the original data into donor and recipient subsets
    donorData = geneticData[rownames(geneticData) %in% pairings$D_Id,]
    recipData = geneticData[rownames(geneticData) %in% pairings$R_Id,]
    #add an id column based on row names for merging
    donorData$D_Id = rownames(donorData)
    recipData$R_Id = rownames(recipData)
    #now merge the genotype data with the ids so that the genotype data has the pair numbers
    donorGenotypes_numbered = merge(donorData, pairings, by.x = "D_Id", by.y = "D_Id")
    recipGenotypes_numbered = merge(recipData, pairings, by.x = "R_Id", by.y = "R_Id")
    #make everything a data frame 
    donorGenotypes_numbered.df = as.data.frame(donorGenotypes_numbered)
    recipGenotypes_numbered.df = as.data.frame(recipGenotypes_numbered)
    #subset to only include Pairs with both D and R info
    donorGenotypes_numbered.df.subset = donorGenotypes_numbered.df[(donorGenotypes_numbered.df$PairNumber %in% recipGenotypes_numbered.df$PairNumber),]
    recipGenotypes_numbered.df.subset = recipGenotypes_numbered.df[(recipGenotypes_numbered.df$PairNumber %in% donorGenotypes_numbered.df$PairNumber),]
    #order by Pair number
    donorGenotypes_numbered.df.ordered = donorGenotypes_numbered.df.subset[order(donorGenotypes_numbered.df.subset$PairNumber),]
    recipGenotypes_numbered.df.ordered = recipGenotypes_numbered.df.subset[order(recipGenotypes_numbered.df.subset$PairNumber),]
    #make rownames into R_Id or D_Id
    rownames(donorGenotypes_numbered.df.ordered) = donorGenotypes_numbered.df.ordered$D_Id
    rownames(recipGenotypes_numbered.df.ordered) = recipGenotypes_numbered.df.ordered$R_Id
    #remove non-SNP columns
    DGenosData = donorGenotypes_numbered.df.ordered[,2:(ncol(donorGenotypes_numbered.df.ordered)-2)]
    RGenosData = recipGenotypes_numbered.df.ordered[,2:(ncol(recipGenotypes_numbered.df.ordered)-2)]
    
    ###############################################
    # Read in the covariate/phenotype file
    ############################################### 
    covariates = read.table(covPhenoFile, header = TRUE)
    #remove nas
    covariates = na.omit(covariates)
    #pull only those Rs from the R geno matrix
    covariates = covariates[covariates[,1] %in% rownames(RGenosData),]
    #set row names to be Ids from 1st row
    rownames(covariates) = covariates[,1]
    #remove the ids column
    covariates = covariates[,-1]
    
    if(dim(RGenosData)[1] > dim(covariates)[1]){
      RGenosData = RGenosData[rownames(RGenosData) %in% rownames(covariates),]
      DGenosData = DGenosData[which(rownames(RGenosData) %in% rownames(covariates)),]
    } else if (dim(RGenosData)[1] < dim(covariates)[1]){
      covariates = covariates[rownames(covariates) %in% rownames(RGenosData),]
    }
    
    #make matrices
    #unlist and make into matrices
    RGenosMat = matrix(unlist(RGenosData), ncol = ncol(RGenosData), byrow = F)
    DGenosMat = matrix(unlist(DGenosData), ncol = ncol(DGenosData), byrow = F)
    #get rid of the Test column if it exists
    if(("Test" %in% colnames(RGenosData)) == TRUE){
      RGenosMat = matrix(RGenosMat[,1], ncol = 1)
      DGenosMat = matrix(DGenosMat[,1], ncol = 1)
    }
    #keep number of pairs for later use
    numPairs = nrow(RGenosMat)
    numSNPs = ncol(RGenosMat)
    
    #make a matrix
    #pull the phenotypes from the covariates file
    phenos = as.matrix(covariates[,1],nrow=numPairs, ncol = 1)
    #remove the phenotypes column from the covariates
    CovData.list = covariates[,-1]
    CovData = matrix(unlist(CovData.list), ncol = ncol(CovData.list), byrow = F)
    
    ###############################################
    # Calculate Gene-based Scores
    ###############################################   
    # First: Calculate single snp scores
    # Need to calculate score based on score variable that was input
    if(score == "IBS"){
      Score.snp = calcIBSMismatch(RGenosMat, DGenosMat)
    } else if(score == "Incomp"){
      Score.snp = calcIncompatibilityScore(RGenosMat, DGenosMat)
    } else if(score == "AMS"){
      Score.snp = calcAMS(RGenosMat, DGenosMat)
    } else if(score == "Bin MM"){
      Score.snp = calcBinaryMM(RGenosMat, DGenosMat)
    } else {
      print("Please enter a valid score type from the following: 'IBS', 'Incomp', 'AMS', 'Bin MM'.")
    }
    
    # Then: Calculate gene based scores
    ##Check to see if weights are used for scores
    if(weightedScores == FALSE){
      #check to see if scores should be standardized (can't be weighted and standardized)
      if(standardizeScores == FALSE){
        Score.gene = calcGeneScore(SingleSNPKernel = Score.snp, standardize = FALSE, useWeights = FALSE)
      } else {
        Score.gene = calcGeneScore(SingleSNPKernel = Score.snp, standardize = TRUE, useWeights = FALSE)
      } 
    } else {
      Score.gene = calcGeneScore(SingleSNPKernel = Score.snp, standardize = FALSE, useWeights = TRUE, scoreWeights)
    }
    
    #combine the R geno and score into matrix of N x (m+1) size
    RGeno.Score.gene = cbind(RGenosData, Score.gene)
    
    ## Generate SKAT Null Models
    # formulas will be Y ~ covariates for continuous and dichotomous Y
    if(phenosBinary == TRUE){
      obj=SKAT_Null_Model(phenos~CovData, out_type="D")
    } else {
      obj=SKAT_Null_Model(phenos~CovData, out_type="C")
    }
    
    #perform SKAT analysis
    RGeno.Score.gene.mat = as.matrix(RGeno.Score.gene)
    Stat_SKAT = SKAT(RGeno.Score.gene.mat, obj, kernel = kernel, is_check_genotype = FALSE)
    
    #pull stat and p value
    statsAndPVals[1,] = c(Stat_SKAT$Q, Stat_SKAT$p.value)
    
    #labeling
    colnames(statsAndPVals) = c("SKAT Q Stat", "p-value")
    rownames(statsAndPVals) = c(paste0(score," Score"))
    
    write.csv(statsAndPVals, file = paste0(outFile))
  }
}

#function to run test for association of gene based score only (not the joint testing!)
realDataPipeline_GeneScoreTest = function(datafile = "genotypeData.csv", matchedPairs = "pairedIDs.txt", covPhenoFile = "covariatesAndPhenotypes.txt", dataPath = path,
                                          standardizeScores = FALSE, weightedScores = FALSE, scoreWeights, phenosBinary = TRUE, score = "All", fitIntercept = TRUE, unpen = c(), outFile = "out.csv"){
  #Inputs: 
  ## datafile: the name of the file containing the genotype data
  ## matchedPairs: a text file containing the information on which donor id is matched to which recipient id
  ## covPhenoFile: a text file containing additional covariate information, the first column of the covariate file should be the Recipient ids,
  ##                the second column should be the phenotype of interest, additional columns are other covariates (all numeric)
  ## dataPath: the path to where the data is located
  ## standardizeScores: TRUE or FALSE for whether scores should be standardized
  ## weightedScores: TRUE or FALSE for whether scores should be weighted
  ## scoreWeights: Weights for each SNP if weightedScores = TRUE
  ## phenosBinary: TRUE or FALSE for whether the phenotype of interest is binary, if FALSE phenotype must be continuous
  ## score: the gene-based score used in the model, either 'IBS', 'Incomp', 'AMS', 'Bin MM', or 'All' if you want to do all 4
  ## fitIntercept: TRUE or FALSE for if you want to force an intercept term for the model, TRUE means the model will be forced to have an intercept term
  ## unpen: A vector of the variables that you want to force into the model, default is to have the model force all covariates to be fit without penalty
  ##        so that only the R genotype SNPs will have penalized model coefficients
  ## outFile: name of file to write output to
  
  #Output: Output file containing the score values and p-values for the R Geno/IBS Score, R Geno/Incomp Score, R Geno/AMS, and R Geno/Bin MM Score
  
  #set working directory to where data is
  setwd(path)
  
  #source the needed functions
  source(paste0(path,"/RealDataSourceFunctions.R"))

  library(stringr)
  
  if(score == "All"){
    #define matrix to hold all Stats and Pvalues
    statsAndPVals = matrix(NA, nrow = 4, ncol = 4)
    
    #Read in the data
    #determine what file type
    fileName = read.table(text = datafile, sep = ".", as.is = TRUE)$V2
    if(fileName == "csv"){ #if it's a csv file, read in with read.csv
      geneticData = read.csv(datafile, header = TRUE, row.names = 1)
    } else { #otherwise should be good with a read.table
      geneticData = read.table(datafile, header = TRUE, row.names = 1)
    }
    
    geneticData = na.omit(geneticData)
    
    #if only 1 SNP, need to add in a column so it doesn't screw everything up
    if(dim(geneticData)[2] == 1){
      geneticData$Test = 0
    }
    ###############################################
    #Determine D/R pairs and match
    ###############################################  
    #read in the file with the donor/recipient pair information
    pairings = read.table(matchedPairs)
    #make sure column names are specified
    colnames(pairings) = c("R_Id", "D_Id")
    #assign each pair a number
    pairings$PairNumber = seq(1:length(pairings$R_Id))
    #split the original data into donor and recipient subsets
    donorData = geneticData[rownames(geneticData) %in% pairings$D_Id,]
    recipData = geneticData[rownames(geneticData) %in% pairings$R_Id,]
    #add an id column based on row names for merging
    donorData$D_Id = rownames(donorData)
    recipData$R_Id = rownames(recipData)
    #now merge the genotype data with the ids so that the genotype data has the pair numbers
    donorGenotypes_numbered = merge(donorData, pairings, by.x = "D_Id", by.y = "D_Id")
    recipGenotypes_numbered = merge(recipData, pairings, by.x = "R_Id", by.y = "R_Id")
    #make everything a data frame 
    donorGenotypes_numbered.df = as.data.frame(donorGenotypes_numbered)
    recipGenotypes_numbered.df = as.data.frame(recipGenotypes_numbered)
    #subset to only include Pairs with both D and R info
    donorGenotypes_numbered.df.subset = donorGenotypes_numbered.df[(donorGenotypes_numbered.df$PairNumber %in% recipGenotypes_numbered.df$PairNumber),]
    recipGenotypes_numbered.df.subset = recipGenotypes_numbered.df[(recipGenotypes_numbered.df$PairNumber %in% donorGenotypes_numbered.df$PairNumber),]
    #order by Pair number
    donorGenotypes_numbered.df.ordered = donorGenotypes_numbered.df.subset[order(donorGenotypes_numbered.df.subset$PairNumber),]
    recipGenotypes_numbered.df.ordered = recipGenotypes_numbered.df.subset[order(recipGenotypes_numbered.df.subset$PairNumber),]
    #make rownames into R_Id or D_Id
    rownames(donorGenotypes_numbered.df.ordered) = donorGenotypes_numbered.df.ordered$D_Id
    rownames(recipGenotypes_numbered.df.ordered) = recipGenotypes_numbered.df.ordered$R_Id
    #remove non-SNP columns
    DGenosData = donorGenotypes_numbered.df.ordered[,2:(ncol(donorGenotypes_numbered.df.ordered)-2)]
    RGenosData = recipGenotypes_numbered.df.ordered[,2:(ncol(recipGenotypes_numbered.df.ordered)-2)]
    
    ###############################################
    # Read in the covariate/phenotype file
    ############################################### 
    covariates = read.table(covPhenoFile, header = TRUE)
    #remove nas
    covariates = na.omit(covariates)
    #pull only those Rs from the R geno matrix
    covariates = covariates[covariates[,1] %in% rownames(RGenosData),]
    #set row names to be Ids from 1st row
    rownames(covariates) = covariates[,1]
    #remove the ids column
    covariates = covariates[,-1]
    
    if(dim(RGenosData)[1] > dim(covariates)[1]){
      RGenosData = RGenosData[rownames(RGenosData) %in% rownames(covariates),]
      DGenosData = DGenosData[which(rownames(RGenosData) %in% rownames(covariates)),]
    } else if (dim(RGenosData)[1] < dim(covariates)[1]){
      covariates = covariates[rownames(covariates) %in% rownames(RGenosData),]
    }
    
    #make matrices
    #unlist and make into matrices
    RGenosMat = matrix(unlist(RGenosData), ncol = ncol(RGenosData), byrow = F)
    DGenosMat = matrix(unlist(DGenosData), ncol = ncol(DGenosData), byrow = F)
    #get rid of the Test column if it exists
    if(("Test" %in% colnames(RGenosData)) == TRUE){
      RGenosMat = matrix(RGenosMat[,1], ncol = 1)
      DGenosMat = matrix(DGenosMat[,1], ncol = 1)
    }
    #keep number of pairs for later use
    numPairs = nrow(RGenosMat)
    numSNPs = ncol(RGenosMat)
    
    #make a matrix
    #pull the phenotypes from the covariates file
    phenos = as.matrix(covariates[,1],nrow=numPairs, ncol = 1)
    #remove the phenotypes column from the covariates
    CovData.list = covariates[,-1]
    CovData = matrix(unlist(CovData.list), ncol = ncol(CovData.list), byrow = F)
    
    ###############################################
    # Calculate Gene-based Scores
    ###############################################   
    # First: Calculate single snp scores
    # Need to calculate score based on score variable that was input
    IBS.snp = calcIBSMismatch(RGenosMat, DGenosMat)
    Incomp.snp = calcIncompatibilityScore(RGenosMat, DGenosMat)
    AMS.snp = calcAMS(RGenosMat, DGenosMat)
    BinMM.snp = calcBinaryMM(RGenosMat, DGenosMat)
    
    # Then: Calculate gene based scores
    #check to see if weights are used for scores
    if(weightedScores == FALSE){
      #check to see if scores should be standardized (can't be weighted and standardized)
      if(standardizeScores == FALSE){
        IBS.gene = calcGeneScore(SingleSNPKernel = IBS.snp, standardize = FALSE, useWeights = FALSE)
        Incomp.gene = calcGeneScore(SingleSNPKernel = Incomp.snp, standardize = FALSE, useWeights = FALSE)
        AMS.gene = calcGeneScore(SingleSNPKernel = AMS.snp, standardize = FALSE, useWeights = FALSE)
        BinMM.gene = calcGeneScore(SingleSNPKernel = BinMM.snp, standardize = FALSE, useWeights = FALSE)
      } else {
        IBS.gene = calcGeneScore(SingleSNPKernel = IBS.snp, standardize = TRUE, useWeights = FALSE)
        Incomp.gene = calcGeneScore(SingleSNPKernel = Incomp.snp, standardize = TRUE, useWeights = FALSE)
        AMS.gene = calcGeneScore(SingleSNPKernel = AMS.snp, standardize = TRUE, useWeights = FALSE)
        BinMM.gene = calcGeneScore(SingleSNPKernel = BinMM.snp, standardize = TRUE, useWeights = FALSE)
      }
    } else {
      IBS.gene = calcGeneScore(SingleSNPKernel = IBS.snp, standardize = FALSE, useWeights = TRUE, scoreWeights)
      Incomp.gene = calcGeneScore(SingleSNPKernel = Incomp.snp, standardize = FALSE, useWeights = TRUE, scoreWeights)
      AMS.gene = calcGeneScore(SingleSNPKernel = AMS.snp, standardize = FALSE, useWeights = TRUE, scoreWeights)
      BinMM.gene = calcGeneScore(SingleSNPKernel = BinMM.snp, standardize = FALSE, useWeights = TRUE, scoreWeights)
    }
    
    # define location of zero components
    # basically, this defines what the null hypothesis is
    # so for gamma = 0 null, we need the non-zero components to be for covariates
    # and the beta, and gamma components to be zero
    numCov = dim(CovData)[2]
    numSNPs = dim(RGenosMat)[2]
    
    #if there is a forced intercept
    if(fitIntercept == TRUE){
      intercept = matrix(1,nrow = dim(RGenosMat)[1], ncol = 1)
      N = c(rep(FALSE, numCov+1+numSNPs), TRUE)
      designMat.IBS = cbind(intercept, CovData, RGenosMat, IBS.gene)
      designMat.Incomp = cbind(intercept, CovData, RGenosMat, Incomp.gene)
      designMat.AMS = cbind(intercept, CovData, RGenosMat, AMS.gene)
      designMat.BinMM = cbind(intercept, CovData, RGenosMat, BinMM.gene)
    } else {
      #otherwise no intercept included
      N = c(rep(FALSE, numCov+numSNPs), TRUE)
      designMat.IBS = cbind(CovData, RGenosMat, IBS.gene)
      designMat.Incomp = cbind(CovData, RGenosMat, Incomp.gene)
      designMat.AMS = cbind(CovData, RGenosMat, AMS.gene)
      designMat.BinMM = cbind(CovData, RGenosMat, BinMM.gene)
    }
    
    #then combine the phenos with the design matrix as a list
    Model.IBS = list(X=designMat.IBS, Y=phenos)
    Model.Incomp = list(X=designMat.Incomp, Y=phenos)
    Model.AMS = list(X=designMat.AMS, Y=phenos)
    Model.BinMM = list(X=designMat.BinMM, Y=phenos)
    
    #determine what to not penalize in the model
    if(length(unpen) == 0){
      #if not specified, do not penalize the covariates
      #number of covariates in the model (with or without intercept) is the length of N, minus the number of SNPs and the gene score
      totalNumCovs =  length(N) - numSNPs - 1
      unpenVals = seq(1,totalNumCovs,1)
    } else {
      unpenVals = unpen
    }
    
    #if outcome is binary
    if(phenosBinary == TRUE){
      #source the needed functions
      source(paste0(path,"/Logistic_ADMM0.r"))
      
      # estimate the uncontrained estimator
      #IBS Score
      beta.unre.IBS <- cv.SCAD_ADMM_unre(X=Model.IBS$X, Y=Model.IBS$Y, N=N, beta0=rep(0, dim(Model.IBS$X)[2]), err=1e-4, tune="cv", unpen=unpenVals)
      indice.unre.IBS <- beta.unre.IBS!=0
      #Incomp Score
      beta.unre.Incomp <- cv.SCAD_ADMM_unre(X=Model.Incomp$X, Y=Model.Incomp$Y, N=N, beta0=rep(0, dim(Model.Incomp$X)[2]), err=1e-4, tune="cv", unpen=unpenVals)
      indice.unre.Incomp <- beta.unre.Incomp!=0
      #AMS Score
      beta.unre.AMS <- cv.SCAD_ADMM_unre(X=Model.AMS$X, Y=Model.AMS$Y, N=N, beta0=rep(0, dim(Model.AMS$X)[2]), err=1e-4, tune="cv", unpen=unpenVals)
      indice.unre.AMS <- beta.unre.AMS!=0
      #Bin MM Score
      beta.unre.BinMM <- cv.SCAD_ADMM_unre(X=Model.BinMM$X, Y=Model.BinMM$Y, N=N, beta0=rep(0, dim(Model.BinMM$X)[2]), err=1e-4, tune="cv", unpen=unpenVals)
      indice.unre.BinMM <- beta.unre.BinMM!=0
      
      # estimate the constrained estimator (only need this for score test)
      #IBS Score
      beta.re.IBS <- cv.SCAD_ADMM_re(X=Model.IBS$X, Y=Model.IBS$Y, N=N, beta0=rep(0, dim(Model.IBS$X)[2]), err=1e-4, tune="cv", unpen=unpenVals)
      indice.re.IBS <- beta.re.IBS!=0
      #Incomp Score
      beta.re.Incomp <- cv.SCAD_ADMM_re(X=Model.Incomp$X, Y=Model.Incomp$Y, N=N, beta0=rep(0, dim(Model.Incomp$X)[2]), err=1e-4, tune="cv", unpen=unpenVals)
      indice.re.Incomp <- beta.re.Incomp!=0
      #AMS Score
      beta.re.AMS <- cv.SCAD_ADMM_re(X=Model.AMS$X, Y=Model.AMS$Y, N=N, beta0=rep(0, dim(Model.AMS$X)[2]), err=1e-4, tune="cv", unpen=unpenVals)
      indice.re.AMS <- beta.re.AMS!=0
      #Bin MM  Score
      beta.re.BinMM <- cv.SCAD_ADMM_re(X=Model.BinMM$X, Y=Model.BinMM$Y, N=N, beta0=rep(0, dim(Model.BinMM$X)[2]), err=1e-4, tune="cv", unpen=unpenVals)
      indice.re.BinMM <- beta.re.BinMM!=0
      
      #IBS probs
      pi.unre.IBS <- logit(Model.IBS$X%*%beta.unre.IBS)
      pi.re.IBS <- logit(Model.IBS$X%*%beta.re.IBS)
      #Incomp probs
      pi.unre.Incomp <- logit(Model.Incomp$X%*%beta.unre.Incomp)
      pi.re.Incomp <- logit(Model.Incomp$X%*%beta.re.Incomp)
      #AMS probs
      pi.unre.AMS <- logit(Model.AMS$X%*%beta.unre.AMS)
      pi.re.AMS <- logit(Model.AMS$X%*%beta.re.AMS)
      #Bin MM probs
      pi.unre.BinMM <- logit(Model.BinMM$X%*%beta.unre.BinMM)
      pi.re.BinMM <- logit(Model.BinMM$X%*%beta.re.BinMM)
      
      # construct the score statistics
      # IBS Score
      eps.IBS <- Model.IBS$Y-pi.re.IBS #Y - e(Y)
      #this is the X^T times Y-E(Y)
      Xeps.IBS <- crossprod(Model.IBS$X[,indice.re.IBS|N], eps.IBS)
      
      C.IBS = crossprod(Model.IBS$X[,indice.re.IBS|N], as.vector(pi.re.IBS*(1-pi.re.IBS))*Model.IBS$X[,indice.re.IBS|N])
      if(rcond(C.IBS) >= 1e-10){
        TS.IBS <- crossprod(Xeps.IBS, solve(C.IBS, Xeps.IBS))
      } else {
        TS.IBS = -1
      }
      #Incomp Score
      eps.Incomp <- Model.Incomp$Y-pi.re.Incomp #Y - e(Y)
      #this is the X^T times Y-E(Y)
      Xeps.Incomp <- crossprod(Model.Incomp$X[,indice.re.Incomp|N], eps.Incomp)
      
      C.Incomp = crossprod(Model.Incomp$X[,indice.re.Incomp|N], as.vector(pi.re.Incomp*(1-pi.re.Incomp))*Model.Incomp$X[,indice.re.Incomp|N])
      if(rcond(C.Incomp) >= 1e-10){
        TS.Incomp <- crossprod(Xeps.Incomp, solve(C.Incomp, Xeps.Incomp))
      } else {
        TS.Incomp = -1
      }
      #AMS Score
      eps.AMS <- Model.AMS$Y-pi.re.AMS #Y - e(Y)
      #this is the X^T times Y-E(Y)
      Xeps.AMS <- crossprod(Model.AMS$X[,indice.re.AMS|N], eps.AMS)
      
      C.AMS = crossprod(Model.AMS$X[,indice.re.AMS|N], as.vector(pi.re.AMS*(1-pi.re.AMS))*Model.AMS$X[,indice.re.AMS|N])
      if(rcond(C.AMS) >= 1e-10){
        TS.AMS <- crossprod(Xeps.AMS, solve(C.AMS, Xeps.AMS))
      } else {
        TS.AMS = -1
      }
      #Bin MM Score
      eps.BinMM <- Model.BinMM$Y-pi.re.BinMM #Y - e(Y)
      #this is the X^T times Y-E(Y)
      Xeps.BinMM <- crossprod(Model.BinMM$X[,indice.re.BinMM|N], eps.BinMM)
      
      C.BinMM = crossprod(Model.BinMM$X[,indice.re.BinMM|N], as.vector(pi.re.BinMM*(1-pi.re.BinMM))*Model.BinMM$X[,indice.re.BinMM|N])
      if(rcond(C.BinMM) >= 1e-10){
        TS.BinMM <- crossprod(Xeps.BinMM, solve(C.BinMM, Xeps.BinMM))
      } else {
        TS.BinMM = -1
      }
      
      #determine if null hyp is rejected or not
      #degrees of freedom is equal to the number of restricted values
      doF = sum(N == TRUE)
      
      #going to try p value using chi squared with doF degrees of freedom (1)
      pv.IBS = 1 - pchisq(TS.IBS, df = doF)
      pv.Incomp = 1 - pchisq(TS.Incomp, df = doF)
      pv.AMS = 1 - pchisq(TS.AMS, df = doF)
      pv.BinMM = 1 - pchisq(TS.BinMM, df = doF)
      
      beta.IBS = c(sum(beta.re.IBS!=0), sum(beta.unre.IBS!=0))
      beta.Incomp = c(sum(beta.re.Incomp!=0), sum(beta.unre.Incomp!=0))
      beta.AMS = c(sum(beta.re.AMS!=0), sum(beta.unre.AMS!=0))
      beta.BinMM = c(sum(beta.re.BinMM!=0), sum(beta.unre.BinMM!=0))
      
      statsAndPVals[1,] = c(TS.IBS, pv.IBS, beta.IBS)
      statsAndPVals[2,] = c(TS.Incomp, pv.Incomp, beta.Incomp)
      statsAndPVals[3,] = c(TS.AMS, pv.AMS, beta.AMS)
      statsAndPVals[4,] = c(TS.BinMM, pv.BinMM, beta.BinMM)
      
      #labeling
      colnames(statsAndPVals) = c("Penalized Score Stat", "p-value", "# Beta - Restricted", "# Betas - Unrestricted")
      rownames(statsAndPVals) = c("IBS Score", "Incomp Score", "AMS Score", "Binary Mismatch Score")
      
      write.csv(statsAndPVals, file = paste0(outFile))
    } else {      
      #if outcome is continuous
      #source the needed functions
      source(paste0(path,"/Linear_ADMM0.r"))
      
      # estimate the uncontrained estimator
      beta.unre.IBS <- cv.SCAD_ADMM_unre(X=Model.IBS$X, Y=Model.IBS$Y, N=N, beta0=rep(0, dim(Model.IBS$X)[2]), err=1e-4, tune="cv", unpen = unpenVals)
      indice.unre.IBS <- beta.unre.IBS!=0
      
      beta.unre.Incomp <- cv.SCAD_ADMM_unre(X=Model.Incomp$X, Y=Model.Incomp$Y, N=N, beta0=rep(0, dim(Model.Incomp$X)[2]), err=1e-4, tune="cv", unpen = unpenVals)
      indice.unre.Incomp <- beta.unre.Incomp!=0
      
      beta.unre.AMS <- cv.SCAD_ADMM_unre(X=Model.AMS$X, Y=Model.AMS$Y, N=N, beta0=rep(0, dim(Model.AMS$X)[2]), err=1e-4, tune="cv", unpen = unpenVals)
      indice.unre.AMS <- beta.unre.AMS!=0
      
      beta.unre.BinMM <- cv.SCAD_ADMM_unre(X=Model.BinMM$X, Y=Model.BinMM$Y, N=N, beta0=rep(0, dim(Model.BinMM$X)[2]), err=1e-4, tune="cv", unpen = unpenVals)
      indice.unre.BinMM <- beta.unre.BinMM!=0
      
      # estimate the constrained estimator
      beta.re.IBS <- cv.SCAD_ADMM_re(X=Model.IBS$X, Y=Model.IBS$Y, N=N, beta0=rep(0, dim(Model.IBS$X)[2]), err=1e-4, tune="cv", unpen = unpenVals)
      indice.re.IBS <- beta.re.IBS!=0
      
      beta.re.Incomp <- cv.SCAD_ADMM_re(X=Model.Incomp$X, Y=Model.Incomp$Y, N=N, beta0=rep(0, dim(Model.Incomp$X)[2]), err=1e-4, tune="cv", unpen = unpenVals)
      indice.re.Incomp <- beta.re.Incomp!=0
      
      beta.re.AMS <- cv.SCAD_ADMM_re(X=Model.AMS$X, Y=Model.AMS$Y, N=N, beta0=rep(0, dim(Model.AMS$X)[2]), err=1e-4, tune="cv", unpen = unpenVals)
      indice.re.AMS <- beta.re.AMS!=0
      
      beta.re.BinMM <- cv.SCAD_ADMM_re(X=Model.BinMM$X, Y=Model.BinMM$Y, N=N, beta0=rep(0, dim(Model.BinMM$X)[2]), err=1e-4, tune="cv", unpen = unpenVals)
      indice.re.BinMM <- beta.re.BinMM!=0
      
      # construct the score statistic
      eps.IBS <- Model.IBS$Y-Model.IBS$X%*%beta.re.IBS
      #this is the X^T times Y-E(Y)
      Xeps.IBS <- crossprod(Model.IBS$X[,indice.re.IBS|N], eps.IBS)
      
      D.IBS = crossprod(Model.IBS$X[,indice.re.IBS|N], Model.IBS$X[,indice.re.IBS|N])
      if(rcond(D.IBS) >= 1e-10){
        TS.IBS <- crossprod(Xeps.IBS, solve(D.IBS, Xeps.IBS))
      } else {
        TS.IBS = -1
      }
      
      eps.Incomp <- Model.Incomp$Y-Model.Incomp$X%*%beta.re.Incomp
      #this is the X^T times Y-E(Y)
      Xeps.Incomp <- crossprod(Model.Incomp$X[,indice.re.Incomp|N], eps.Incomp)
      
      D.Incomp = crossprod(Model.Incomp$X[,indice.re.Incomp|N], Model.Incomp$X[,indice.re.Incomp|N])
      if(rcond(D.Incomp) >= 1e-10){
        TS.Incomp <- crossprod(Xeps.Incomp, solve(D.Incomp, Xeps.Incomp))
      } else {
        TS.Incomp = -1
      }
      
      eps.AMS <- Model.AMS$Y-Model.AMS$X%*%beta.re.AMS
      #this is the X^T times Y-E(Y)
      Xeps.AMS <- crossprod(Model.AMS$X[,indice.re.AMS|N], eps.AMS)
      
      D.AMS = crossprod(Model.AMS$X[,indice.re.AMS|N], Model.AMS$X[,indice.re.AMS|N])
      if(rcond(D.AMS) >= 1e-10){
        TS.AMS <- crossprod(Xeps.AMS, solve(D.AMS, Xeps.AMS))
      } else {
        TS.AMS = -1
      }
      
      eps.BinMM <- Model.BinMM$Y-Model.BinMM$X%*%beta.re.BinMM
      #this is the X^T times Y-E(Y)
      Xeps.BinMM <- crossprod(Model.BinMM$X[,indice.re.BinMM|N], eps.BinMM)
      
      D.BinMM = crossprod(Model.BinMM$X[,indice.re.BinMM|N], Model.BinMM$X[,indice.re.BinMM|N])
      if(rcond(D.BinMM) >= 1e-10){
        TS.BinMM <- crossprod(Xeps.BinMM, solve(D.BinMM, Xeps.BinMM))
      } else {
        TS.BinMM = -1
      }
      
      #determine if null hyp is rejected or not
      #degrees of freedom is equal to the number of restricted values
      doF = sum(N == TRUE)
      
      #going to try p value using chi squared with doF degrees of freedom (1)
      pv.IBS = 1 - pchisq(TS.IBS, df = doF)
      pv.Incomp = 1 - pchisq(TS.Incomp, df = doF)
      pv.AMS = 1 - pchisq(TS.AMS, df = doF)
      pv.BinMM = 1 - pchisq(TS.BinMM, df = doF)
      
      beta.IBS = c(sum(beta.re.IBS!=0), sum(beta.unre.IBS!=0))
      beta.Incomp = c(sum(beta.re.Incomp!=0), sum(beta.unre.Incomp!=0))
      beta.AMS = c(sum(beta.re.AMS!=0), sum(beta.unre.AMS!=0))
      beta.BinMM = c(sum(beta.re.BinMM!=0), sum(beta.unre.BinMM!=0))
      
      statsAndPVals[1,] = c(TS.IBS, pv.IBS, beta.IBS)
      statsAndPVals[2,] = c(TS.Incomp, pv.Incomp, beta.Incomp)
      statsAndPVals[3,] = c(TS.AMS, pv.AMS, beta.AMS)
      statsAndPVals[4,] = c(TS.BinMM, pv.BinMM, beta.BinMM)
      
      #labeling
      colnames(statsAndPVals) = c("Penalized Score Stat", "p-value", "# Beta - Restricted", "# Betas - Unrestricted")
      rownames(statsAndPVals) = c("IBS Score", "Incomp Score", "AMS Score", "Binary Mismatch Score")
      
      write.csv(statsAndPVals, file = paste0(outFile))
    }
  } else {
    #define matrix to hold all Stats and Pvalues
    statsAndPVals = matrix(NA, nrow = 1, ncol = 4)
    
    #Read in the data
    #determine what file type
    fileName = read.table(text = datafile, sep = ".", as.is = TRUE)$V2
    if(fileName == "csv"){ #if it's a csv file, read in with read.csv
      geneticData = read.csv(datafile, header = TRUE, row.names = 1)
    } else { #otherwise should be good with a read.table
      geneticData = read.table(datafile, header = TRUE, row.names = 1)
    }
    
    geneticData = na.omit(geneticData)
    
    #if only 1 SNP, need to add in a column so it doesn't screw everything up
    if(dim(geneticData)[2] == 1){
      geneticData$Test = 0
    }
    ###############################################
    #Determine D/R pairs and match
    ###############################################  
    #read in the file with the donor/recipient pair information
    pairings = read.table(matchedPairs)
    #make sure column names are specified
    colnames(pairings) = c("R_Id", "D_Id")
    #assign each pair a number
    pairings$PairNumber = seq(1:length(pairings$R_Id))
    #split the original data into donor and recipient subsets
    donorData = geneticData[rownames(geneticData) %in% pairings$D_Id,]
    recipData = geneticData[rownames(geneticData) %in% pairings$R_Id,]
    #add an id column based on row names for merging
    donorData$D_Id = rownames(donorData)
    recipData$R_Id = rownames(recipData)
    #now merge the genotype data with the ids so that the genotype data has the pair numbers
    donorGenotypes_numbered = merge(donorData, pairings, by.x = "D_Id", by.y = "D_Id")
    recipGenotypes_numbered = merge(recipData, pairings, by.x = "R_Id", by.y = "R_Id")
    #make everything a data frame 
    donorGenotypes_numbered.df = as.data.frame(donorGenotypes_numbered)
    recipGenotypes_numbered.df = as.data.frame(recipGenotypes_numbered)
    #subset to only include Pairs with both D and R info
    donorGenotypes_numbered.df.subset = donorGenotypes_numbered.df[(donorGenotypes_numbered.df$PairNumber %in% recipGenotypes_numbered.df$PairNumber),]
    recipGenotypes_numbered.df.subset = recipGenotypes_numbered.df[(recipGenotypes_numbered.df$PairNumber %in% donorGenotypes_numbered.df$PairNumber),]
    #order by Pair number
    donorGenotypes_numbered.df.ordered = donorGenotypes_numbered.df.subset[order(donorGenotypes_numbered.df.subset$PairNumber),]
    recipGenotypes_numbered.df.ordered = recipGenotypes_numbered.df.subset[order(recipGenotypes_numbered.df.subset$PairNumber),]
    #make rownames into R_Id or D_Id
    rownames(donorGenotypes_numbered.df.ordered) = donorGenotypes_numbered.df.ordered$D_Id
    rownames(recipGenotypes_numbered.df.ordered) = recipGenotypes_numbered.df.ordered$R_Id
    #remove non-SNP columns
    DGenosData = donorGenotypes_numbered.df.ordered[,2:(ncol(donorGenotypes_numbered.df.ordered)-2)]
    RGenosData = recipGenotypes_numbered.df.ordered[,2:(ncol(recipGenotypes_numbered.df.ordered)-2)]
    
    ###############################################
    # Read in the covariate/phenotype file
    ############################################### 
    covariates = read.table(covPhenoFile, header = TRUE)
    #remove nas
    covariates = na.omit(covariates)
    #pull only those Rs from the R geno matrix
    covariates = covariates[covariates[,1] %in% rownames(RGenosData),]
    #set row names to be Ids from 1st row
    rownames(covariates) = covariates[,1]
    #remove the ids column
    covariates = covariates[,-1]
    
    if(dim(RGenosData)[1] > dim(covariates)[1]){
      RGenosData = RGenosData[rownames(RGenosData) %in% rownames(covariates),]
      DGenosData = DGenosData[which(rownames(RGenosData) %in% rownames(covariates)),]
    } else if (dim(RGenosData)[1] < dim(covariates)[1]){
      covariates = covariates[rownames(covariates) %in% rownames(RGenosData),]
    }
    
    #make matrices
    #unlist and make into matrices
    RGenosMat = matrix(unlist(RGenosData), ncol = ncol(RGenosData), byrow = F)
    DGenosMat = matrix(unlist(DGenosData), ncol = ncol(DGenosData), byrow = F)
    #get rid of the Test column if it exists
    if(("Test" %in% colnames(RGenosData)) == TRUE){
      RGenosMat = matrix(RGenosMat[,1], ncol = 1)
      DGenosMat = matrix(DGenosMat[,1], ncol = 1)
    }
    #keep number of pairs for later use
    numPairs = nrow(RGenosMat)
    numSNPs = ncol(RGenosMat)
  
    #make a matrix
    #pull the phenotypes from the covariates file
    phenos = as.matrix(covariates[,1],nrow=numPairs, ncol = 1)
    #remove the phenotypes column from the covariates
    CovData.list = covariates[,-1]
    CovData = matrix(unlist(CovData.list), ncol = ncol(CovData.list), byrow = F)
    
    #calculate single snp scores
    if(score == "IBS"){
      Score.snp = calcIBSMismatch(RGenosMat, DGenosMat)
    } else if(score == "Incomp"){
      Score.snp = calcIncompatibilityScore(RGenosMat, DGenosMat)
    } else if(score == "AMS"){
      Score.snp = calcAMS(RGenosMat, DGenosMat)
    } else{
      Score.snp = calcBinaryMM(RGenosMat, DGenosMat)
    }
    
    #calculate gene based scores
    #check to see if weights are used for scores
    if(weightedScores == FALSE){
      #check to see if scores should be standardized (can't be weighted and standardized)
      if(standardizeScores == FALSE){
        Score.gene = calcGeneScore(SingleSNPKernel = Score.snp, standardize = FALSE, useWeights = FALSE)
      } else {
        Score.gene = calcGeneScore(SingleSNPKernel = Score.snp, standardize = TRUE, useWeights = FALSE)
      }
    } else {
      Score.gene = calcGeneScore(SingleSNPKernel = Score.snp, standardize = FALSE, useWeights = TRUE, scoreWeights)
    }
    
    # define location of zero components
    # basically, this defines what the null hypothesis is
    # so for gamma = 0 null, we need the non-zero components to be for covariates
    # and the beta, and gamma components to be zero
    numCov = dim(CovData)[2]
    numSNPs = dim(RGenosMat)[2]
    
    #if there is a forced intercept
    if(fitIntercept == TRUE){
      intercept = matrix(1,nrow = dim(RGenosMat)[1], ncol = 1)
      N = c(rep(FALSE, numCov+1+numSNPs), TRUE)
      designMat = cbind(intercept, CovData, RGenosMat, Score.gene)
    } else {
      #otherwise no intercept included
      N = c(rep(FALSE, numCov+numSNPs), TRUE)
      designMat = cbind(CovData, RGenosMat, Score.gene)
    }
    
    #then combine the phenos with the design matrix as a list
    Model = list(X=designMat, Y=phenos)
    
    #determine what to not penalize in the model
    if(length(unpen) == 0){
      #if not specified, do not penalize the covariates
      #number of covariates in the model (with or without intercept) is the length of N, minus the number of SNPs and the gene score
      totalNumCovs =  length(N) - numSNPs - 1
      unpenVals = seq(1,totalNumCovs,1)
    } else {
      unpenVals = unpen
    }
    
    #if outcome is binary
    if(phenosBinary == TRUE){
      #source the needed functions
      source(paste0(path,"/Logistic_ADMM0.r"))
      
      # estimate the uncontrained estimator
      beta.unre <- cv.SCAD_ADMM_unre(X=Model$X, Y=Model$Y, N=N, beta0=rep(0, dim(Model$X)[2]), err=1e-4, tune="cv", unpen=unpenVals)
      indice.unre <- beta.unre!=0
      
      # estimate the constrained estimator (only need this for score test)
      beta.re <- cv.SCAD_ADMM_re(X=Model$X, Y=Model$Y, N=N, beta0=rep(0, dim(Model$X)[2]), err=1e-4, tune="cv", unpen=unpenVals)
      indice.re <- beta.re!=0
      
      pi.unre <- logit(Model$X%*%beta.unre)
      pi.re <- logit(Model$X%*%beta.re)
      
      # construct the score statistics
      eps <- Model$Y-pi.re #Y - e(Y)
      #this is the X^T times Y-E(Y)
      Xeps <- crossprod(Model$X[,indice.re|N], eps)
      
      C = crossprod(Model$X[,indice.re|N], as.vector(pi.re*(1-pi.re))*Model$X[,indice.re|N])
      if(rcond(C) >= 1e-10){
        TS <- crossprod(Xeps, solve(C, Xeps))
      } else {
        TS = -1
      }
      
      #determine if null hyp is rejected or not
      #degrees of freedom is equal to the number of restricted values
      doF = sum(N == TRUE)
      
      #going to try p value using chi squared with doF degrees of freedom (1)
      pv = 1 - pchisq(TS, df = doF)
      
      beta = c(sum(beta.re!=0), sum(beta.unre!=0))
      
      statsAndPVals[1,] = c(TS, pv, beta)
      
      #labeling
      colnames(statsAndPVals) = c("Penalized Score Stat", "p-value", "# Beta - Restricted", "# Betas - Unrestricted")
      rownames(statsAndPVals) = score
      
      write.csv(statsAndPVals, file = paste0(outFile))
    } else {      
      #if outcome is continuous
      #source the needed functions
      source(paste0(path,"/Linear_ADMM0.r"))
      
      # estimate the uncontrained estimator
      beta.unre <- cv.SCAD_ADMM_unre(X=Model$X, Y=Model$Y, N=N, beta0=rep(0, dim(Model$X)[2]), err=1e-4, tune="cv", unpen = unpenVals)
      indice.unre <- beta.unre!=0
      
      # estimate the constrained estimator
      beta.re <- cv.SCAD_ADMM_re(X=Model$X, Y=Model$Y, N=N, beta0=rep(0, dim(Model$X)[2]), err=1e-4, tune="cv", unpen = unpenVals)
      indice.re <- beta.re!=0
      
      # construct the score statistic
      eps <- Model$Y-Model$X%*%beta.re
      #this is the X^T times Y-E(Y)
      Xeps <- crossprod(Model$X[,indice.re|N], eps)
      
      D = crossprod(Model$X[,indice.re|N], Model$X[,indice.re|N])
      if(rcond(D) >= 1e-10){
        TS <- crossprod(Xeps, solve(D, Xeps))
      } else {
        TS = -1
      }
      
      #determine if null hyp is rejected or not
      #degrees of freedom is equal to the number of restricted values
      doF = sum(N == TRUE)
      
      #going to try p value using chi squared with doF degrees of freedom (1)
      pv = 1 - pchisq(TS, df = doF)
      
      beta = c(sum(beta.re!=0), sum(beta.unre!=0))
      
      statsAndPVals[1,] = c(TS, pv, beta)
      
      #labeling
      colnames(statsAndPVals) = c("Penalized Score Stat", "p-value", "# Beta - Restricted", "# Betas - Unrestricted")
      rownames(statsAndPVals) = score
        
      write.csv(statsAndPVals, file = paste0(outFile))
    }
  }
}
