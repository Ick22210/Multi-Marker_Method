#############################################
## File to source for whole pipelines
## Includes TIE and Power pipelines
## Using long form eqn and Ghat for calculations
#############################################
#edit 5/21/20 added in ability to standardize scores and weight scores to ghat pipelines
#edited version for HapGen Data, Only using Ghat versions
#edit 6/18/2020 added in code for running SKAT for comparison
#edited 9/24/2020 to run with ProjectIISourceFunctions_v2.R so that s is equal to the number of PCs, not the PVE (this was changed back!)

#Relies on ProjectIISourceFunctions_v2.R

##################################################
## Original Pipelines
##################################################
# TIE Pipeline
#TIE Pipeline with Ghat
RunTIEPipelineLocalGHat = function(chr, gene, numPairs, YPrev, s, standardizeScores = FALSE, weightedScores = FALSE, scoreWeights){
  #function to run whole TIE pipeline locally, uses Ghat instead of long form eqn for calculations
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #s = PVE for choosing the number of PCs
  #standardizeScores = T or F whether the scores should be standardized based on maximum score value
  #weightedScores = T or F whether the scores will be weighted
  #scoreWeights = m x 1 vector of weights, one weight for each SNP
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes TIE values to csv files
  
  library(parallel)
  
  #always the same
  numSims = 5000  
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_II_Sims/HapGen_Files/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/ProjectIISourceFunctions_v2.R")
  
  myList = lapply(1:numSims, rep, times = 1)
  
  statsAndPVals = mclapply(myList, function(ii){    
    #define matrix to hold all Stats and Pvalues
    statsAndPVals = matrix(NA, nrow = numSims, ncol = 16)
    
    #pull recipient and donor genotypes
    RGenos = obtainRGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    DGenos = obtainDGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    
    #calculate single snp scores
    IBS.snp = calcIBSMismatch(RGenosMat = RGenos, DGenosMat = DGenos)
    Incomp.snp = calcIncompatibilityScore(RGenosMat = RGenos, DGenosMat = DGenos)
    AMS.snp = calcAMS(RGenosMat = RGenos, DGenosMat = DGenos)
    BinMM.snp = calcBinaryMM(RGenosMat = RGenos, DGenosMat = DGenos)
    
    #calculate gene based scores
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
    
    #generate covariates
    #for now, a single binary and a single continous covariate
    CovData = GenCovData(SampleSize = numPairs, BinaryValues = 1, ContinuousValues = 1)
    
    #generate phenotypes, both continuous and binary
    CatPhenos = GenNullPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = TRUE, YPrev = YPrev,  Covariates = CovData)
    ContPhenos = GenNullPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = FALSE,  Covariates = CovData)
    
    #Need to calculate UR and US scores
    #will have then for continuous and binary phenos
    #need separate US for each score type
    #Recipient genotype UR values:
    UR_CatPhenos_Ghat = CalcUScoreGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = TRUE, RGenoData = RGenos, Phenos = CatPhenos, BinPhenos = TRUE)
    UR_ContPhenos_Ghat = CalcUScoreGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = TRUE, RGenoData = RGenos, Phenos = ContPhenos, BinPhenos = FALSE)
    #US for IBS Score
    US_IBS_CatPhenos_Ghat = CalcUScoreGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = IBS.gene, Phenos = CatPhenos, BinPhenos = TRUE)
    US_IBS_ContPhenos_Ghat = CalcUScoreGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = IBS.gene, Phenos = ContPhenos, BinPhenos = FALSE)
    #US for Incomp Score
    US_Incomp_CatPhenos_Ghat = CalcUScoreGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = Incomp.gene, Phenos = CatPhenos, BinPhenos = TRUE)
    US_Incomp_ContPhenos_Ghat = CalcUScoreGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = Incomp.gene, Phenos = ContPhenos, BinPhenos = FALSE)
    #US for AMS Score
    US_AMS_CatPhenos_Ghat = CalcUScoreGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = AMS.gene, Phenos = CatPhenos, BinPhenos = TRUE)
    US_AMS_ContPhenos_Ghat = CalcUScoreGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = AMS.gene, Phenos = ContPhenos, BinPhenos = FALSE)
    #US for Bin MM Score
    US_BinMM_CatPhenos_Ghat = CalcUScoreGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = BinMM.gene, Phenos = CatPhenos, BinPhenos = TRUE)
    US_BinMM_ContPhenos_Ghat = CalcUScoreGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = BinMM.gene, Phenos = ContPhenos, BinPhenos = FALSE)
    
    #Need to calculate Q values
    #R geno QR values
    QR_CatPhenos_Ghat = CalcQValuesGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = TRUE, RGenoData = RGenos, Phenos = CatPhenos, BinPhenos = TRUE)
    QR_ContPhenos_Ghat = CalcQValuesGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = TRUE, RGenoData = RGenos, Phenos = ContPhenos, BinPhenos = FALSE)
    #QS for IBS Score
    QS_IBS_CatPhenos_Ghat = CalcQValuesGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = IBS.gene, Phenos = CatPhenos, BinPhenos = TRUE)
    QS_IBS_ContPhenos_Ghat = CalcQValuesGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = IBS.gene, Phenos = ContPhenos, BinPhenos = FALSE)
    #QS for Incomp Score
    QS_Incomp_CatPhenos_Ghat = CalcQValuesGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = Incomp.gene, Phenos = CatPhenos, BinPhenos = TRUE)
    QS_Incomp_ContPhenos_Ghat = CalcQValuesGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = Incomp.gene, Phenos = ContPhenos, BinPhenos = FALSE)
    #QS for AMS score
    QS_AMS_CatPhenos_Ghat = CalcQValuesGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = AMS.gene, Phenos = CatPhenos, BinPhenos = TRUE)
    QS_AMS_ContPhenos_Ghat = CalcQValuesGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = AMS.gene, Phenos = ContPhenos, BinPhenos = FALSE)
    #QS for Bin MM Score
    QS_BinMM_CatPhenos_Ghat = CalcQValuesGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = BinMM.gene, Phenos = CatPhenos, BinPhenos = TRUE)
    QS_BinMM_ContPhenos_Ghat = CalcQValuesGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = BinMM.gene, Phenos = ContPhenos, BinPhenos = FALSE)
    
    #Create combined Qs (QR, QS)
    #R geno and IBS score
    Q_IBS_CatPhenos_Ghat = cbind(QR_CatPhenos_Ghat, QS_IBS_CatPhenos_Ghat)
    Q_IBS_ContPhenos_Ghat = cbind(QR_ContPhenos_Ghat, QS_IBS_ContPhenos_Ghat) 
    #R geno and Incomp Score
    Q_Incomp_CatPhenos_Ghat = cbind(QR_CatPhenos_Ghat, QS_Incomp_CatPhenos_Ghat)
    Q_Incomp_ContPhenos_Ghat = cbind(QR_ContPhenos_Ghat, QS_Incomp_ContPhenos_Ghat)
    #R geno and AMS score
    Q_AMS_CatPhenos_Ghat = cbind(QR_CatPhenos_Ghat, QS_AMS_CatPhenos_Ghat)
    Q_AMS_ContPhenos_Ghat = cbind(QR_ContPhenos_Ghat, QS_AMS_ContPhenos_Ghat)
    #R geno and Bin MM score
    Q_BinMM_CatPhenos_Ghat = cbind(QR_CatPhenos_Ghat, QS_BinMM_CatPhenos_Ghat)
    Q_BinMM_ContPhenos_Ghat = cbind(QR_ContPhenos_Ghat, QS_BinMM_ContPhenos_Ghat)
    
    #calculate full variance
    #each set of Qs will have a variance calculation
    #R genos + IBS Score
    OrgVar_Q_IBS_CatPhenos_Ghat = CalcVariance(SampleSize = numPairs, QValues = Q_IBS_CatPhenos_Ghat)
    OrgVar_Q_IBS_ContPhenos_Ghat = CalcVariance(SampleSize = numPairs, QValues = Q_IBS_ContPhenos_Ghat)
    #R genos + Incomp Score
    OrgVar_Q_Incomp_CatPhenos_Ghat = CalcVariance(SampleSize = numPairs, QValues = Q_Incomp_CatPhenos_Ghat)
    OrgVar_Q_Incomp_ContPhenos_Ghat = CalcVariance(SampleSize = numPairs, QValues = Q_Incomp_ContPhenos_Ghat)
    #R genos + AMS Score
    OrgVar_Q_AMS_CatPhenos_Ghat = CalcVariance(SampleSize = numPairs, QValues = Q_AMS_CatPhenos_Ghat)
    OrgVar_Q_AMS_ContPhenos_Ghat = CalcVariance(SampleSize = numPairs, QValues = Q_AMS_ContPhenos_Ghat)
    #R genos + Bin MM Score
    OrgVar_Q_BinMM_CatPhenos_Ghat = CalcVariance(SampleSize = numPairs, QValues = Q_BinMM_CatPhenos_Ghat)
    OrgVar_Q_BinMM_ContPhenos_Ghat = CalcVariance(SampleSize = numPairs, QValues = Q_BinMM_ContPhenos_Ghat)
    
    #calculate final stat and p value
    #R geno and IBS score
    Stat_IBS_CatPhenos = CalcStatisticPVal(SampleSize = numPairs, Variance = OrgVar_Q_IBS_CatPhenos_Ghat, UscoresR = UR_CatPhenos_Ghat, UscoreS = US_IBS_CatPhenos_Ghat, s = s)
    Stat_IBS_ContPhenos = CalcStatisticPVal(SampleSize = numPairs, Variance = OrgVar_Q_IBS_ContPhenos_Ghat, UscoresR = UR_ContPhenos_Ghat, UscoreS = US_IBS_ContPhenos_Ghat, s = s)
    #R geno and Incomp score
    Stat_Incomp_CatPhenos = CalcStatisticPVal(SampleSize = numPairs, Variance = OrgVar_Q_Incomp_CatPhenos_Ghat, UscoresR = UR_CatPhenos_Ghat, UscoreS = US_Incomp_CatPhenos_Ghat, s = s)
    Stat_Incomp_ContPhenos = CalcStatisticPVal(SampleSize = numPairs, Variance = OrgVar_Q_Incomp_ContPhenos_Ghat, UscoresR = UR_ContPhenos_Ghat, UscoreS = US_Incomp_ContPhenos_Ghat, s = s)
    #R geno and AMS Score
    Stat_AMS_CatPhenos = CalcStatisticPVal(SampleSize = numPairs, Variance = OrgVar_Q_AMS_CatPhenos_Ghat, UscoresR = UR_CatPhenos_Ghat, UscoreS = US_AMS_CatPhenos_Ghat, s = s)
    Stat_AMS_ContPhenos = CalcStatisticPVal(SampleSize = numPairs, Variance = OrgVar_Q_AMS_ContPhenos_Ghat, UscoresR = UR_ContPhenos_Ghat, UscoreS = US_AMS_ContPhenos_Ghat, s = s)
    #R geno and Bin MM score
    Stat_BinMM_CatPhenos = CalcStatisticPVal(SampleSize = numPairs, Variance = OrgVar_Q_BinMM_CatPhenos_Ghat, UscoresR = UR_CatPhenos_Ghat, UscoreS = US_BinMM_CatPhenos_Ghat, s = s)
    Stat_BinMM_ContPhenos = CalcStatisticPVal(SampleSize = numPairs, Variance = OrgVar_Q_BinMM_ContPhenos_Ghat, UscoresR = UR_ContPhenos_Ghat, UscoreS = US_BinMM_ContPhenos_Ghat, s = s)
    
    #fill columns in order
    ##Binary Phenos
    ##IBS Cat, Incomp Cat, AMS Cat, Bin MM Cat, 
    statsAndPValsCatPhenos_IBS = Stat_IBS_CatPhenos[1:2]
    statsAndPValsCatPhenos_Incomp = Stat_Incomp_CatPhenos[1:2]
    statsAndPValsCatPhenos_AMS = Stat_AMS_CatPhenos[1:2]
    statsAndPValsCatPhenos_BinMM = Stat_BinMM_CatPhenos[1:2]
    ##Cont Phenos
    ##IBS Cont, Incomp Cont, AMS Cont, Bin MM Cont
    statsAndPValsContPhenos_IBS = Stat_IBS_ContPhenos[1:2]
    statsAndPValsContPhenos_Incomp = Stat_Incomp_ContPhenos[1:2]
    statsAndPValsContPhenos_AMS = Stat_AMS_ContPhenos[1:2]
    statsAndPValsContPhenos_BinMM = Stat_BinMM_ContPhenos[1:2]
    
    IBSCat.mat = as.matrix(statsAndPValsCatPhenos_IBS)
    IncompCat.mat = as.matrix(statsAndPValsCatPhenos_Incomp)
    AMSCat.mat = as.matrix(statsAndPValsCatPhenos_AMS)
    BinMMCat.mat = as.matrix(statsAndPValsCatPhenos_BinMM)
    
    IBSCont.mat = as.matrix(statsAndPValsContPhenos_IBS)
    IncompCont.mat = as.matrix(statsAndPValsContPhenos_Incomp)
    AMSCont.mat = as.matrix(statsAndPValsContPhenos_AMS)
    BinMMCont.mat = as.matrix(statsAndPValsContPhenos_BinMM)
    
    statsAndPVals = cbind(t(IBSCat.mat),t(IncompCat.mat),t(AMSCat.mat),t(BinMMCat.mat),t(IBSCont.mat),t(IncompCont.mat),t(AMSCont.mat),t(BinMMCont.mat))
    
    print(paste0("Simulation ",ii," is complete."))
    
    statsAndPVals
  }, mc.cores=4)
  
  statsAndPValsAll = matrix(unlist(statsAndPVals),nrow = numSims, ncol = 16, byrow = TRUE)
  statsAndPValsCatPhenos = statsAndPValsAll[,1:8]
  statsAndPValsContPhenos = statsAndPValsAll[,9:16]
  
  #write out the Stats and p values
  if(weightedScores == FALSE){
    if(standardizeScores == FALSE){
      write.csv(statsAndPValsCatPhenos, file = paste0(path,"/TIE_CatPhenos_Prev",YPrev*100,"_Ghat_StatsAndPValues_s",s*100,".csv"))
      write.csv(statsAndPValsContPhenos, file = paste0(path,"/TIE_ContPhenos_Ghat_StatsAndPValues_s",s*100,".csv"))
    } else { #scores are standardized
      write.csv(statsAndPValsCatPhenos, file = paste0(path,"/TIE_CatPhenos_Prev",YPrev*100,"_Ghat_StandardizedScores_StatsAndPValues_s",s*100,".csv"))
      write.csv(statsAndPValsContPhenos, file = paste0(path,"/TIE_ContPhenos_Ghat_StandardizedScores_StatsAndPValues_s",s*100,".csv"))
    }
  } else { #scores are weighted
    write.csv(statsAndPValsCatPhenos, file = paste0(path,"/TIE_CatPhenos_Prev",YPrev*100,"_Ghat_WeightedScores_StatsAndPValues_s",s*100,".csv"))
    write.csv(statsAndPValsContPhenos, file = paste0(path,"/TIE_ContPhenos_Ghat_WeightedScores_StatsAndPValues_s",s*100,".csv"))
  }
  
  #calculate the type I errors for each combination of score and r geno
  #first column for cat phenos, second for cont phenos
  #rows go IBS, Incomp, AMS, Bin MM
  ties = matrix(nrow = 4, ncol = 2)
  
  #pull only the pvalues
  PValsCatPhenos = statsAndPValsCatPhenos[,c(2,4,6,8)]
  PValsContPhenos = statsAndPValsContPhenos[,c(2,4,6,8)]
  
  #calc type I error
  for(jj in 1:4){
    ties[jj,1] = sum(PValsCatPhenos[,jj] <= 0.05)/numSims
    ties[jj,2] = sum(PValsContPhenos[,jj] <= 0.05)/numSims
  }
  
  #write out type I error results
  if(weightedScores == FALSE){
    if(standardizeScores == FALSE){
      write.csv(ties, file = paste0(path,"/TIE_Results_Ghat_CatAndContPhenos_Prev",YPrev*100,"_s",s*100,".csv"))
    } else { #scores are standardized
      write.csv(ties, file = paste0(path,"/TIE_Results_Ghat_StandardizedScores_CatAndContPhenos_Prev",YPrev*100,"_s",s*100,".csv"))
    }
  } else { #scores are weighted 
    write.csv(ties, file = paste0(path,"/TIE_Results_Ghat_WeightedScores_CatAndContPhenos_Prev",YPrev*100,"_s",s*100,".csv"))
  }
}

#Power Pipeline with Ghat for Scores being true
RunPowerPipelineLocalGhat_Scores = function(chr, gene, numPairs, YPrev, s, Gamma, TrueScore, ORSize, standardizeScores = FALSE, weightedScores = FALSE, scoreWeights, percentageAssoc, LowLD){
  #function to run whole power pipeline locally, uses Ghat instead of long form eqn for calculations
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #s = PVE for choosing the number of PCs
  #Gamma = effect size for score, length 1, can be 0
  #TrueScore = IBS.gene, Incomp.gene, AMS.gene, or BinMM.gene
  #ORSize = Small, Medium, or Large for what OR was used for the associated SNP/score
  #standardizeScores = T or F whether the scores should be standardized based on maximum score value
  #weightedScores = T or F whether the scores will be weighted
  #scoreWeights = m x 1 vector of weights, one weight for each SNP
  #percentageAssoc = percentage of SNPs associated with outcome (either 5, 25, 50, 75, or 100) 
  #LowLD = True or FALSE whether the associated SNPs are in low LD or high LD
  #### This value shows which score is potentially being used to create phenotypes
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes power values to csv files
  
  library(parallel)
  
  #need to define this for naming at the end
  snpOrScore = TrueScore
  
  #number of sims always the same
  numSims = 5000
  
  #define path to data
  #for sampled haplotype data
  # path = paste0("/home/vlynn/Paper_II_Sims/",gene,"_Results_",numPairs,"Pairs")
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_II_Sims/HapGen_Files/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/ProjectIISourceFunctions_v2.R")
  
  #define matrix to hold all Stats and Pvalues
  myList = lapply(1:numSims, rep, times = 1)
  statsAndPValsMat = matrix(NA,ncol=24)
  
  statsAndPValsMat = mclapply(myList, function(ii){
    statsAndPValsMat = matrix(NA,ncol=24)
    
    #pull recipient and donor genotypes
    RGenos = obtainRGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    DGenos = obtainDGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    
    #calculate single snp scores
    IBS.snp = calcIBSMismatch(RGenosMat = RGenos, DGenosMat = DGenos)
    Incomp.snp = calcIncompatibilityScore(RGenosMat = RGenos, DGenosMat = DGenos)
    AMS.snp = calcAMS(RGenosMat = RGenos, DGenosMat = DGenos)
    BinMM.snp = calcBinaryMM(RGenosMat = RGenos, DGenosMat = DGenos)
    
    #calculate gene based scores
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
    
    #also need to calculate gene based scores if not all SNPs are associated
    IBS.gene.PercentOfSNPs = calcGeneScorePercentOfSNPs(SingleSNPKernel = IBS.snp, gene = gene, percentageAssoc = percentageAssoc, LowLD = LowLD, standardize = FALSE, useWeights = FALSE)
    Incomp.gene.PercentOfSNPs = calcGeneScorePercentOfSNPs(SingleSNPKernel = Incomp.snp, gene =  gene, percentageAssoc = percentageAssoc, LowLD = LowLD, standardize = FALSE, useWeights = FALSE)
    AMS.gene.PercentOfSNPs = calcGeneScorePercentOfSNPs(SingleSNPKernel = AMS.snp, gene =  gene, percentageAssoc = percentageAssoc, LowLD = LowLD, standardize = FALSE, useWeights = FALSE)
    BinMM.gene.PercentOfSNPs = calcGeneScorePercentOfSNPs(SingleSNPKernel = BinMM.snp, gene =  gene, percentageAssoc = percentageAssoc, LowLD = LowLD, standardize = FALSE, useWeights = FALSE)
    
    #need to use TrueScore to pull gene based scores matrix for generating phenotypes
    if(TrueScore == "IBS.gene"){
      PhenoScore = IBS.gene.PercentOfSNPs
    } else if(TrueScore == "Incomp.gene"){
      PhenoScore = Incomp.gene.PercentOfSNPs
    } else if(TrueScore == "AMS.gene"){
      PhenoScore = AMS.gene.PercentOfSNPs
    } else {
      PhenoScore = BinMM.gene.PercentOfSNPs
    }
    
    #generate covariates
    #for now, a single binary and a single continous covariate
    CovData = GenCovData(SampleSize = numPairs, BinaryValues = 1, ContinuousValues = 1)
    
    #need to define null Betas for phenotype generation
    nSNP = ncol(RGenos) #this should be the number of SNPs
    Betas = rep(0,nSNP) #generate null beta values
    Betas = as.matrix(Betas, ncol = 1)
    
    #generate phenotypes, both continuous and binary
    #Based on single true score
    CatPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = TRUE, YPrev = YPrev,  Covariates = CovData, RGenoData = RGenos, ScoreData = PhenoScore, Betas = Betas, Gamma = Gamma)
    ContPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = FALSE,  Covariates = CovData, RGenoData = RGenos, ScoreData = PhenoScore, Betas = Betas, Gamma = Gamma)
    
    #Need to calculate UR and US scores
    #will have then for continuous and binary phenos
    #need separate US for each score type
    #Recipient genotype UR values:
    UR_CatPhenos_Ghat = CalcUScoreGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = TRUE, RGenoData = RGenos, Phenos = CatPhenos, BinPhenos = TRUE)
    UR_ContPhenos_Ghat = CalcUScoreGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = TRUE, RGenoData = RGenos, Phenos = ContPhenos, BinPhenos = FALSE)
    #US for IBS Score
    US_IBS_CatPhenos_Ghat = CalcUScoreGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = IBS.gene, Phenos = CatPhenos, BinPhenos = TRUE)
    US_IBS_ContPhenos_Ghat = CalcUScoreGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = IBS.gene, Phenos = ContPhenos, BinPhenos = FALSE)
    #US for Incomp Score
    US_Incomp_CatPhenos_Ghat = CalcUScoreGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = Incomp.gene, Phenos = CatPhenos, BinPhenos = TRUE)
    US_Incomp_ContPhenos_Ghat = CalcUScoreGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = Incomp.gene, Phenos = ContPhenos, BinPhenos = FALSE)
    #US for AMS Score
    US_AMS_CatPhenos_Ghat = CalcUScoreGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = AMS.gene, Phenos = CatPhenos, BinPhenos = TRUE)
    US_AMS_ContPhenos_Ghat = CalcUScoreGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = AMS.gene, Phenos = ContPhenos, BinPhenos = FALSE)
    #US for Bin MM Score
    US_BinMM_CatPhenos_Ghat = CalcUScoreGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = BinMM.gene, Phenos = CatPhenos, BinPhenos = TRUE)
    US_BinMM_ContPhenos_Ghat = CalcUScoreGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = BinMM.gene, Phenos = ContPhenos, BinPhenos = FALSE)
    
    #Need to calculate Q values
    #R geno QR values
    QR_CatPhenos_Ghat = CalcQValuesGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = TRUE, RGenoData = RGenos, Phenos = CatPhenos, BinPhenos = TRUE)
    QR_ContPhenos_Ghat = CalcQValuesGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = TRUE, RGenoData = RGenos, Phenos = ContPhenos, BinPhenos = FALSE)
    #QS for IBS Score
    QS_IBS_CatPhenos_Ghat = CalcQValuesGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = IBS.gene, Phenos = CatPhenos, BinPhenos = TRUE)
    QS_IBS_ContPhenos_Ghat = CalcQValuesGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = IBS.gene, Phenos = ContPhenos, BinPhenos = FALSE)
    #QS for Incomp Score
    QS_Incomp_CatPhenos_Ghat = CalcQValuesGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = Incomp.gene, Phenos = CatPhenos, BinPhenos = TRUE)
    QS_Incomp_ContPhenos_Ghat = CalcQValuesGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = Incomp.gene, Phenos = ContPhenos, BinPhenos = FALSE)
    #QS for AMS score
    QS_AMS_CatPhenos_Ghat = CalcQValuesGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = AMS.gene, Phenos = CatPhenos, BinPhenos = TRUE)
    QS_AMS_ContPhenos_Ghat = CalcQValuesGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = AMS.gene, Phenos = ContPhenos, BinPhenos = FALSE)
    #QS for Bin MM Score
    QS_BinMM_CatPhenos_Ghat = CalcQValuesGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = BinMM.gene, Phenos = CatPhenos, BinPhenos = TRUE)
    QS_BinMM_ContPhenos_Ghat = CalcQValuesGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = BinMM.gene, Phenos = ContPhenos, BinPhenos = FALSE)
    
    #Create combined Qs (QR, QS)
    #R geno and IBS score
    Q_IBS_CatPhenos_Ghat = cbind(QR_CatPhenos_Ghat, QS_IBS_CatPhenos_Ghat)
    Q_IBS_ContPhenos_Ghat = cbind(QR_ContPhenos_Ghat, QS_IBS_ContPhenos_Ghat) 
    #R geno and Incomp Score
    Q_Incomp_CatPhenos_Ghat = cbind(QR_CatPhenos_Ghat, QS_Incomp_CatPhenos_Ghat)
    Q_Incomp_ContPhenos_Ghat = cbind(QR_ContPhenos_Ghat, QS_Incomp_ContPhenos_Ghat)
    #R geno and AMS score
    Q_AMS_CatPhenos_Ghat = cbind(QR_CatPhenos_Ghat, QS_AMS_CatPhenos_Ghat)
    Q_AMS_ContPhenos_Ghat = cbind(QR_ContPhenos_Ghat, QS_AMS_ContPhenos_Ghat)
    #R geno and Bin MM score
    Q_BinMM_CatPhenos_Ghat = cbind(QR_CatPhenos_Ghat, QS_BinMM_CatPhenos_Ghat)
    Q_BinMM_ContPhenos_Ghat = cbind(QR_ContPhenos_Ghat, QS_BinMM_ContPhenos_Ghat)
    
    #calculate full variance
    #each set of Qs will have a variance calculation
    #R genos + IBS Score
    OrgVar_Q_IBS_CatPhenos_Ghat = CalcVariance(SampleSize = numPairs, QValues = Q_IBS_CatPhenos_Ghat)
    OrgVar_Q_IBS_ContPhenos_Ghat = CalcVariance(SampleSize = numPairs, QValues = Q_IBS_ContPhenos_Ghat)
    #R genos + Incomp Score
    OrgVar_Q_Incomp_CatPhenos_Ghat = CalcVariance(SampleSize = numPairs, QValues = Q_Incomp_CatPhenos_Ghat)
    OrgVar_Q_Incomp_ContPhenos_Ghat = CalcVariance(SampleSize = numPairs, QValues = Q_Incomp_ContPhenos_Ghat)
    #R genos + AMS Score
    OrgVar_Q_AMS_CatPhenos_Ghat = CalcVariance(SampleSize = numPairs, QValues = Q_AMS_CatPhenos_Ghat)
    OrgVar_Q_AMS_ContPhenos_Ghat = CalcVariance(SampleSize = numPairs, QValues = Q_AMS_ContPhenos_Ghat)
    #R genos + Bin MM Score
    OrgVar_Q_BinMM_CatPhenos_Ghat = CalcVariance(SampleSize = numPairs, QValues = Q_BinMM_CatPhenos_Ghat)
    OrgVar_Q_BinMM_ContPhenos_Ghat = CalcVariance(SampleSize = numPairs, QValues = Q_BinMM_ContPhenos_Ghat)
    
    #calculate final stat and p value
    #R geno and IBS score
    Stat_IBS_CatPhenos = CalcStatisticPVal(SampleSize = numPairs, Variance = OrgVar_Q_IBS_CatPhenos_Ghat, UscoresR = UR_CatPhenos_Ghat, UscoreS = US_IBS_CatPhenos_Ghat, s = s)
    Stat_IBS_ContPhenos = CalcStatisticPVal(SampleSize = numPairs, Variance = OrgVar_Q_IBS_ContPhenos_Ghat, UscoresR = UR_ContPhenos_Ghat, UscoreS = US_IBS_ContPhenos_Ghat, s = s)
    #R geno and Incomp score
    Stat_Incomp_CatPhenos = CalcStatisticPVal(SampleSize = numPairs, Variance = OrgVar_Q_Incomp_CatPhenos_Ghat, UscoresR = UR_CatPhenos_Ghat, UscoreS = US_Incomp_CatPhenos_Ghat, s = s)
    Stat_Incomp_ContPhenos = CalcStatisticPVal(SampleSize = numPairs, Variance = OrgVar_Q_Incomp_ContPhenos_Ghat, UscoresR = UR_ContPhenos_Ghat, UscoreS = US_Incomp_ContPhenos_Ghat, s = s)
    #R geno and AMS Score
    Stat_AMS_CatPhenos = CalcStatisticPVal(SampleSize = numPairs, Variance = OrgVar_Q_AMS_CatPhenos_Ghat, UscoresR = UR_CatPhenos_Ghat, UscoreS = US_AMS_CatPhenos_Ghat, s = s)
    Stat_AMS_ContPhenos = CalcStatisticPVal(SampleSize = numPairs, Variance = OrgVar_Q_AMS_ContPhenos_Ghat, UscoresR = UR_ContPhenos_Ghat, UscoreS = US_AMS_ContPhenos_Ghat, s = s)
    #R geno and Bin MM score
    Stat_BinMM_CatPhenos = CalcStatisticPVal(SampleSize = numPairs, Variance = OrgVar_Q_BinMM_CatPhenos_Ghat, UscoresR = UR_CatPhenos_Ghat, UscoreS = US_BinMM_CatPhenos_Ghat, s = s)
    Stat_BinMM_ContPhenos = CalcStatisticPVal(SampleSize = numPairs, Variance = OrgVar_Q_BinMM_ContPhenos_Ghat, UscoresR = UR_ContPhenos_Ghat, UscoreS = US_BinMM_ContPhenos_Ghat, s = s)
    
    #fill columns in order
    ##Binary Phenos
    ##IBS Cat, Incomp Cat, AMS Cat, Bin MM Cat, 
    statsAndPValsMat[,1:3] = Stat_IBS_CatPhenos
    statsAndPValsMat[,4:6] = Stat_Incomp_CatPhenos
    statsAndPValsMat[,7:9] = Stat_AMS_CatPhenos
    statsAndPValsMat[,10:12] = Stat_BinMM_CatPhenos
    ##Cont Phenos
    ##IBS Cont, Incomp Cont, AMS Cont, Bin MM Cont
    statsAndPValsMat[,13:15] = Stat_IBS_ContPhenos
    statsAndPValsMat[,16:18] = Stat_Incomp_ContPhenos
    statsAndPValsMat[,19:21] = Stat_AMS_ContPhenos
    statsAndPValsMat[,22:24] = Stat_BinMM_ContPhenos
    
    print(paste0("Simulation ",ii," is complete."))
    
    statsAndPValsMat
  }, mc.cores=4)
  
  statsAndPValsAll = matrix(unlist(statsAndPValsMat),nrow = numSims, ncol = 24, byrow = TRUE)
  statsAndPValsCatPhenos = statsAndPValsAll[,1:12]
  statsAndPValsContPhenos = statsAndPValsAll[,13:24]
  
  
  #define values for naming conventions
  if(LowLD == TRUE){
    ld = "LowLD"
  } else {
    ld = "HighLD"
  }
  
  if(weightedScores == FALSE){
    weighted = ""
  } else{
    weighted = "WeightedScores"
  }
  
  if(standardizeScores == FALSE){
    standardized = ""
  } else {
    standardized = "StandardizedScores"
  }
  
  #write out the Stats and p values
  write.csv(statsAndPValsCatPhenos, file = paste0(path,"/Power_CatPhenos_",weighted,"_",standardized,"_Prev",YPrev*100,"_Ghat_StatsAndPValues_",percentageAssoc,"SNPsAssociated_",ld,"_s",s,"_TrueScore",TrueScore,"_",ORSize,"OR_assocSNPOrScore",snpOrScore,".csv"))
  write.csv(statsAndPValsContPhenos, file = paste0(path,"/Power_ContPhenos__",weighted,"_",standardized,"_Ghat_StatsAndPValues_",percentageAssoc,"SNPsAssociated_",ld,"_s",s,"_TrueScore",TrueScore,"_",ORSize,"OR_assocSNPOrScore",snpOrScore,".csv"))
  
  #calculate the power for each combination of score and r geno
  #first column for cat phenos, second for cont phenos
  #rows go IBS, Incomp, AMS, Bin MM
  power = matrix(nrow = 4, ncol = 2)
  
  #pull only the pvalues
  PValsCatPhenos = statsAndPValsCatPhenos[,c(2,5,8,11)]
  PValsContPhenos = statsAndPValsContPhenos[,c(2,5,8,11)]
  
  #calc power
  for(jj in 1:4){
    power[jj,1] = sum(PValsCatPhenos[,jj] <= 0.05)/numSims
    power[jj,2] = sum(PValsContPhenos[,jj] <= 0.05)/numSims
  }
  
  #write out power results
  write.csv(power, file = paste0(path,"/Power_Results_",percentageAssoc,"SNPsAssociated_",ld,"_",weighted,"_",standardized,"_Ghat_CatAndContPhenos_Prev",YPrev*100,"_s",s,"_TrueScore",TrueScore,"_",ORSize,"OR_assocSNPOrScore",snpOrScore,".csv"))
  write.csv(power, file = paste0(path,"/Power_Results_",percentageAssoc,"SNPsAssociated_",ld,"_",weighted,"_",standardized,"_Ghat_CatAndContPhenos_Prev",YPrev*100,"_s",s,"_TrueScore",TrueScore,"_",ORSize,"OR_assocSNPOrScore",snpOrScore,".csv"))
}

#Power Pipeline with Ghat for RSNPs being true
RunPowerPipelineLocalGhat_RSNPs = function(chr, gene, numPairs, YPrev, s, Gamma, TrueScore, ORSize, standardizeScores = FALSE, weightedScores = FALSE, scoreWeights, percentageAssoc, LowLD){
  #function to run whole power pipeline locally, uses Ghat instead of long form eqn for calculations
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #s = PVE for choosing the number of PCs
  #Gamma = effect size for score, length 1, can be 0
  #TrueScore = IBS.gene, Incomp.gene, AMS.gene, or BinMM.gene
  #ORSize = small, medium or large for effect size of associated R geno SNP
  #standardizeScores = T or F whether the scores should be standardized based on maximum score value
  #weightedScores = T or F whether the scores will be weighted
  #scoreWeights = m x 1 vector of weights, one weight for each SNP
  #### This value shows which score is potentially being used to create phenotypes
  #percentageAssoc = percentage of SNPs associated with outcome (either 5, 25, 50, 75, or 100) 
  #LowLD = True or FALSE whether the associated SNPs are in low LD or high LD
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes power values to csv files
  
  library(parallel)
  
  #need to define this for naming at the end
  snpOrScore = "RSNP"
  
  #define effect based on OR size
  if(ORSize == "small"){
    effect = 0.14
  } else if(ORSize == "medium"){
    effect = 0.41
  } else {
    effect = 0.69
  }
  
  #number of sims always the same
  numSims = 5000
  
  #define path to data
  #for sampled haplotype data
  # path = paste0("/home/vlynn/Paper_II_Sims/",gene,"_Results_",numPairs,"Pairs")
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_II_Sims/HapGen_Files/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/ProjectIISourceFunctions.R")
  
  #determine which SNPs to actually set as assoc. based on gene
  assocSNPs = DetermineAssocRSNPs(gene = gene, LowLD = LowLD, percentageAssoc = percentageAssoc)
  
  #turn output matrices into output lists
  myList = lapply(1:numSims, rep, times = 1)
  statsAndPValsMat = matrix(NA,ncol=24)
  
  statsAndPValsMat = mclapply(myList, function(ii){
    statsAndPValsMat = matrix(NA,ncol=24)
    
    #pull recipient and donor genotypes
    RGenos = obtainRGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    DGenos = obtainDGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    
    #calculate single snp scores
    IBS.snp = calcIBSMismatch(RGenosMat = RGenos, DGenosMat = DGenos)
    Incomp.snp = calcIncompatibilityScore(RGenosMat = RGenos, DGenosMat = DGenos)
    AMS.snp = calcAMS(RGenosMat = RGenos, DGenosMat = DGenos)
    BinMM.snp = calcBinaryMM(RGenosMat = RGenos, DGenosMat = DGenos)
    
    #calculate gene based scores
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
    
    #need to use TrueScore to pull gene based scores matrix for generating phenotypes
    if(TrueScore == "IBS.gene"){
      PhenoScore = IBS.gene
    } else if(TrueScore == "Incomp.gene"){
      PhenoScore = Incomp.gene
    } else if(TrueScore == "AMS.gene"){
      PhenoScore = AMS.gene
    } else {
      PhenoScore = BinMM.gene
    }
    
    #generate covariates
    #for now, a single binary and a single continous covariate
    CovData = GenCovData(SampleSize = numPairs, BinaryValues = 1, ContinuousValues = 1)
    
    #need to define null Betas for phenotype generation
    nSNP = ncol(RGenos) #this should be the number of SNPs
    nullBetas = rep(0,nSNP) #generate null beta values
    Betas = nullBetas
    
    #set assoc Betas
    #all betas have same effect for now
    for(jj in assocSNPs){
      Betas[jj] = effect
      Betas = as.matrix(Betas, ncol = 1)
    }
    
    #generate phenotypes, both continuous and binary
    #Based on single true score
    CatPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = TRUE, YPrev = YPrev,  Covariates = CovData, RGenoData = RGenos, ScoreData = PhenoScore, Betas = Betas, Gamma = Gamma)
    ContPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = FALSE,  Covariates = CovData, RGenoData = RGenos, ScoreData = PhenoScore, Betas = Betas, Gamma = Gamma)
    
    #Need to calculate UR and US scores
    #will have then for continuous and binary phenos
    #need separate US for each score type
    #Recipient genotype UR values:
    UR_CatPhenos_Ghat = CalcUScoreGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = TRUE, RGenoData = RGenos, Phenos = CatPhenos, BinPhenos = TRUE)
    UR_ContPhenos_Ghat = CalcUScoreGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = TRUE, RGenoData = RGenos, Phenos = ContPhenos, BinPhenos = FALSE)
    #US for IBS Score
    US_IBS_CatPhenos_Ghat = CalcUScoreGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = IBS.gene, Phenos = CatPhenos, BinPhenos = TRUE)
    US_IBS_ContPhenos_Ghat = CalcUScoreGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = IBS.gene, Phenos = ContPhenos, BinPhenos = FALSE)
    #US for Incomp Score
    US_Incomp_CatPhenos_Ghat = CalcUScoreGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = Incomp.gene, Phenos = CatPhenos, BinPhenos = TRUE)
    US_Incomp_ContPhenos_Ghat = CalcUScoreGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = Incomp.gene, Phenos = ContPhenos, BinPhenos = FALSE)
    #US for AMS Score
    US_AMS_CatPhenos_Ghat = CalcUScoreGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = AMS.gene, Phenos = CatPhenos, BinPhenos = TRUE)
    US_AMS_ContPhenos_Ghat = CalcUScoreGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = AMS.gene, Phenos = ContPhenos, BinPhenos = FALSE)
    #US for Bin MM Score
    US_BinMM_CatPhenos_Ghat = CalcUScoreGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = BinMM.gene, Phenos = CatPhenos, BinPhenos = TRUE)
    US_BinMM_ContPhenos_Ghat = CalcUScoreGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = BinMM.gene, Phenos = ContPhenos, BinPhenos = FALSE)
    
    #Need to calculate Q values
    #R geno QR values
    QR_CatPhenos_Ghat = CalcQValuesGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = TRUE, RGenoData = RGenos, Phenos = CatPhenos, BinPhenos = TRUE)
    QR_ContPhenos_Ghat = CalcQValuesGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = TRUE, RGenoData = RGenos, Phenos = ContPhenos, BinPhenos = FALSE)
    #QS for IBS Score
    QS_IBS_CatPhenos_Ghat = CalcQValuesGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = IBS.gene, Phenos = CatPhenos, BinPhenos = TRUE)
    QS_IBS_ContPhenos_Ghat = CalcQValuesGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = IBS.gene, Phenos = ContPhenos, BinPhenos = FALSE)
    #QS for Incomp Score
    QS_Incomp_CatPhenos_Ghat = CalcQValuesGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = Incomp.gene, Phenos = CatPhenos, BinPhenos = TRUE)
    QS_Incomp_ContPhenos_Ghat = CalcQValuesGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = Incomp.gene, Phenos = ContPhenos, BinPhenos = FALSE)
    #QS for AMS score
    QS_AMS_CatPhenos_Ghat = CalcQValuesGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = AMS.gene, Phenos = CatPhenos, BinPhenos = TRUE)
    QS_AMS_ContPhenos_Ghat = CalcQValuesGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = AMS.gene, Phenos = ContPhenos, BinPhenos = FALSE)
    #QS for Bin MM Score
    QS_BinMM_CatPhenos_Ghat = CalcQValuesGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = BinMM.gene, Phenos = CatPhenos, BinPhenos = TRUE)
    QS_BinMM_ContPhenos_Ghat = CalcQValuesGhat(SampleSize = numPairs, includeCov = TRUE, CovData = CovData, CalcUR = FALSE, ScoreData = BinMM.gene, Phenos = ContPhenos, BinPhenos = FALSE)
    
    #Create combined Qs (QR, QS)
    #R geno and IBS score
    Q_IBS_CatPhenos_Ghat = cbind(QR_CatPhenos_Ghat, QS_IBS_CatPhenos_Ghat)
    Q_IBS_ContPhenos_Ghat = cbind(QR_ContPhenos_Ghat, QS_IBS_ContPhenos_Ghat) 
    #R geno and Incomp Score
    Q_Incomp_CatPhenos_Ghat = cbind(QR_CatPhenos_Ghat, QS_Incomp_CatPhenos_Ghat)
    Q_Incomp_ContPhenos_Ghat = cbind(QR_ContPhenos_Ghat, QS_Incomp_ContPhenos_Ghat)
    #R geno and AMS score
    Q_AMS_CatPhenos_Ghat = cbind(QR_CatPhenos_Ghat, QS_AMS_CatPhenos_Ghat)
    Q_AMS_ContPhenos_Ghat = cbind(QR_ContPhenos_Ghat, QS_AMS_ContPhenos_Ghat)
    #R geno and Bin MM score
    Q_BinMM_CatPhenos_Ghat = cbind(QR_CatPhenos_Ghat, QS_BinMM_CatPhenos_Ghat)
    Q_BinMM_ContPhenos_Ghat = cbind(QR_ContPhenos_Ghat, QS_BinMM_ContPhenos_Ghat)
    
    #calculate full variance
    #each set of Qs will have a variance calculation
    #R genos + IBS Score
    OrgVar_Q_IBS_CatPhenos_Ghat = CalcVariance(SampleSize = numPairs, QValues = Q_IBS_CatPhenos_Ghat)
    OrgVar_Q_IBS_ContPhenos_Ghat = CalcVariance(SampleSize = numPairs, QValues = Q_IBS_ContPhenos_Ghat)
    #R genos + Incomp Score
    OrgVar_Q_Incomp_CatPhenos_Ghat = CalcVariance(SampleSize = numPairs, QValues = Q_Incomp_CatPhenos_Ghat)
    OrgVar_Q_Incomp_ContPhenos_Ghat = CalcVariance(SampleSize = numPairs, QValues = Q_Incomp_ContPhenos_Ghat)
    #R genos + AMS Score
    OrgVar_Q_AMS_CatPhenos_Ghat = CalcVariance(SampleSize = numPairs, QValues = Q_AMS_CatPhenos_Ghat)
    OrgVar_Q_AMS_ContPhenos_Ghat = CalcVariance(SampleSize = numPairs, QValues = Q_AMS_ContPhenos_Ghat)
    #R genos + Bin MM Score
    OrgVar_Q_BinMM_CatPhenos_Ghat = CalcVariance(SampleSize = numPairs, QValues = Q_BinMM_CatPhenos_Ghat)
    OrgVar_Q_BinMM_ContPhenos_Ghat = CalcVariance(SampleSize = numPairs, QValues = Q_BinMM_ContPhenos_Ghat)
    
    #calculate final stat and p value
    #R geno and IBS score
    Stat_IBS_CatPhenos = CalcStatisticPVal(SampleSize = numPairs, Variance = OrgVar_Q_IBS_CatPhenos_Ghat, UscoresR = UR_CatPhenos_Ghat, UscoreS = US_IBS_CatPhenos_Ghat, s = s)
    Stat_IBS_ContPhenos = CalcStatisticPVal(SampleSize = numPairs, Variance = OrgVar_Q_IBS_ContPhenos_Ghat, UscoresR = UR_ContPhenos_Ghat, UscoreS = US_IBS_ContPhenos_Ghat, s = s)
    #R geno and Incomp score
    Stat_Incomp_CatPhenos = CalcStatisticPVal(SampleSize = numPairs, Variance = OrgVar_Q_Incomp_CatPhenos_Ghat, UscoresR = UR_CatPhenos_Ghat, UscoreS = US_Incomp_CatPhenos_Ghat, s = s)
    Stat_Incomp_ContPhenos = CalcStatisticPVal(SampleSize = numPairs, Variance = OrgVar_Q_Incomp_ContPhenos_Ghat, UscoresR = UR_ContPhenos_Ghat, UscoreS = US_Incomp_ContPhenos_Ghat, s = s)
    #R geno and AMS Score
    Stat_AMS_CatPhenos = CalcStatisticPVal(SampleSize = numPairs, Variance = OrgVar_Q_AMS_CatPhenos_Ghat, UscoresR = UR_CatPhenos_Ghat, UscoreS = US_AMS_CatPhenos_Ghat, s = s)
    Stat_AMS_ContPhenos = CalcStatisticPVal(SampleSize = numPairs, Variance = OrgVar_Q_AMS_ContPhenos_Ghat, UscoresR = UR_ContPhenos_Ghat, UscoreS = US_AMS_ContPhenos_Ghat, s = s)
    #R geno and Bin MM score
    Stat_BinMM_CatPhenos = CalcStatisticPVal(SampleSize = numPairs, Variance = OrgVar_Q_BinMM_CatPhenos_Ghat, UscoresR = UR_CatPhenos_Ghat, UscoreS = US_BinMM_CatPhenos_Ghat, s = s)
    Stat_BinMM_ContPhenos = CalcStatisticPVal(SampleSize = numPairs, Variance = OrgVar_Q_BinMM_ContPhenos_Ghat, UscoresR = UR_ContPhenos_Ghat, UscoreS = US_BinMM_ContPhenos_Ghat, s = s)
    
    #fill columns in order
    ##Binary Phenos
    ##IBS Cat, Incomp Cat, AMS Cat, Bin MM Cat, 
    statsAndPValsMat[,1:3] = Stat_IBS_CatPhenos
    statsAndPValsMat[,4:6] = Stat_Incomp_CatPhenos
    statsAndPValsMat[,7:9] = Stat_AMS_CatPhenos
    statsAndPValsMat[,10:12] = Stat_BinMM_CatPhenos
    ##Cont Phenos
    ##IBS Cont, Incomp Cont, AMS Cont, Bin MM Cont
    statsAndPValsMat[,13:15] = Stat_IBS_ContPhenos
    statsAndPValsMat[,16:18] = Stat_Incomp_ContPhenos
    statsAndPValsMat[,19:21] = Stat_AMS_ContPhenos
    statsAndPValsMat[,22:24] = Stat_BinMM_ContPhenos
    
    print(paste0("Simulation ",ii," is complete."))
    
    #reset Betas
    Betas = nullBetas
    
    statsAndPValsMat
  }, mc.cores = 4)
  
  statsAndPValsAll = matrix(unlist(statsAndPValsMat),nrow = numSims, ncol = 24, byrow = TRUE)
  statsAndPValsCatPhenos = statsAndPValsAll[,1:12]
  statsAndPValsContPhenos = statsAndPValsAll[,13:24]
  
  #define values for naming conventions
  if(LowLD == TRUE){
    ld = "LowLD"
  } else {
    ld = "HighLD"
  }
  
  if(weightedScores == FALSE){
    weighted = ""
  } else{
    weighted = "WeightedScores"
  }
  
  if(standardizeScores == FALSE){
    standardized = ""
  } else {
    standardized = "StandardizedScores"
  }
  
  #write out the Stats and p values
  write.csv(statsAndPValsCatPhenos, file = paste0(path,"/Power_CatPhenos_",weighted,"_",standardized,"_Prev",YPrev*100,"_Ghat_StatsAndPValues_",percentageAssoc,"SNPsAssociated_",ld,"_s",s,"_TrueScore",TrueScore,"_",ORSize,"OR_assocSNPOrScore",snpOrScore,".csv"))
  write.csv(statsAndPValsContPhenos, file = paste0(path,"/Power_ContPhenos__",weighted,"_",standardized,"_Ghat_StatsAndPValues_",percentageAssoc,"SNPsAssociated_",ld,"_s",s,"_TrueScore",TrueScore,"_",ORSize,"OR_assocSNPOrScore",snpOrScore,".csv"))
  
  #calculate the power for each combination of score and r geno
  #first column for cat phenos, second for cont phenos
  #rows go IBS, Incomp, AMS, Bin MM
  power = matrix(nrow = 4, ncol = 2)
  
  #pull only the pvalues
  PValsCatPhenos = statsAndPValsCatPhenos[,c(2,5,8,11)]
  PValsContPhenos = statsAndPValsContPhenos[,c(2,5,8,11)]
  
  #calc power
  for(jj in 1:4){
    power[jj,1] = sum(PValsCatPhenos[,jj] <= 0.05)/numSims
    power[jj,2] = sum(PValsContPhenos[,jj] <= 0.05)/numSims
  }
  
  #write out power results
  write.csv(power, file = paste0(path,"/Power_Results_",percentageAssoc,"SNPsAssociated_",ld,"_",weighted,"_",standardized,"_Ghat_CatAndContPhenos_Prev",YPrev*100,"_s",s,"_TrueScore",TrueScore,"_",ORSize,"OR_assocSNPOrScore",snpOrScore,".csv"))
  write.csv(power, file = paste0(path,"/Power_Results_",percentageAssoc,"SNPsAssociated_",ld,"_",weighted,"_",standardized,"_Ghat_CatAndContPhenos_Prev",YPrev*100,"_s",s,"_TrueScore",TrueScore,"_",ORSize,"OR_assocSNPOrScore",snpOrScore,".csv"))
}

##################################################
## SKAT Pipelines
##################################################
#TIE Pipeline with SKAT
#N x 2m matrix instead of N x (m+1)

RunTIEPipelineSKAT = function(chr, gene, numPairs, YPrev, kernel, kernelWeights=c(), standardizeScores = FALSE, weightedScores = FALSE, scoreWeights){
  #function to run whole TIE pipeline, Calculates SKAT stat and pvalue for comparison
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #kernel = which kernel to use for SKAT (linear, linear.weighted, IBS, IBS.weighted)
  #kernelWeights = numeric vector of weights for weighted kernels
  #standardizeScores = T or F whether the scores should be standardized based on maximum score value
  #weightedScores = T or F whether the scores will be weighted
  #scoreWeights = m x 1 vector of weights, one weight for each SNP
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes TIE values to csv files
  
  #load SKAT library
  library(SKAT)
  library(parallel)
  
  #always the same
  numSims = 5000  
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_II_Sims/HapGen_Files/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/ProjectIISourceFunctions_v2.R")
  
  myList = lapply(1:numSims, rep, times = 1)
  
  statsAndPVals = mclapply(myList, function(ii){    
    #define matrix to hold all Stats and Pvalues
    statsAndPVals = matrix(NA, nrow = numSims, ncol = 16)
    
    #pull recipient and donor genotypes
    RGenos = obtainRGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    DGenos = obtainDGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    
    #calculate single snp scores
    IBS.snp = calcIBSMismatch(RGenosMat = RGenos, DGenosMat = DGenos)
    Incomp.snp = calcIncompatibilityScore(RGenosMat = RGenos, DGenosMat = DGenos)
    AMS.snp = calcAMS(RGenosMat = RGenos, DGenosMat = DGenos)
    BinMM.snp = calcBinaryMM(RGenosMat = RGenos, DGenosMat = DGenos)
    
    #calculate gene based scores
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
    
    #generate covariates
    #for now, a single binary and a single continous covariate
    CovData = GenCovData(SampleSize = numPairs, BinaryValues = 1, ContinuousValues = 1)
    
    #generate phenotypes, both continuous and binary
    CatPhenos = GenNullPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = TRUE, YPrev = YPrev,  Covariates = CovData)
    ContPhenos = GenNullPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = FALSE,  Covariates = CovData)
    
    #Combine R geno and Scores into 4 separate datasets, size: N x (m+1)
    RGeno.IBS.gene = cbind(RGenos, IBS.gene)
    RGeno.Incomp.gene = cbind(RGenos, Incomp.gene)
    RGeno.AMS.gene = cbind(RGenos, AMS.gene)
    RGeno.BinMM.gene = cbind(RGenos, BinMM.gene)
    
    ## Generate SKAT Null Models
    # formulas will be Y ~ covariates for continuous and dichotomous Y
    obj_dich=SKAT_Null_Model(CatPhenos~CovData, out_type="D")
    obj_cont=SKAT_Null_Model(ContPhenos~CovData, out_type="C")
    
    #Perform SKAT for all 8 combos of score and cont/dich outcome
    #unweighted SKAT
    if(length(kernelWeights) == 0){
      Stat_IBS_CatPhenos = SKAT(RGeno.IBS.gene, obj_dich, kernel = kernel, is_check_genotype = FALSE)
      Stat_Incomp_CatPhenos = SKAT(RGeno.Incomp.gene, obj_dich, kernel = kernel, is_check_genotype = FALSE)
      Stat_AMS_CatPhenos = SKAT(RGeno.AMS.gene, obj_dich, kernel = kernel, is_check_genotype = FALSE)
      Stat_BinMM_CatPhenos = SKAT(RGeno.BinMM.gene, obj_dich, kernel = kernel, is_check_genotype = FALSE)
      
      Stat_IBS_ContPhenos = SKAT(RGeno.IBS.gene, obj_cont, kernel = kernel, is_check_genotype = FALSE)
      Stat_Incomp_ContPhenos = SKAT(RGeno.Incomp.gene, obj_cont, kernel = kernel, is_check_genotype = FALSE)
      Stat_AMS_ContPhenos = SKAT(RGeno.AMS.gene, obj_cont, kernel = kernel, is_check_genotype = FALSE)
      Stat_BinMM_ContPhenos = SKAT(RGeno.BinMM.gene, obj_cont, kernel = kernel, is_check_genotype = FALSE)
    } else {#weighted SKAT
      Stat_IBS_CatPhenos = SKAT(RGeno.IBS.gene, obj_dich, kernel = kernel, weights = kernelWeights, is_check_genotype = FALSE)
      Stat_Incomp_CatPhenos = SKAT(RGeno.Incomp.gene, obj_dich, kernel = kernel, weights = kernelWeights, is_check_genotype = FALSE)
      Stat_AMS_CatPhenos = SKAT(RGeno.AMS.gene, obj_dich, kernel = kernel, weights = kernelWeights, is_check_genotype = FALSE)
      Stat_BinMM_CatPhenos = SKAT(RGeno.BinMM.gene, obj_dich, kernel = kernel, weights = kernelWeights, is_check_genotype = FALSE)
      
      Stat_IBS_ContPhenos = SKAT(RGeno.IBS.gene, obj_cont, kernel = kernel, weights = kernelWeights, is_check_genotype = FALSE)
      Stat_Incomp_ContPhenos = SKAT(RGeno.Incomp.gene, obj_cont, kernel = kernel, weights = kernelWeights, is_check_genotype = FALSE)
      Stat_AMS_ContPhenos = SKAT(RGeno.AMS.gene, obj_cont, kernel = kernel, weights = kernelWeights, is_check_genotype = FALSE)
      Stat_BinMM_ContPhenos = SKAT(RGeno.BinMM.gene, obj_cont, kernel = kernel, weights = kernelWeights, is_check_genotype = FALSE)
    }	
    
    #fill columns in order
    ##Binary Phenos first, then cont
    ##IBS, Incomp, AMS, Bin MM
    IBSCat = c(Stat_IBS_CatPhenos$Q, Stat_IBS_CatPhenos$p.value)
    IncompCat = c(Stat_Incomp_CatPhenos$Q, Stat_Incomp_CatPhenos$p.value)
    AMSCat = c(Stat_AMS_CatPhenos$Q, Stat_AMS_CatPhenos$p.value)
    BinMMCat = c(Stat_BinMM_CatPhenos$Q, Stat_BinMM_CatPhenos$p.value)
    ##Cont Phenos
    ##IBS Cont, Incomp Cont, AMS Cont, Bin MM Cont
    IBSCont = c(Stat_IBS_ContPhenos$Q, Stat_IBS_ContPhenos$p.value)
    IncompCont = c(Stat_Incomp_ContPhenos$Q, Stat_Incomp_ContPhenos$p.value)
    AMSCont = c(Stat_AMS_ContPhenos$Q, Stat_AMS_ContPhenos$p.value)
    BinMMCont = c(Stat_BinMM_ContPhenos$Q, Stat_BinMM_ContPhenos$p.value)
    
    IBSCat.mat = as.matrix(IBSCat)
    IncompCat.mat = as.matrix(IncompCat)
    AMSCat.mat = as.matrix(AMSCat)
    BinMMCat.mat = as.matrix(BinMMCat)
    
    IBSCont.mat = as.matrix(IBSCont)
    IncompCont.mat = as.matrix(IncompCont)
    AMSCont.mat = as.matrix(AMSCont)
    BinMMCont.mat = as.matrix(BinMMCont)
    
    statsAndPVals = cbind(t(IBSCat.mat),t(IncompCat.mat),t(AMSCat.mat),t(BinMMCat.mat),t(IBSCont.mat),t(IncompCont.mat),t(AMSCont.mat),t(BinMMCont.mat))
    
    print(paste0("Simulation ",ii," is complete."))
    
    statsAndPVals
  }, mc.cores=4)
  
  statsAndPValsAll = matrix(unlist(statsAndPVals),nrow = numSims, ncol = 16, byrow = TRUE)
  
  #separate into cat and cont phenos
  statsAndPValsCatPhenos = statsAndPValsAll[,1:8]
  statsAndPValsContPhenos = statsAndPValsAll[,9:16]
  
  #write out the Stats and p values
  if(weightedScores == FALSE){
    if(standardizeScores == FALSE){
      write.csv(statsAndPValsCatPhenos, file = paste0(path,"/TIE_mPlus1_CatPhenos_Prev",YPrev*100,"_SKAT_kernel",kernel,"_StatsAndPValues.csv"))
      write.csv(statsAndPValsContPhenos, file = paste0(path,"/TIE_mPlus1_ContPhenos_SKAT_kernel",kernel,"_StatsAndPValues.csv"))
    } else { #scores are standardized
      write.csv(statsAndPValsCatPhenos, file = paste0(path,"/TIE_mPlus1_CatPhenos_Prev",YPrev*100,"_SKAT_kernel",kernel,"_StandardizedScores_StatsAndPValues.csv"))
      write.csv(statsAndPValsContPhenos, file = paste0(path,"/TIE_mPlus1_ContPhenos_SKAT_kernel",kernel,"_StandardizedScores_StatsAndPValues.csv"))
    }
  } else { #scores are weighted
    write.csv(statsAndPValsCatPhenos, file = paste0(path,"/TIE_mPlus1_CatPhenos_Prev",YPrev*100,"_SKAT_kernel",kernel,"_WeightedScores_StatsAndPValues.csv"))
    write.csv(statsAndPValsContPhenos, file = paste0(path,"/TIE_mPlus1_ContPhenos_SKAT_kernel",kernel,"_WeightedScores_StatsAndPValues.csv"))
  }
  
  #calculate the type I errors for each combination of score and r geno
  #first column for cat phenos, second for cont phenos
  #rows go IBS, Incomp, AMS, Bin MM
  ties = matrix(nrow = 4, ncol = 2)
  
  #pull only the pvalues
  PValsCatPhenos = statsAndPValsCatPhenos[,c(2,4,6,8)]
  PValsContPhenos = statsAndPValsContPhenos[,c(2,4,6,8)]
  
  #calc type I error
  for(jj in 1:4){
    ties[jj,1] = sum(PValsCatPhenos[,jj] <= 0.05)/numSims
    ties[jj,2] = sum(PValsContPhenos[,jj] <= 0.05)/numSims
  }
  
  #write out type I error results
  if(weightedScores == FALSE){
    if(standardizeScores == FALSE){
      write.csv(ties, file = paste0(path,"/TIE_Results_mPlus1_SKAT_kernel",kernel,"_CatAndContPhenos_Prev",YPrev*100,".csv"))
    } else { #scores are standardized
      write.csv(ties, file = paste0(path,"/TIE_Results_mPlus1_SKAT_kernel",kernel,"_StandardizedScores_CatAndContPhenos_Prev",YPrev*100,".csv"))
    }
  } else { #scores are weighted 
    write.csv(ties, file = paste0(path,"/TIE_Results_mPlus1_SKAT_kernel",kernel,"_WeightedScores_CatAndContPhenos_Prev",YPrev*100,".csv"))
  }
}

#Power Pipeline with SKAT for one of the Scores being true
# N x (m+1) instead of N x 2m
RunPowerPipelineSKAT_Scores = function(chr, gene, numPairs, YPrev, kernel, kernelWeights=c(), Gamma, TrueScore, ORSize, standardizeScores = FALSE, weightedScores = FALSE, scoreWeights, percentageAssoc, LowLD){
  #function to run whole power pipeline , Calculates SKAT stat and pvalue for comparison
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #kernel = which kernel to use for SKAT (linear, linear.weighted, IBS, IBS.weighted)
  #kernelWeights = numeric vector of weights for weighted kernels
  #Gamma = effect size for score, length 1, can be 0
  #TrueScore = IBS.gene, Incomp.gene, AMS.gene, or BinMM.gene
  #ORSize = Small, Medium, or Large for what OR was used for the associated SNP/score
  #standardizeScores = T or F whether the scores should be standardized based on maximum score value
  #weightedScores = T or F whether the scores will be weighted
  #scoreWeights = m x 1 vector of weights, one weight for each SNP
  #percentageAssoc = percentage of SNPs associated with outcome (either 5, 25, 50, 75, or 100) 
  #LowLD = TRUE or FALSE whether the associated SNPs are in low LD or high LD
  #### This value shows which score is potentially being used to create phenotypes
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes power values to csv files
  
  #load SKAT library
  library(SKAT)
  library(parallel)
  
  #need to define this for naming at the end
  snpOrScore = TrueScore
  
  #number of sims always the same
  numSims = 5000
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_II_Sims/HapGen_Files/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/ProjectIISourceFunctions_v2.R")
  
  #turn output matrices into output lists
  myList = lapply(1:numSims, rep, times = 1)
  statsAndPValsMat = matrix(NA,ncol=16)
  
  statsAndPValsMat = mclapply(myList, function(ii){
    statsAndPValsMat = matrix(NA,ncol=16)
    
    #pull recipient and donor genotypes
    RGenos = obtainRGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    DGenos = obtainDGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    
    #calculate single snp scores
    IBS.snp = calcIBSMismatch(RGenosMat = RGenos, DGenosMat = DGenos)
    Incomp.snp = calcIncompatibilityScore(RGenosMat = RGenos, DGenosMat = DGenos)
    AMS.snp = calcAMS(RGenosMat = RGenos, DGenosMat = DGenos)
    BinMM.snp = calcBinaryMM(RGenosMat = RGenos, DGenosMat = DGenos)
    
    #calc gene based score for all SNPs
    IBS.gene = calcGeneScore(SingleSNPKernel = IBS.snp, standardize = standardizeScores, useWeights = FALSE)
    Incomp.gene = calcGeneScore(SingleSNPKernel = Incomp.snp, standardize = standardizeScores, useWeights = FALSE)
    AMS.gene = calcGeneScore(SingleSNPKernel = AMS.snp, standardize = standardizeScores, useWeights = FALSE)
    BinMM.gene = calcGeneScore(SingleSNPKernel = BinMM.snp, standardize = standardizeScores, useWeights = FALSE)
    
    #also need to calculate gene based scores if not all SNPs are associated
    IBS.gene.PercentOfSNPs = calcGeneScorePercentOfSNPs(SingleSNPKernel = IBS.snp, gene = gene, percentageAssoc = percentageAssoc, LowLD = LowLD, standardize = FALSE, useWeights = FALSE)
    Incomp.gene.PercentOfSNPs = calcGeneScorePercentOfSNPs(SingleSNPKernel = Incomp.snp, gene =  gene, percentageAssoc = percentageAssoc, LowLD = LowLD, standardize = FALSE, useWeights = FALSE)
    AMS.gene.PercentOfSNPs = calcGeneScorePercentOfSNPs(SingleSNPKernel = AMS.snp, gene =  gene, percentageAssoc = percentageAssoc, LowLD = LowLD, standardize = FALSE, useWeights = FALSE)
    BinMM.gene.PercentOfSNPs = calcGeneScorePercentOfSNPs(SingleSNPKernel = BinMM.snp, gene =  gene, percentageAssoc = percentageAssoc, LowLD = LowLD, standardize = FALSE, useWeights = FALSE)
    
    #need to use TrueScore to pull gene based scores matrix for generating phenotypes
    if(TrueScore == "IBS.gene"){
      PhenoScore = IBS.gene.PercentOfSNPs
    } else if(TrueScore == "Incomp.gene"){
      PhenoScore = Incomp.gene.PercentOfSNPs
    } else if(TrueScore == "AMS.gene"){
      PhenoScore = AMS.gene.PercentOfSNPs
    } else {
      PhenoScore = BinMM.gene.PercentOfSNPs
    }
    
    #generate covariates
    #for now, a single binary and a single continous covariate
    CovData = GenCovData(SampleSize = numPairs, BinaryValues = 1, ContinuousValues = 1)
    
    #need to define null Betas for phenotype generation
    nSNP = ncol(RGenos) #this should be the number of SNPs
    Betas = rep(0,nSNP) #generate null beta values
    Betas = as.matrix(Betas, ncol = 1)
    
    #generate phenotypes, both continuous and binary
    #Based on single true score
    CatPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = TRUE, YPrev = YPrev,  Covariates = CovData, RGenoData = RGenos, ScoreData = PhenoScore, Betas = Betas, Gamma = Gamma)
    ContPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = FALSE,  Covariates = CovData, RGenoData = RGenos, ScoreData = PhenoScore, Betas = Betas, Gamma = Gamma)
    
    #Combine R geno and Scores into 4 separate datasets, size: N x (m+1)
    RGeno.IBS.snp = cbind(RGenos, IBS.gene)
    RGeno.Incomp.snp = cbind(RGenos, Incomp.gene)
    RGeno.AMS.snp = cbind(RGenos, AMS.gene)
    RGeno.BinMM.snp = cbind(RGenos, BinMM.gene)
    
    ## Generate SKAT Null Models
    # formulas will be Y ~ covariates for continuous and dichotomous Y
    obj_dich=SKAT_Null_Model(CatPhenos~CovData, out_type="D", Adjustment = FALSE)
    obj_cont=SKAT_Null_Model(ContPhenos~CovData, out_type="C", Adjustment = FALSE)
    
    #Perform SKAT for all 8 combos of score and cont/dich outcome
    #unweighted SKAT
    if(length(kernelWeights) == 0){
      Stat_IBS_CatPhenos = SKAT(RGeno.IBS.snp, obj_dich, kernel = kernel, is_check_genotype = FALSE)
      Stat_Incomp_CatPhenos = SKAT(RGeno.Incomp.snp, obj_dich, kernel = kernel, is_check_genotype = FALSE)
      Stat_AMS_CatPhenos = SKAT(RGeno.AMS.snp, obj_dich, kernel = kernel, is_check_genotype = FALSE)
      Stat_BinMM_CatPhenos = SKAT(RGeno.BinMM.snp, obj_dich, kernel = kernel, is_check_genotype = FALSE)
      
      Stat_IBS_ContPhenos = SKAT(RGeno.IBS.snp, obj_cont, kernel = kernel, is_check_genotype = FALSE)
      Stat_Incomp_ContPhenos = SKAT(RGeno.Incomp.snp, obj_cont, kernel = kernel, is_check_genotype = FALSE)
      Stat_AMS_ContPhenos = SKAT(RGeno.AMS.snp, obj_cont, kernel = kernel, is_check_genotype = FALSE)
      Stat_BinMM_ContPhenos = SKAT(RGeno.BinMM.snp, obj_cont, kernel = kernel, is_check_genotype = FALSE)
    } else {
      Stat_IBS_CatPhenos = SKAT(RGeno.IBS.snp, obj_dich, kernel = kernel, weights = kernelWeights, is_check_genotype = FALSE)
      Stat_Incomp_CatPhenos = SKAT(RGeno.Incomp.snp, obj_dich, kernel = kernel, weights = kernelWeights, is_check_genotype = FALSE)
      Stat_AMS_CatPhenos = SKAT(RGeno.AMS.snp, obj_dich, kernel = kernel, weights = kernelWeights, is_check_genotype = FALSE)
      Stat_BinMM_CatPhenos = SKAT(RGeno.BinMM.snp, obj_dich, kernel = kernel, weights = kernelWeights, is_check_genotype = FALSE)
      
      Stat_IBS_ContPhenos = SKAT(RGeno.IBS.snp, obj_cont, kernel = kernel, weights = kernelWeights, is_check_genotype = FALSE)
      Stat_Incomp_ContPhenos = SKAT(RGeno.Incomp.snp, obj_cont, kernel = kernel, weights = kernelWeights, is_check_genotype = FALSE)
      Stat_AMS_ContPhenos = SKAT(RGeno.AMS.snp, obj_cont, kernel = kernel, weights = kernelWeights, is_check_genotype = FALSE)
      Stat_BinMM_ContPhenos = SKAT(RGeno.BinMM.snp, obj_cont, kernel = kernel, weights = kernelWeights, is_check_genotype = FALSE)
    }
    #fill columns in order
    ##Binary Phenos
    ##IBS Cat, Incomp Cat, AMS Cat, Bin MM Cat, 
    IBSCat = c(Stat_IBS_CatPhenos$Q, Stat_IBS_CatPhenos$p.value)
    IncompCat = c(Stat_Incomp_CatPhenos$Q, Stat_Incomp_CatPhenos$p.value)
    AMSCat = c(Stat_AMS_CatPhenos$Q, Stat_AMS_CatPhenos$p.value)
    BinMMCat = c(Stat_BinMM_CatPhenos$Q, Stat_BinMM_CatPhenos$p.value)
    ##Cont Phenos
    ##IBS Cont, Incomp Cont, AMS Cont, Bin MM Cont
    IBSCont = c(Stat_IBS_ContPhenos$Q, Stat_IBS_ContPhenos$p.value)
    IncompCont = c(Stat_Incomp_ContPhenos$Q, Stat_Incomp_ContPhenos$p.value)
    AMSCont = c(Stat_AMS_ContPhenos$Q, Stat_AMS_ContPhenos$p.value)
    BinMMCont = c(Stat_BinMM_ContPhenos$Q, Stat_BinMM_ContPhenos$p.value)
    
    statsAndPValsMat = cbind(IBSCat,IncompCat,AMSCat,BinMMCat,IBSCont,IncompCont,AMSCont,BinMMCont)
    
    print(paste0("Simulation ",ii," is complete."))
    
    statsAndPValsMat
  }, mc.cores = 4)
  
  statsAndPValsAll = matrix(unlist(statsAndPValsMat),nrow = numSims, ncol = 16, byrow = TRUE)
  
  #separate into cat and cont phenos
  statsAndPValsCatPhenos = statsAndPValsAll[,1:8]
  statsAndPValsContPhenos = statsAndPValsAll[,9:16] 
  
  #define values for naming conventions
  if(LowLD == TRUE){
    ld = "LowLD"
  } else {
    ld = "HighLD"
  }
  
  if(weightedScores == FALSE){
    weighted = ""
  } else{
    weighted = "WeightedScores"
  }
  
  if(standardizeScores == FALSE){
    standardized = ""
  } else {
    standardized = "StandardizedScores"
  }
  
  #write out the Stats and p values
  write.csv(statsAndPValsCatPhenos, file = paste0(path,"/Power_CatPhenos_",weighted,"_",standardized,"_Prev",YPrev*100,"_mPlus1_SKAT_kernel",kernel,"_StatsAndPValues_",percentageAssoc,"SNPsAssociated_",ld,"_TrueScore",TrueScore,"_",ORSize,"OR_assocSNPOrScore",snpOrScore,".csv"))
  write.csv(statsAndPValsContPhenos, file = paste0(path,"/Power_ContPhenos_",weighted,"_",standardized,"_mPlus1_SKAT_kernel",kernel,"_StatsAndPValues_",percentageAssoc,"SNPsAssociated_",ld,"_TrueScore",TrueScore,"_",ORSize,"OR_assocSNPOrScore",snpOrScore,".csv"))
  
  #calculate the power for each combination of score and r geno
  #first column for cat phenos, second for cont phenos
  #rows go IBS, Incomp, AMS, Bin MM
  power = matrix(nrow = 4, ncol = 2)
  
  #pull only the pvalues
  PValsCatPhenos = statsAndPValsCatPhenos[,c(2,4,6,8)]
  PValsContPhenos = statsAndPValsContPhenos[,c(2,4,6,8)]
  
  #calc power
  for(jj in 1:4){
    power[jj,1] = sum(PValsCatPhenos[,jj] <= 0.05)/numSims
    power[jj,2] = sum(PValsContPhenos[,jj] <= 0.05)/numSims
  }
  
  #write out power results
  write.csv(power, file = paste0(path,"/Power_Results_",percentageAssoc,"SNPsAssociated_",ld,"_",weighted,"_",standardized,"_mPlus1_SKAT_kernel",kernel,"_CatAndContPhenos_Prev",YPrev*100,"_TrueScore",TrueScore,"_",ORSize,"OR_assocSNPOrScore",snpOrScore,".csv"))
}

#Power Pipeline with SKAT for one of the RSNPs being true
# N x (m+1) instead of N x 2m
RunPowerPipelineSKAT_RSNPs = function(chr, gene, numPairs, YPrev, kernel, kernelWeights=c(), Gamma, TrueScore, ORSize, standardizeScores = FALSE, weightedScores = FALSE, scoreWeights, percentageAssoc, LowLD){
  #function to run whole power pipeline , Calculates SKAT stat and pvalue for comparison
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #kernel = which kernel to use for SKAT (linear, linear.weighted, IBS, IBS.weighted)
  #kernelWeights = numeric vector of weights for weighted kernels
  #Gamma = effect size for score, length 1, can be 0
  #TrueScore = IBS.gene, Incomp.gene, AMS.gene, or BinMM.gene
  #ORSize = Small, Medium, or Large for what OR was used for the associated SNP/score
  #standardizeScores = T or F whether the scores should be standardized based on maximum score value
  #weightedScores = T or F whether the scores will be weighted
  #scoreWeights = m x 1 vector of weights, one weight for each SNP
  #percentageAssoc = percentage of SNPs associated with outcome (either 5, 25, 50, 75, or 100) 
  #LowLD = TRUE or FALSE whether the associated SNPs are in low LD or high LD
  #### This value shows which score is potentially being used to create phenotypes
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes power values to csv files
  
  #load libraries
  library(SKAT)
  library(parallel)
  
  #need to define this for naming at the end
  snpOrScore = "RSNP"
  
  #define effect based on OR size
  if(ORSize == "small"){
    effect = 0.14
  } else if(ORSize == "medium"){
    effect = 0.41
  } else {
    effect = 0.69
  }
  
  #number of sims always the same
  numSims = 5000
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_II_Sims/HapGen_Files/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/ProjectIISourceFunctions.R")
  
  #determine which SNPs to actually set as assoc. based on gene
  assocSNPs = DetermineAssocRSNPs(gene = gene, LowLD = LowLD, percentageAssoc = percentageAssoc)
  
  #turn output matrices into output lists
  myList = lapply(1:numSims, rep, times = 1)
  statsAndPValsMat = matrix(NA,ncol=16)
  
  statsAndPValsMat = mclapply(myList, function(ii){
    statsAndPValsMat = matrix(NA,ncol=16)
    
    #pull recipient and donor genotypes
    RGenos = obtainRGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    DGenos = obtainDGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    
    #calculate single snp scores
    IBS.snp = calcIBSMismatch(RGenosMat = RGenos, DGenosMat = DGenos)
    Incomp.snp = calcIncompatibilityScore(RGenosMat = RGenos, DGenosMat = DGenos)
    AMS.snp = calcAMS(RGenosMat = RGenos, DGenosMat = DGenos)
    BinMM.snp = calcBinaryMM(RGenosMat = RGenos, DGenosMat = DGenos)
    
    #calc gene based score for all SNPs
    IBS.gene = calcGeneScore(SingleSNPKernel = IBS.snp, standardize = standardizeScores, useWeights = FALSE)
    Incomp.gene = calcGeneScore(SingleSNPKernel = Incomp.snp, standardize = standardizeScores, useWeights = FALSE)
    AMS.gene = calcGeneScore(SingleSNPKernel = AMS.snp, standardize = standardizeScores, useWeights = FALSE)
    BinMM.gene = calcGeneScore(SingleSNPKernel = BinMM.snp, standardize = standardizeScores, useWeights = FALSE)
    
    #need to use TrueScore to pull gene based scores matrix for generating phenotypes
    if(TrueScore == "IBS.gene"){
      PhenoScore = IBS.gene
    } else if(TrueScore == "Incomp.gene"){
      PhenoScore = Incomp.gene
    } else if(TrueScore == "AMS.gene"){
      PhenoScore = AMS.gene
    } else {
      PhenoScore = BinMM.gene
    }
    
    #generate covariates
    #for now, a single binary and a single continous covariate
    CovData = GenCovData(SampleSize = numPairs, BinaryValues = 1, ContinuousValues = 1)
    
    #need to define null Betas for phenotype generation
    nSNP = ncol(RGenos) #this should be the number of SNPs
    nullBetas = rep(0,nSNP) #generate null beta values
    Betas = nullBetas
    
    #set assoc Betas
    #all betas have same effect for now
    for(jj in assocSNPs){
      Betas[jj] = effect
    }
    Betas = as.matrix(Betas, ncol = 1)
    
    #generate phenotypes, both continuous and binary
    #Based on single true score
    CatPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = TRUE, YPrev = YPrev,  Covariates = CovData, RGenoData = RGenos, ScoreData = PhenoScore, Betas = Betas, Gamma = Gamma)
    ContPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = FALSE,  Covariates = CovData, RGenoData = RGenos, ScoreData = PhenoScore, Betas = Betas, Gamma = Gamma)
    
    #Combine R geno and Scores into 4 separate datasets, size: N x (m+1)
    RGeno.IBS.snp = cbind(RGenos, IBS.gene)
    RGeno.Incomp.snp = cbind(RGenos, Incomp.gene)
    RGeno.AMS.snp = cbind(RGenos, AMS.gene)
    RGeno.BinMM.snp = cbind(RGenos, BinMM.gene)
    
    ## Generate SKAT Null Models
    # formulas will be Y ~ covariates for continuous and dichotomous Y
    obj_dich=SKAT_Null_Model(CatPhenos~CovData, out_type="D", Adjustment = FALSE)
    obj_cont=SKAT_Null_Model(ContPhenos~CovData, out_type="C", Adjustment = FALSE)
    
    #Perform SKAT for all 8 combos of score and cont/dich outcome
    #unweighted SKAT
    if(length(kernelWeights) == 0){
      Stat_IBS_CatPhenos = SKAT(RGeno.IBS.snp, obj_dich, kernel = kernel, is_check_genotype = FALSE)
      Stat_Incomp_CatPhenos = SKAT(RGeno.Incomp.snp, obj_dich, kernel = kernel, is_check_genotype = FALSE)
      Stat_AMS_CatPhenos = SKAT(RGeno.AMS.snp, obj_dich, kernel = kernel, is_check_genotype = FALSE)
      Stat_BinMM_CatPhenos = SKAT(RGeno.BinMM.snp, obj_dich, kernel = kernel, is_check_genotype = FALSE)
      
      Stat_IBS_ContPhenos = SKAT(RGeno.IBS.snp, obj_cont, kernel = kernel, is_check_genotype = FALSE)
      Stat_Incomp_ContPhenos = SKAT(RGeno.Incomp.snp, obj_cont, kernel = kernel, is_check_genotype = FALSE)
      Stat_AMS_ContPhenos = SKAT(RGeno.AMS.snp, obj_cont, kernel = kernel, is_check_genotype = FALSE)
      Stat_BinMM_ContPhenos = SKAT(RGeno.BinMM.snp, obj_cont, kernel = kernel, is_check_genotype = FALSE)
    } else {
      Stat_IBS_CatPhenos = SKAT(RGeno.IBS.snp, obj_dich, kernel = kernel, weights = kernelWeights, is_check_genotype = FALSE)
      Stat_Incomp_CatPhenos = SKAT(RGeno.Incomp.snp, obj_dich, kernel = kernel, weights = kernelWeights, is_check_genotype = FALSE)
      Stat_AMS_CatPhenos = SKAT(RGeno.AMS.snp, obj_dich, kernel = kernel, weights = kernelWeights, is_check_genotype = FALSE)
      Stat_BinMM_CatPhenos = SKAT(RGeno.BinMM.snp, obj_dich, kernel = kernel, weights = kernelWeights, is_check_genotype = FALSE)
      
      Stat_IBS_ContPhenos = SKAT(RGeno.IBS.snp, obj_cont, kernel = kernel, weights = kernelWeights, is_check_genotype = FALSE)
      Stat_Incomp_ContPhenos = SKAT(RGeno.Incomp.snp, obj_cont, kernel = kernel, weights = kernelWeights, is_check_genotype = FALSE)
      Stat_AMS_ContPhenos = SKAT(RGeno.AMS.snp, obj_cont, kernel = kernel, weights = kernelWeights, is_check_genotype = FALSE)
      Stat_BinMM_ContPhenos = SKAT(RGeno.BinMM.snp, obj_cont, kernel = kernel, weights = kernelWeights, is_check_genotype = FALSE)
    }
    
    #fill columns in order
    ##Binary Phenos
    ##IBS Cat, Incomp Cat, AMS Cat, Bin MM Cat, 
    IBSCat = c(Stat_IBS_CatPhenos$Q, Stat_IBS_CatPhenos$p.value)
    IncompCat = c(Stat_Incomp_CatPhenos$Q, Stat_Incomp_CatPhenos$p.value)
    AMSCat = c(Stat_AMS_CatPhenos$Q, Stat_AMS_CatPhenos$p.value)
    BinMMCat = c(Stat_BinMM_CatPhenos$Q, Stat_BinMM_CatPhenos$p.value)
    ##Cont Phenos
    ##IBS Cont, Incomp Cont, AMS Cont, Bin MM Cont
    IBSCont = c(Stat_IBS_ContPhenos$Q, Stat_IBS_ContPhenos$p.value)
    IncompCont = c(Stat_Incomp_ContPhenos$Q, Stat_Incomp_ContPhenos$p.value)
    AMSCont = c(Stat_AMS_ContPhenos$Q, Stat_AMS_ContPhenos$p.value)
    BinMMCont = c(Stat_BinMM_ContPhenos$Q, Stat_BinMM_ContPhenos$p.value)
    
    statsAndPValsMat = cbind(IBSCat,IncompCat,AMSCat,BinMMCat,IBSCont,IncompCont,AMSCont,BinMMCont)
    
    print(paste0("Simulation ",ii," is complete."))
    
    #reset Betas
    Betas = nullBetas
    
    statsAndPValsMat
  }, mc.cores = 4)
  
  statsAndPValsAll = matrix(unlist(statsAndPValsMat),nrow = numSims, ncol = 16, byrow = TRUE)
  
  #separate into cat and cont phenos
  statsAndPValsCatPhenos = statsAndPValsAll[,1:8]
  statsAndPValsContPhenos = statsAndPValsAll[,9:16] 
  
  #define values for naming conventions
  if(LowLD == TRUE){
    ld = "LowLD"
  } else {
    ld = "HighLD"
  }
  
  if(weightedScores == FALSE){
    weighted = ""
  } else{
    weighted = "WeightedScores"
  }
  
  if(standardizeScores == FALSE){
    standardized = ""
  } else {
    standardized = "StandardizedScores"
  }
  
  #write out the Stats and p values
  write.csv(statsAndPValsCatPhenos, file = paste0(path,"/Power_CatPhenos_",weighted,"_",standardized,"_Prev",YPrev*100,"_mPlus1_SKAT_kernel",kernel,"_StatsAndPValues_",percentageAssoc,"SNPsAssociated_",ld,"_TrueScore",TrueScore,"_",ORSize,"OR_assocSNPOrScore",snpOrScore,".csv"))
  write.csv(statsAndPValsContPhenos, file = paste0(path,"/Power_ContPhenos_",weighted,"_",standardized,"_mPlus1_SKAT_kernel",kernel,"_StatsAndPValues_",percentageAssoc,"SNPsAssociated_",ld,"_TrueScore",TrueScore,"_",ORSize,"OR_assocSNPOrScore",snpOrScore,".csv"))
  
  #calculate the power for each combination of score and r geno
  #first column for cat phenos, second for cont phenos
  #rows go IBS, Incomp, AMS, Bin MM
  power = matrix(nrow = 4, ncol = 2)
  
  #pull only the pvalues
  PValsCatPhenos = statsAndPValsCatPhenos[,c(2,4,6,8)]
  PValsContPhenos = statsAndPValsContPhenos[,c(2,4,6,8)]
  
  #calc power
  for(jj in 1:4){
    power[jj,1] = sum(PValsCatPhenos[,jj] <= 0.05)/numSims
    power[jj,2] = sum(PValsContPhenos[,jj] <= 0.05)/numSims
  }
  
  #write out power results
  write.csv(power, file = paste0(path,"/Power_Results_",percentageAssoc,"SNPsAssociated_",ld,"_",weighted,"_",standardized,"_mPlus1_SKAT_kernel",kernel,"_CatAndContPhenos_Prev",YPrev*100,"_TrueScore",TrueScore,"_",ORSize,"OR_assocSNPOrScore",snpOrScore,".csv"))
}

#########################################################################################################################################################
## These two functions were for testing whether standardizing the gene score to between 0-1 would improve power for the SKAT with linear kernel
RunPowerPipelineSKAT_RSNPs_standardized = function(chr, gene, numPairs, YPrev, kernel, kernelWeights=c(), Gamma, TrueScore, ORSize, standardizeScores = FALSE, weightedScores = FALSE, scoreWeights, percentageAssoc, LowLD){
  #function to run whole power pipeline , Calculates SKAT stat and pvalue for comparison
  #Standardizes the gene score to between 0-1 
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #kernel = which kernel to use for SKAT (linear, linear.weighted, IBS, IBS.weighted)
  #kernelWeights = numeric vector of weights for weighted kernels
  #Gamma = effect size for score, length 1, can be 0
  #TrueScore = IBS.gene, Incomp.gene, AMS.gene, or BinMM.gene
  #ORSize = Small, Medium, or Large for what OR was used for the associated SNP/score
  #standardizeScores = T or F whether the scores should be standardized based on maximum score value
  #weightedScores = T or F whether the scores will be weighted
  #scoreWeights = m x 1 vector of weights, one weight for each SNP
  #percentageAssoc = percentage of SNPs associated with outcome (either 5, 25, 50, 75, or 100) 
  #LowLD = TRUE or FALSE whether the associated SNPs are in low LD or high LD
  #### This value shows which score is potentially being used to create phenotypes
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes power values to csv files
  
  #load libraries
  library(SKAT)
  library(parallel)
  
  #need to define this for naming at the end
  snpOrScore = "RSNP"
  
  #define effect based on OR size
  if(ORSize == "small"){
    effect = 0.14
  } else if(ORSize == "medium"){
    effect = 0.41
  } else {
    effect = 0.69
  }
  
  #number of sims always the same
  numSims = 5000
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_II_Sims/HapGen_Files/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/ProjectIISourceFunctions.R")
  
  #determine which SNPs to actually set as assoc. based on gene
  assocSNPs = DetermineAssocRSNPs(gene = gene, LowLD = LowLD, percentageAssoc = percentageAssoc)
  
  #turn output matrices into output lists
  myList = lapply(1:numSims, rep, times = 1)
  statsAndPValsMat = matrix(NA,ncol=16)
  
  statsAndPValsMat = mclapply(myList, function(ii){
    statsAndPValsMat = matrix(NA,ncol=16)
    
    #pull recipient and donor genotypes
    RGenos = obtainRGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    DGenos = obtainDGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    
    #calculate single snp scores
    IBS.snp = calcIBSMismatch(RGenosMat = RGenos, DGenosMat = DGenos)
    Incomp.snp = calcIncompatibilityScore(RGenosMat = RGenos, DGenosMat = DGenos)
    AMS.snp = calcAMS(RGenosMat = RGenos, DGenosMat = DGenos)
    BinMM.snp = calcBinaryMM(RGenosMat = RGenos, DGenosMat = DGenos)
    
    #calc gene based score for all SNPs
    IBS.gene = calcGeneScore(SingleSNPKernel = IBS.snp, standardize = FALSE, useWeights = FALSE)
    Incomp.gene = calcGeneScore(SingleSNPKernel = Incomp.snp, standardize = FALSE, useWeights = FALSE)
    AMS.gene = calcGeneScore(SingleSNPKernel = AMS.snp, standardize = FALSE, useWeights = FALSE)
    BinMM.gene = calcGeneScore(SingleSNPKernel = BinMM.snp, standardize = FALSE, useWeights = FALSE)
    
    #need to use TrueScore to pull gene based scores matrix for generating phenotypes
    if(TrueScore == "IBS.gene"){
      PhenoScore = IBS.gene
    } else if(TrueScore == "Incomp.gene"){
      PhenoScore = Incomp.gene
    } else if(TrueScore == "AMS.gene"){
      PhenoScore = AMS.gene
    } else {
      PhenoScore = BinMM.gene
    }
    
    #generate covariates
    #for now, a single binary and a single continous covariate
    CovData = GenCovData(SampleSize = numPairs, BinaryValues = 1, ContinuousValues = 1)
    
    #need to define null Betas for phenotype generation
    nSNP = ncol(RGenos) #this should be the number of SNPs
    nullBetas = rep(0,nSNP) #generate null beta values
    Betas = nullBetas
    
    #set assoc Betas
    #all betas have same effect for now
    for(jj in assocSNPs){
      Betas[jj] = effect
    }
    Betas = as.matrix(Betas, ncol = 1)
    
    #generate phenotypes, both continuous and binary
    #Based on single true score
    CatPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = TRUE, YPrev = YPrev,  Covariates = CovData, RGenoData = RGenos, ScoreData = PhenoScore, Betas = Betas, Gamma = Gamma)
    ContPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = FALSE,  Covariates = CovData, RGenoData = RGenos, ScoreData = PhenoScore, Betas = Betas, Gamma = Gamma)
    
    #standardize to range from 0-1
    IBS.stand = IBS.gene/(nSNP*2)
    Incomp.stand = Incomp.gene/(nSNP)
    AMS.stand = AMS.gene/(nSNP*2)
    BinMM.stand = BinMM.gene/(nSNP)
    #Combine R geno and Scores into 4 separate datasets, size: N x (m+1)
    RGeno.IBS.snp = cbind(RGenos, IBS.stand)
    RGeno.Incomp.snp = cbind(RGenos, Incomp.stand)
    RGeno.AMS.snp = cbind(RGenos, AMS.stand)
    RGeno.BinMM.snp = cbind(RGenos, BinMM.stand)
    
    ## Generate SKAT Null Models
    # formulas will be Y ~ covariates for continuous and dichotomous Y
    obj_dich=SKAT_Null_Model(CatPhenos~CovData, out_type="D", Adjustment = FALSE)
    obj_cont=SKAT_Null_Model(ContPhenos~CovData, out_type="C", Adjustment = FALSE)
    
    #Perform SKAT for all 8 combos of score and cont/dich outcome
    #unweighted SKAT
    if(length(kernelWeights) == 0){
      Stat_IBS_CatPhenos = SKAT(RGeno.IBS.snp, obj_dich, kernel = kernel, is_check_genotype = FALSE)
      Stat_Incomp_CatPhenos = SKAT(RGeno.Incomp.snp, obj_dich, kernel = kernel, is_check_genotype = FALSE)
      Stat_AMS_CatPhenos = SKAT(RGeno.AMS.snp, obj_dich, kernel = kernel, is_check_genotype = FALSE)
      Stat_BinMM_CatPhenos = SKAT(RGeno.BinMM.snp, obj_dich, kernel = kernel, is_check_genotype = FALSE)
      
      Stat_IBS_ContPhenos = SKAT(RGeno.IBS.snp, obj_cont, kernel = kernel, is_check_genotype = FALSE)
      Stat_Incomp_ContPhenos = SKAT(RGeno.Incomp.snp, obj_cont, kernel = kernel, is_check_genotype = FALSE)
      Stat_AMS_ContPhenos = SKAT(RGeno.AMS.snp, obj_cont, kernel = kernel, is_check_genotype = FALSE)
      Stat_BinMM_ContPhenos = SKAT(RGeno.BinMM.snp, obj_cont, kernel = kernel, is_check_genotype = FALSE)
    } else {
      Stat_IBS_CatPhenos = SKAT(RGeno.IBS.snp, obj_dich, kernel = kernel, weights = kernelWeights, is_check_genotype = FALSE)
      Stat_Incomp_CatPhenos = SKAT(RGeno.Incomp.snp, obj_dich, kernel = kernel, weights = kernelWeights, is_check_genotype = FALSE)
      Stat_AMS_CatPhenos = SKAT(RGeno.AMS.snp, obj_dich, kernel = kernel, weights = kernelWeights, is_check_genotype = FALSE)
      Stat_BinMM_CatPhenos = SKAT(RGeno.BinMM.snp, obj_dich, kernel = kernel, weights = kernelWeights, is_check_genotype = FALSE)
      
      Stat_IBS_ContPhenos = SKAT(RGeno.IBS.snp, obj_cont, kernel = kernel, weights = kernelWeights, is_check_genotype = FALSE)
      Stat_Incomp_ContPhenos = SKAT(RGeno.Incomp.snp, obj_cont, kernel = kernel, weights = kernelWeights, is_check_genotype = FALSE)
      Stat_AMS_ContPhenos = SKAT(RGeno.AMS.snp, obj_cont, kernel = kernel, weights = kernelWeights, is_check_genotype = FALSE)
      Stat_BinMM_ContPhenos = SKAT(RGeno.BinMM.snp, obj_cont, kernel = kernel, weights = kernelWeights, is_check_genotype = FALSE)
    }
    
    #fill columns in order
    ##Binary Phenos
    ##IBS Cat, Incomp Cat, AMS Cat, Bin MM Cat, 
    IBSCat = c(Stat_IBS_CatPhenos$Q, Stat_IBS_CatPhenos$p.value)
    IncompCat = c(Stat_Incomp_CatPhenos$Q, Stat_Incomp_CatPhenos$p.value)
    AMSCat = c(Stat_AMS_CatPhenos$Q, Stat_AMS_CatPhenos$p.value)
    BinMMCat = c(Stat_BinMM_CatPhenos$Q, Stat_BinMM_CatPhenos$p.value)
    ##Cont Phenos
    ##IBS Cont, Incomp Cont, AMS Cont, Bin MM Cont
    IBSCont = c(Stat_IBS_ContPhenos$Q, Stat_IBS_ContPhenos$p.value)
    IncompCont = c(Stat_Incomp_ContPhenos$Q, Stat_Incomp_ContPhenos$p.value)
    AMSCont = c(Stat_AMS_ContPhenos$Q, Stat_AMS_ContPhenos$p.value)
    BinMMCont = c(Stat_BinMM_ContPhenos$Q, Stat_BinMM_ContPhenos$p.value)
    
    statsAndPValsMat = cbind(IBSCat,IncompCat,AMSCat,BinMMCat,IBSCont,IncompCont,AMSCont,BinMMCont)
    
    print(paste0("Simulation ",ii," is complete."))
    
    #reset Betas
    Betas = nullBetas
    
    statsAndPValsMat
  }, mc.cores = 4)
  
  statsAndPValsAll = matrix(unlist(statsAndPValsMat),nrow = numSims, ncol = 16, byrow = TRUE)
  
  #separate into cat and cont phenos
  statsAndPValsCatPhenos = statsAndPValsAll[,1:8]
  statsAndPValsContPhenos = statsAndPValsAll[,9:16] 
  
  #define values for naming conventions
  if(LowLD == TRUE){
    ld = "LowLD"
  } else {
    ld = "HighLD"
  }
  
  if(weightedScores == FALSE){
    weighted = ""
  } else{
    weighted = "WeightedScores"
  }
  
  if(standardizeScores == FALSE){
    standardized = ""
  } else {
    standardized = "StandardizedScores"
  }
  
  #write out the Stats and p values
  write.csv(statsAndPValsCatPhenos, file = paste0(path,"/Power_CatPhenos_",weighted,"_",standardized,"_Prev",YPrev*100,"_mPlus1_SKAT_kernel",kernel,"standardized_StatsAndPValues_",percentageAssoc,"SNPsAssociated_",ld,"_TrueScore",TrueScore,"_",ORSize,"OR_assocSNPOrScore",snpOrScore,".csv"))
  write.csv(statsAndPValsContPhenos, file = paste0(path,"/Power_ContPhenos_",weighted,"_",standardized,"_mPlus1_SKAT_kernel",kernel,"standardized_StatsAndPValues_",percentageAssoc,"SNPsAssociated_",ld,"_TrueScore",TrueScore,"_",ORSize,"OR_assocSNPOrScore",snpOrScore,".csv"))
  
  #calculate the power for each combination of score and r geno
  #first column for cat phenos, second for cont phenos
  #rows go IBS, Incomp, AMS, Bin MM
  power = matrix(nrow = 4, ncol = 2)
  
  #pull only the pvalues
  PValsCatPhenos = statsAndPValsCatPhenos[,c(2,4,6,8)]
  PValsContPhenos = statsAndPValsContPhenos[,c(2,4,6,8)]
  
  #calc power
  for(jj in 1:4){
    power[jj,1] = sum(PValsCatPhenos[,jj] <= 0.05)/numSims
    power[jj,2] = sum(PValsContPhenos[,jj] <= 0.05)/numSims
  }
  
  #write out power results
  write.csv(power, file = paste0(path,"/Power_Results_",percentageAssoc,"SNPsAssociated_",ld,"_",weighted,"_",standardized,"_mPlus1_SKAT_kernel",kernel,"standardized_CatAndContPhenos_Prev",YPrev*100,"_TrueScore",TrueScore,"_",ORSize,"OR_assocSNPOrScore",snpOrScore,".csv"))
}

RunPowerPipelineSKAT_RSNPs_WFunctions = function(chr, gene, numPairs, YPrev, kernel, kernelWeights=c(), Gamma, TrueScore, ORSize, standardizeScores = FALSE, weightedScores = FALSE, scoreWeights, percentageAssoc, LowLD){
  #function to run whole power pipeline , Calculates SKAT stat and pvalue for comparison
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #kernel = which kernel to use for SKAT (linear, linear.weighted, IBS, IBS.weighted)
  #kernelWeights = numeric vector of weights for weighted kernels
  #Gamma = effect size for score, length 1, can be 0
  #TrueScore = IBS.gene, Incomp.gene, AMS.gene, or BinMM.gene
  #ORSize = Small, Medium, or Large for what OR was used for the associated SNP/score
  #standardizeScores = T or F whether the scores should be standardized based on maximum score value
  #weightedScores = T or F whether the scores will be weighted
  #scoreWeights = m x 1 vector of weights, one weight for each SNP
  #percentageAssoc = percentage of SNPs associated with outcome (either 5, 25, 50, 75, or 100) 
  #LowLD = TRUE or FALSE whether the associated SNPs are in low LD or high LD
  #### This value shows which score is potentially being used to create phenotypes
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes power values to csv files
  
  #load libraries
  library(SKAT)
  library(parallel)
  
  #need to define this for naming at the end
  snpOrScore = "RSNP"
  
  #number of sims always the same
  numSims = 100
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_II_Sims/HapGen_Files/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/ProjectIISourceFunctions_v2.R")
  
  #turn output matrices into output lists
  myList = lapply(1:numSims, rep, times = 1)
  statsAndPValsMat = matrix(NA,ncol=16)
  
  statsAndPValsMat = mclapply(myList, function(ii){
    statsAndPValsMat = matrix(NA,ncol=16)
    
    #generate alternate phenotypes with R SNPs associated
    PhenosList = CalcAltPhenotypeData_RSNPs(chr = chr, numPairs = numPairs, simNum = ii, YPrev = YPrev, gene = gene, path = path, ORSize = ORSize, LowLD = LowLD, percentAssoc = percentAssoc, TrueScore = TrueScore)
    
    #run SKAT analysis
    # statsAndPValsMat = RunSKATAnalysis(PhenoList = PhenosList, kernel = kernel)
    statsAndPValsMat = RunSKATBinaryAnalysis(PhenoList = PhenosList, kernel = kernel)
    
    print(paste0("Simulation ",ii," is complete."))
    
    statsAndPValsMat
  }, mc.cores = 4)
  
  statsAndPValsAll = matrix(unlist(statsAndPValsMat),nrow = numSims, ncol = 16, byrow = TRUE)
  
  #separate into cat and cont phenos
  statsAndPValsCatPhenos = statsAndPValsAll[,1:8]
  statsAndPValsContPhenos = statsAndPValsAll[,9:16] 
  
  #define values for naming conventions
  if(LowLD == TRUE){
    ld = "LowLD"
  } else {
    ld = "HighLD"
  }
  
  if(weightedScores == FALSE){
    weighted = ""
  } else{
    weighted = "WeightedScores"
  }
  
  if(standardizeScores == FALSE){
    standardized = ""
  } else {
    standardized = "StandardizedScores"
  }
  
  #write out the Stats and p values
  write.csv(statsAndPValsCatPhenos, file = paste0(path,"/Power_CatPhenos_",weighted,"_",standardized,"_Prev",YPrev*100,"_mPlus1_SKAT_kernel",kernel,"_StatsAndPValues_",percentageAssoc,"SNPsAssociated_",ld,"_TrueScore",TrueScore,"_",ORSize,"OR_assocSNPOrScore",snpOrScore,".csv"))
  write.csv(statsAndPValsContPhenos, file = paste0(path,"/Power_ContPhenos_",weighted,"_",standardized,"_mPlus1_SKAT_kernel",kernel,"_StatsAndPValues_",percentageAssoc,"SNPsAssociated_",ld,"_TrueScore",TrueScore,"_",ORSize,"OR_assocSNPOrScore",snpOrScore,".csv"))
  
  #calculate the power for each combination of score and r geno
  #first column for cat phenos, second for cont phenos
  #rows go IBS, Incomp, AMS, Bin MM
  power = matrix(nrow = 4, ncol = 2)
  
  #pull only the pvalues
  PValsCatPhenos = statsAndPValsCatPhenos[,c(2,4,6,8)]
  PValsContPhenos = statsAndPValsContPhenos[,c(2,4,6,8)]
  
  #calc power
  for(jj in 1:4){
    power[jj,1] = sum(PValsCatPhenos[,jj] <= 0.05)/numSims
    power[jj,2] = sum(PValsContPhenos[,jj] <= 0.05)/numSims
  }
  
  #write out power results
  write.csv(power, file = paste0(path,"/Power_Results_",percentageAssoc,"SNPsAssociated_",ld,"_",weighted,"_",standardized,"_mPlus1_SKAT_kernel",kernel,"_CatAndContPhenos_Prev",YPrev*100,"_TrueScore",TrueScore,"_",ORSize,"OR_assocSNPOrScore",snpOrScore,".csv"))
}
