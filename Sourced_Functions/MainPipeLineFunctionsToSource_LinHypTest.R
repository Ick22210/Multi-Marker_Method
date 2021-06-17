##########################################################################################################################################
## Shi, Song, Chen, and Li Paper Pipelines
##########################################################################################################################################
# rely on code from 2019 Annals of Statistics Paper: Linear Hyp Testing for High Dim GLMs
##########################################################################################################################################
#################          Joint Testing         #########################################################################################
##########################################################################################################################################   
## TIE Pipelines:
#categorical outcome - standard GLM
RunTIEPipelineGLMCat = function(chr, gene, numPairs, YPrev, standardizeScores = FALSE, weightedScores = FALSE, scoreWeights, score, start, stop){
  #function to run whole TIE pipeline, Calculates LRT test stat for Lin Hyp test method and whether p-value <0.05
  # Only use when Y is binary
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #standardizeScores = T or F whether the scores should be standardized based on maximum score value
  #weightedScores = T or F whether the scores will be weighted
  #scoreWeights = m x 1 vector of weights, one weight for each SNP
  # score - which score we are using in the main JST model
  # start - simulation number to start at
  # stop - simulation number to stop at
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes TIE values to csv files
  
  library(parallel)
  suppressMessages(library(epicalc))
  
  #always the same
  numSims = stop  
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_II_Sims/HapGen_Files/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/ProjectIISourceFunctions_v2.R")
  
  myList = lapply(start:numSims, rep, times = 1)
  
  statsAndPVals = mclapply(myList, function(ii){
    #define matrix to hold all Stats and Pvalues
    statsAndPVals = matrix(NA, nrow = numSims, ncol = 3)
    
    #pull recipient and donor genotypes
    RGenos = obtainRGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    DGenos = obtainDGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    
    #calculate single snp scores
    if(score == "IBS"){
      Score.snp = calcIBSMismatch(RGenosMat = RGenos, DGenosMat = DGenos)
    } else if(score == "Incomp"){
      Score.snp = calcIncompatibilityScore(RGenosMat = RGenos, DGenosMat = DGenos)
    } else if(score == "AMS"){
      Score.snp = calcAMS(RGenosMat = RGenos, DGenosMat = DGenos)
    } else{
      Score.snp = calcBinaryMM(RGenosMat = RGenos, DGenosMat = DGenos)
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
    
    #generate covariates
    #for now, a single binary and a single continous covariate
    CovData = GenCovData(SampleSize = numPairs, BinaryValues = 1, ContinuousValues = 1)
    
    #generate phenotypes, both continuous and binary
    CatPhenos = GenNullPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = TRUE, YPrev = YPrev,  Covariates = CovData)
    
    allData = cbind(CovData, RGenos, Score.gene, CatPhenos)
    
    #fit the null and alternative models
    fitNull = glm(CatPhenos~CovData, family = binomial)
    fitAlt = glm(CatPhenos~CovData+RGenos+Score.gene, family = binomial, epsilon = 1e-6)	
    
    if(fitNull$deviance - fitAlt$deviance >= 0){
      # construct the likelihood ratio statistics
      TL = lrtest(fitNull, fitAlt)
      
      pv = TL$p.value
      stat = TL$Chisquared
      doF = TL$df
    } else {
      pv = 1
      stat = 0
      doF = (dim(RGenos)[2]+1)
    }
    
    print(paste0("Simulation ",ii," is complete."))
    
    statsAndPVals = c(pv=pv, TL=stat, dof=doF)
    statsAndPVals
  }, mc.cores=4)
  
  statsAndPValsAll = matrix(unlist(statsAndPVals),nrow = numSims, ncol = 3, byrow = TRUE)
  
  #rename 
  statsAndPValsCatPhenos = statsAndPValsAll
  
  #write out the Stats and p values
  if(weightedScores == FALSE){
    if(standardizeScores == FALSE){
      write.csv(statsAndPValsCatPhenos, file = paste0(path,"/TIE_CatPhenos_Prev",YPrev*100,"_LinHypTest_Score",score,"_Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    } else { #scores are standardized
      write.csv(statsAndPValsCatPhenos, file = paste0(path,"/TIE_CatPhenos_Prev",YPrev*100,"_LinHypTest_Score",score,"_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
    }
  } else { #scores are weighted
    write.csv(statsAndPValsCatPhenos, file = paste0(path,"/TIE_CatPhenos_Prev",YPrev*100,"_LinHypTest_Score",score,"_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
  }
}

#continuous outcome - standard GLM
RunTIEPipelineGLMCont = function(chr, gene, numPairs, standardizeScores = FALSE, weightedScores = FALSE, scoreWeights, score, start, stop){
  #function to run whole TIE pipeline, Calculates LRT test stat for Lin Hyp test method and whether p-value <0.05
  # only use when Y is continuous (Normal)
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #standardizeScores = T or F whether the scores should be standardized based on maximum score value
  #weightedScores = T or F whether the scores will be weighted
  #scoreWeights = m x 1 vector of weights, one weight for each SNP
  # score - which score we are using in the main JST model
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes TIE values to csv files
  
  library(parallel)
  suppressMessages(library(epicalc))
  
  #always the same
  numSims = stop  
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_II_Sims/HapGen_Files/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/ProjectIISourceFunctions_v2.R")
  #source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/Linear_ADMM0.r")
  
  myList = lapply(start:numSims, rep, times = 1)
  
  statsAndPVals = mclapply(myList, function(ii){
    #define matrix to hold all Stats and Pvalues
    statsAndPVals = matrix(NA, nrow = numSims, ncol = 3)
    
    #pull recipient and donor genotypes
    RGenos = obtainRGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    DGenos = obtainDGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    
    #calculate single snp scores
    if(score == "IBS"){
      Score.snp = calcIBSMismatch(RGenosMat = RGenos, DGenosMat = DGenos)
    } else if(score == "Incomp"){
      Score.snp = calcIncompatibilityScore(RGenosMat = RGenos, DGenosMat = DGenos)
    } else if(score == "AMS"){
      Score.snp = calcAMS(RGenosMat = RGenos, DGenosMat = DGenos)
    } else{
      Score.snp = calcBinaryMM(RGenosMat = RGenos, DGenosMat = DGenos)
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
    
    #generate covariates
    #for now, a single binary and a single continous covariate
    CovData = GenCovData(SampleSize = numPairs, BinaryValues = 1, ContinuousValues = 1)
    
    #generate phenotypes, both continuous and binary
    ContPhenos = GenNullPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = FALSE, Covariates = CovData)
    
    #fit the null and alternative models
    fitNull = glm(ContPhenos~CovData, family = gaussian)
    fitAlt = glm(ContPhenos~CovData+RGenos+Score.gene, family = gaussian, epsilon = 1e-6)	
    
    if(fitNull$deviance - fitAlt$deviance >= 0){
      # construct the likelihood ratio statistics
      TL = lrtest(fitNull, fitAlt)
      
      pv = TL$p.value
      stat = TL$Chisquared
      doF = TL$df
    } else {
      pv = 1
      stat = 0
      doF = (dim(RGenos)[2]+1)
    }
    
    print(paste0("Simulation ",ii," is complete."))
    
    statsAndPVals = c(pv=pv, TL=stat, dof=doF)
    statsAndPVals
  }, mc.cores=1)
  
  statsAndPValsAll = matrix(unlist(statsAndPVals),nrow = numSims, ncol = 3, byrow = TRUE)
  
  #separate into cat and cont phenos
  statsAndPValsContPhenos = statsAndPValsAll
  
  #write out the Stats and p values
  if(weightedScores == FALSE){
    if(standardizeScores == FALSE){
      write.csv(statsAndPValsContPhenos, file = paste0(path,"/TIE_ContPhenos_LinHypTest_Score",score,"_Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    } else { #scores are standardized
      write.csv(statsAndPValsContPhenos, file = paste0(path,"/TIE_ContPhenos_LinHypTest_Score",score,"_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
    }
  } else { #scores are weighted
    write.csv(statsAndPValsContPhenos, file = paste0(path,"/TIE_ContPhenos_LinHypTest_Score",score,"_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
  }
}

#categorical outcome - Linear Hypothesis Testing
RunTIEPipelineLinHypTestCat = function(chr, gene, numPairs, YPrev, standardizeScores = FALSE, weightedScores = FALSE, scoreWeights, score, start, stop){
  #function to run whole TIE pipeline, Calculates LRT test stat for Lin Hyp test method and whether p-value <0.05
  # Only use when Y is binary
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #standardizeScores = T or F whether the scores should be standardized based on maximum score value
  #weightedScores = T or F whether the scores will be weighted
  #scoreWeights = m x 1 vector of weights, one weight for each SNP
  # score - which score we are using in the main JST model
  # start - simulation number to start at
  # stop - simulation number to stop at
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes TIE values to csv files
  
  library(parallel)

  #always the same
  numSims = stop  
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_II_Sims/HapGen_Files/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/ProjectIISourceFunctions_v2.R")
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/Logistic_ADMM0.r")
  
  myList = lapply(start:numSims, rep, times = 1)
  # p-value initialization
  pv <- rep(0, 3)
  Tall <- matrix(0, 1, 3)
  beta.al <- matrix(0, 1, 2)
  
  statsAndPVals = mclapply(myList, function(ii){
    #define matrix to hold all Stats and Pvalues
    statsAndPVals = matrix(NA, nrow = numSims, ncol = 3)
    
    #pull recipient and donor genotypes
    RGenos = obtainRGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    DGenos = obtainDGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    
    #calculate single snp scores
    if(score == "IBS"){
      Score.snp = calcIBSMismatch(RGenosMat = RGenos, DGenosMat = DGenos)
    } else if(score == "Incomp"){
      Score.snp = calcIncompatibilityScore(RGenosMat = RGenos, DGenosMat = DGenos)
    } else if(score == "AMS"){
      Score.snp = calcAMS(RGenosMat = RGenos, DGenosMat = DGenos)
    } else{
      Score.snp = calcBinaryMM(RGenosMat = RGenos, DGenosMat = DGenos)
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
    
    #generate covariates
    #for now, a single binary and a single continous covariate
    CovData = GenCovData(SampleSize = numPairs, BinaryValues = 1, ContinuousValues = 1)
    
    #generate phenotypes, both continuous and binary
    CatPhenos = GenNullPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = TRUE, YPrev = YPrev,  Covariates = CovData)
    
    # define location of zero components
    # basically, this defines what the null hypothesis is, I think
    # so for joint null, we need the non-zero components to be for covariates
    # and the beta and gamma components to be zero
    numCov = dim(CovData)[2]
    numSNPs = dim(RGenos)[2]
    N = c(rep(FALSE, numCov+1), rep(TRUE,1+numSNPs))
    
    #need to combine the covariates, r genos, and score matrices together
    #need to include column of 1s for intercept?
    intercept = matrix(1,nrow = numPairs, ncol = 1)
    designMat = cbind(intercept, CovData, RGenos, Score.gene)
    
    #then combine the phenos with the design matrix as a list
    Model = list(X=designMat, Y=CatPhenos)
    
    # estimate the uncontrained estimator
    beta.unre <- cv.SCAD_ADMM_unre(X=Model$X, Y=Model$Y, N=N, beta0=rep(0, dim(Model$X)[2]), err=1e-4, tune="cv", unpen=c(1,2,3))
    indice.unre <- beta.unre!=0
    # estimate the constrained estimator (only need this for score test)
    beta.re <- cv.SCAD_ADMM_re(X=Model$X, Y=Model$Y, N=N, beta0=rep(0, dim(Model$X)[2]), err=1e-4, tune="cv", unpen=c(1,2,3))
    indice.re <- beta.re!=0
    
    pi.unre <- logit(Model$X%*%beta.unre)
    pi.re <- logit(Model$X%*%beta.re)
    
    # construct the likelihood ratio statistics
    TL <- 2*(sum(log(1+exp(Model$X%*%beta.re))-Model$Y*(Model$X%*%beta.re))-
               sum(log(1+exp(Model$X%*%beta.unre))-Model$Y*(Model$X%*%beta.unre)))
    #if LRT stat is less than 0, there is error, so set to -1
    if(TL < 0){
      TL = -1
    }
    
    # construct the Wald statistics
    #B_0 should give you Omega_a hat
    A = crossprod(Model$X[,indice.unre|N], as.vector(pi.unre*(1-pi.unre))*Model$X[,indice.unre|N])
    if(rcond(A) >= 1e-10){
      B_0 <- solve(A)
      #d_0 should determine which rows to subset
      #gives the first value that is restricted
      d_0 <- sum(indice.re == TRUE) + 1
      #in case the inverse doesn't exist, need another if else
      B = B_0[d_0:ncol(B_0),d_0:ncol(B_0)]
      if(rcond(B) >= 1e-10){
        #so the B_0 needs to be subset to m rows and columns that are restricted based on H0
        TW <- crossprod(beta.unre[N], solve(B, beta.unre[N]))
      } else {
        TW = -2
      }
    } else {
      TW = -1
    }
    
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
    
    #set pv to 1 if reject, 0 otherwise
    #LRT
    if (TL>=qchisq(0.95, df = doF)){
      pv[1] <- pv[1]+1/5000
    }
    #Wald
    if (TW>=qchisq(0.95, df = doF)){
      pv[2] <- pv[2]+1/5000
    }
    #Score
    if (TS>=qchisq(0.95, df = doF)){
      pv[3] <- pv[3]+1/5000
    }
    
    print(paste0("Simulation ",ii," is complete."))
    
    Tall[1, ] <- c(TL, TW, TS)
    beta.al[1, ] <- c(sum(beta.re!=0), sum(beta.unre!=0))
    
    statsAndPVals = list(pv=pv, TScores=Tall, beta=beta.al)
    statsAndPVals
  }, mc.cores=4)
  
  n = length(statsAndPVals)
  
  pValsCV = matrix(nrow = n, ncol = 3)
  ScoresCV = matrix(nrow = n, ncol = 3)
  BetasCV = matrix(nrow = n, ncol = 2)
  
  for(ll in 1:n){
    pValsCV[ll,] = unlist(statsAndPVals[[ll]]$pv)
    ScoresCV[ll,] = unlist(statsAndPVals[[ll]]$TScores)
    BetasCV[ll,] = unlist(statsAndPVals[[ll]]$beta)
  }
  
  #write out the Stats and p values
  if(weightedScores == FALSE){
    if(standardizeScores == FALSE){
      write.csv(pValsCV, file = paste0(path,"/TIE_CatY_Prev",YPrev*100,"_LinHypTest_Score",score,"_JointTesting_Sim",start,"to",numSims,"_PValuesCV_ForceCovFit.csv"))
      write.csv(ScoresCV, file = paste0(path,"/TIE_CatY_Prev",YPrev*100,"_LinHypTest_Score",score,"_JointTesting_Sim",start,"to",numSims,"_StatsCV_ForceCovFit.csv"))
      write.csv(BetasCV, file = paste0(path,"/TIE_CatY_Prev",YPrev*100,"_LinHypTest_Score",score,"_JointTesting_Sim",start,"to",numSims,"_BetasCV_ForceCovFit.csv"))
    } else { #scores are standardized
      write.csv(statsAndPValsCatPhenos, file = paste0(path,"/TIE_CatPhenos_Prev",YPrev*100,"_LinHypTest_Score",score,"_JointTesting_NoCovs_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
    }
  } else { #scores are weighted
    write.csv(statsAndPValsCatPhenos, file = paste0(path,"/TIE_CatPhenos_Prev",YPrev*100,"_LinHypTest_Score",score,"_JointTesting_NoCovs_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
  }
}

#continuous outcome - Linear Hypothesis Testing
RunTIEPipelineLinHypTestCont = function(chr, gene, numPairs, standardizeScores = FALSE, weightedScores = FALSE, scoreWeights, score, start, stop){
  #function to run whole TIE pipeline, Calculates LRT test stat for Lin Hyp test method and whether p-value <0.05
  # only use when Y is continuous (Normal)
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #standardizeScores = T or F whether the scores should be standardized based on maximum score value
  #weightedScores = T or F whether the scores will be weighted
  #scoreWeights = m x 1 vector of weights, one weight for each SNP
  # score - which score we are using in the main JST model
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes TIE values to csv files
  
  library(parallel)

  #always the same
  numSims = stop  
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_II_Sims/HapGen_Files/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/ProjectIISourceFunctions_v2.R")
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/Linear_ADMM0.r")
  
  myList = lapply(start:numSims, rep, times = 1)
  # p-value initialization
  pv <- rep(0, 6)
  Tall <- matrix(0, 1, 3)
  beta.al <- matrix(0, 1, 2)
  
  statsAndPVals = mclapply(myList, function(ii){
    #define matrix to hold all Stats and Pvalues
    statsAndPVals = matrix(NA, nrow = numSims, ncol = 3)
    
    #pull recipient and donor genotypes
    RGenos = obtainRGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    DGenos = obtainDGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    
    #calculate single snp scores
    if(score == "IBS"){
      Score.snp = calcIBSMismatch(RGenosMat = RGenos, DGenosMat = DGenos)
    } else if(score == "Incomp"){
      Score.snp = calcIncompatibilityScore(RGenosMat = RGenos, DGenosMat = DGenos)
    } else if(score == "AMS"){
      Score.snp = calcAMS(RGenosMat = RGenos, DGenosMat = DGenos)
    } else{
      Score.snp = calcBinaryMM(RGenosMat = RGenos, DGenosMat = DGenos)
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
    
    #generate covariates
    #for now, a single binary and a single continous covariate
    CovData = GenCovData(SampleSize = numPairs, BinaryValues = 1, ContinuousValues = 1)
    
    #generate phenotypes, both continuous and binary
    ContPhenos = GenNullPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = FALSE, Covariates = CovData)
    
    # define location of zero components
    # basically, this defines what the null hypothesis is, I think
    # so for joint null, we need the non-zero components to be for covariates
    # and the beta and gamma components to be zero
    numCov = dim(CovData)[2]
    numSNPs = dim(RGenos)[2]
    N = c(rep(FALSE, numCov+1), rep(TRUE,1+numSNPs))
    
    #need to combine the covariates, r genos, and score matrices together
    #need to include column of 1s for intercept?
    intercept = matrix(1,nrow = numPairs, ncol = 1)
    designMat = cbind(intercept, CovData, RGenos, Score.gene)
    
    #then combine the phenos with the design matrix as a list
    Model = list(X=designMat, Y=ContPhenos)
    
    # estimate the uncontrained estimator
    beta.unre <- cv.SCAD_ADMM_unre(X=Model$X, Y=Model$Y, N=N, beta0=rep(0, dim(Model$X)[2]), err=1e-4, tune="cv", unpen = c(1,2,3))
    indice.unre <- beta.unre!=0
    
    # estimate the constrained estimator
    beta.re <- cv.SCAD_ADMM_re(X=Model$X, Y=Model$Y, N=N, beta0=rep(0, dim(Model$X)[2]), err=1e-4, tune="cv", unpen = c(1,2,3))
    indice.re <- beta.re!=0
    
    # estimate the conditional variance
    #    sig2 <- mean((Model$Y-Model$X%*%beta.unre)^2)
    n = numSims
    sig2 <- mean((Model$Y-Model$X%*%beta.unre)^2)*n/(n-sum(beta.unre!=0))
    
    # construct the likelihood ratio statistic
    TL <- sum((Model$X%*%beta.re-Model$Y)^2)-sum((Model$X%*%beta.unre-Model$Y)^2)
    #if LRT stat is less than 0, there is error, so set to -1
    if(TL < 0){
      TL = -1
    }
    
    # construct the Wald statistic
    #B_0 should give you Omega_a hat
    B = crossprod(Model$X[,indice.unre|N], Model$X[,indice.unre|N])
    if(rcond(B) >= 1e-10){
      B_0 <- solve(B)
      #d_0 should determine which rows to subset
      #gives the first value that is restricted
      d_0 <- sum(indice.re == TRUE) + 1
      #in case the inverse doesn't exist, need another if else
      A = B_0[d_0:ncol(B_0),d_0:ncol(B_0)]
      if(rcond(A) >= 1e-10){
        TW <- crossprod(beta.unre[N], solve(A, beta.unre[N]))
      } else {
        TW = -2
      }
    } else {
      TW = -1
    }
    
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
    
    #set pv to 1 if reject, 0 otherwise
    #LRT
    if (TL>=sig2*qchisq(0.95, df = doF)){
      pv[1] <- pv[1]+1/5000
    }
    if (TL>=qchisq(0.95, df = doF)){
      pv[2] <- pv[2]+1/5000
    }
    #Wald
    if (TW>=sig2*qchisq(0.95, df = doF)){
      pv[3] <- pv[3]+1/5000
    }
    if (TW>=qchisq(0.95, df = doF)){
      pv[4] <- pv[4]+1/5000
    }
    #Score
    if (TS>=sig2*qchisq(0.95, df = doF)){
      pv[5] <- pv[5]+1/5000
    }
    if (TS>=qchisq(0.95, df = doF)){
      pv[6] <- pv[6]+1/5000
    }
    
    print(paste0("Simulation ",ii," is complete."))
    
    Tall[1, ] <- c(TL, TW, TS)
    beta.al[1, ] <- c(sum(beta.re!=0), sum(beta.unre!=0))
    
    statsAndPVals = list(pv=pv, TScores=Tall, beta=beta.al)
    statsAndPVals
  }, mc.cores=1)
  
  n = length(statsAndPVals)
  
  pVals = matrix(nrow = n, ncol = 6)
  Scores = matrix(nrow = n, ncol = 3)
  Betas = matrix(nrow = n, ncol = 2)
  
  for(ll in 1:n){
    pVals[ll,] = unlist(statsAndPVals[[ll]]$pv)
    Scores[ll,] = unlist(statsAndPVals[[ll]]$TScores)
    Betas[ll,] = unlist(statsAndPVals[[ll]]$beta)
  }
  
  #write out the Stats and p values
  if(weightedScores == FALSE){
    if(standardizeScores == FALSE){
      write.csv(pVals, file = paste0(path,"/TIE_ContY_LinHypTest_Score",score,"_JointTesting_Sim",start,"to",numSims,"_PValuesCV_ForceCovs.csv"))
      write.csv(Scores, file = paste0(path,"/TIE_ContY_LinHypTest_Score",score,"_JointTesting_Sim",start,"to",numSims,"_StatsCV_ForceCovs.csv"))
      write.csv(Betas, file = paste0(path,"/TIE_ContY_LinHypTest_Score",score,"_JointTesting_Sim",start,"to",numSims,"_BetasCV_ForceCovs.csv"))
    } else { #scores are standardized
      write.csv(statsAndPValsContPhenos, file = paste0(path,"/TIE_ContPhenos_LinHypTest_Score",score,"_GammaOnly_NoCovs_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
    }
  } else { #scores are weighted
    write.csv(statsAndPValsContPhenos, file = paste0(path,"/TIE_ContPhenos_LinHypTest_Score",score,"_GammaOnly_NoCovs_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
  }
}

#################           
## Power Pipelines:
# Score is associated, R SNP is not
# binary outcome - standard GLM
RunPowerPipelineGLMCat_Score = function(chr, gene, numPairs, YPrev, Gamma, TrueScore, ORSize, standardizeScores = FALSE, weightedScores = FALSE, scoreWeights, start, stop, percentageAssoc, LowLD){
  # Function to determine power of standard GLM method when score is associated, joint testing
  # Only use when Y is binary
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #Gamma = effect size for score, length 1, can be 0
  #TrueScore = IBS.gene, Incomp.gene, AMS.gene, or BinMM.gene
  #ORSize = Small, Medium, or Large for what OR was used for the associated SNP/score
  #standardizeScores = T or F whether the scores should be standardized based on maximum score value
  #weightedScores = T or F whether the scores will be weighted
  #scoreWeights = m x 1 vector of weights, one weight for each SNP
  # start - simulation number to start at
  # stop - simulation number to stop at
  #percentageAssoc = percentage of SNPs associated with outcome (either 5, 25, 50, 75, or 100) 
  #LowLD = True or FALSE whether the associated SNPs are in low LD or high LD
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes power values to csv files
  
  library(parallel)
  suppressMessages(library(lmtest))
  
  #need to define this for naming at the end
  snpOrScore = TrueScore
  
  #always the same
  numSims = stop  
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_II_Sims/HapGen_Files/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/ProjectIISourceFunctions_v2.R")
  
  myList = lapply(start:numSims, rep, times = 1)
  
  statsAndPVals = mclapply(myList, function(ii){
    #define matrix to hold all Stats and Pvalues
    statsAndPVals = matrix(NA, nrow = numSims, ncol = 3)
    
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
    CatPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = TRUE, YPrev = YPrev,  Covariates = CovData, RGenoData = RGenos, ScoreData = PhenoScore, Betas = Betas, Gamma = Gamma)
    
    #fit the null and alternative models
    fitNull = glm(CatPhenos~CovData, family = binomial)
    fitAltIBS = glm(CatPhenos~CovData+RGenos+IBS.gene, family = binomial, epsilon = 1e-6)	
    fitAltIncomp = glm(CatPhenos~CovData+RGenos+Incomp.gene, family = binomial, epsilon = 1e-6)	
    fitAltAMS = glm(CatPhenos~CovData+RGenos+AMS.gene, family = binomial, epsilon = 1e-6)	
    fitAltBinMM = glm(CatPhenos~CovData+RGenos+BinMM.gene, family = binomial, epsilon = 1e-6)	
    
    if(fitNull$deviance - fitAltIBS$deviance >= 0){
      # construct the likelihood ratio statistics
      TL.IBS = lrtest(fitNull, fitAltIBS)
      
      pv.IBS = TL.IBS$p.value
      stat.IBS = TL.IBS$Chisquared
      doF.IBS = TL.IBS$df
    } else {
      pv.IBS = 1
      stat.IBS = 0
      doF.IBS = (dim(RGenos)[2]+1)
    }
    
    if(fitNull$deviance - fitAltIncomp$deviance >= 0){
      # construct the likelihood ratio statistics
      TL.Incomp = lrtest(fitNull, fitAltIncomp)
      
      pv.Incomp = TL.Incomp$p.value
      stat.Incomp = TL.Incomp$Chisquared
      doF.Incomp = TL.Incomp$df
    } else {
      pv.Incomp = 1
      stat.Incomp = 0
      doF.Incomp = (dim(RGenos)[2]+1)
    }
    
    if(fitNull$deviance - fitAltAMS$deviance >= 0){
      # construct the likelihood ratio statistics
      TL.AMS = lrtest(fitNull, fitAltAMS)
      
      pv.AMS = TL.AMS$p.value
      stat.AMS = TL.AMS$Chisquared
      doF.AMS = TL.AMS$df
    } else {
      pv.AMS = 1
      stat.AMS = 0
      doF.AMS = (dim(RGenos)[2]+1)
    }
    
    if(fitNull$deviance - fitAltBinMM$deviance >= 0){
      # construct the likelihood ratio statistics
      TL.BinMM = lrtest(fitNull, fitAltBinMM)
      
      pv.BinMM = TL.BinMM$p.value
      stat.BinMM = TL.BinMM$Chisquared
      doF.BinMM = TL.BinMM$df
    } else {
      pv.BinMM = 1
      stat.BinMM = 0
      doF.BinMM = (dim(RGenos)[2]+1)
    }
    
    print(paste0("Simulation ",ii," is complete."))
    
    pv = c(pv.IBS, pv.Incomp, pv.AMS, pv.BinMM)
    stat = c(stat.IBS, stat.Incomp, stat.AMS, stat.BinMM)
    doF = c(doF.IBS, doF.Incomp, doF.AMS, doF.BinMM)
    
    statsAndPVals = list(pv=pv, TL=stat, dof=doF)
    statsAndPVals
  }, mc.cores=4)
  
  n = length(statsAndPVals)
  
  pVals = matrix(nrow = n, ncol = 4)
  Scores = matrix(nrow = n, ncol = 4)
  doF = matrix(nrow = n, ncol = 4)
  
  for(ll in 1:n){
    pVals[ll,] = unlist(statsAndPVals[[ll]]$pv)
    Scores[ll,] = unlist(statsAndPVals[[ll]]$TL)
    doF[ll,] = unlist(statsAndPVals[[ll]]$dof)
  }
  
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
  write.csv(pVals, file = paste0(path,"/Power_CatPhenos_Prev",YPrev*100,"_GLM_Pvals_",percentageAssoc,"SNPsAssoc_",ld,"_TrueScore",TrueScore,"_",ORSize,"OR_assocScore",snpOrScore,"_Sim",start,"to",numSims,".csv"))
  write.csv(Scores, file = paste0(path,"/Power_CatPhenos_Prev",YPrev*100,"_GLM_Scores_",percentageAssoc,"SNPsAssoc_",ld,"_TrueScore",TrueScore,"_",ORSize,"OR_assocScore",snpOrScore,"_Sim",start,"to",numSims,".csv"))
  write.csv(doF, file = paste0(path,"/Power_CatPhenos_Prev",YPrev*100,"_GLM_dof_",percentageAssoc,"SNPsAssoc_",ld,"_TrueScore",TrueScore,"_",ORSize,"OR_assocScore",snpOrScore,"_Sim",start,"to",numSims,".csv"))
}

# continuous outcome - standard GLM
RunPowerPipelineGLMCont_Score = function(chr, gene, numPairs, Gamma, TrueScore, ORSize, standardizeScores = FALSE, weightedScores = FALSE, scoreWeights, start, stop, percentageAssoc, LowLD){
  # Function to determine power of standard GLM method when score is associated, joint testing
  # only use when Y is continuous (Normal)
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #Gamma = effect size for score, length 1, can be 0
  #TrueScore = IBS.gene, Incomp.gene, AMS.gene, or BinMM.gene
  #ORSize = Small, Medium, or Large for what OR was used for the associated SNP/score
  #standardizeScores = T or F whether the scores should be standardized based on maximum score value
  #weightedScores = T or F whether the scores will be weighted
  #scoreWeights = m x 1 vector of weights, one weight for each SNP
  #percentageAssoc = percentage of SNPs associated with outcome (either 5, 25, 50, 75, or 100) 
  #LowLD = True or FALSE whether the associated SNPs are in low LD or high LD
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes power values to csv files
  
  library(parallel)
  suppressMessages(library(lmtest))
  
  #need to define this for naming at the end
  snpOrScore = TrueScore
  
  #always the same
  numSims = stop  
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_II_Sims/HapGen_Files/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/ProjectIISourceFunctions_v2.R")

  myList = lapply(start:numSims, rep, times = 1)
  
  statsAndPVals = mclapply(myList, function(ii){
    #define matrix to hold all Stats and Pvalues
    statsAndPVals = matrix(NA, nrow = numSims, ncol = 3)
    
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
    ContPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = FALSE,  Covariates = CovData, RGenoData = RGenos, ScoreData = PhenoScore, Betas = Betas, Gamma = Gamma)
    
    #fit the null and alternative models
    fitNull = glm(ContPhenos~CovData, family = gaussian)
    fitAlt.IBS = glm(ContPhenos~CovData+RGenos+IBS.gene, family = gaussian, epsilon = 1e-6)	
    fitAlt.Incomp = glm(ContPhenos~CovData+RGenos+Incomp.gene, family = gaussian, epsilon = 1e-6)	
    fitAlt.AMS = glm(ContPhenos~CovData+RGenos+AMS.gene, family = gaussian, epsilon = 1e-6)	
    fitAlt.BinMM = glm(ContPhenos~CovData+RGenos+BinMM.gene, family = gaussian, epsilon = 1e-6)	
    
    if(fitNull$deviance - fitAlt.IBS$deviance >= 0){
      # construct the likelihood ratio statistics
      TL.IBS = lrtest(fitNull, fitAlt.IBS)
      
      pv.IBS = TL.IBS[2,5]
      stat.IBS = TL.IBS[2,4]
      doF.IBS = TL.IBS[2,1]
    } else {
      pv.IBS = 1
      stat.IBS = 0
      doF.IBS = (dim(RGenos)[2]+1)
    }
    if(fitNull$deviance - fitAlt.Incomp$deviance >= 0){
      # construct the likelihood ratio statistics
      TL.Incomp = lrtest(fitNull, fitAlt.Incomp)
      
      pv.Incomp = TL.Incomp[2,5]
      stat.Incomp = TL.Incomp[2,4]
      doF.Incomp = TL.Incomp[2,1]
    } else {
      pv.Incomp = 1
      stat.Incomp = 0
      doF.Incomp = (dim(RGenos)[2]+1)
    }
    if(fitNull$deviance - fitAlt.AMS$deviance >= 0){
      # construct the likelihood ratio statistics
      TL.AMS = lrtest(fitNull, fitAlt.AMS)
      
      pv.AMS = TL.AMS[2,5]
      stat.AMS = TL.AMS[2,4]
      doF.AMS = TL.AMS[2,1]
    } else {
      pv.AMS = 1
      stat.AMS = 0
      doF.AMS = (dim(RGenos)[2]+1)
    }
    if(fitNull$deviance - fitAlt.BinMM$deviance >= 0){
      # construct the likelihood ratio statistics
      TL.BinMM = lrtest(fitNull, fitAlt.BinMM)
      
      pv.BinMM = TL.BinMM[2,5]
      stat.BinMM = TL.BinMM[2,4]
      doF.BinMM = TL.BinMM[2,1]
    } else {
      pv.BinMM = 1
      stat.BinMM = 0
      doF.BinMM = (dim(RGenos)[2]+1)
    }
    
    print(paste0("Simulation ",ii," is complete."))
    
    pv = c(pv.IBS, pv.Incomp, pv.AMS, pv.BinMM)
    stat = c(stat.IBS, stat.Incomp, stat.AMS, stat.BinMM)
    doF = c(doF.IBS, doF.Incomp, doF.AMS, doF.BinMM)
    
    statsAndPVals = list(pv=pv, TL=stat, dof=doF)
    statsAndPVals
  }, mc.cores=4)
  
  n = length(statsAndPVals)
  
  pVals = matrix(nrow = n, ncol = 4)
  Scores = matrix(nrow = n, ncol = 4)
  doF = matrix(nrow = n, ncol = 4)
  
  for(ll in 1:n){
    pVals[ll,] = unlist(statsAndPVals[[ll]]$pv)
    Scores[ll,] = unlist(statsAndPVals[[ll]]$TL)
    doF[ll,] = unlist(statsAndPVals[[ll]]$dof)
  }
  
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
  write.csv(pVals, file = paste0(path,"/Power_ContPhenos_GLM_Pvals_",percentageAssoc,"SNPsAssoc_",ld,"_TrueScore",TrueScore,"_",ORSize,"OR_assocScore",snpOrScore,"_Sim",start,"to",numSims,".csv"))
  write.csv(Scores, file = paste0(path,"/Power_ContPhenos_GLM_Scores_",percentageAssoc,"SNPsAssoc_",ld,"_TrueScore",TrueScore,"_",ORSize,"OR_assocScore",snpOrScore,"_Sim",start,"to",numSims,".csv"))
  write.csv(doF, file = paste0(path,"/Power_ContPhenos_GLM_dof_",percentageAssoc,"SNPsAssoc_",ld,"_TrueScore",TrueScore,"_",ORSize,"OR_assocScore",snpOrScore,"_Sim",start,"to",numSims,".csv"))
}

# binary outcome - linear hypothesis testing
RunPowerPipelineLinHypTestCat_Score = function(chr, gene, numPairs, YPrev, Gamma, TrueScore, ORSize, standardizeScores = FALSE, weightedScores = FALSE, scoreWeights, start, stop, percentageAssoc, LowLD){
  # Function to determine power of linear hyp testing method when score is associated, joint testing
  # Only use when Y is binary
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #Gamma = effect size for score, length 1, can be 0
  #TrueScore = IBS.gene, Incomp.gene, AMS.gene, or BinMM.gene
  #ORSize = Small, Medium, or Large for what OR was used for the associated SNP/score
  #standardizeScores = T or F whether the scores should be standardized based on maximum score value
  #weightedScores = T or F whether the scores will be weighted
  #scoreWeights = m x 1 vector of weights, one weight for each SNP
  # start - simulation number to start at
  # stop - simulation number to stop at
  #percentageAssoc = percentage of SNPs associated with outcome (either 5, 25, 50, 75, or 100) 
  #LowLD = True or FALSE whether the associated SNPs are in low LD or high LD
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes power values to csv files
  
  library(parallel)
  
  #need to define this for naming at the end
  snpOrScore = TrueScore
  
  #always the same
  numSims = stop  
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_II_Sims/HapGen_Files/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/ProjectIISourceFunctions_v2.R")
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/Logistic_ADMM0.r")
  
  myList = lapply(start:numSims, rep, times = 1)
  # p-value initialization
  pv.IBS <- rep(0, 3)
  pv.Incomp <- rep(0, 3)
  pv.AMS <- rep(0, 3)
  pv.BinMM <- rep(0, 3)
  
  Tall.IBS <- matrix(0, 1, 3)
  Tall.Incomp <- matrix(0, 1, 3)
  Tall.AMS <- matrix(0, 1, 3)
  Tall.BinMM <- matrix(0, 1, 3)
  
  beta.al.IBS <- matrix(0, 1, 2)
  beta.al.Incomp <- matrix(0, 1, 2)
  beta.al.AMS <- matrix(0, 1, 2)
  beta.al.BinMM <- matrix(0, 1, 2)
  
  statsAndPVals = mclapply(myList, function(ii){
    #define matrix to hold all Stats and Pvalues
    statsAndPVals = matrix(NA, nrow = numSims, ncol = 3)
    
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
    CatPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = TRUE, YPrev = YPrev,  Covariates = CovData, RGenoData = RGenos, ScoreData = PhenoScore, Betas = Betas, Gamma = Gamma)
    
    # define location of zero components
    # basically, this defines what the null hypothesis is, I think
    # so for joint null, we need the non-zero components to be for covariates
    # and the beta and gamma components to be zero
    numCov = dim(CovData)[2]
    numSNPs = dim(RGenos)[2]
    N = c(rep(FALSE, numCov+1), rep(TRUE,1+numSNPs))
    
    #need to combine the covariates, r genos, and score matrices together
    #need to include column of 1s for intercept?
    intercept = matrix(1,nrow = numPairs, ncol = 1)
    designMat.IBS = cbind(intercept, CovData, RGenos, IBS.gene)
    designMat.Incomp = cbind(intercept, CovData, RGenos, Incomp.gene)
    designMat.AMS = cbind(intercept, CovData, RGenos, AMS.gene)
    designMat.BinMM = cbind(intercept, CovData, RGenos, BinMM.gene)
    
    #then combine the phenos with the design matrix as a list
    Model.IBS = list(X=designMat.IBS, Y=CatPhenos)
    Model.Incomp = list(X=designMat.Incomp, Y=CatPhenos)
    Model.AMS = list(X=designMat.AMS, Y=CatPhenos)
    Model.BinMM = list(X=designMat.BinMM, Y=CatPhenos)
    
    # estimate the uncontrained estimator
    beta.unre.IBS <- cv.SCAD_ADMM_unre(X=Model.IBS$X, Y=Model.IBS$Y, N=N, beta0=rep(0, dim(Model.IBS$X)[2]), err=1e-4, tune="cv", unpen=c(1,2,3))
    indice.unre.IBS <- beta.unre.IBS!=0
    
    beta.unre.Incomp <- cv.SCAD_ADMM_unre(X=Model.Incomp$X, Y=Model.Incomp$Y, N=N, beta0=rep(0, dim(Model.Incomp$X)[2]), err=1e-4, tune="cv", unpen=c(1,2,3))
    indice.unre.Incomp <- beta.unre.Incomp!=0
    
    beta.unre.AMS <- cv.SCAD_ADMM_unre(X=Model.AMS$X, Y=Model.AMS$Y, N=N, beta0=rep(0, dim(Model.AMS$X)[2]), err=1e-4, tune="cv", unpen=c(1,2,3))
    indice.unre.AMS <- beta.unre.AMS!=0
    
    beta.unre.BinMM <- cv.SCAD_ADMM_unre(X=Model.BinMM$X, Y=Model.BinMM$Y, N=N, beta0=rep(0, dim(Model.BinMM$X)[2]), err=1e-4, tune="cv", unpen=c(1,2,3))
    indice.unre.BinMM <- beta.unre.BinMM!=0
    
    # estimate the constrained estimator (only need this for score test)
    beta.re.IBS <- cv.SCAD_ADMM_re(X=Model.IBS$X, Y=Model.IBS$Y, N=N, beta0=rep(0, dim(Model.IBS$X)[2]), err=1e-4, tune="cv", unpen=c(1,2,3))
    indice.re.IBS <- beta.re.IBS!=0
    
    beta.re.Incomp <- cv.SCAD_ADMM_re(X=Model.Incomp$X, Y=Model.Incomp$Y, N=N, beta0=rep(0, dim(Model.Incomp$X)[2]), err=1e-4, tune="cv", unpen=c(1,2,3))
    indice.re.Incomp <- beta.re.Incomp!=0
    
    beta.re.AMS <- cv.SCAD_ADMM_re(X=Model.AMS$X, Y=Model.AMS$Y, N=N, beta0=rep(0, dim(Model.AMS$X)[2]), err=1e-4, tune="cv", unpen=c(1,2,3))
    indice.re.AMS <- beta.re.AMS!=0
    
    beta.re.BinMM <- cv.SCAD_ADMM_re(X=Model.BinMM$X, Y=Model.BinMM$Y, N=N, beta0=rep(0, dim(Model.BinMM$X)[2]), err=1e-4, tune="cv", unpen=c(1,2,3))
    indice.re.BinMM <- beta.re.BinMM!=0
    
    pi.unre.AMS <- logit(Model.AMS$X%*%beta.unre.AMS)
    pi.re.AMS <- logit(Model.AMS$X%*%beta.re.AMS)

    pi.unre.BinMM <- logit(Model.BinMM$X%*%beta.unre.BinMM)
    pi.re.BinMM <- logit(Model.BinMM$X%*%beta.re.BinMM)

    pi.unre.IBS <- logit(Model.IBS$X%*%beta.unre.IBS)
    pi.re.IBS <- logit(Model.IBS$X%*%beta.re.IBS)

    pi.unre.Incomp <- logit(Model.Incomp$X%*%beta.unre.Incomp)
    pi.re.Incomp <- logit(Model.Incomp$X%*%beta.re.Incomp)
    
    # construct the likelihood ratio statistics
    TL.IBS <- 2*(sum(log(1+exp(Model.IBS$X%*%beta.re.IBS))-Model.IBS$Y*(Model.IBS$X%*%beta.re.IBS))-
               sum(log(1+exp(Model.IBS$X%*%beta.unre.IBS))-Model.IBS$Y*(Model.IBS$X%*%beta.unre.IBS)))
    #if LRT stat is less than 0, there is error, so set to -1
    if(TL.IBS < 0){
      TL.IBS = -1
    }
    
    TL.Incomp <- 2*(sum(log(1+exp(Model.Incomp$X%*%beta.re.Incomp))-Model.Incomp$Y*(Model.Incomp$X%*%beta.re.Incomp))-
               sum(log(1+exp(Model.Incomp$X%*%beta.unre.Incomp))-Model.Incomp$Y*(Model.Incomp$X%*%beta.unre.Incomp)))
    #if LRT stat is less than 0, there is error, so set to -1
    if(TL.Incomp < 0){
      TL.Incomp = -1
    }
    
    TL.AMS <- 2*(sum(log(1+exp(Model.AMS$X%*%beta.re.AMS))-Model.AMS$Y*(Model.AMS$X%*%beta.re.AMS))-
               sum(log(1+exp(Model.AMS$X%*%beta.unre.AMS))-Model.AMS$Y*(Model.AMS$X%*%beta.unre.AMS)))
    #if LRT stat is less than 0, there is error, so set to -1
    if(TL.AMS < 0){
      TL.AMS = -1
    }
    
    TL.BinMM <- 2*(sum(log(1+exp(Model.BinMM$X%*%beta.re.BinMM))-Model.BinMM$Y*(Model.BinMM$X%*%beta.re.BinMM))-
               sum(log(1+exp(Model.BinMM$X%*%beta.unre.BinMM))-Model.BinMM$Y*(Model.BinMM$X%*%beta.unre.BinMM)))
    #if LRT stat is less than 0, there is error, so set to -1
    if(TL.BinMM < 0){
      TL.BinMM = -1
    }
    
    # construct the Wald statistics
    #B_0 should give you Omega_a hat
    A.IBS = crossprod(Model.IBS$X[,indice.unre.IBS|N], as.vector(pi.unre.IBS*(1-pi.unre.IBS))*Model.IBS$X[,indice.unre.IBS|N])
    if(rcond(A.IBS) >= 1e-10){
      B_0.IBS <- solve(A.IBS)
      #d_0 should determine which rows to subset
      #gives the first value that is restricted
      d_0.IBS <- sum(indice.re.IBS == TRUE) + 1
      #in case the inverse doesn't exist, need another if else
      B.IBS = B_0.IBS[d_0.IBS:ncol(B_0.IBS),d_0.IBS:ncol(B_0.IBS)]
      if(rcond(B.IBS) >= 1e-10){
        #so the B_0 needs to be subset to m rows and columns that are restricted based on H0
        TW.IBS <- crossprod(beta.unre.IBS[N], solve(B.IBS, beta.unre.IBS[N]))
      } else {
        TW.IBS = -2
      }
    } else {
      TW.IBS = -1
    }
    
    A.Incomp = crossprod(Model.Incomp$X[,indice.unre.Incomp|N], as.vector(pi.unre.Incomp*(1-pi.unre.Incomp))*Model.Incomp$X[,indice.unre.Incomp|N])
    if(rcond(A.Incomp) >= 1e-10){
      B_0.Incomp <- solve(A.Incomp)
      #d_0 should determine which rows to subset
      #gives the first value that is restricted
      d_0.Incomp <- sum(indice.re.Incomp == TRUE) + 1
      #in case the inverse doesn't exist, need another if else
      B.Incomp = B_0.Incomp[d_0.Incomp:ncol(B_0.Incomp),d_0.Incomp:ncol(B_0.Incomp)]
      if(rcond(B.Incomp) >= 1e-10){
        #so the B_0 needs to be subset to m rows and columns that are restricted based on H0
        TW.Incomp <- crossprod(beta.unre.Incomp[N], solve(B.Incomp, beta.unre.Incomp[N]))
      } else {
        TW.Incomp = -2
      }
    } else {
      TW.Incomp = -1
    }
    
    A.AMS = crossprod(Model.AMS$X[,indice.unre.AMS|N], as.vector(pi.unre.AMS*(1-pi.unre.AMS))*Model.AMS$X[,indice.unre.AMS|N])
    if(rcond(A.AMS) >= 1e-10){
      B_0.AMS <- solve(A.AMS)
      #d_0 should determine which rows to subset
      #gives the first value that is restricted
      d_0.AMS <- sum(indice.re.AMS == TRUE) + 1
      #in case the inverse doesn't exist, need another if else
      B.AMS = B_0.AMS[d_0.AMS:ncol(B_0.AMS),d_0.AMS:ncol(B_0.AMS)]
      if(rcond(B.AMS) >= 1e-10){
        #so the B_0 needs to be subset to m rows and columns that are restricted based on H0
        TW.AMS <- crossprod(beta.unre.AMS[N], solve(B.AMS, beta.unre.AMS[N]))
      } else {
        TW.AMS = -2
      }
    } else {
      TW.AMS = -1
    }
    
    A.BinMM = crossprod(Model.BinMM$X[,indice.unre.BinMM|N], as.vector(pi.unre.BinMM*(1-pi.unre.BinMM))*Model.BinMM$X[,indice.unre.BinMM|N])
    if(rcond(A.BinMM) >= 1e-10){
      B_0.BinMM <- solve(A.BinMM)
      #d_0 should determine which rows to subset
      #gives the first value that is restricted
      d_0.BinMM <- sum(indice.re.BinMM == TRUE) + 1
      #in case the inverse doesn't exist, need another if else
      B.BinMM = B_0.BinMM[d_0.BinMM:ncol(B_0.BinMM),d_0.BinMM:ncol(B_0.BinMM)]
      if(rcond(B.BinMM) >= 1e-10){
        #so the B_0 needs to be subset to m rows and columns that are restricted based on H0
        TW.BinMM <- crossprod(beta.unre.BinMM[N], solve(B.BinMM, beta.unre.BinMM[N]))
      } else {
        TW.BinMM = -2
      }
    } else {
      TW.BinMM = -1
    }
    
    # construct the score statistics
    eps.IBS <- Model.IBS$Y-pi.re.IBS #Y - e(Y)
    #this is the X^T times Y-E(Y)
    Xeps.IBS <- crossprod(Model.IBS$X[,indice.re.IBS|N], eps.IBS)
    
    C.IBS = crossprod(Model.IBS$X[,indice.re.IBS|N], as.vector(pi.re.IBS*(1-pi.re.IBS))*Model.IBS$X[,indice.re.IBS|N])
    if(rcond(C.IBS) >= 1e-10){
      TS.IBS <- crossprod(Xeps.IBS, solve(C.IBS, Xeps.IBS))
    } else {
      TS.IBS = -1
    }
    
    eps.Incomp <- Model.Incomp$Y-pi.re.Incomp #Y - e(Y)
    #this is the X^T times Y-E(Y)
    Xeps.Incomp <- crossprod(Model.Incomp$X[,indice.re.Incomp|N], eps.Incomp)
    
    C.Incomp = crossprod(Model.Incomp$X[,indice.re.Incomp|N], as.vector(pi.re.Incomp*(1-pi.re.Incomp))*Model.Incomp$X[,indice.re.Incomp|N])
    if(rcond(C.Incomp) >= 1e-10){
      TS.Incomp <- crossprod(Xeps.Incomp, solve(C.Incomp, Xeps.Incomp))
    } else {
      TS.Incomp = -1
    }
    
    eps.AMS <- Model.AMS$Y-pi.re.AMS #Y - e(Y)
    #this is the X^T times Y-E(Y)
    Xeps.AMS <- crossprod(Model.AMS$X[,indice.re.AMS|N], eps.AMS)
    
    C.AMS = crossprod(Model.AMS$X[,indice.re.AMS|N], as.vector(pi.re.AMS*(1-pi.re.AMS))*Model.AMS$X[,indice.re.AMS|N])
    if(rcond(C.AMS) >= 1e-10){
      TS.AMS <- crossprod(Xeps.AMS, solve(C.AMS, Xeps.AMS))
    } else {
      TS.AMS = -1
    }
    
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
    
    #set pv to 1 if reject, 0 otherwise
    #LRT
    if (TL.IBS>=qchisq(0.95, df = doF)){
      pv.IBS[1] <- pv.IBS[1]+1/5000
    }
    if (TL.Incomp>=qchisq(0.95, df = doF)){
      pv.Incomp[1] <- pv.Incomp[1]+1/5000
    }
    if (TL.AMS>=qchisq(0.95, df = doF)){
      pv.AMS[1] <- pv.AMS[1]+1/5000
    }
    if (TL.BinMM>=qchisq(0.95, df = doF)){
      pv.BinMM[1] <- pv.BinMM[1]+1/5000
    }
    #Wald
    if (TW.IBS>=qchisq(0.95, df = doF)){
      pv.IBS[2] <- pv.IBS[2]+1/5000
    }
    if (TW.Incomp>=qchisq(0.95, df = doF)){
      pv.Incomp[2] <- pv.Incomp[2]+1/5000
    }
    if (TW.AMS>=qchisq(0.95, df = doF)){
      pv.AMS[2] <- pv.AMS[2]+1/5000
    }
    if (TW.BinMM>=qchisq(0.95, df = doF)){
      pv.BinMM[2] <- pv.BinMM[2]+1/5000
    }
    #Score
    if (TS.IBS>=qchisq(0.95, df = doF)){
      pv.IBS[3] <- pv.IBS[3]+1/5000
    }
    if (TS.Incomp>=qchisq(0.95, df = doF)){
      pv.Incomp[3] <- pv.Incomp[3]+1/5000
    }
    if (TS.AMS>=qchisq(0.95, df = doF)){
      pv.AMS[3] <- pv.AMS[3]+1/5000
    }
    if (TS.BinMM>=qchisq(0.95, df = doF)){
      pv.BinMM[3] <- pv.BinMM[3]+1/5000
    }
    
    print(paste0("Simulation ",ii," is complete."))
    
    Tall.IBS[1, ] <- c(TL.IBS, TW.IBS, TS.IBS)
    beta.al.IBS[1, ] <- c(sum(beta.re.IBS!=0), sum(beta.unre.IBS!=0))
    
    Tall.Incomp[1, ] <- c(TL.Incomp, TW.Incomp, TS.Incomp)
    beta.al.Incomp[1, ] <- c(sum(beta.re.Incomp!=0), sum(beta.unre.Incomp!=0))
    
    Tall.AMS[1, ] <- c(TL.AMS, TW.AMS, TS.AMS)
    beta.al.AMS[1, ] <- c(sum(beta.re.AMS!=0), sum(beta.unre.AMS!=0))
    
    Tall.BinMM[1, ] <- c(TL.BinMM, TW.BinMM, TS.BinMM)
    beta.al.BinMM[1, ] <- c(sum(beta.re.BinMM!=0), sum(beta.unre.BinMM!=0))
    
    statsAndPVals = list(pv.IBS=pv.IBS, TScores.IBS=Tall.IBS, beta.IBS=beta.al.IBS,
                         pv.Incomp=pv.Incomp, TScores.Incomp=Tall.Incomp, beta.Incomp=beta.al.Incomp,
                         pv.AMS=pv.AMS, TScores.AMS=Tall.AMS, beta.AMS=beta.al.AMS,
                         pv.BinMM=pv.BinMM, TScores.BinMM=Tall.BinMM, beta.BinMM=beta.al.BinMM)
    statsAndPVals
  }, mc.cores=4)
  
  n = length(statsAndPVals)
  
  pVals.IBS = pVals.Incomp = pVals.AMS = pVals.BinMM = matrix(nrow = n, ncol = 3)
  Scores.IBS = Scores.Incomp = Scores.AMS = Scores.BinMM = matrix(nrow = n, ncol = 3)
  Betas.IBS = Betas.Incomp = Betas.AMS = Betas.BinMM = matrix(nrow = n, ncol = 2)
  
  for(ll in 1:n){
    pVals.IBS[ll,] = unlist(statsAndPVals[[ll]]$pv.IBS)
    Scores.IBS[ll,] = unlist(statsAndPVals[[ll]]$TScores.IBS)
    Betas.IBS[ll,] = unlist(statsAndPVals[[ll]]$beta.IBS)
    
    pVals.Incomp[ll,] = unlist(statsAndPVals[[ll]]$pv.Incomp)
    Scores.Incomp[ll,] = unlist(statsAndPVals[[ll]]$TScores.Incomp)
    Betas.Incomp[ll,] = unlist(statsAndPVals[[ll]]$beta.Incomp)
    
    pVals.AMS[ll,] = unlist(statsAndPVals[[ll]]$pv.AMS)
    Scores.AMS[ll,] = unlist(statsAndPVals[[ll]]$TScores.AMS)
    Betas.AMS[ll,] = unlist(statsAndPVals[[ll]]$beta.AMS)
    
    pVals.BinMM[ll,] = unlist(statsAndPVals[[ll]]$pv.BinMM)
    Scores.BinMM[ll,] = unlist(statsAndPVals[[ll]]$TScores.BinMM)
    Betas.BinMM[ll,] = unlist(statsAndPVals[[ll]]$beta.BinMM)
  }
  
  pVals = cbind(pVals.IBS, pVals.Incomp, pVals.AMS, pVals.BinMM)
  Scores = cbind(Scores.IBS, Scores.Incomp, Scores.AMS, Scores.BinMM)
  Betas = cbind(Betas.IBS, Betas.Incomp, Betas.AMS, Betas.BinMM)
  
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
  write.csv(pVals, file = paste0(path,"/Power_CatY_Prev",YPrev*100,"_LinHypTest_Score",snpOrScore,"_",percentageAssoc,"SNPsAssoc_JointTesting_OR",ORSize,"_",ld,"_Sim",start,"to",numSims,"_PValues_ForceCovFit.csv"))
  write.csv(Scores, file = paste0(path,"/Power_CatY_Prev",YPrev*100,"_LinHypTest_Score",snpOrScore,"_",percentageAssoc,"SNPsAssoc_JointTesting_OR",ORSize,"_",ld,"_Sim",start,"to",numSims,"_Stats_ForceCovFit.csv"))
  write.csv(Betas, file = paste0(path,"/Power_CatY_Prev",YPrev*100,"_LinHypTest_Score",snpOrScore,"_",percentageAssoc,"SNPsAssoc_JointTesting_OR",ORSize,"_",ld,"_Sim",start,"to",numSims,"_Betas_ForceCovFit.csv"))
  
}

# continuous outcome - linear hypothesis testing
RunPowerPipelineLinHypTestCont_Score = function(chr, gene, numPairs, Gamma, TrueScore, ORSize, standardizeScores = FALSE, weightedScores = FALSE, scoreWeights, start, stop, percentageAssoc, LowLD){
  # Function to determine power of linear hyp testing method when score is associated, joint testing
  # Only use when Y is continuous
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #Gamma = effect size for score, length 1, can be 0
  #TrueScore = IBS.gene, Incomp.gene, AMS.gene, or BinMM.gene
  #ORSize = Small, Medium, or Large for what OR was used for the associated SNP/score
  #standardizeScores = T or F whether the scores should be standardized based on maximum score value
  #weightedScores = T or F whether the scores will be weighted
  #scoreWeights = m x 1 vector of weights, one weight for each SNP
  # start - simulation number to start at
  # stop - simulation number to stop at
  #percentageAssoc = percentage of SNPs associated with outcome (either 5, 25, 50, 75, or 100) 
  #LowLD = True or FALSE whether the associated SNPs are in low LD or high LD
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes power values to csv files
  
  library(parallel)
  
  #need to define this for naming at the end
  snpOrScore = TrueScore
  
  #always the same
  numSims = stop  
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_II_Sims/HapGen_Files/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/ProjectIISourceFunctions_v2.R")
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/Linear_ADMM0.r")
  
  myList = lapply(start:numSims, rep, times = 1)
  # p-value initialization
  pv.IBS <- rep(0, 6)
  pv.Incomp <- rep(0, 6)
  pv.AMS <- rep(0, 6)
  pv.BinMM <- rep(0, 6)
  
  Tall.IBS <- matrix(0, 1, 3)
  Tall.Incomp <- matrix(0, 1, 3)
  Tall.AMS <- matrix(0, 1, 3)
  Tall.BinMM <- matrix(0, 1, 3)
  
  beta.al.IBS <- matrix(0, 1, 2)
  beta.al.Incomp <- matrix(0, 1, 2)
  beta.al.AMS <- matrix(0, 1, 2)
  beta.al.BinMM <- matrix(0, 1, 2)
  
  statsAndPVals = mclapply(myList, function(ii){
    #define matrix to hold all Stats and Pvalues
    statsAndPVals = matrix(NA, nrow = numSims, ncol = 3)
    
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
    ContPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = FALSE,  Covariates = CovData, RGenoData = RGenos, ScoreData = PhenoScore, Betas = Betas, Gamma = Gamma)
    
    # define location of zero components
    # basically, this defines what the null hypothesis is, I think
    # so for joint null, we need the non-zero components to be for covariates
    # and the beta and gamma components to be zero
    numCov = dim(CovData)[2]
    numSNPs = dim(RGenos)[2]
    N = c(rep(FALSE, numCov+1), rep(TRUE,1+numSNPs))
    
    #need to combine the covariates, r genos, and score matrices together
    #need to include column of 1s for intercept?
    intercept = matrix(1,nrow = numPairs, ncol = 1)
    designMat.IBS = cbind(intercept, CovData, RGenos, IBS.gene)
    designMat.Incomp = cbind(intercept, CovData, RGenos, Incomp.gene)
    designMat.AMS = cbind(intercept, CovData, RGenos, AMS.gene)
    designMat.BinMM = cbind(intercept, CovData, RGenos, BinMM.gene)
    
    #then combine the phenos with the design matrix as a list
    Model.IBS = list(X=designMat.IBS, Y=ContPhenos)
    Model.Incomp = list(X=designMat.Incomp, Y=ContPhenos)
    Model.AMS = list(X=designMat.AMS, Y=ContPhenos)
    Model.BinMM = list(X=designMat.BinMM, Y=ContPhenos)
    
    # estimate the uncontrained estimator
    beta.unre.IBS <- cv.SCAD_ADMM_unre(X=Model.IBS$X, Y=Model.IBS$Y, N=N, beta0=rep(0, dim(Model.IBS$X)[2]), err=1e-4, tune="cv", unpen = c(1,2,3))
    indice.unre.IBS <- beta.unre.IBS!=0
    
    beta.unre.Incomp <- cv.SCAD_ADMM_unre(X=Model.Incomp$X, Y=Model.Incomp$Y, N=N, beta0=rep(0, dim(Model.Incomp$X)[2]), err=1e-4, tune="cv", unpen = c(1,2,3))
    indice.unre.Incomp <- beta.unre.Incomp!=0
    
    beta.unre.AMS <- cv.SCAD_ADMM_unre(X=Model.AMS$X, Y=Model.AMS$Y, N=N, beta0=rep(0, dim(Model.AMS$X)[2]), err=1e-4, tune="cv", unpen = c(1,2,3))
    indice.unre.AMS <- beta.unre.AMS!=0
    
    beta.unre.BinMM <- cv.SCAD_ADMM_unre(X=Model.BinMM$X, Y=Model.BinMM$Y, N=N, beta0=rep(0, dim(Model.BinMM$X)[2]), err=1e-4, tune="cv", unpen = c(1,2,3))
    indice.unre.BinMM <- beta.unre.BinMM!=0
    
    # estimate the constrained estimator
    beta.re.IBS <- cv.SCAD_ADMM_re(X=Model.IBS$X, Y=Model.IBS$Y, N=N, beta0=rep(0, dim(Model.IBS$X)[2]), err=1e-4, tune="cv", unpen = c(1,2,3))
    indice.re.IBS <- beta.re.IBS!=0
    
    beta.re.Incomp <- cv.SCAD_ADMM_re(X=Model.Incomp$X, Y=Model.Incomp$Y, N=N, beta0=rep(0, dim(Model.Incomp$X)[2]), err=1e-4, tune="cv", unpen = c(1,2,3))
    indice.re.Incomp <- beta.re.Incomp!=0
    
    beta.re.AMS <- cv.SCAD_ADMM_re(X=Model.AMS$X, Y=Model.AMS$Y, N=N, beta0=rep(0, dim(Model.AMS$X)[2]), err=1e-4, tune="cv", unpen = c(1,2,3))
    indice.re.AMS <- beta.re.AMS!=0
    
    beta.re.BinMM <- cv.SCAD_ADMM_re(X=Model.BinMM$X, Y=Model.BinMM$Y, N=N, beta0=rep(0, dim(Model.BinMM$X)[2]), err=1e-4, tune="cv", unpen = c(1,2,3))
    indice.re.BinMM <- beta.re.BinMM!=0
    
    # estimate the conditional variance
    #    sig2 <- mean((Model$Y-Model$X%*%beta.unre)^2)
    n = numSims
    sig2.IBS <- mean((Model.IBS$Y-Model.IBS$X%*%beta.unre.IBS)^2)*n/(n-sum(beta.unre.IBS!=0))
    sig2.Incomp <- mean((Model.Incomp$Y-Model.Incomp$X%*%beta.unre.Incomp)^2)*n/(n-sum(beta.unre.Incomp!=0))
    sig2.AMS <- mean((Model.AMS$Y-Model.AMS$X%*%beta.unre.AMS)^2)*n/(n-sum(beta.unre.AMS!=0))
    sig2.BinMM <- mean((Model.BinMM$Y-Model.BinMM$X%*%beta.unre.BinMM)^2)*n/(n-sum(beta.unre.BinMM!=0))
    
    # construct the likelihood ratio statistic
    TL.IBS <- sum((Model.IBS$X%*%beta.re.IBS-Model.IBS$Y)^2)-sum((Model.IBS$X%*%beta.unre.IBS-Model.IBS$Y)^2)
    #if LRT stat is less than 0, there is error, so set to -1
    if(TL.IBS < 0){
      TL.IBS = -1
    }
    TL.Incomp <- sum((Model.Incomp$X%*%beta.re.Incomp-Model.Incomp$Y)^2)-sum((Model.Incomp$X%*%beta.unre.Incomp-Model.Incomp$Y)^2)
    #if LRT stat is less than 0, there is error, so set to -1
    if(TL.Incomp < 0){
      TL.Incomp = -1
    }
    TL.AMS <- sum((Model.AMS$X%*%beta.re.AMS-Model.AMS$Y)^2)-sum((Model.AMS$X%*%beta.unre.AMS-Model.AMS$Y)^2)
    #if LRT stat is less than 0, there is error, so set to -1
    if(TL.AMS < 0){
      TL.AMS = -1
    }
    TL.BinMM <- sum((Model.BinMM$X%*%beta.re.BinMM-Model.BinMM$Y)^2)-sum((Model.BinMM$X%*%beta.unre.BinMM-Model.BinMM$Y)^2)
    #if LRT stat is less than 0, there is error, so set to -1
    if(TL.BinMM < 0){
      TL.BinMM = -1
    }
    
    # construct the Wald statistic
    #B_0 should give you Omega_a hat
    B.IBS = crossprod(Model.IBS$X[,indice.unre.IBS|N], Model.IBS$X[,indice.unre.IBS|N])
    if(rcond(B.IBS) >= 1e-10){
      B_0.IBS <- solve(B.IBS)
      #d_0.IBS should determine which rows to subset
      #gives the first value that is restricted
      d_0.IBS <- sum(indice.re.IBS == TRUE) + 1
      #in case the inverse doesn't exist, need another if else
      A.IBS = B_0.IBS[d_0.IBS:ncol(B_0.IBS),d_0.IBS:ncol(B_0.IBS)]
      if(rcond(A.IBS) >= 1e-10){
        TW.IBS <- crossprod(beta.unre.IBS[N], solve(A.IBS, beta.unre.IBS[N]))
      } else {
        TW.IBS = -2
      }
    } else {
      TW.IBS = -1
    }
    
    B.Incomp = crossprod(Model.Incomp$X[,indice.unre.Incomp|N], Model.Incomp$X[,indice.unre.Incomp|N])
    if(rcond(B.Incomp) >= 1e-10){
      B_0.Incomp <- solve(B.Incomp)
      #d_0.Incomp should determine which rows to subset
      #gives the first value that is restricted
      d_0.Incomp <- sum(indice.re.Incomp == TRUE) + 1
      #in case the inverse doesn't exist, need another if else
      A.Incomp = B_0.Incomp[d_0.Incomp:ncol(B_0.Incomp),d_0.Incomp:ncol(B_0.Incomp)]
      if(rcond(A.Incomp) >= 1e-10){
        TW.Incomp <- crossprod(beta.unre.Incomp[N], solve(A.Incomp, beta.unre.Incomp[N]))
      } else {
        TW.Incomp = -2
      }
    } else {
      TW.Incomp = -1
    }
    
    B.AMS = crossprod(Model.AMS$X[,indice.unre.AMS|N], Model.AMS$X[,indice.unre.AMS|N])
    if(rcond(B.AMS) >= 1e-10){
      B_0.AMS <- solve(B.AMS)
      #d_0.AMS should determine which rows to subset
      #gives the first value that is restricted
      d_0.AMS <- sum(indice.re.AMS == TRUE) + 1
      #in case the inverse doesn't exist, need another if else
      A.AMS = B_0.AMS[d_0.AMS:ncol(B_0.AMS),d_0.AMS:ncol(B_0.AMS)]
      if(rcond(A.AMS) >= 1e-10){
        TW.AMS <- crossprod(beta.unre.AMS[N], solve(A.AMS, beta.unre.AMS[N]))
      } else {
        TW.AMS = -2
      }
    } else {
      TW.AMS = -1
    }
    
    B.BinMM = crossprod(Model.BinMM$X[,indice.unre.BinMM|N], Model.BinMM$X[,indice.unre.BinMM|N])
    if(rcond(B.BinMM) >= 1e-10){
      B_0.BinMM <- solve(B.BinMM)
      #d_0.BinMM should determine which rows to subset
      #gives the first value that is restricted
      d_0.BinMM <- sum(indice.re.BinMM == TRUE) + 1
      #in case the inverse doesn't exist, need another if else
      A.BinMM = B_0.BinMM[d_0.BinMM:ncol(B_0.BinMM),d_0.BinMM:ncol(B_0.BinMM)]
      if(rcond(A.BinMM) >= 1e-10){
        TW.BinMM <- crossprod(beta.unre.BinMM[N], solve(A.BinMM, beta.unre.BinMM[N]))
      } else {
        TW.BinMM = -2
      }
    } else {
      TW.BinMM = -1
    }
    
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
    
    #LRT
    if (TL.IBS>=qchisq(0.95, df = doF)){
      pv.IBS[1] <- pv.IBS[1]+1/5000
    }
    if (TL.Incomp>=qchisq(0.95, df = doF)){
      pv.Incomp[1] <- pv.Incomp[1]+1/5000
    }
    if (TL.AMS>=qchisq(0.95, df = doF)){
      pv.AMS[1] <- pv.AMS[1]+1/5000
    }
    if (TL.BinMM>=qchisq(0.95, df = doF)){
      pv.BinMM[1] <- pv.BinMM[1]+1/5000
    }
    #Wald
    if (TW.IBS>=qchisq(0.95, df = doF)){
      pv.IBS[2] <- pv.IBS[2]+1/5000
    }
    if (TW.Incomp>=qchisq(0.95, df = doF)){
      pv.Incomp[2] <- pv.Incomp[2]+1/5000
    }
    if (TW.AMS>=qchisq(0.95, df = doF)){
      pv.AMS[2] <- pv.AMS[2]+1/5000
    }
    if (TW.BinMM>=qchisq(0.95, df = doF)){
      pv.BinMM[2] <- pv.BinMM[2]+1/5000
    }
    #Score
    if (TS.IBS>=qchisq(0.95, df = doF)){
      pv.IBS[3] <- pv.IBS[3]+1/5000
    }
    if (TS.Incomp>=qchisq(0.95, df = doF)){
      pv.Incomp[3] <- pv.Incomp[3]+1/5000
    }
    if (TS.AMS>=qchisq(0.95, df = doF)){
      pv.AMS[3] <- pv.AMS[3]+1/5000
    }
    if (TS.BinMM>=qchisq(0.95, df = doF)){
      pv.BinMM[3] <- pv.BinMM[3]+1/5000
    }
    
    print(paste0("Simulation ",ii," is complete."))
    
    Tall.IBS[1, ] <- c(TL.IBS, TW.IBS, TS.IBS)
    beta.al.IBS[1, ] <- c(sum(beta.re.IBS!=0), sum(beta.unre.IBS!=0))
    
    Tall.Incomp[1, ] <- c(TL.Incomp, TW.Incomp, TS.Incomp)
    beta.al.Incomp[1, ] <- c(sum(beta.re.Incomp!=0), sum(beta.unre.Incomp!=0))
    
    Tall.AMS[1, ] <- c(TL.AMS, TW.AMS, TS.AMS)
    beta.al.AMS[1, ] <- c(sum(beta.re.AMS!=0), sum(beta.unre.AMS!=0))
    
    Tall.BinMM[1, ] <- c(TL.BinMM, TW.BinMM, TS.BinMM)
    beta.al.BinMM[1, ] <- c(sum(beta.re.BinMM!=0), sum(beta.unre.BinMM!=0))
    
    statsAndPVals = list(pv.IBS=pv.IBS, TScores.IBS=Tall.IBS, beta.IBS=beta.al.IBS,
                         pv.Incomp=pv.Incomp, TScores.Incomp=Tall.Incomp, beta.Incomp=beta.al.Incomp,
                         pv.AMS=pv.AMS, TScores.AMS=Tall.AMS, beta.AMS=beta.al.AMS,
                         pv.BinMM=pv.BinMM, TScores.BinMM=Tall.BinMM, beta.BinMM=beta.al.BinMM)
    statsAndPVals
  }, mc.cores=4)
  
  n = length(statsAndPVals)
  
  pVals.IBS = pVals.Incomp = pVals.AMS = pVals.BinMM = matrix(nrow = n, ncol = 6)
  Scores.IBS = Scores.Incomp = Scores.AMS = Scores.BinMM = matrix(nrow = n, ncol = 3)
  Betas.IBS = Betas.Incomp = Betas.AMS = Betas.BinMM = matrix(nrow = n, ncol = 2)
  
  for(ll in 1:n){
    pVals.IBS[ll,] = unlist(statsAndPVals[[ll]]$pv.IBS)
    Scores.IBS[ll,] = unlist(statsAndPVals[[ll]]$TScores.IBS)
    Betas.IBS[ll,] = unlist(statsAndPVals[[ll]]$beta.IBS)
    
    pVals.Incomp[ll,] = unlist(statsAndPVals[[ll]]$pv.Incomp)
    Scores.Incomp[ll,] = unlist(statsAndPVals[[ll]]$TScores.Incomp)
    Betas.Incomp[ll,] = unlist(statsAndPVals[[ll]]$beta.Incomp)
    
    pVals.AMS[ll,] = unlist(statsAndPVals[[ll]]$pv.AMS)
    Scores.AMS[ll,] = unlist(statsAndPVals[[ll]]$TScores.AMS)
    Betas.AMS[ll,] = unlist(statsAndPVals[[ll]]$beta.AMS)
    
    pVals.BinMM[ll,] = unlist(statsAndPVals[[ll]]$pv.BinMM)
    Scores.BinMM[ll,] = unlist(statsAndPVals[[ll]]$TScores.BinMM)
    Betas.BinMM[ll,] = unlist(statsAndPVals[[ll]]$beta.BinMM)
  }
  
  pVals = cbind(pVals.IBS, pVals.Incomp, pVals.AMS, pVals.BinMM)
  Scores = cbind(Scores.IBS, Scores.Incomp, Scores.AMS, Scores.BinMM)
  Betas = cbind(Betas.IBS, Betas.Incomp, Betas.AMS, Betas.BinMM)
  
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
  write.csv(pVals, file = paste0(path,"/Power_ContY_LinHypTest_Score",snpOrScore,"_",percentageAssoc,"SNPsAssoc_JointTesting_OR",ORSize,"_",ld,"_Sim",start,"to",numSims,"_PValues_ForceCovFit.csv"))
  write.csv(Scores, file = paste0(path,"/Power_ContY_LinHypTest_Score",snpOrScore,"_",percentageAssoc,"SNPsAssoc_JointTesting_OR",ORSize,"_",ld,"_Sim",start,"to",numSims,"_Stats_ForceCovFit.csv"))
  write.csv(Betas, file = paste0(path,"/Power_ContY_LinHypTest_Score",snpOrScore,"_",percentageAssoc,"SNPsAssoc_JointTesting_OR",ORSize,"_",ld,"_Sim",start,"to",numSims,"_Betas_ForceCovFit.csv"))
}

#RSNP is associated, Score is not
# binary outcome - standard GLM
RunPowerPipelineGLMCat_RSNP = function(chr, gene, numPairs, YPrev, Gamma, TrueScore, ORSize, standardizeScores = FALSE, weightedScores = FALSE, scoreWeights, start, stop, percentageAssoc, LowLD){
  # Function to determine power of standard GLM method when score is associated, joint testing
  # Only use when Y is binary
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #Gamma = effect size for score, length 1, can be 0
  #TrueScore = IBS.gene, Incomp.gene, AMS.gene, or BinMM.gene
  #ORSize = Small, Medium, or Large for what OR was used for the associated SNP/score
  #standardizeScores = T or F whether the scores should be standardized based on maximum score value
  #weightedScores = T or F whether the scores will be weighted
  #scoreWeights = m x 1 vector of weights, one weight for each SNP
  # start - simulation number to start at
  # stop - simulation number to stop at
  #percentageAssoc = percentage of SNPs associated with outcome (either 5, 25, 50, 75, or 100) 
  #LowLD = True or FALSE whether the associated SNPs are in low LD or high LD
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes power values to csv files
  
  library(parallel)
  suppressMessages(library(epicalc))
  
  #need to define this for naming at the end
  snpOrScore = "RSNP"
  
  #always the same
  numSims = stop  
  
  #define effect based on OR size
  if(ORSize == "small"){
    effect = 0.14
  } else if(ORSize == "medium"){
    effect = 0.41
  } else {
    effect = 0.69
  }
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_II_Sims/HapGen_Files/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/ProjectIISourceFunctions_v2.R")
  
  #determine which SNPs to actually set as assoc. based on gene
  #check gene
  assocSNPs = DetermineAssocRSNPs(gene = gene, LowLD = LowLD, percentageAssoc = percentageAssoc)
  
  #turn output matrices into output lists
  myList = lapply(start:numSims, rep, times = 1)
  
  statsAndPVals = mclapply(myList, function(ii){
    #define matrix to hold all Stats and Pvalues
    statsAndPVals = matrix(NA, nrow = numSims, ncol = 3)
    
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
    CatPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = TRUE, YPrev = YPrev,  Covariates = CovData, RGenoData = RGenos, ScoreData = PhenoScore, Betas = Betas, Gamma = Gamma)
    
    #fit the null and alternative models
    fitNull = glm(CatPhenos~CovData, family = binomial)
    fitAltIBS = glm(CatPhenos~CovData+RGenos+IBS.gene, family = binomial, epsilon = 1e-6)	
    fitAltIncomp = glm(CatPhenos~CovData+RGenos+Incomp.gene, family = binomial, epsilon = 1e-6)	
    fitAltAMS = glm(CatPhenos~CovData+RGenos+AMS.gene, family = binomial, epsilon = 1e-6)	
    fitAltBinMM = glm(CatPhenos~CovData+RGenos+BinMM.gene, family = binomial, epsilon = 1e-6)	
    
    if(fitNull$deviance - fitAltIBS$deviance >= 0){
      # construct the likelihood ratio statistics
      TL.IBS = lrtest(fitNull, fitAltIBS)
      
      pv.IBS = TL.IBS$p.value
      stat.IBS = TL.IBS$Chisquared
      doF.IBS = TL.IBS$df
    } else {
      pv.IBS = 1
      stat.IBS = 0
      doF.IBS = (dim(RGenos)[2]+1)
    }
    
    if(fitNull$deviance - fitAltIncomp$deviance >= 0){
      # construct the likelihood ratio statistics
      TL.Incomp = lrtest(fitNull, fitAltIncomp)
      
      pv.Incomp = TL.Incomp$p.value
      stat.Incomp = TL.Incomp$Chisquared
      doF.Incomp = TL.Incomp$df
    } else {
      pv.Incomp = 1
      stat.Incomp = 0
      doF.Incomp = (dim(RGenos)[2]+1)
    }
    
    if(fitNull$deviance - fitAltAMS$deviance >= 0){
      # construct the likelihood ratio statistics
      TL.AMS = lrtest(fitNull, fitAltAMS)
      
      pv.AMS = TL.AMS$p.value
      stat.AMS = TL.AMS$Chisquared
      doF.AMS = TL.AMS$df
    } else {
      pv.AMS = 1
      stat.AMS = 0
      doF.AMS = (dim(RGenos)[2]+1)
    }
    
    if(fitNull$deviance - fitAltBinMM$deviance >= 0){
      # construct the likelihood ratio statistics
      TL.BinMM = lrtest(fitNull, fitAltBinMM)
      
      pv.BinMM = TL.BinMM$p.value
      stat.BinMM = TL.BinMM$Chisquared
      doF.BinMM = TL.BinMM$df
    } else {
      pv.BinMM = 1
      stat.BinMM = 0
      doF.BinMM = (dim(RGenos)[2]+1)
    }
    
    print(paste0("Simulation ",ii," is complete."))
    
    pv = c(pv.IBS, pv.Incomp, pv.AMS, pv.BinMM)
    stat = c(stat.IBS, stat.Incomp, stat.AMS, stat.BinMM)
    doF = c(doF.IBS, doF.Incomp, doF.AMS, doF.BinMM)
    
    statsAndPVals = list(pv=pv, TL=stat, dof=doF)
    statsAndPVals
  }, mc.cores=4)
  
  n = length(statsAndPVals)
  
  pVals = matrix(nrow = n, ncol = 4)
  Scores = matrix(nrow = n, ncol = 4)
  doF = matrix(nrow = n, ncol = 4)
  
  for(ll in 1:n){
    pVals[ll,] = unlist(statsAndPVals[[ll]]$pv)
    Scores[ll,] = unlist(statsAndPVals[[ll]]$TL)
    doF[ll,] = unlist(statsAndPVals[[ll]]$dof)
  }
  
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
  write.csv(pVals, file = paste0(path,"/Power_CatPhenos_Prev",YPrev*100,"_GLM_Pvals_",percentageAssoc,"SNPsAssoc_",ld,"_TrueScore",TrueScore,"_",ORSize,"OR_assocScore",snpOrScore,"_Sim",start,"to",numSims,".csv"))
  write.csv(Scores, file = paste0(path,"/Power_CatPhenos_Prev",YPrev*100,"_GLM_Scores_",percentageAssoc,"SNPsAssoc_",ld,"_TrueScore",TrueScore,"_",ORSize,"OR_assocScore",snpOrScore,"_Sim",start,"to",numSims,".csv"))
  write.csv(doF, file = paste0(path,"/Power_CatPhenos_Prev",YPrev*100,"_GLM_dof_",percentageAssoc,"SNPsAssoc_",ld,"_TrueScore",TrueScore,"_",ORSize,"OR_assocScore",snpOrScore,"_Sim",start,"to",numSims,".csv"))
}

# continuous outcome - standard GLM
RunPowerPipelineGLMCont_RSNP = function(chr, gene, numPairs, Gamma, TrueScore, ORSize, standardizeScores = FALSE, weightedScores = FALSE, scoreWeights, start, stop, percentageAssoc, LowLD){
  # Function to determine power of standard GLM method when score is associated, joint testing
  # only use when Y is continuous (Normal)
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #Gamma = effect size for score, length 1, can be 0
  #TrueScore = IBS.gene, Incomp.gene, AMS.gene, or BinMM.gene
  #ORSize = Small, Medium, or Large for what OR was used for the associated SNP/score
  #standardizeScores = T or F whether the scores should be standardized based on maximum score value
  #weightedScores = T or F whether the scores will be weighted
  #scoreWeights = m x 1 vector of weights, one weight for each SNP
  #percentageAssoc = percentage of SNPs associated with outcome (either 5, 25, 50, 75, or 100) 
  #LowLD = True or FALSE whether the associated SNPs are in low LD or high LD
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes power values to csv files
  
  library(parallel)
  suppressMessages(library(lmtest))
  
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
  
  #always the same
  numSims = stop  
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_II_Sims/HapGen_Files/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/ProjectIISourceFunctions_v2.R")
  
  #determine which SNPs to actually set as assoc. based on gene
  #check gene
  assocSNPs = DetermineAssocRSNPs(gene = gene, LowLD = LowLD, percentageAssoc = percentageAssoc)
  
  myList = lapply(start:numSims, rep, times = 1)
  
  statsAndPVals = mclapply(myList, function(ii){
    #define matrix to hold all Stats and Pvalues
    statsAndPVals = matrix(NA, nrow = numSims, ncol = 3)
    
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
    ContPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = FALSE,  Covariates = CovData, RGenoData = RGenos, ScoreData = PhenoScore, Betas = Betas, Gamma = Gamma)
    
    #fit the null and alternative models
    fitNull = glm(ContPhenos~CovData, family = gaussian)
    fitAlt.IBS = glm(ContPhenos~CovData+RGenos+IBS.gene, family = gaussian, epsilon = 1e-6)	
    fitAlt.Incomp = glm(ContPhenos~CovData+RGenos+Incomp.gene, family = gaussian, epsilon = 1e-6)	
    fitAlt.AMS = glm(ContPhenos~CovData+RGenos+AMS.gene, family = gaussian, epsilon = 1e-6)	
    fitAlt.BinMM = glm(ContPhenos~CovData+RGenos+BinMM.gene, family = gaussian, epsilon = 1e-6)	
    
    if(fitNull$deviance - fitAlt.IBS$deviance >= 0){
      # construct the likelihood ratio statistics
      TL.IBS = lrtest(fitNull, fitAlt.IBS)
      
      pv.IBS = TL.IBS[2,5]
      stat.IBS = TL.IBS[2,4]
      doF.IBS = TL.IBS[2,1]
    } else {
      pv.IBS = 1
      stat.IBS = 0
      doF.IBS = (dim(RGenos)[2]+1)
    }
    if(fitNull$deviance - fitAlt.Incomp$deviance >= 0){
      # construct the likelihood ratio statistics
      TL.Incomp = lrtest(fitNull, fitAlt.Incomp)
      
      pv.Incomp = TL.Incomp[2,5]
      stat.Incomp = TL.Incomp[2,4]
      doF.Incomp = TL.Incomp[2,1]
    } else {
      pv.Incomp = 1
      stat.Incomp = 0
      doF.Incomp = (dim(RGenos)[2]+1)
    }
    if(fitNull$deviance - fitAlt.AMS$deviance >= 0){
      # construct the likelihood ratio statistics
      TL.AMS = lrtest(fitNull, fitAlt.AMS)
      
      pv.AMS = TL.AMS[2,5]
      stat.AMS = TL.AMS[2,4]
      doF.AMS = TL.AMS[2,1]
    } else {
      pv.AMS = 1
      stat.AMS = 0
      doF.AMS = (dim(RGenos)[2]+1)
    }
    if(fitNull$deviance - fitAlt.BinMM$deviance >= 0){
      # construct the likelihood ratio statistics
      TL.BinMM = lrtest(fitNull, fitAlt.BinMM)
      
      pv.BinMM = TL.BinMM[2,5]
      stat.BinMM = TL.BinMM[2,4]
      doF.BinMM = TL.BinMM[2,1]
    } else {
      pv.BinMM = 1
      stat.BinMM = 0
      doF.BinMM = (dim(RGenos)[2]+1)
    }
    
    print(paste0("Simulation ",ii," is complete."))
    
    pv = c(pv.IBS, pv.Incomp, pv.AMS, pv.BinMM)
    stat = c(stat.IBS, stat.Incomp, stat.AMS, stat.BinMM)
    doF = c(doF.IBS, doF.Incomp, doF.AMS, doF.BinMM)
    
    statsAndPVals = list(pv=pv, TL=stat, dof=doF)
    statsAndPVals
  }, mc.cores=4)
  
  n = length(statsAndPVals)
  
  pVals = matrix(nrow = n, ncol = 4)
  Scores = matrix(nrow = n, ncol = 4)
  doF = matrix(nrow = n, ncol = 4)
  
  for(ll in 1:n){
    pVals[ll,] = unlist(statsAndPVals[[ll]]$pv)
    Scores[ll,] = unlist(statsAndPVals[[ll]]$TL)
    doF[ll,] = unlist(statsAndPVals[[ll]]$dof)
  }
  
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
  write.csv(pVals, file = paste0(path,"/Power_ContPhenos_GLM_Pvals_",percentageAssoc,"SNPsAssoc_",ld,"_TrueScore",TrueScore,"_",ORSize,"OR_assocScore",snpOrScore,"_Sim",start,"to",numSims,".csv"))
  write.csv(Scores, file = paste0(path,"/Power_ContPhenos_GLM_Scores_",percentageAssoc,"SNPsAssoc_",ld,"_TrueScore",TrueScore,"_",ORSize,"OR_assocScore",snpOrScore,"_Sim",start,"to",numSims,".csv"))
  write.csv(doF, file = paste0(path,"/Power_ContPhenos_GLM_dof_",percentageAssoc,"SNPsAssoc_",ld,"_TrueScore",TrueScore,"_",ORSize,"OR_assocScore",snpOrScore,"_Sim",start,"to",numSims,".csv"))
}

# binary outcome - linear hypothesis testing
RunPowerPipelineLinHypTestCat_RSNP = function(chr, gene, numPairs, YPrev, Gamma, TrueScore, ORSize, standardizeScores = FALSE, weightedScores = FALSE, scoreWeights, start, stop, percentageAssoc, LowLD){
  # Function to determine power of linear hyp testing method when score is associated, joint testing
  # Only use when Y is binary
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #Gamma = effect size for score, length 1, can be 0
  #TrueScore = IBS.gene, Incomp.gene, AMS.gene, or BinMM.gene
  #ORSize = Small, Medium, or Large for what OR was used for the associated SNP/score
  #standardizeScores = T or F whether the scores should be standardized based on maximum score value
  #weightedScores = T or F whether the scores will be weighted
  #scoreWeights = m x 1 vector of weights, one weight for each SNP
  # start - simulation number to start at
  # stop - simulation number to stop at
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
  
  #always the same
  numSims = stop  
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_II_Sims/HapGen_Files/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/ProjectIISourceFunctions_v2.R")
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/Logistic_ADMM0.r")
  
  #determine which SNPs to actually set as assoc. based on gene
  #check gene
  assocSNPs = DetermineAssocRSNPs(gene = gene, LowLD = LowLD, percentageAssoc = percentageAssoc)
  
  myList = lapply(start:numSims, rep, times = 1)
  
  # p-value initialization
  pv.IBS <- rep(0, 3)
  pv.Incomp <- rep(0, 3)
  pv.AMS <- rep(0, 3)
  pv.BinMM <- rep(0, 3)
  
  Tall.IBS <- matrix(0, 1, 3)
  Tall.Incomp <- matrix(0, 1, 3)
  Tall.AMS <- matrix(0, 1, 3)
  Tall.BinMM <- matrix(0, 1, 3)
  
  beta.al.IBS <- matrix(0, 1, 2)
  beta.al.Incomp <- matrix(0, 1, 2)
  beta.al.AMS <- matrix(0, 1, 2)
  beta.al.BinMM <- matrix(0, 1, 2)
  
  statsAndPVals = mclapply(myList, function(ii){
    #define matrix to hold all Stats and Pvalues
    statsAndPVals = matrix(NA, nrow = numSims, ncol = 3)
    
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
    CatPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = TRUE, YPrev = YPrev,  Covariates = CovData, RGenoData = RGenos, ScoreData = PhenoScore, Betas = Betas, Gamma = Gamma)
    
    # define location of zero components
    # basically, this defines what the null hypothesis is, I think
    # so for joint null, we need the non-zero components to be for covariates
    # and the beta and gamma components to be zero
    numCov = dim(CovData)[2]
    numSNPs = dim(RGenos)[2]
    N = c(rep(FALSE, numCov+1), rep(TRUE,1+numSNPs))
    
    #need to combine the covariates, r genos, and score matrices together
    #need to include column of 1s for intercept?
    intercept = matrix(1,nrow = numPairs, ncol = 1)
    designMat.IBS = cbind(intercept, CovData, RGenos, IBS.gene)
    designMat.Incomp = cbind(intercept, CovData, RGenos, Incomp.gene)
    designMat.AMS = cbind(intercept, CovData, RGenos, AMS.gene)
    designMat.BinMM = cbind(intercept, CovData, RGenos, BinMM.gene)
    
    #then combine the phenos with the design matrix as a list
    Model.IBS = list(X=designMat.IBS, Y=CatPhenos)
    Model.Incomp = list(X=designMat.Incomp, Y=CatPhenos)
    Model.AMS = list(X=designMat.AMS, Y=CatPhenos)
    Model.BinMM = list(X=designMat.BinMM, Y=CatPhenos)
    
    # estimate the uncontrained estimator
    beta.unre.IBS <- cv.SCAD_ADMM_unre(X=Model.IBS$X, Y=Model.IBS$Y, N=N, beta0=rep(0, dim(Model.IBS$X)[2]), err=1e-4, tune="cv", unpen=c(1,2,3))
    indice.unre.IBS <- beta.unre.IBS!=0
    
    beta.unre.Incomp <- cv.SCAD_ADMM_unre(X=Model.Incomp$X, Y=Model.Incomp$Y, N=N, beta0=rep(0, dim(Model.Incomp$X)[2]), err=1e-4, tune="cv", unpen=c(1,2,3))
    indice.unre.Incomp <- beta.unre.Incomp!=0
    
    beta.unre.AMS <- cv.SCAD_ADMM_unre(X=Model.AMS$X, Y=Model.AMS$Y, N=N, beta0=rep(0, dim(Model.AMS$X)[2]), err=1e-4, tune="cv", unpen=c(1,2,3))
    indice.unre.AMS <- beta.unre.AMS!=0
    
    beta.unre.BinMM <- cv.SCAD_ADMM_unre(X=Model.BinMM$X, Y=Model.BinMM$Y, N=N, beta0=rep(0, dim(Model.BinMM$X)[2]), err=1e-4, tune="cv", unpen=c(1,2,3))
    indice.unre.BinMM <- beta.unre.BinMM!=0
    
    # estimate the constrained estimator (only need this for score test)
    beta.re.IBS <- cv.SCAD_ADMM_re(X=Model.IBS$X, Y=Model.IBS$Y, N=N, beta0=rep(0, dim(Model.IBS$X)[2]), err=1e-4, tune="cv", unpen=c(1,2,3))
    indice.re.IBS <- beta.re.IBS!=0
    
    beta.re.Incomp <- cv.SCAD_ADMM_re(X=Model.Incomp$X, Y=Model.Incomp$Y, N=N, beta0=rep(0, dim(Model.Incomp$X)[2]), err=1e-4, tune="cv", unpen=c(1,2,3))
    indice.re.Incomp <- beta.re.Incomp!=0
    
    beta.re.AMS <- cv.SCAD_ADMM_re(X=Model.AMS$X, Y=Model.AMS$Y, N=N, beta0=rep(0, dim(Model.AMS$X)[2]), err=1e-4, tune="cv", unpen=c(1,2,3))
    indice.re.AMS <- beta.re.AMS!=0
    
    beta.re.BinMM <- cv.SCAD_ADMM_re(X=Model.BinMM$X, Y=Model.BinMM$Y, N=N, beta0=rep(0, dim(Model.BinMM$X)[2]), err=1e-4, tune="cv", unpen=c(1,2,3))
    indice.re.BinMM <- beta.re.BinMM!=0
    
    pi.unre.AMS <- logit(Model.AMS$X%*%beta.unre.AMS)
    pi.re.AMS <- logit(Model.AMS$X%*%beta.re.AMS)
    
    pi.unre.BinMM <- logit(Model.BinMM$X%*%beta.unre.BinMM)
    pi.re.BinMM <- logit(Model.BinMM$X%*%beta.re.BinMM)
    
    pi.unre.IBS <- logit(Model.IBS$X%*%beta.unre.IBS)
    pi.re.IBS <- logit(Model.IBS$X%*%beta.re.IBS)
    
    pi.unre.Incomp <- logit(Model.Incomp$X%*%beta.unre.Incomp)
    pi.re.Incomp <- logit(Model.Incomp$X%*%beta.re.Incomp)
    
    # construct the likelihood ratio statistics
    TL.IBS <- 2*(sum(log(1+exp(Model.IBS$X%*%beta.re.IBS))-Model.IBS$Y*(Model.IBS$X%*%beta.re.IBS))-
                   sum(log(1+exp(Model.IBS$X%*%beta.unre.IBS))-Model.IBS$Y*(Model.IBS$X%*%beta.unre.IBS)))
    #if LRT stat is less than 0, there is error, so set to -1
    if(TL.IBS < 0){
      TL.IBS = -1
    }
    
    TL.Incomp <- 2*(sum(log(1+exp(Model.Incomp$X%*%beta.re.Incomp))-Model.Incomp$Y*(Model.Incomp$X%*%beta.re.Incomp))-
                      sum(log(1+exp(Model.Incomp$X%*%beta.unre.Incomp))-Model.Incomp$Y*(Model.Incomp$X%*%beta.unre.Incomp)))
    #if LRT stat is less than 0, there is error, so set to -1
    if(TL.Incomp < 0){
      TL.Incomp = -1
    }
    
    TL.AMS <- 2*(sum(log(1+exp(Model.AMS$X%*%beta.re.AMS))-Model.AMS$Y*(Model.AMS$X%*%beta.re.AMS))-
                   sum(log(1+exp(Model.AMS$X%*%beta.unre.AMS))-Model.AMS$Y*(Model.AMS$X%*%beta.unre.AMS)))
    #if LRT stat is less than 0, there is error, so set to -1
    if(TL.AMS < 0){
      TL.AMS = -1
    }
    
    TL.BinMM <- 2*(sum(log(1+exp(Model.BinMM$X%*%beta.re.BinMM))-Model.BinMM$Y*(Model.BinMM$X%*%beta.re.BinMM))-
                     sum(log(1+exp(Model.BinMM$X%*%beta.unre.BinMM))-Model.BinMM$Y*(Model.BinMM$X%*%beta.unre.BinMM)))
    #if LRT stat is less than 0, there is error, so set to -1
    if(TL.BinMM < 0){
      TL.BinMM = -1
    }
    
    # construct the Wald statistics
    #B_0 should give you Omega_a hat
    A.IBS = crossprod(Model.IBS$X[,indice.unre.IBS|N], as.vector(pi.unre.IBS*(1-pi.unre.IBS))*Model.IBS$X[,indice.unre.IBS|N])
    if(rcond(A.IBS) >= 1e-10){
      B_0.IBS <- solve(A.IBS)
      #d_0 should determine which rows to subset
      #gives the first value that is restricted
      d_0.IBS <- sum(indice.re.IBS == TRUE) + 1
      #in case the inverse doesn't exist, need another if else
      B.IBS = B_0.IBS[d_0.IBS:ncol(B_0.IBS),d_0.IBS:ncol(B_0.IBS)]
      if(rcond(B.IBS) >= 1e-10){
        #so the B_0 needs to be subset to m rows and columns that are restricted based on H0
        TW.IBS <- crossprod(beta.unre.IBS[N], solve(B.IBS, beta.unre.IBS[N]))
      } else {
        TW.IBS = -2
      }
    } else {
      TW.IBS = -1
    }
    
    A.Incomp = crossprod(Model.Incomp$X[,indice.unre.Incomp|N], as.vector(pi.unre.Incomp*(1-pi.unre.Incomp))*Model.Incomp$X[,indice.unre.Incomp|N])
    if(rcond(A.Incomp) >= 1e-10){
      B_0.Incomp <- solve(A.Incomp)
      #d_0 should determine which rows to subset
      #gives the first value that is restricted
      d_0.Incomp <- sum(indice.re.Incomp == TRUE) + 1
      #in case the inverse doesn't exist, need another if else
      B.Incomp = B_0.Incomp[d_0.Incomp:ncol(B_0.Incomp),d_0.Incomp:ncol(B_0.Incomp)]
      if(rcond(B.Incomp) >= 1e-10){
        #so the B_0 needs to be subset to m rows and columns that are restricted based on H0
        TW.Incomp <- crossprod(beta.unre.Incomp[N], solve(B.Incomp, beta.unre.Incomp[N]))
      } else {
        TW.Incomp = -2
      }
    } else {
      TW.Incomp = -1
    }
    
    A.AMS = crossprod(Model.AMS$X[,indice.unre.AMS|N], as.vector(pi.unre.AMS*(1-pi.unre.AMS))*Model.AMS$X[,indice.unre.AMS|N])
    if(rcond(A.AMS) >= 1e-10){
      B_0.AMS <- solve(A.AMS)
      #d_0 should determine which rows to subset
      #gives the first value that is restricted
      d_0.AMS <- sum(indice.re.AMS == TRUE) + 1
      #in case the inverse doesn't exist, need another if else
      B.AMS = B_0.AMS[d_0.AMS:ncol(B_0.AMS),d_0.AMS:ncol(B_0.AMS)]
      if(rcond(B.AMS) >= 1e-10){
        #so the B_0 needs to be subset to m rows and columns that are restricted based on H0
        TW.AMS <- crossprod(beta.unre.AMS[N], solve(B.AMS, beta.unre.AMS[N]))
      } else {
        TW.AMS = -2
      }
    } else {
      TW.AMS = -1
    }
    
    A.BinMM = crossprod(Model.BinMM$X[,indice.unre.BinMM|N], as.vector(pi.unre.BinMM*(1-pi.unre.BinMM))*Model.BinMM$X[,indice.unre.BinMM|N])
    if(rcond(A.BinMM) >= 1e-10){
      B_0.BinMM <- solve(A.BinMM)
      #d_0 should determine which rows to subset
      #gives the first value that is restricted
      d_0.BinMM <- sum(indice.re.BinMM == TRUE) + 1
      #in case the inverse doesn't exist, need another if else
      B.BinMM = B_0.BinMM[d_0.BinMM:ncol(B_0.BinMM),d_0.BinMM:ncol(B_0.BinMM)]
      if(rcond(B.BinMM) >= 1e-10){
        #so the B_0 needs to be subset to m rows and columns that are restricted based on H0
        TW.BinMM <- crossprod(beta.unre.BinMM[N], solve(B.BinMM, beta.unre.BinMM[N]))
      } else {
        TW.BinMM = -2
      }
    } else {
      TW.BinMM = -1
    }
    
    # construct the score statistics
    eps.IBS <- Model.IBS$Y-pi.re.IBS #Y - e(Y)
    #this is the X^T times Y-E(Y)
    Xeps.IBS <- crossprod(Model.IBS$X[,indice.re.IBS|N], eps.IBS)
    
    C.IBS = crossprod(Model.IBS$X[,indice.re.IBS|N], as.vector(pi.re.IBS*(1-pi.re.IBS))*Model.IBS$X[,indice.re.IBS|N])
    if(rcond(C.IBS) >= 1e-10){
      TS.IBS <- crossprod(Xeps.IBS, solve(C.IBS, Xeps.IBS))
    } else {
      TS.IBS = -1
    }
    
    eps.Incomp <- Model.Incomp$Y-pi.re.Incomp #Y - e(Y)
    #this is the X^T times Y-E(Y)
    Xeps.Incomp <- crossprod(Model.Incomp$X[,indice.re.Incomp|N], eps.Incomp)
    
    C.Incomp = crossprod(Model.Incomp$X[,indice.re.Incomp|N], as.vector(pi.re.Incomp*(1-pi.re.Incomp))*Model.Incomp$X[,indice.re.Incomp|N])
    if(rcond(C.Incomp) >= 1e-10){
      TS.Incomp <- crossprod(Xeps.Incomp, solve(C.Incomp, Xeps.Incomp))
    } else {
      TS.Incomp = -1
    }
    
    eps.AMS <- Model.AMS$Y-pi.re.AMS #Y - e(Y)
    #this is the X^T times Y-E(Y)
    Xeps.AMS <- crossprod(Model.AMS$X[,indice.re.AMS|N], eps.AMS)
    
    C.AMS = crossprod(Model.AMS$X[,indice.re.AMS|N], as.vector(pi.re.AMS*(1-pi.re.AMS))*Model.AMS$X[,indice.re.AMS|N])
    if(rcond(C.AMS) >= 1e-10){
      TS.AMS <- crossprod(Xeps.AMS, solve(C.AMS, Xeps.AMS))
    } else {
      TS.AMS = -1
    }
    
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
    
    #set pv to 1 if reject, 0 otherwise
    #LRT
    if (TL.IBS>=qchisq(0.95, df = doF)){
      pv.IBS[1] <- pv.IBS[1]+1/5000
    }
    if (TL.Incomp>=qchisq(0.95, df = doF)){
      pv.Incomp[1] <- pv.Incomp[1]+1/5000
    }
    if (TL.AMS>=qchisq(0.95, df = doF)){
      pv.AMS[1] <- pv.AMS[1]+1/5000
    }
    if (TL.BinMM>=qchisq(0.95, df = doF)){
      pv.BinMM[1] <- pv.BinMM[1]+1/5000
    }
    #Wald
    if (TW.IBS>=qchisq(0.95, df = doF)){
      pv.IBS[2] <- pv.IBS[2]+1/5000
    }
    if (TW.Incomp>=qchisq(0.95, df = doF)){
      pv.Incomp[2] <- pv.Incomp[2]+1/5000
    }
    if (TW.AMS>=qchisq(0.95, df = doF)){
      pv.AMS[2] <- pv.AMS[2]+1/5000
    }
    if (TW.BinMM>=qchisq(0.95, df = doF)){
      pv.BinMM[2] <- pv.BinMM[2]+1/5000
    }
    #Score
    if (TS.IBS>=qchisq(0.95, df = doF)){
      pv.IBS[3] <- pv.IBS[3]+1/5000
    }
    if (TS.Incomp>=qchisq(0.95, df = doF)){
      pv.Incomp[3] <- pv.Incomp[3]+1/5000
    }
    if (TS.AMS>=qchisq(0.95, df = doF)){
      pv.AMS[3] <- pv.AMS[3]+1/5000
    }
    if (TS.BinMM>=qchisq(0.95, df = doF)){
      pv.BinMM[3] <- pv.BinMM[3]+1/5000
    }
    
    print(paste0("Simulation ",ii," is complete."))
    
    Tall.IBS[1, ] <- c(TL.IBS, TW.IBS, TS.IBS)
    beta.al.IBS[1, ] <- c(sum(beta.re.IBS!=0), sum(beta.unre.IBS!=0))
    
    Tall.Incomp[1, ] <- c(TL.Incomp, TW.Incomp, TS.Incomp)
    beta.al.Incomp[1, ] <- c(sum(beta.re.Incomp!=0), sum(beta.unre.Incomp!=0))
    
    Tall.AMS[1, ] <- c(TL.AMS, TW.AMS, TS.AMS)
    beta.al.AMS[1, ] <- c(sum(beta.re.AMS!=0), sum(beta.unre.AMS!=0))
    
    Tall.BinMM[1, ] <- c(TL.BinMM, TW.BinMM, TS.BinMM)
    beta.al.BinMM[1, ] <- c(sum(beta.re.BinMM!=0), sum(beta.unre.BinMM!=0))
    
    statsAndPVals = list(pv.IBS=pv.IBS, TScores.IBS=Tall.IBS, beta.IBS=beta.al.IBS,
                         pv.Incomp=pv.Incomp, TScores.Incomp=Tall.Incomp, beta.Incomp=beta.al.Incomp,
                         pv.AMS=pv.AMS, TScores.AMS=Tall.AMS, beta.AMS=beta.al.AMS,
                         pv.BinMM=pv.BinMM, TScores.BinMM=Tall.BinMM, beta.BinMM=beta.al.BinMM)
    statsAndPVals
  }, mc.cores=4)
  
  n = length(statsAndPVals)
  
  pVals.IBS = pVals.Incomp = pVals.AMS = pVals.BinMM = matrix(nrow = n, ncol = 3)
  Scores.IBS = Scores.Incomp = Scores.AMS = Scores.BinMM = matrix(nrow = n, ncol = 3)
  Betas.IBS = Betas.Incomp = Betas.AMS = Betas.BinMM = matrix(nrow = n, ncol = 2)
  
  for(ll in 1:n){
    pVals.IBS[ll,] = unlist(statsAndPVals[[ll]]$pv.IBS)
    Scores.IBS[ll,] = unlist(statsAndPVals[[ll]]$TScores.IBS)
    Betas.IBS[ll,] = unlist(statsAndPVals[[ll]]$beta.IBS)
    
    pVals.Incomp[ll,] = unlist(statsAndPVals[[ll]]$pv.Incomp)
    Scores.Incomp[ll,] = unlist(statsAndPVals[[ll]]$TScores.Incomp)
    Betas.Incomp[ll,] = unlist(statsAndPVals[[ll]]$beta.Incomp)
    
    pVals.AMS[ll,] = unlist(statsAndPVals[[ll]]$pv.AMS)
    Scores.AMS[ll,] = unlist(statsAndPVals[[ll]]$TScores.AMS)
    Betas.AMS[ll,] = unlist(statsAndPVals[[ll]]$beta.AMS)
    
    pVals.BinMM[ll,] = unlist(statsAndPVals[[ll]]$pv.BinMM)
    Scores.BinMM[ll,] = unlist(statsAndPVals[[ll]]$TScores.BinMM)
    Betas.BinMM[ll,] = unlist(statsAndPVals[[ll]]$beta.BinMM)
  }
  
  pVals = cbind(pVals.IBS, pVals.Incomp, pVals.AMS, pVals.BinMM)
  Scores = cbind(Scores.IBS, Scores.Incomp, Scores.AMS, Scores.BinMM)
  Betas = cbind(Betas.IBS, Betas.Incomp, Betas.AMS, Betas.BinMM)
  
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
  write.csv(pVals, file = paste0(path,"/Power_CatY_Prev",YPrev*100,"_LinHypTest_Score",snpOrScore,"_",percentageAssoc,"SNPsAssoc_JointTesting_OR",ORSize,"_",ld,"_Sim",start,"to",numSims,"_PValues_ForceCovFit.csv"))
  write.csv(Scores, file = paste0(path,"/Power_CatY_Prev",YPrev*100,"_LinHypTest_Score",snpOrScore,"_",percentageAssoc,"SNPsAssoc_JointTesting_OR",ORSize,"_",ld,"_Sim",start,"to",numSims,"_Stats_ForceCovFit.csv"))
  write.csv(Betas, file = paste0(path,"/Power_CatY_Prev",YPrev*100,"_LinHypTest_Score",snpOrScore,"_",percentageAssoc,"SNPsAssoc_JointTesting_OR",ORSize,"_",ld,"_Sim",start,"to",numSims,"_Betas_ForceCovFit.csv"))
  
}

# continuous outcome - linear hypothesis testing
RunPowerPipelineLinHypTestCont_RSNP = function(chr, gene, numPairs, Gamma, TrueScore, ORSize, standardizeScores = FALSE, weightedScores = FALSE, scoreWeights, start, stop, percentageAssoc, LowLD){
  # Function to determine power of linear hyp testing method when score is associated, joint testing
  # Only use when Y is continuous
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #Gamma = effect size for score, length 1, can be 0
  #TrueScore = IBS.gene, Incomp.gene, AMS.gene, or BinMM.gene
  #ORSize = Small, Medium, or Large for what OR was used for the associated SNP/score
  #standardizeScores = T or F whether the scores should be standardized based on maximum score value
  #weightedScores = T or F whether the scores will be weighted
  #scoreWeights = m x 1 vector of weights, one weight for each SNP
  # start - simulation number to start at
  # stop - simulation number to stop at
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
  
  #always the same
  numSims = stop  
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_II_Sims/HapGen_Files/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/ProjectIISourceFunctions_v2.R")
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/Linear_ADMM0.r")
  
  #determine which SNPs to actually set as assoc. based on gene
  #check gene
  assocSNPs = DetermineAssocRSNPs(gene = gene, LowLD = LowLD, percentageAssoc = percentageAssoc)
  
  myList = lapply(start:numSims, rep, times = 1)
  
  # p-value initialization
  pv.IBS <- rep(0, 6)
  pv.Incomp <- rep(0, 6)
  pv.AMS <- rep(0, 6)
  pv.BinMM <- rep(0, 6)
  
  Tall.IBS <- matrix(0, 1, 3)
  Tall.Incomp <- matrix(0, 1, 3)
  Tall.AMS <- matrix(0, 1, 3)
  Tall.BinMM <- matrix(0, 1, 3)
  
  beta.al.IBS <- matrix(0, 1, 2)
  beta.al.Incomp <- matrix(0, 1, 2)
  beta.al.AMS <- matrix(0, 1, 2)
  beta.al.BinMM <- matrix(0, 1, 2)
  
  statsAndPVals = mclapply(myList, function(ii){
    #define matrix to hold all Stats and Pvalues
    statsAndPVals = matrix(NA, nrow = numSims, ncol = 3)
    
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
    ContPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = FALSE,  Covariates = CovData, RGenoData = RGenos, ScoreData = PhenoScore, Betas = Betas, Gamma = Gamma)
    
    # define location of zero components
    # basically, this defines what the null hypothesis is, I think
    # so for joint null, we need the non-zero components to be for covariates
    # and the beta and gamma components to be zero
    numCov = dim(CovData)[2]
    numSNPs = dim(RGenos)[2]
    N = c(rep(FALSE, numCov+1), rep(TRUE,1+numSNPs))
    
    #need to combine the covariates, r genos, and score matrices together
    #need to include column of 1s for intercept?
    intercept = matrix(1,nrow = numPairs, ncol = 1)
    designMat.IBS = cbind(intercept, CovData, RGenos, IBS.gene)
    designMat.Incomp = cbind(intercept, CovData, RGenos, Incomp.gene)
    designMat.AMS = cbind(intercept, CovData, RGenos, AMS.gene)
    designMat.BinMM = cbind(intercept, CovData, RGenos, BinMM.gene)
    
    #then combine the phenos with the design matrix as a list
    Model.IBS = list(X=designMat.IBS, Y=ContPhenos)
    Model.Incomp = list(X=designMat.Incomp, Y=ContPhenos)
    Model.AMS = list(X=designMat.AMS, Y=ContPhenos)
    Model.BinMM = list(X=designMat.BinMM, Y=ContPhenos)
    
    # estimate the uncontrained estimator
    beta.unre.IBS <- cv.SCAD_ADMM_unre(X=Model.IBS$X, Y=Model.IBS$Y, N=N, beta0=rep(0, dim(Model.IBS$X)[2]), err=1e-4, tune="cv", unpen = c(1,2,3))
    indice.unre.IBS <- beta.unre.IBS!=0
    
    beta.unre.Incomp <- cv.SCAD_ADMM_unre(X=Model.Incomp$X, Y=Model.Incomp$Y, N=N, beta0=rep(0, dim(Model.Incomp$X)[2]), err=1e-4, tune="cv", unpen = c(1,2,3))
    indice.unre.Incomp <- beta.unre.Incomp!=0
    
    beta.unre.AMS <- cv.SCAD_ADMM_unre(X=Model.AMS$X, Y=Model.AMS$Y, N=N, beta0=rep(0, dim(Model.AMS$X)[2]), err=1e-4, tune="cv", unpen = c(1,2,3))
    indice.unre.AMS <- beta.unre.AMS!=0
    
    beta.unre.BinMM <- cv.SCAD_ADMM_unre(X=Model.BinMM$X, Y=Model.BinMM$Y, N=N, beta0=rep(0, dim(Model.BinMM$X)[2]), err=1e-4, tune="cv", unpen = c(1,2,3))
    indice.unre.BinMM <- beta.unre.BinMM!=0
    
    # estimate the constrained estimator
    beta.re.IBS <- cv.SCAD_ADMM_re(X=Model.IBS$X, Y=Model.IBS$Y, N=N, beta0=rep(0, dim(Model.IBS$X)[2]), err=1e-4, tune="cv", unpen = c(1,2,3))
    indice.re.IBS <- beta.re.IBS!=0
    
    beta.re.Incomp <- cv.SCAD_ADMM_re(X=Model.Incomp$X, Y=Model.Incomp$Y, N=N, beta0=rep(0, dim(Model.Incomp$X)[2]), err=1e-4, tune="cv", unpen = c(1,2,3))
    indice.re.Incomp <- beta.re.Incomp!=0
    
    beta.re.AMS <- cv.SCAD_ADMM_re(X=Model.AMS$X, Y=Model.AMS$Y, N=N, beta0=rep(0, dim(Model.AMS$X)[2]), err=1e-4, tune="cv", unpen = c(1,2,3))
    indice.re.AMS <- beta.re.AMS!=0
    
    beta.re.BinMM <- cv.SCAD_ADMM_re(X=Model.BinMM$X, Y=Model.BinMM$Y, N=N, beta0=rep(0, dim(Model.BinMM$X)[2]), err=1e-4, tune="cv", unpen = c(1,2,3))
    indice.re.BinMM <- beta.re.BinMM!=0
    
    # estimate the conditional variance
    #    sig2 <- mean((Model$Y-Model$X%*%beta.unre)^2)
    n = numSims
    sig2.IBS <- mean((Model.IBS$Y-Model.IBS$X%*%beta.unre.IBS)^2)*n/(n-sum(beta.unre.IBS!=0))
    sig2.Incomp <- mean((Model.Incomp$Y-Model.Incomp$X%*%beta.unre.Incomp)^2)*n/(n-sum(beta.unre.Incomp!=0))
    sig2.AMS <- mean((Model.AMS$Y-Model.AMS$X%*%beta.unre.AMS)^2)*n/(n-sum(beta.unre.AMS!=0))
    sig2.BinMM <- mean((Model.BinMM$Y-Model.BinMM$X%*%beta.unre.BinMM)^2)*n/(n-sum(beta.unre.BinMM!=0))
    
    # construct the likelihood ratio statistic
    TL.IBS <- sum((Model.IBS$X%*%beta.re.IBS-Model.IBS$Y)^2)-sum((Model.IBS$X%*%beta.unre.IBS-Model.IBS$Y)^2)
    #if LRT stat is less than 0, there is error, so set to -1
    if(TL.IBS < 0){
      TL.IBS = -1
    }
    TL.Incomp <- sum((Model.Incomp$X%*%beta.re.Incomp-Model.Incomp$Y)^2)-sum((Model.Incomp$X%*%beta.unre.Incomp-Model.Incomp$Y)^2)
    #if LRT stat is less than 0, there is error, so set to -1
    if(TL.Incomp < 0){
      TL.Incomp = -1
    }
    TL.AMS <- sum((Model.AMS$X%*%beta.re.AMS-Model.AMS$Y)^2)-sum((Model.AMS$X%*%beta.unre.AMS-Model.AMS$Y)^2)
    #if LRT stat is less than 0, there is error, so set to -1
    if(TL.AMS < 0){
      TL.AMS = -1
    }
    TL.BinMM <- sum((Model.BinMM$X%*%beta.re.BinMM-Model.BinMM$Y)^2)-sum((Model.BinMM$X%*%beta.unre.BinMM-Model.BinMM$Y)^2)
    #if LRT stat is less than 0, there is error, so set to -1
    if(TL.BinMM < 0){
      TL.BinMM = -1
    }
    
    # construct the Wald statistic
    #B_0 should give you Omega_a hat
    B.IBS = crossprod(Model.IBS$X[,indice.unre.IBS|N], Model.IBS$X[,indice.unre.IBS|N])
    if(rcond(B.IBS) >= 1e-10){
      B_0.IBS <- solve(B.IBS)
      #d_0.IBS should determine which rows to subset
      #gives the first value that is restricted
      d_0.IBS <- sum(indice.re.IBS == TRUE) + 1
      #in case the inverse doesn't exist, need another if else
      A.IBS = B_0.IBS[d_0.IBS:ncol(B_0.IBS),d_0.IBS:ncol(B_0.IBS)]
      if(rcond(A.IBS) >= 1e-10){
        TW.IBS <- crossprod(beta.unre.IBS[N], solve(A.IBS, beta.unre.IBS[N]))
      } else {
        TW.IBS = -2
      }
    } else {
      TW.IBS = -1
    }
    
    B.Incomp = crossprod(Model.Incomp$X[,indice.unre.Incomp|N], Model.Incomp$X[,indice.unre.Incomp|N])
    if(rcond(B.Incomp) >= 1e-10){
      B_0.Incomp <- solve(B.Incomp)
      #d_0.Incomp should determine which rows to subset
      #gives the first value that is restricted
      d_0.Incomp <- sum(indice.re.Incomp == TRUE) + 1
      #in case the inverse doesn't exist, need another if else
      A.Incomp = B_0.Incomp[d_0.Incomp:ncol(B_0.Incomp),d_0.Incomp:ncol(B_0.Incomp)]
      if(rcond(A.Incomp) >= 1e-10){
        TW.Incomp <- crossprod(beta.unre.Incomp[N], solve(A.Incomp, beta.unre.Incomp[N]))
      } else {
        TW.Incomp = -2
      }
    } else {
      TW.Incomp = -1
    }
    
    B.AMS = crossprod(Model.AMS$X[,indice.unre.AMS|N], Model.AMS$X[,indice.unre.AMS|N])
    if(rcond(B.AMS) >= 1e-10){
      B_0.AMS <- solve(B.AMS)
      #d_0.AMS should determine which rows to subset
      #gives the first value that is restricted
      d_0.AMS <- sum(indice.re.AMS == TRUE) + 1
      #in case the inverse doesn't exist, need another if else
      A.AMS = B_0.AMS[d_0.AMS:ncol(B_0.AMS),d_0.AMS:ncol(B_0.AMS)]
      if(rcond(A.AMS) >= 1e-10){
        TW.AMS <- crossprod(beta.unre.AMS[N], solve(A.AMS, beta.unre.AMS[N]))
      } else {
        TW.AMS = -2
      }
    } else {
      TW.AMS = -1
    }
    
    B.BinMM = crossprod(Model.BinMM$X[,indice.unre.BinMM|N], Model.BinMM$X[,indice.unre.BinMM|N])
    if(rcond(B.BinMM) >= 1e-10){
      B_0.BinMM <- solve(B.BinMM)
      #d_0.BinMM should determine which rows to subset
      #gives the first value that is restricted
      d_0.BinMM <- sum(indice.re.BinMM == TRUE) + 1
      #in case the inverse doesn't exist, need another if else
      A.BinMM = B_0.BinMM[d_0.BinMM:ncol(B_0.BinMM),d_0.BinMM:ncol(B_0.BinMM)]
      if(rcond(A.BinMM) >= 1e-10){
        TW.BinMM <- crossprod(beta.unre.BinMM[N], solve(A.BinMM, beta.unre.BinMM[N]))
      } else {
        TW.BinMM = -2
      }
    } else {
      TW.BinMM = -1
    }
    
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
    
    #LRT
    if (TL.IBS>=qchisq(0.95, df = doF)){
      pv.IBS[1] <- pv.IBS[1]+1/5000
    }
    if (TL.Incomp>=qchisq(0.95, df = doF)){
      pv.Incomp[1] <- pv.Incomp[1]+1/5000
    }
    if (TL.AMS>=qchisq(0.95, df = doF)){
      pv.AMS[1] <- pv.AMS[1]+1/5000
    }
    if (TL.BinMM>=qchisq(0.95, df = doF)){
      pv.BinMM[1] <- pv.BinMM[1]+1/5000
    }
    #Wald
    if (TW.IBS>=qchisq(0.95, df = doF)){
      pv.IBS[2] <- pv.IBS[2]+1/5000
    }
    if (TW.Incomp>=qchisq(0.95, df = doF)){
      pv.Incomp[2] <- pv.Incomp[2]+1/5000
    }
    if (TW.AMS>=qchisq(0.95, df = doF)){
      pv.AMS[2] <- pv.AMS[2]+1/5000
    }
    if (TW.BinMM>=qchisq(0.95, df = doF)){
      pv.BinMM[2] <- pv.BinMM[2]+1/5000
    }
    #Score
    if (TS.IBS>=qchisq(0.95, df = doF)){
      pv.IBS[3] <- pv.IBS[3]+1/5000
    }
    if (TS.Incomp>=qchisq(0.95, df = doF)){
      pv.Incomp[3] <- pv.Incomp[3]+1/5000
    }
    if (TS.AMS>=qchisq(0.95, df = doF)){
      pv.AMS[3] <- pv.AMS[3]+1/5000
    }
    if (TS.BinMM>=qchisq(0.95, df = doF)){
      pv.BinMM[3] <- pv.BinMM[3]+1/5000
    }
    
    print(paste0("Simulation ",ii," is complete."))
    
    Tall.IBS[1, ] <- c(TL.IBS, TW.IBS, TS.IBS)
    beta.al.IBS[1, ] <- c(sum(beta.re.IBS!=0), sum(beta.unre.IBS!=0))
    
    Tall.Incomp[1, ] <- c(TL.Incomp, TW.Incomp, TS.Incomp)
    beta.al.Incomp[1, ] <- c(sum(beta.re.Incomp!=0), sum(beta.unre.Incomp!=0))
    
    Tall.AMS[1, ] <- c(TL.AMS, TW.AMS, TS.AMS)
    beta.al.AMS[1, ] <- c(sum(beta.re.AMS!=0), sum(beta.unre.AMS!=0))
    
    Tall.BinMM[1, ] <- c(TL.BinMM, TW.BinMM, TS.BinMM)
    beta.al.BinMM[1, ] <- c(sum(beta.re.BinMM!=0), sum(beta.unre.BinMM!=0))
    
    statsAndPVals = list(pv.IBS=pv.IBS, TScores.IBS=Tall.IBS, beta.IBS=beta.al.IBS,
                         pv.Incomp=pv.Incomp, TScores.Incomp=Tall.Incomp, beta.Incomp=beta.al.Incomp,
                         pv.AMS=pv.AMS, TScores.AMS=Tall.AMS, beta.AMS=beta.al.AMS,
                         pv.BinMM=pv.BinMM, TScores.BinMM=Tall.BinMM, beta.BinMM=beta.al.BinMM)
    statsAndPVals
  }, mc.cores=4)
  
  n = length(statsAndPVals)
  
  pVals.IBS = pVals.Incomp = pVals.AMS = pVals.BinMM = matrix(nrow = n, ncol = 6)
  Scores.IBS = Scores.Incomp = Scores.AMS = Scores.BinMM = matrix(nrow = n, ncol = 3)
  Betas.IBS = Betas.Incomp = Betas.AMS = Betas.BinMM = matrix(nrow = n, ncol = 2)
  
  for(ll in 1:n){
    pVals.IBS[ll,] = unlist(statsAndPVals[[ll]]$pv.IBS)
    Scores.IBS[ll,] = unlist(statsAndPVals[[ll]]$TScores.IBS)
    Betas.IBS[ll,] = unlist(statsAndPVals[[ll]]$beta.IBS)
    
    pVals.Incomp[ll,] = unlist(statsAndPVals[[ll]]$pv.Incomp)
    Scores.Incomp[ll,] = unlist(statsAndPVals[[ll]]$TScores.Incomp)
    Betas.Incomp[ll,] = unlist(statsAndPVals[[ll]]$beta.Incomp)
    
    pVals.AMS[ll,] = unlist(statsAndPVals[[ll]]$pv.AMS)
    Scores.AMS[ll,] = unlist(statsAndPVals[[ll]]$TScores.AMS)
    Betas.AMS[ll,] = unlist(statsAndPVals[[ll]]$beta.AMS)
    
    pVals.BinMM[ll,] = unlist(statsAndPVals[[ll]]$pv.BinMM)
    Scores.BinMM[ll,] = unlist(statsAndPVals[[ll]]$TScores.BinMM)
    Betas.BinMM[ll,] = unlist(statsAndPVals[[ll]]$beta.BinMM)
  }
  
  pVals = cbind(pVals.IBS, pVals.Incomp, pVals.AMS, pVals.BinMM)
  Scores = cbind(Scores.IBS, Scores.Incomp, Scores.AMS, Scores.BinMM)
  Betas = cbind(Betas.IBS, Betas.Incomp, Betas.AMS, Betas.BinMM)
  
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
  write.csv(pVals, file = paste0(path,"/Power_ContY_LinHypTest_Score",snpOrScore,"_",percentageAssoc,"SNPsAssoc_JointTesting_OR",ORSize,"_",ld,"_Sim",start,"to",numSims,"_PValues_ForceCovFit.csv"))
  write.csv(Scores, file = paste0(path,"/Power_ContY_LinHypTest_Score",snpOrScore,"_",percentageAssoc,"SNPsAssoc_JointTesting_OR",ORSize,"_",ld,"_Sim",start,"to",numSims,"_Stats_ForceCovFit.csv"))
  write.csv(Betas, file = paste0(path,"/Power_ContY_LinHypTest_Score",snpOrScore,"_",percentageAssoc,"SNPsAssoc_JointTesting_OR",ORSize,"_",ld,"_Sim",start,"to",numSims,"_Betas_ForceCovFit.csv"))
}

##########################################################################################################################################
#################           Gamma = 0 Testing Only         ###############################################################################
#################### TIE Pipeline with AnnalsOfStats Paper
#################### Case 1: Testing Gamma = 0 when Beta = 0

#categorical outcome
RunTIEPipelineLinHypTestCat_GammaTest_GLM = function(chr, gene, numPairs, YPrev, standardizeScores = FALSE, weightedScores = FALSE, scoreWeights, score, start, stop){
  #function to run whole TIE pipeline, Calculates LRT test stat for Lin Hyp test method and whether p-value <0.05
  # Only use when Y is binary
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #standardizeScores = T or F whether the scores should be standardized based on maximum score value
  #weightedScores = T or F whether the scores will be weighted
  #scoreWeights = m x 1 vector of weights, one weight for each SNP
  # score - which score we are using in the main JST model
  # start - simulation number to start at
  # stop - simulation number to stop at
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes TIE values to csv files

  library(parallel)
  suppressMessages(library(epicalc))

  #always the same
  numSims = stop

  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_II_Sims/HapGen_Files/",gene,"_Results_",numPairs,"Pairs")

  #source the needed functions
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/ProjectIISourceFunctions_v2.R")

  myList = lapply(start:numSims, rep, times = 1)
  nullValues =  altValues = statsAndPValsmat = finalOutput = list()

  statsAndPVals = mclapply(myList, function(ii){
    #pull recipient and donor genotypes
    RGenos = obtainRGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    DGenos = obtainDGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)

    #calculate single snp scores
    if(score == "IBS"){
      Score.snp = calcIBSMismatch(RGenosMat = RGenos, DGenosMat = DGenos)
    } else if(score == "Incomp"){
      Score.snp = calcIncompatibilityScore(RGenosMat = RGenos, DGenosMat = DGenos)
    } else if(score == "AMS"){
      Score.snp = calcAMS(RGenosMat = RGenos, DGenosMat = DGenos)
    } else{
      Score.snp = calcBinaryMM(RGenosMat = RGenos, DGenosMat = DGenos)
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

    #generate covariates
    #for now, a single binary and a single continous covariate
    # CovData = GenCovData(SampleSize = numPairs, BinaryValues = 1, ContinuousValues = 1)

    #generate phenotypes, both continuous and binary
    # CatPhenos = GenNullPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = TRUE, YPrev = YPrev,  Covariates = CovData)
    CatPhenos = GenNullPhenos(SampleSize = numPairs, includeCov = FALSE, YCat = TRUE, YPrev = YPrev)

    # allData = cbind(CovData, RGenos, Score.gene, CatPhenos)
    allData = cbind(RGenos, Score.gene, CatPhenos)

    #fit the null and alternative models
    fitNull = glm(CatPhenos~RGenos, family = binomial)
    fitAlt = glm(CatPhenos~RGenos+Score.gene, family = binomial, epsilon = 1e-6)

    nullValues = summary(fitNull)$coefficients
    altValues = summary(fitAlt)$coefficients

    if(fitNull$deviance - fitAlt$deviance >= 0){
      # construct the likelihood ratio statistics
      TL = lrtest(fitNull, fitAlt)

      pv = TL$p.value
      stat = TL$Chisquared
      doF = TL$df
    } else {
      pv = 1
      stat = 0
      doF = (dim(RGenos)[2]+1)
    }

    print(paste0("Simulation ",ii," is complete."))

    statsAndPValsmat = c(pv, stat)

    finalOutput[[ii]] = list(statsAndPValsmat, nullValues, altValues)
  }, mc.cores=4)

  statsAndPValsAll = list()
  nulValuesAll = list()
  altValuesAll = list()

  for(ii in 1:length(statsAndPVals)){
    statsAndPValsAll[[ii]] = statsAndPVals[[ii]][1]
    nulValuesAll[[ii]] =  statsAndPVals[[ii]][2]
    altValuesAll[[ii]] = statsAndPVals[[ii]][3]
  }

  listLength = length(statsAndPValsAll)
  statsAndPValsAll.mat = matrix(unlist(statsAndPValsAll),nrow = listLength, ncol = 2, byrow = TRUE)
  dim1 = dim(nulValuesAll[[1]][[1]])
  dim2 = dim(altValuesAll[[1]][[1]])
  rownames1 = c()
  rownames2 = c()
  for(jj in 1:length(statsAndPVals)){
    rownames1 = c(rownames1, rownames(nulValuesAll[[jj]][[1]]))
    rownames2 = c(rownames2, rownames(altValuesAll[[jj]][[1]]))
  }
  nulValuesAll.mat = matrix(unlist(nulValuesAll), ncol = dim1[2], byrow = TRUE, dimnames = list(rownames1, colnames(nulValuesAll[[1]][[1]])))
  altValuesAll.mat = matrix(unlist(altValuesAll), ncol = dim2[2], byrow = TRUE, dimnames = list(rownames2, colnames(altValuesAll[[1]][[1]])))

  #rename
  statsAndPValsCatPhenos = statsAndPValsAll.mat

  #write out the Stats and p values, and Summary stats
  if(weightedScores == FALSE){
    if(standardizeScores == FALSE){
      write.csv(statsAndPValsCatPhenos, file = paste0(path,"/TIE_CatPhenos_Prev",YPrev*100,"_GLM_Score",score,"_Sim",start,"to",numSims,"_StatsAndPValues.csv"))
      write.csv(nulValuesAll.mat, file = paste0(path,"/TIE_CatNullBetasSummary_Prev",YPrev*100,"_GLM_Score",score,"_Sim",start,"to",numSims,"_StatsAndPValues.csv"))
      write.csv(altValuesAll.mat, file = paste0(path,"/TIE_CatAltBetasSummary_Prev",YPrev*100,"_GLM_Score",score,"_Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    } else { #scores are standardized
      write.csv(statsAndPValsCatPhenos, file = paste0(path,"/TIE_CatPhenos_Prev",YPrev*100,"_GLM_Score",score,"_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
      write.csv(nulValuesAll.mat, file = paste0(path,"/TIE_CatBetasSummary_Prev",YPrev*100,"_GLM_Score",score,"_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
      write.csv(altValuesAll.mat, file = paste0(path,"/TIE_CatAltBetasSummary_Prev",YPrev*100,"_GLM_Score",score,"_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
    }
  } else { #scores are weighted
    write.csv(statsAndPValsCatPhenos, file = paste0(path,"/TIE_CatPhenos_Prev",YPrev*100,"_GLM_Score",score,"_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
    write.csv(nulValuesAll.mat, file = paste0(path,"/TIE_CatBetasSummary_Prev",YPrev*100,"_GLM_Score",score,"_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
    write.csv(altValuesAll.mat, file = paste0(path,"/TIE_CatAltBetasSummary_Prev",YPrev*100,"_GLM_Score",score,"_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
  }
}

#continuous outcome
RunTIEPipelineLinHypTestCont_GammaTest_GLM = function(chr, gene, numPairs, standardizeScores = FALSE, weightedScores = FALSE, scoreWeights, score, start, stop){
  #function to run whole TIE pipeline, Calculates LRT test stat for Lin Hyp test method and whether p-value <0.05
  # only use when Y is continuous (Normal)
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #standardizeScores = T or F whether the scores should be standardized based on maximum score value
  #weightedScores = T or F whether the scores will be weighted
  #scoreWeights = m x 1 vector of weights, one weight for each SNP
  # score - which score we are using in the main JST model
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes TIE values to csv files
  
  library(parallel)
  suppressMessages(library(epicalc))
  
  #always the same
  numSims = stop  
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_II_Sims/HapGen_Files/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/ProjectIISourceFunctions_v2.R")
  
  myList = lapply(start:numSims, rep, times = 1)
  nullValues =  altValues = statsAndPValsmat = finalOutput = list()
  
  statsAndPVals = mclapply(myList, function(ii){
    #pull recipient and donor genotypes
    RGenos = obtainRGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    DGenos = obtainDGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    
    #calculate single snp scores
    if(score == "IBS"){
      Score.snp = calcIBSMismatch(RGenosMat = RGenos, DGenosMat = DGenos)
    } else if(score == "Incomp"){
      Score.snp = calcIncompatibilityScore(RGenosMat = RGenos, DGenosMat = DGenos)
    } else if(score == "AMS"){
      Score.snp = calcAMS(RGenosMat = RGenos, DGenosMat = DGenos)
    } else{
      Score.snp = calcBinaryMM(RGenosMat = RGenos, DGenosMat = DGenos)
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
    
    #generate covariates
    #for now, a single binary and a single continous covariate
    # CovData = GenCovData(SampleSize = numPairs, BinaryValues = 1, ContinuousValues = 1)
    
    #generate phenotypes, both continuous and binary
    # ContPhenos = GenNullPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = FALSE, Covariates = CovData)
    ContPhenos = GenNullPhenos(SampleSize = numPairs, includeCov = FALSE, YCat = FALSE)
    
    #fit the null and alternative models
    fitNull = glm(ContPhenos~RGenos, family = gaussian)
    fitAlt = glm(ContPhenos~RGenos+Score.gene, family = gaussian, epsilon = 1e-6)	
    
    nullValues = summary(fitNull)$coefficients
    altValues = summary(fitAlt)$coefficients
    
    if(fitNull$deviance - fitAlt$deviance >= 0){
      # construct the likelihood ratio statistics
      TL = lrtest(fitNull, fitAlt)
      
      pv = TL$p.value
      stat = TL$Chisquared
      doF = TL$df
    } else {
      pv = 1
      stat = 0
      doF = (dim(RGenos)[2]+1)
    }
    
    print(paste0("Simulation ",ii," is complete."))
    
    statsAndPValsmat = c(pv, stat)
    
    finalOutput[[ii]] = list(statsAndPValsmat, nullValues, altValues)
  }, mc.cores=1)
  
  
  statsAndPValsAll = list()
  nulValuesAll = list()
  altValuesAll = list()
  
  for(ii in 1:length(statsAndPVals)){
    statsAndPValsAll[[ii]] = statsAndPVals[[ii]][1]
    nulValuesAll[[ii]] =  statsAndPVals[[ii]][2]
    altValuesAll[[ii]] = statsAndPVals[[ii]][3]
  }
  
  listLength = length(statsAndPValsAll)  
  statsAndPValsAll.mat = matrix(unlist(statsAndPValsAll),nrow = listLength, ncol = 2, byrow = TRUE)
  dim1 = dim(nulValuesAll[[1]][[1]])
  dim2 = dim(altValuesAll[[1]][[1]])
  rownames1 = c()
  rownames2 = c()
  for(jj in 1:length(statsAndPVals)){
    rownames1 = c(rownames1, rownames(nulValuesAll[[jj]][[1]]))
    rownames2 = c(rownames2, rownames(altValuesAll[[jj]][[1]]))
  }
  nulValuesAll.mat = matrix(unlist(nulValuesAll), ncol = dim1[2], byrow = TRUE, dimnames = list(rownames1, colnames(nulValuesAll[[1]][[1]])))
  altValuesAll.mat = matrix(unlist(altValuesAll), ncol = dim2[2], byrow = TRUE, dimnames = list(rownames2, colnames(altValuesAll[[1]][[1]])))
  
  #separate into cat and cont phenos
  statsAndPValsContPhenos = statsAndPValsAll.mat
  
  #write out the Stats and p values
  if(weightedScores == FALSE){
    if(standardizeScores == FALSE){
      write.csv(statsAndPValsContPhenos, file = paste0(path,"/TIE_ContPhenos_LinHypTest_Score",score,"_Sim",start,"to",numSims,"_StatsAndPValues.csv"))
      write.csv(nulValuesAll.mat, file = paste0(path,"/TIE_ContNullBetasSummary_LinHypTest_Score",score,"_Sim",start,"to",numSims,"_StatsAndPValues.csv"))
      write.csv(altValuesAll.mat, file = paste0(path,"/TIE_ContAltBetasSummary_LinHypTest_Score",score,"_Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    } else { #scores are standardized
      write.csv(statsAndPValsContPhenos, file = paste0(path,"/TIE_ContPhenos_LinHypTest_Score",score,"_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
      write.csv(nulValuesAll.mat, file = paste0(path,"/TIE_ContNullBetasSummary_LinHypTest_Score",score,"_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
      write.csv(altValuesAll.mat, file = paste0(path,"/TIE_ContAltBetasSummary_LinHypTest_Score",score,"_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
    }
  } else { #scores are weighted
    write.csv(statsAndPValsContPhenos, file = paste0(path,"/TIE_ContPhenos_LinHypTest_Score",score,"_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
    write.csv(nulValuesAll.mat, file = paste0(path,"/TIE_ContNullBetasSummary_LinHypTest_Score",score,"_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
    write.csv(altValuesAll.mat, file = paste0(path,"/TIE_ContAltBetasSummary_LinHypTest_Score",score,"_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
  }
}

#categorical outcome
RunTIEPipelineLinHypTestCat_GammaTest_Lasso = function(chr, gene, numPairs, YPrev, standardizeScores = FALSE, weightedScores = FALSE, scoreWeights, score, start, stop){
  #function to run whole TIE pipeline, Calculates LRT test stat for Lin Hyp test method and whether p-value <0.05
  # Only use when Y is binary
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #standardizeScores = T or F whether the scores should be standardized based on maximum score value
  #weightedScores = T or F whether the scores will be weighted
  #scoreWeights = m x 1 vector of weights, one weight for each SNP
  # score - which score we are using in the main JST model
  # start - simulation number to start at
  # stop - simulation number to stop at
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes TIE values to csv files
  
  library(parallel)
  suppressMessages(library(epicalc))
  suppressMessages(library(glmnet))
  
  #always the same
  numSims = stop  
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_II_Sims/HapGen_Files/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/ProjectIISourceFunctions_v2.R")
  
  myList = lapply(start:numSims, rep, times = 1)
  nullValues =  altValues = statsAndPValsmat = finalOutput = list()
  
  statsAndPVals = mclapply(myList, function(ii){
    #pull recipient and donor genotypes
    RGenos = obtainRGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    DGenos = obtainDGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    
    #calculate single snp scores
    if(score == "IBS"){
      Score.snp = calcIBSMismatch(RGenosMat = RGenos, DGenosMat = DGenos)
    } else if(score == "Incomp"){
      Score.snp = calcIncompatibilityScore(RGenosMat = RGenos, DGenosMat = DGenos)
    } else if(score == "AMS"){
      Score.snp = calcAMS(RGenosMat = RGenos, DGenosMat = DGenos)
    } else{
      Score.snp = calcBinaryMM(RGenosMat = RGenos, DGenosMat = DGenos)
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
    
    #generate covariates
    #for now, a single binary and a single continous covariate
    # CovData = GenCovData(SampleSize = numPairs, BinaryValues = 1, ContinuousValues = 1)
    
    #generate phenotypes, both continuous and binary
    # CatPhenos = GenNullPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = TRUE, YPrev = YPrev,  Covariates = CovData)
    CatPhenos = GenNullPhenos(SampleSize = numPairs, includeCov = FALSE, YCat = TRUE, YPrev = YPrev)
    
    #combine the data for x design matrix
    allData = cbind(RGenos, Score.gene)
    
    #split into training and test data
    #train with 80% of data
    train.null.index = sample(1:nrow(RGenos), 5*(nrow(RGenos)/10))
    train.all.index = sample(1:nrow(allData), 5*(nrow(allData)/10))
    #pull the index values for the test data (20%)
    test.null.index = -train.null.index
    test.all.index = -train.all.index
    #separate the data
    x.train.null = RGenos[train.null.index,]
    x.train.all = allData[train.all.index,]
    x.test.null = RGenos[test.null.index,]
    x.test.all = allData[test.all.index,]
    
    #convert y phenos outcome to numerical variable
    CatPhenos.fac = ifelse(CatPhenos == 1,1,0)
    
    #pull the y phenos for the 
    y.train.null = CatPhenos.fac[train.null.index]
    y.train.all = CatPhenos.fac[train.all.index]
    y.test.null = CatPhenos.fac[test.null.index]
    y.test.all = CatPhenos.fac[test.all.index]
    
    #use cv to find best lambda value
    cv_output_null = cv.glmnet(x.train.null, y.train.null, family = "binomial", alpha=1, intercept=FALSE)
    cv_output_alt = cv.glmnet(x.train.all, y.train.all, family = "binomial", alpha=1, intercept=FALSE)
    cv_output_null.int = cv.glmnet(x.train.null, y.train.null, family = "binomial", alpha=1, intercept=TRUE)
    cv_output_alt.int = cv.glmnet(x.train.all, y.train.all, family = "binomial", alpha=1, intercept=TRUE)
    
    #identifying the best lambda
    best_lambda_null = cv_output_null$lambda.min
    best_lambda_alt = cv_output_alt$lambda.min
    best_lambda_null.int = cv_output_null.int$lambda.min
    best_lambda_alt.int = cv_output_alt.int$lambda.min
    
    #fit the null and alternative models with the best lambda values
    fitNull = glmnet(RGenos, CatPhenos.fac, family = "binomial", alpha=1, lambda = best_lambda_null, intercept=FALSE)
    fitAlt = glmnet(allData, CatPhenos.fac, family = "binomial", alpha=1, lambda = best_lambda_alt, intercept=FALSE)	
    fitNull.int = glmnet(RGenos, CatPhenos.fac, family = "binomial", alpha=1, lambda = best_lambda_null.int, intercept=TRUE)
    fitAlt.int = glmnet(allData, CatPhenos.fac, family = "binomial", alpha=1, lambda = best_lambda_alt.int, intercept=TRUE)	
    
    nullValues = coef(fitNull)
    altValues = coef(fitAlt)
    nullValues.int = coef(fitNull.int)
    altValues.int = coef(fitAlt.int)
    
    nullValues.lambda = fitNull$lambda
    altValues.lambda = fitAlt$lambda
    nullValues.int.lambda = fitNull.int$lambda
    altValues.int.lambda = fitAlt.int$lambda
    
    lambdas = c(nullValues.lambda, altValues.lambda, nullValues.int.lambda, altValues.int.lambda)
    
    print(paste0("Simulation ",ii," is complete."))
    
    finalOutput[[ii]] = list(nullValues, altValues, nullValues.int, altValues.int, lambdas)
  }, mc.cores=4)
  
  betasNull = list()
  betasAll = list()
  betasNull.int = list()
  betasAll.int = list()
  lambdas = list()
  
  for(ii in 1:length(statsAndPVals)){
    betasNull[[ii]] = statsAndPVals[[ii]][1]
    betasAll[[ii]] =  statsAndPVals[[ii]][2]
    betasNull.int[[ii]] = statsAndPVals[[ii]][3]
    betasAll.int[[ii]] = statsAndPVals[[ii]][4]
    lambdas[[ii]] = statsAndPVals[[ii]][5]
  }
  
  betasNull.mat =  do.call(rbind, unlist(betasNull, recursive = FALSE))
  betasAll.mat = do.call(rbind, unlist(betasAll, recursive = FALSE))
  betasNull.int.mat = do.call(rbind, unlist(betasNull.int, recursive = FALSE))
  betasAll.int.mat = do.call(rbind, unlist(betasAll.int, recursive = FALSE))
  lambdas.mat = matrix(unlist(lambdas), ncol = 4, byrow = TRUE)
  
  betasNullmat = as.matrix(betasNull.mat)
  betasAllmat = as.matrix(betasAll.mat)
  betasNull.intmat = as.matrix(betasNull.int.mat)
  betasAll.intmat = as.matrix(betasAll.int.mat)
  
  #write out the Stats and p values, and Summary stats
  if(weightedScores == FALSE){
    if(standardizeScores == FALSE){
      write.csv(betasNullmat, file = paste0(path,"/TIE_CatNullBetas_Prev",YPrev*100,"_GLM_Score",score,"_Sim",start,"to",numSims,"_StatsAndPValues.csv"))
      write.csv(betasAllmat, file = paste0(path,"/TIE_CatAllBetas_Prev",YPrev*100,"_GLM_Score",score,"_Sim",start,"to",numSims,"_StatsAndPValues.csv"))
      write.csv(betasNull.intmat, file = paste0(path,"/TIE_CatNullBetasWIntercept_Prev",YPrev*100,"_GLM_Score",score,"_Sim",start,"to",numSims,"_StatsAndPValues.csv"))
      write.csv(betasAll.intmat, file = paste0(path,"/TIE_CatAllBetasWIntercept_Prev",YPrev*100,"_GLM_Score",score,"_Sim",start,"to",numSims,"_StatsAndPValues.csv"))
      write.csv(lambdas.mat, file = paste0(path,"/TIE_CatLambdas_Prev",YPrev*100,"_GLM_Score",score,"_Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    } else { #scores are standardized
      write.csv(betasNullmat, file = paste0(path,"/TIE_CatNullBetas_Prev",YPrev*100,"_GLM_Score",score,"_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
      write.csv(betasAllmat, file = paste0(path,"/TIE_CatAllBetas_Prev",YPrev*100,"_GLM_Score",score,"_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
      write.csv(betasNull.intmat, file = paste0(path,"/TIE_CatNullBetasWIntercept_Prev",YPrev*100,"_GLM_Score",score,"_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
      write.csv(betasAll.intmat, file = paste0(path,"/TIE_CatAllBetasWIntercept_Prev",YPrev*100,"_GLM_Score",score,"_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
      write.csv(lambdas.mat, file = paste0(path,"/TIE_CatLambdas_Prev",YPrev*100,"_GLM_Score",score,"_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
    }
  } else { #scores are weighted
    write.csv(betasNullmat, file = paste0(path,"/TIE_CatNullBetas_Prev",YPrev*100,"_GLM_Score",score,"_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
    write.csv(betasAllmat, file = paste0(path,"/TIE_CatAllBetas_Prev",YPrev*100,"_GLM_Score",score,"_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
    write.csv(betasNull.intmat, file = paste0(path,"/TIE_CatNullBetasWIntercept_Prev",YPrev*100,"_GLM_Score",score,"_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
    write.csv(betasAll.intmat, file = paste0(path,"/TIE_CatAllBetasWIntercept_Prev",YPrev*100,"_GLM_Score",score,"_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
    write.csv(lambdas.mat, file = paste0(path,"/TIE_CatLambdas_Prev",YPrev*100,"_GLM_Score",score,"_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
  }
}

#continuous outcome
RunTIEPipelineLinHypTestCont_GammaTest_Lasso = function(chr, gene, numPairs, standardizeScores = FALSE, weightedScores = FALSE, scoreWeights, score, start, stop){
  #function to run whole TIE pipeline, Calculates LRT test stat for Lin Hyp test method and whether p-value <0.05
  # only use when Y is continuous (Normal)
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #standardizeScores = T or F whether the scores should be standardized based on maximum score value
  #weightedScores = T or F whether the scores will be weighted
  #scoreWeights = m x 1 vector of weights, one weight for each SNP
  # score - which score we are using in the main JST model
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes TIE values to csv files
  
  library(parallel)
  suppressMessages(library(epicalc))
  suppressMessages(library(glmnet))
  
  #always the same
  numSims = stop  
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_II_Sims/HapGen_Files/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/ProjectIISourceFunctions_v2.R")
  
  myList = lapply(start:numSims, rep, times = 1)
  nullValues =  altValues = statsAndPValsmat = finalOutput = list()
  
  statsAndPVals = mclapply(myList, function(ii){
    #pull recipient and donor genotypes
    RGenos = obtainRGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    DGenos = obtainDGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    
    #calculate single snp scores
    if(score == "IBS"){
      Score.snp = calcIBSMismatch(RGenosMat = RGenos, DGenosMat = DGenos)
    } else if(score == "Incomp"){
      Score.snp = calcIncompatibilityScore(RGenosMat = RGenos, DGenosMat = DGenos)
    } else if(score == "AMS"){
      Score.snp = calcAMS(RGenosMat = RGenos, DGenosMat = DGenos)
    } else{
      Score.snp = calcBinaryMM(RGenosMat = RGenos, DGenosMat = DGenos)
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
    
    #generate covariates
    #for now, a single binary and a single continous covariate
    # CovData = GenCovData(SampleSize = numPairs, BinaryValues = 1, ContinuousValues = 1)
    
    #generate phenotypes, both continuous and binary
    # ContPhenos = GenNullPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = FALSE, Covariates = CovData)
    ContPhenos = GenNullPhenos(SampleSize = numPairs, includeCov = FALSE, YCat = FALSE)
    
    #combine the data for x design matrix
    allData = cbind(RGenos, Score.gene)
    
    #split into training and test data
    #train with 80% of data
    train.null.index = sample(1:nrow(RGenos), 5*(nrow(RGenos)/10))
    train.all.index = sample(1:nrow(allData), 5*(nrow(allData)/10))
    #pull the index values for the test data (20%)
    test.null.index = -train.null.index
    test.all.index = -train.all.index
    #separate the data
    x.train.null = RGenos[train.null.index,]
    x.train.all = allData[train.all.index,]
    x.test.null = RGenos[test.null.index,]
    x.test.all = allData[test.all.index,]
    
    #pull the y phenos for the 
    y.train.null = ContPhenos[train.null.index]
    y.train.all = ContPhenos[train.all.index]
    y.test.null = ContPhenos[test.null.index]
    y.test.all = ContPhenos[test.all.index]
    
    #use cv to find best lambda value
    cv_output_null = cv.glmnet(x.train.null, y.train.null, family = "gaussian", alpha=1, intercept=FALSE)
    cv_output_alt = cv.glmnet(x.train.all, y.train.all, family = "gaussian", alpha=1, intercept=FALSE)
    cv_output_null.int = cv.glmnet(x.train.null, y.train.null, family = "gaussian", alpha=1, intercept=TRUE)
    cv_output_alt.int = cv.glmnet(x.train.all, y.train.all, family = "gaussian", alpha=1, intercept=TRUE)
    
    #identifying the best lambda
    best_lambda_null = cv_output_null$lambda.min
    best_lambda_alt = cv_output_alt$lambda.min
    best_lambda_null.int = cv_output_null.int$lambda.min
    best_lambda_alt.int = cv_output_alt.int$lambda.min
    
    #fit the null and alternative models with the best lambda values
    fitNull = glmnet(RGenos, ContPhenos, family = "gaussian", alpha=1, lambda = best_lambda_null, intercept=FALSE)
    fitAlt = glmnet(allData, ContPhenos, family = "gaussian", alpha=1, lambda = best_lambda_alt, intercept=FALSE)	
    fitNull.int = glmnet(RGenos, ContPhenos, family = "gaussian", alpha=1, lambda = best_lambda_null.int, intercept=TRUE)
    fitAlt.int = glmnet(allData, ContPhenos, family = "gaussian", alpha=1, lambda = best_lambda_alt.int, intercept=TRUE)	
    
    nullValues = coef(fitNull)
    altValues = coef(fitAlt)
    nullValues.int = coef(fitNull.int)
    altValues.int = coef(fitAlt.int)
    
    nullValues.lambda = fitNull$lambda
    altValues.lambda = fitAlt$lambda
    nullValues.int.lambda = fitNull.int$lambda
    altValues.int.lambda = fitAlt.int$lambda
    
    lambdas = c(nullValues.lambda, altValues.lambda, nullValues.int.lambda, altValues.int.lambda)
    
    print(paste0("Simulation ",ii," is complete."))
    
    finalOutput[[ii]] = list(nullValues, altValues, nullValues.int, altValues.int, lambdas)
  }, mc.cores=1)
  
  betasNull = list()
  betasAll = list()
  betasNull.int = list()
  betasAll.int = list()
  lambdas = list()
  
  for(ii in 1:length(statsAndPVals)){
    betasNull[[ii]] = statsAndPVals[[ii]][1]
    betasAll[[ii]] =  statsAndPVals[[ii]][2]
    betasNull.int[[ii]] = statsAndPVals[[ii]][3]
    betasAll.int[[ii]] = statsAndPVals[[ii]][4]
    lambdas[[ii]] = statsAndPVals[[ii]][5]
  }
  
  betasNull.mat =  do.call(rbind, unlist(betasNull, recursive = FALSE))
  betasAll.mat = do.call(rbind, unlist(betasAll, recursive = FALSE))
  betasNull.int.mat = do.call(rbind, unlist(betasNull.int, recursive = FALSE))
  betasAll.int.mat = do.call(rbind, unlist(betasAll.int, recursive = FALSE))
  lambdas.mat = matrix(unlist(lambdas), ncol = 4, byrow = TRUE)
  
  betasNullmat = as.matrix(betasNull.mat)
  betasAllmat = as.matrix(betasAll.mat)
  betasNull.intmat = as.matrix(betasNull.int.mat)
  betasAll.intmat = as.matrix(betasAll.int.mat)
  
  #write out the Stats and p values
  if(weightedScores == FALSE){
    if(standardizeScores == FALSE){
      write.csv(betasNullmat, file = paste0(path,"/TIE_ContNullBetas_LinHypTest_Score",score,"_Sim",start,"to",numSims,"_StatsAndPValues.csv"))
      write.csv(betasAllmat, file = paste0(path,"/TIE_ContAltBetas_LinHypTest_Score",score,"_Sim",start,"to",numSims,"_StatsAndPValues.csv"))
      write.csv(betasNull.intmat, file = paste0(path,"/TIE_ContNullBetasIntercept_LinHypTest_Score",score,"_Sim",start,"to",numSims,"_StatsAndPValues.csv"))
      write.csv(betasAll.intmat, file = paste0(path,"/TIE_ContAltBetasIntercept_LinHypTest_Score",score,"_Sim",start,"to",numSims,"_StatsAndPValues.csv"))
      write.csv(lambdas.mat, file = paste0(path,"/TIE_ContLamdbas_LinHypTest_Score",score,"_Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    } else { #scores are standardized
      write.csv(betasNullmat, file = paste0(path,"/TIE_ContNullBetas_LinHypTest_Score",score,"_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
      write.csv(betasAllmat, file = paste0(path,"/TIE_ContAltBetas_LinHypTest_Score",score,"_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
      write.csv(betasNull.intmat, file = paste0(path,"/TIE_ContNullBetasIntercept_LinHypTest_Score",score,"_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
      write.csv(betasAll.intmat, file = paste0(path,"/TIE_ContAltBetasIntercept_LinHypTest_Score",score,"_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
      write.csv(lambdas.mat, file = paste0(path,"/TIE_ContLamdbas_LinHypTest_Score",score,"_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
    }
  } else { #scores are weighted
    write.csv(betasNullmat, file = paste0(path,"/TIE_ContNullBetas_LinHypTest_Score",score,"_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
    write.csv(betasAllmat, file = paste0(path,"/TIE_ContAltBetas_LinHypTest_Score",score,"_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
    write.csv(betasNull.intmat, file = paste0(path,"/TIE_ContNullBetasIntercept_LinHypTest_Score",score,"_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
    write.csv(betasAll.intmat, file = paste0(path,"/TIE_ContAltBetasIntercept_LinHypTest_Score",score,"_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
    write.csv(lambdas.mat, file = paste0(path,"/TIE_ContLamdbas_LinHypTest_Score",score,"_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
  }
}

#categorical outcome
RunTIEPipelineLinHypTestCat_GammaTest_BICvsCV = function(chr, gene, numPairs, YPrev, standardizeScores = FALSE, weightedScores = FALSE, scoreWeights, score, start, stop){
  #function to run whole TIE pipeline, Calculates Score test stat for Lin Hyp test method and whether p-value <0.05
  # Only use when Y is binary
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #standardizeScores = T or F whether the scores should be standardized based on maximum score value
  #weightedScores = T or F whether the scores will be weighted
  #scoreWeights = m x 1 vector of weights, one weight for each SNP
  # score - which score we are using in the main JST model
  # start - simulation number to start at
  # stop - simulation number to stop at
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes TIE values to csv files

  library(parallel)

  #always the same
  numSims = stop

  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_II_Sims/HapGen_Files/",gene,"_Results_",numPairs,"Pairs")

  #source the needed functions
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/ProjectIISourceFunctions_v2.R")
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/Logistic_ADMM0.r")
  
  myList = lapply(start:numSims, rep, times = 1)
  # p-value initialization
  pv1 <- rep(0, 3)
  pv2 <- rep(0, 3)
  Tall1 <- matrix(0, 1, 3)
  Tall2 <- matrix(0, 1, 3)
  beta.al1 <- matrix(0, 1, 2)
  beta.al2 <- matrix(0, 1, 2)

  statsAndPVals = mclapply(myList, function(ii){
    #define matrix to hold all Stats and Pvalues
    statsAndPVals = matrix(NA, nrow = numSims, ncol = 4)

    #pull recipient and donor genotypes
    RGenos = obtainRGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    DGenos = obtainDGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)

    #calculate single snp scores
    if(score == "IBS"){
      Score.snp = calcIBSMismatch(RGenosMat = RGenos, DGenosMat = DGenos)
    } else if(score == "Incomp"){
      Score.snp = calcIncompatibilityScore(RGenosMat = RGenos, DGenosMat = DGenos)
    } else if(score == "AMS"){
      Score.snp = calcAMS(RGenosMat = RGenos, DGenosMat = DGenos)
    } else{
      Score.snp = calcBinaryMM(RGenosMat = RGenos, DGenosMat = DGenos)
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

    #generate covariates
    #for now, a single binary and a single continous covariate
    # CovData = GenCovData(SampleSize = numPairs, BinaryValues = 1, ContinuousValues = 1)

    #generate phenotypes, both continuous and binary
    CatPhenos = GenNullPhenos(SampleSize = numPairs, includeCov = FALSE, YCat = TRUE, YPrev = YPrev)
    # CatPhenos = GenNullPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = TRUE, YPrev = YPrev,  Covariates = CovData)

    # define location of zero components
    # basically, this defines what the null hypothesis is, I think
    # so for joint null, we need the non-zero components to be for covariates
    # and the beta and gamma components to be zero
    # numCov = dim(CovData)[2]
    numSNPs = dim(RGenos)[2]
    # N = c(rep(FALSE, numCov+1+numSNPs), TRUE)
    N = c(rep(FALSE, numSNPs+1), TRUE)

    #need to combine the covariates, r genos, and score matrices together
    #need to include column of 1s for intercept?
    intercept = matrix(1,nrow = numPairs, ncol = 1)
    designMat = cbind(intercept, RGenos, Score.gene)
    # designMat = cbind(intercept, CovData, RGenos, Score.gene)

    #then combine the phenos with the design matrix as a list
    Model = list(X=designMat, Y=CatPhenos)

    # estimate the uncontrained estimator
    beta.unre1 <- cv.SCAD_ADMM_unre(X=Model$X, Y=Model$Y, N=N, beta0=rep(0, dim(Model$X)[2]), err=1e-4, tune="BIC", unpen=c())
    indice.unre1 <- beta.unre1!=0

    beta.unre2 <- cv.SCAD_ADMM_unre(X=Model$X, Y=Model$Y, N=N, beta0=rep(0, dim(Model$X)[2]), err=1e-4, tune="cv", unpen=c())
    indice.unre2 <- beta.unre2!=0
    # estimate the constrained estimator (only need this for score test)
    beta.re1 <- cv.SCAD_ADMM_re(X=Model$X, Y=Model$Y, N=N, beta0=rep(0, dim(Model$X)[2]), err=1e-4, tune="BIC", unpen=c())
    indice.re1 <- beta.re1!=0

    beta.re2 <- cv.SCAD_ADMM_re(X=Model$X, Y=Model$Y, N=N, beta0=rep(0, dim(Model$X)[2]), err=1e-4, tune="cv", unpen=c())
    indice.re2 <- beta.re2!=0
    
    pi.unre1 <- logit(Model$X%*%beta.unre1)
    pi.re1 <- logit(Model$X%*%beta.re1)

    pi.unre2 <- logit(Model$X%*%beta.unre2)
    pi.re2 <- logit(Model$X%*%beta.re2)
    
    # construct the likelihood ratio statistics
    TL1 <- 2*(sum(log(1+exp(Model$X%*%beta.re1))-Model$Y*(Model$X%*%beta.re1))-
               sum(log(1+exp(Model$X%*%beta.unre1))-Model$Y*(Model$X%*%beta.unre1)))
    #if LRT stat is less than 0, there is error, so set to -1
    if(TL1 < 0){
      TL1 = -1
    }
    
    TL2 <- 2*(sum(log(1+exp(Model$X%*%beta.re2))-Model$Y*(Model$X%*%beta.re2))-
               sum(log(1+exp(Model$X%*%beta.unre2))-Model$Y*(Model$X%*%beta.unre2)))
    if(TL2 < 0){
      TL2 = -1
    }
    
    # construct the Wald statistics
    #B_0 should give you Omega_a hat
    A = crossprod(Model$X[,indice.unre1|N], as.vector(pi.unre1*(1-pi.unre1))*Model$X[,indice.unre1|N])
    if(rcond(A) >= 1e-10){
      B_01 <- solve(A)
      #d_0 should determine which rows to subset
      #gives the first value that is restricted
      d_01 <- sum(which(indice.unre1|N)<=length(N))
      #so the B_0 needs to be subset to m rows and columns that are restricted based on H0
      TW1 <- crossprod(beta.unre1[N], solve(B_01[d_01,d_01], beta.unre1[N]))
      
    } else {
      TW1 = -1
    }
    
    B = crossprod(Model$X[,indice.unre2|N], as.vector(pi.unre2*(1-pi.unre2))*Model$X[,indice.unre2|N])
    if(rcond(B) >= 1e-10){
      B_02 <- solve(B)
      #d_0 should determine which rows to subset
      #gives the first value that is restricted
      d_02 <- sum(which(indice.unre2|N)<=length(N))
      #so the B_0 needs to be subset to m rows and columns that are restricted based on H0
      TW2 <- crossprod(beta.unre2[N], solve(B_02[d_02,d_02], beta.unre2[N]))
    } else {
      TW2 = -1
    }
    
    # construct the score statistics
    eps1 <- Model$Y-pi.re1 #Y - e(Y)
    eps2 <- Model$Y-pi.re2 #Y - e(Y)
    #this is the X^T times Y-E(Y)
    Xeps1 <- crossprod(Model$X[,indice.re1|N], eps1)
    Xeps2 <- crossprod(Model$X[,indice.re2|N], eps2)
    
    C = crossprod(Model$X[,indice.re1|N], as.vector(pi.re1*(1-pi.re1))*Model$X[,indice.re1|N])
    if(rcond(C) >= 1e-10){
      TS1 <- crossprod(Xeps1, solve(C, Xeps1))
    } else {
      TS1 = -1
    }
    
    D = crossprod(Model$X[,indice.re2|N], as.vector(pi.re2*(1-pi.re2))*Model$X[,indice.re2|N])
    if(rcond(D) >= 1e-10){
      TS2 <- crossprod(Xeps2, solve(D, Xeps2))
    } else {
      TS2 = -1
    }
    
    #determine if null hyp is rejected or not
    #degrees of freedom is equal to the number of restricted values
    doF = sum(N == TRUE)

    #set pv to 1 if reject, 0 otherwise
    #LRT
    if (TL1>=qchisq(0.95, df = doF)){
      pv1[1] <- pv1[1]+1/5000
    }
    if (TL2>=qchisq(0.95, df = doF)){
      pv2[1] <- pv2[1]+1/5000
    }
    #Wald
    if (TW1>=qchisq(0.95, df = doF)){
      pv1[2] <- pv1[2]+1/5000
    }
    if (TW2>=qchisq(0.95, df = doF)){
      pv2[2] <- pv2[2]+1/5000
    }
    #Score
    if (TS1>=qchisq(0.95, df = doF)){
      pv1[3] <- pv1[3]+1/5000
    }
    if (TS2>=qchisq(0.95, df = doF)){
      pv2[3] <- pv2[3]+1/5000
    }
    
    print(paste0("Simulation ",ii," is complete."))

    Tall1[1, ] <- c(TL1, TW1, TS1)
    beta.al1[1, ] <- c(sum(beta.re1!=0), sum(beta.unre1!=0))
    
    Tall2[1, ] <- c(TL2, TW2, TS2)
    beta.al2[1, ] <- c(sum(beta.re2!=0), sum(beta.unre2!=0))
    
    statsAndPVals = list(pv1=pv1, TS1=Tall1, beta1=beta.al1, 
                         pv2=pv2, TS2=Tall2, beta2=beta.al2)
    statsAndPVals
  }, mc.cores=4)

  n = length(statsAndPVals)
  
  pValsBIC = matrix(nrow = n, ncol = 3)
  pValsCV = matrix(nrow = n, ncol = 3)
  ScoresBIC = matrix(nrow = n, ncol = 3)
  ScoresCV = matrix(nrow = n, ncol = 3)
  BetasBIC = matrix(nrow = n, ncol = 2)
  BetasCV = matrix(nrow = n, ncol = 2)

  for(ll in 1:n){
    pValsBIC[ll,] = unlist(statsAndPVals[[ll]]$pv1)
    pValsCV[ll,] = unlist(statsAndPVals[[ll]]$pv2)
    ScoresBIC[ll,] = unlist(statsAndPVals[[ll]]$TS1)
    ScoresCV[ll,] = unlist(statsAndPVals[[ll]]$TS2)
    BetasBIC[ll,] = unlist(statsAndPVals[[ll]]$beta1)
    BetasCV[ll,] = unlist(statsAndPVals[[ll]]$beta2)
  }
  
  #write out the Stats and p values
  if(weightedScores == FALSE){
    if(standardizeScores == FALSE){
      write.csv(pValsBIC, file = paste0(path,"/TIE_CatY_Prev",YPrev*100,"_LinHypTest_Score",score,"_GammaOnly_NoCovs_Sim",start,"to",numSims,"_PValuesBIC.csv"))
      write.csv(pValsCV, file = paste0(path,"/TIE_CatY_Prev",YPrev*100,"_LinHypTest_Score",score,"_GammaOnly_NoCovs_Sim",start,"to",numSims,"_PValuesCV.csv"))
      write.csv(ScoresBIC, file = paste0(path,"/TIE_CatY_Prev",YPrev*100,"_LinHypTest_Score",score,"_GammaOnly_NoCovs_Sim",start,"to",numSims,"_StatsBIC.csv"))
      write.csv(ScoresCV, file = paste0(path,"/TIE_CatY_Prev",YPrev*100,"_LinHypTest_Score",score,"_GammaOnly_NoCovs_Sim",start,"to",numSims,"_StatsCV.csv"))
      write.csv(BetasBIC, file = paste0(path,"/TIE_CatY_Prev",YPrev*100,"_LinHypTest_Score",score,"_GammaOnly_NoCovs_Sim",start,"to",numSims,"_BetasBIC.csv"))
      write.csv(BetasCV, file = paste0(path,"/TIE_CatY_Prev",YPrev*100,"_LinHypTest_Score",score,"_GammaOnly_NoCovs_Sim",start,"to",numSims,"_BetasCV.csv"))
    } else { #scores are standardized
      write.csv(statsAndPValsCatPhenos, file = paste0(path,"/TIE_CatPhenos_Prev",YPrev*100,"_LinHypTest_Score",score,"_GammaOnly_NoCovs_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
    }
  } else { #scores are weighted
    write.csv(statsAndPValsCatPhenos, file = paste0(path,"/TIE_CatPhenos_Prev",YPrev*100,"_LinHypTest_Score",score,"_GammaOnly_NoCovs_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
  }
}

#continuous outcome
RunTIEPipelineLinHypTestCont_GammaTest_BICvsCV = function(chr, gene, numPairs, standardizeScores = FALSE, weightedScores = FALSE, scoreWeights, score, start, stop){
  #function to run whole TIE pipeline, Calculates LRT test stat for Lin Hyp test method and whether p-value <0.05
  # only use when Y is continuous (Normal)
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #standardizeScores = T or F whether the scores should be standardized based on maximum score value
  #weightedScores = T or F whether the scores will be weighted
  #scoreWeights = m x 1 vector of weights, one weight for each SNP
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes TIE values to csv files

  library(parallel)

  #always the same
  numSims = stop

  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_II_Sims/HapGen_Files/",gene,"_Results_",numPairs,"Pairs")

  #source the needed functions
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/ProjectIISourceFunctions_v2.R")
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/Linear_ADMM0.r")

  myList = lapply(start:numSims, rep, times = 1)
  # p-value initialization
  pv1 <- rep(0, 6)
  pv2 <- rep(0, 6)
  Tall1 <- matrix(0, 1, 3)
  Tall2 <- matrix(0, 1, 3)
  beta.al1 <- matrix(0, 1, 2)
  beta.al2 <- matrix(0, 1, 2)

  statsAndPVals = mclapply(myList, function(ii){
    #define matrix to hold all Stats and Pvalues
    statsAndPVals = matrix(NA, nrow = numSims, ncol = 6)

    #pull recipient and donor genotypes
    RGenos = obtainRGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    DGenos = obtainDGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)

    #calculate single snp scores
    if(score == "IBS"){
      Score.snp = calcIBSMismatch(RGenosMat = RGenos, DGenosMat = DGenos)
    } else if(score == "Incomp"){
      Score.snp = calcIncompatibilityScore(RGenosMat = RGenos, DGenosMat = DGenos)
    } else if(score == "AMS"){
      Score.snp = calcAMS(RGenosMat = RGenos, DGenosMat = DGenos)
    } else{
      Score.snp = calcBinaryMM(RGenosMat = RGenos, DGenosMat = DGenos)
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

    #generate covariates
    #for now, a single binary and a single continous covariate
    # CovData = GenCovData(SampleSize = numPairs, BinaryValues = 1, ContinuousValues = 1)

    #generate phenotypes, both continuous and binary
    ContPhenos = GenNullPhenos(SampleSize = numPairs, includeCov = FALSE, YCat = FALSE)
    # ContPhenos = GenNullPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = FALSE,  Covariates = CovData)

    # define location of zero components
    # basically, this defines what the null hypothesis is, I think
    # so for joint null, we need the non-zero components to be for covariates
    # and the beta and gamma components to be zero
    # numCov = dim(CovData)[2]
    numSNPs = dim(RGenos)[2]
    N = c(rep(FALSE, numSNPs+1), TRUE)
    # N = c(rep(FALSE, numCov+1+numSNPs), TRUE)

    #need to combine the covariates, r genos, and score matrices together
    #need to include column of 1s for intercept?
    intercept = matrix(1,nrow = numPairs, ncol = 1)
    designMat = cbind(intercept, RGenos, Score.gene)
    # designMat = cbind(intercept, CovData, RGenos, Score.gene)

    #then combine the phenos with the design matrix as a list
    Model = list(X=designMat, Y=ContPhenos)

    # estimate the uncontrained estimator
    beta.unre1 <- cv.SCAD_ADMM_unre(X=Model$X, Y=Model$Y, N=N, beta0=rep(0, dim(Model$X)[2]), err=1e-4, tune="BIC")
    indice.unre1 <- beta.unre1!=0
    
    beta.unre2 <- cv.SCAD_ADMM_unre(X=Model$X, Y=Model$Y, N=N, beta0=rep(0, dim(Model$X)[2]), err=1e-4, tune="cv")
    indice.unre2 <- beta.unre2!=0

    # estimate the constrained estimator
    beta.re1 <- cv.SCAD_ADMM_re(X=Model$X, Y=Model$Y, N=N, beta0=rep(0, dim(Model$X)[2]), err=1e-4, tune="BIC")
    indice.re1 <- beta.re1!=0
    
    beta.re2 <- cv.SCAD_ADMM_re(X=Model$X, Y=Model$Y, N=N, beta0=rep(0, dim(Model$X)[2]), err=1e-4, tune="cv")
    indice.re2 <- beta.re2!=0
    
    # estimate the conditional variance
    #    sig2 <- mean((Model$Y-Model$X%*%beta.unre)^2)
    n = numSims
    sig21 <- mean((Model$Y-Model$X%*%beta.unre1)^2)*n/(n-sum(beta.unre1!=0))
    sig22 <- mean((Model$Y-Model$X%*%beta.unre2)^2)*n/(n-sum(beta.unre2!=0))
    
    # construct the likelihood ratio statistic
    TL1 <- sum((Model$X%*%beta.re1-Model$Y)^2)-sum((Model$X%*%beta.unre1-Model$Y)^2)
    if(TL1 < 0){
      TL1 = -1
    }
    
    TL2 <- sum((Model$X%*%beta.re2-Model$Y)^2)-sum((Model$X%*%beta.unre2-Model$Y)^2)
    if(TL2 < 0){
      TL2 = -1
    }
    
    # construct the Wald statistic
    A = crossprod(Model$X[,indice.unre1|N], Model$X[,indice.unre1|N])
    if(rcond(A) >= 1e-10){
      B_01 <- solve(A)
      d_01 <- sum(which(indice.unre1|N)<=length(N))
      TW1 <- crossprod(beta.unre1[N], solve(B_01[d_01,d_01], beta.unre1[N]))
    } else {
      TW1 = -1
    }

    B = crossprod(Model$X[,indice.unre2|N], Model$X[,indice.unre2|N])
    if(rcond(B) >= 1e-10){
      B_02 <- solve(B)
      d_02 <- sum(which(indice.unre2|N)<=length(N))
      TW2 <- crossprod(beta.unre2[N], solve(B_02[d_02,d_02], beta.unre2[N]))
    } else {
      TW2 = -1
    }
    
    # construct the score statistic
    eps1 <- Model$Y-Model$X%*%beta.re1
    eps2 <- Model$Y-Model$X%*%beta.re2

    Xeps1 <- crossprod(Model$X[,indice.re1|N], eps1)
    Xeps2 <- crossprod(Model$X[,indice.re2|N], eps2)

    C = crossprod(Model$X[,indice.re1|N], Model$X[,indice.re1|N])
    if(rcond(C) >= 1e-10){
      TS1 <- crossprod(Xeps1, solve(C, Xeps1))
    } else {
      TS1 = -1
    }
    
    D = crossprod(Model$X[,indice.re2|N], Model$X[,indice.re2|N])
    if(rcond(D) >= 1e-10){
      TS2 <- crossprod(Xeps2, solve(D, Xeps2))
    } else {
      TS2 = -1
    }
    
    #determine if null hyp is rejected or not
    #degrees of freedom is equal to the number of restricted values
    doF = sum(N == TRUE)

    #set pv to 1 if reject, 0 otherwise
    #LRT
    if (TL1>=sig21*qchisq(0.95, df = doF)){
      pv1[1] <- pv1[1]+1/5000
    }
    if (TL2>=sig22*qchisq(0.95, df = doF)){
      pv2[1] <- pv2[1]+1/5000
    }
    if (TL1>=qchisq(0.95, df = doF)){
      pv1[2] <- pv1[2]+1/5000
    }
    if (TL2>=qchisq(0.95, df = doF)){
      pv2[2] <- pv2[2]+1/5000
    }
    #Wald
    if (TW1>=sig21*qchisq(0.95, df = doF)){
      pv1[3] <- pv1[3]+1/5000
    }
    if (TW2>=sig22*qchisq(0.95, df = doF)){
      pv2[3] <- pv2[3]+1/5000
    }
    if (TW1>=qchisq(0.95, df = doF)){
      pv1[4] <- pv1[4]+1/5000
    }
    if (TW2>=qchisq(0.95, df = doF)){
      pv2[4] <- pv2[4]+1/5000
    }
    #Score
    if (TS1>=sig21*qchisq(0.95, df = doF)){
      pv1[5] <- pv1[5]+1/5000
    }
    if (TS2>=sig22*qchisq(0.95, df = doF)){
      pv2[5] <- pv2[5]+1/5000
    }
    if (TS1>=qchisq(0.95, df = doF)){
      pv1[6] <- pv1[6]+1/5000
    }
    if (TS2>=qchisq(0.95, df = doF)){
      pv2[6] <- pv2[6]+1/5000
    }

    print(paste0("Simulation ",ii," is complete."))

    Tall1[1, ] <- c(TL1, TW1, TS1)
    beta.al1[1, ] <- c(sum(beta.re1!=0), sum(beta.unre1!=0))
    
    Tall2[1, ] <- c(TL2, TW2, TS2)
    beta.al2[1, ] <- c(sum(beta.re2!=0), sum(beta.unre2!=0))
    
    statsAndPVals = list(pv1=pv1, TS1=Tall1, beta1=beta.al1, 
                         pv2=pv2, TS2=Tall2, beta2=beta.al2)
    statsAndPVals
  }, mc.cores=4)
  
  n = length(statsAndPVals)
  
  pValsBIC = matrix(nrow = n, ncol = 6)
  pValsCV = matrix(nrow = n, ncol = 6)
  ScoresBIC = matrix(nrow = n, ncol = 3)
  ScoresCV = matrix(nrow = n, ncol = 3)
  BetasBIC = matrix(nrow = n, ncol = 2)
  BetasCV = matrix(nrow = n, ncol = 2)
  
  for(ll in 1:n){
    pValsBIC[ll,] = unlist(statsAndPVals[[ll]]$pv1)
    pValsCV[ll,] = unlist(statsAndPVals[[ll]]$pv2)
    ScoresBIC[ll,] = unlist(statsAndPVals[[ll]]$TS1)
    ScoresCV[ll,] = unlist(statsAndPVals[[ll]]$TS2)
    BetasBIC[ll,] = unlist(statsAndPVals[[ll]]$beta1)
    BetasCV[ll,] = unlist(statsAndPVals[[ll]]$beta2)
  }
  
  
  #write out the Stats and p values
  if(weightedScores == FALSE){
    if(standardizeScores == FALSE){
      write.csv(pValsBIC, file = paste0(path,"/TIE_ContY_LinHypTest_Score",score,"_GammaOnly_NoCovs_Sim",start,"to",numSims,"_PValuesBIC.csv"))
      write.csv(pValsCV, file = paste0(path,"/TIE_ContY_LinHypTest_Score",score,"_GammaOnly_NoCovs_Sim",start,"to",numSims,"_PValuesCV.csv"))
      write.csv(ScoresBIC, file = paste0(path,"/TIE_ContY_LinHypTest_Score",score,"_GammaOnly_NoCovs_Sim",start,"to",numSims,"_StatsBIC.csv"))
      write.csv(ScoresCV, file = paste0(path,"/TIE_ContY_LinHypTest_Score",score,"_GammaOnly_NoCovs_Sim",start,"to",numSims,"_StatsCV.csv"))
      write.csv(BetasBIC, file = paste0(path,"/TIE_ContY_LinHypTest_Score",score,"_GammaOnly_NoCovs_Sim",start,"to",numSims,"_BetasBIC.csv"))
      write.csv(BetasCV, file = paste0(path,"/TIE_ContY_LinHypTest_Score",score,"_GammaOnly_NoCovs_Sim",start,"to",numSims,"_BetasCV.csv"))
    } else { #scores are standardized
      write.csv(statsAndPValsContPhenos, file = paste0(path,"/TIE_ContPhenos_LinHypTest_Score",score,"_GammaOnly_NoCovs_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
    }
  } else { #scores are weighted
    write.csv(statsAndPValsContPhenos, file = paste0(path,"/TIE_ContPhenos_LinHypTest_Score",score,"_GammaOnly_NoCovs_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
  }
}

#categorical outcome
RunTIEPipelineLinHypTestCat_GammaTest = function(chr, gene, numPairs, YPrev, standardizeScores = FALSE, weightedScores = FALSE, scoreWeights, score, start, stop){
  #function to run whole TIE pipeline, Calculates Score test stat for Lin Hyp test method and whether p-value <0.05
  # Only use when Y is binary
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #standardizeScores = T or F whether the scores should be standardized based on maximum score value
  #weightedScores = T or F whether the scores will be weighted
  #scoreWeights = m x 1 vector of weights, one weight for each SNP
  # score - which score we are using in the main JST model
  # start - simulation number to start at
  # stop - simulation number to stop at
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes TIE values to csv files
  
  library(parallel)
  
  #always the same
  numSims = stop
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_II_Sims/HapGen_Files/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/ProjectIISourceFunctions_v2.R")
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/Logistic_ADMM0.r")
  
  myList = lapply(start:numSims, rep, times = 1)
  # p-value initialization
  pv <- rep(0, 3)
  Tall <- matrix(0, 1, 3)
  beta.al <- matrix(0, 1, 2)

  statsAndPVals = mclapply(myList, function(ii){
    #define matrix to hold all Stats and Pvalues
    statsAndPVals = matrix(NA, nrow = numSims, ncol = 4)
    
    #pull recipient and donor genotypes
    RGenos = obtainRGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    DGenos = obtainDGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    
    #calculate single snp scores
    if(score == "IBS"){
      Score.snp = calcIBSMismatch(RGenosMat = RGenos, DGenosMat = DGenos)
    } else if(score == "Incomp"){
      Score.snp = calcIncompatibilityScore(RGenosMat = RGenos, DGenosMat = DGenos)
    } else if(score == "AMS"){
      Score.snp = calcAMS(RGenosMat = RGenos, DGenosMat = DGenos)
    } else{
      Score.snp = calcBinaryMM(RGenosMat = RGenos, DGenosMat = DGenos)
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
    
    #generate covariates
    #for now, a single binary and a single continous covariate
    CovData = GenCovData(SampleSize = numPairs, BinaryValues = 1, ContinuousValues = 1)
    
    #generate phenotypes, both continuous and binary
    # CatPhenos = GenNullPhenos(SampleSize = numPairs, includeCov = FALSE, YCat = TRUE, YPrev = YPrev)
    CatPhenos = GenNullPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = TRUE, YPrev = YPrev,  Covariates = CovData)
    
    # define location of zero components
    # basically, this defines what the null hypothesis is, I think
    # so for joint null, we need the non-zero components to be for covariates
    # and the beta and gamma components to be zero
    numCov = dim(CovData)[2]
    numSNPs = dim(RGenos)[2]
    N = c(rep(FALSE, numCov+1+numSNPs), TRUE)
    # N = c(rep(FALSE, numSNPs+1), TRUE)
    
    #need to combine the covariates, r genos, and score matrices together
    #need to include column of 1s for intercept?
    intercept = matrix(1,nrow = numPairs, ncol = 1)
    # designMat = cbind(intercept, RGenos, Score.gene)
    designMat = cbind(intercept, CovData, RGenos, Score.gene)
    
    #then combine the phenos with the design matrix as a list
    Model = list(X=designMat, Y=CatPhenos)
    
    # estimate the uncontrained estimator
    beta.unre <- cv.SCAD_ADMM_unre(X=Model$X, Y=Model$Y, N=N, beta0=rep(0, dim(Model$X)[2]), err=1e-4, tune="cv", unpen=c(1,2,3))
    indice.unre <- beta.unre!=0
    # estimate the constrained estimator (only need this for score test)
    beta.re <- cv.SCAD_ADMM_re(X=Model$X, Y=Model$Y, N=N, beta0=rep(0, dim(Model$X)[2]), err=1e-4, tune="cv", unpen=c(1,2,3))
    indice.re <- beta.re!=0
    
    pi.unre <- logit(Model$X%*%beta.unre)
    pi.re <- logit(Model$X%*%beta.re)
    
    # construct the likelihood ratio statistics
    TL <- 2*(sum(log(1+exp(Model$X%*%beta.re))-Model$Y*(Model$X%*%beta.re))-
                sum(log(1+exp(Model$X%*%beta.unre))-Model$Y*(Model$X%*%beta.unre)))
    #if LRT stat is less than 0, there is error, so set to -1
    if(TL < 0){
      TL = -1
    }
    
    # construct the Wald statistics
    #B_0 should give you Omega_a hat
    A = crossprod(Model$X[,indice.unre|N], as.vector(pi.unre*(1-pi.unre))*Model$X[,indice.unre|N])
    if(rcond(A) >= 1e-10){
      B_01 <- solve(A)
      #d_0 should determine which rows to subset
      #gives the first value that is restricted
      d_01 <- sum(which(indice.unre|N)<=length(N))
      #so the B_0 needs to be subset to m rows and columns that are restricted based on H0
      TW <- crossprod(beta.unre[N], solve(B_01[d_01,d_01], beta.unre[N]))
    } else {
      TW = -1
    }
    
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
    
    #set pv to 1 if reject, 0 otherwise
    #LRT
    if (TL>=qchisq(0.95, df = doF)){
      pv[1] <- pv[1]+1/5000
    }
    #Wald
    if (TW>=qchisq(0.95, df = doF)){
      pv[2] <- pv[2]+1/5000
    }
    #Score
    if (TS>=qchisq(0.95, df = doF)){
      pv[3] <- pv[3]+1/5000
    }
    
    print(paste0("Simulation ",ii," is complete."))
    
    Tall[1, ] <- c(TL, TW, TS)
    beta.al[1, ] <- c(sum(beta.re!=0), sum(beta.unre!=0))
    
    statsAndPVals = list(pv=pv, TScores=Tall, beta=beta.al)
    statsAndPVals
  }, mc.cores=4)
  
  n = length(statsAndPVals)
  
  pValsCV = matrix(nrow = n, ncol = 3)
  ScoresCV = matrix(nrow = n, ncol = 3)
  BetasCV = matrix(nrow = n, ncol = 2)
  
  for(ll in 1:n){
    pValsCV[ll,] = unlist(statsAndPVals[[ll]]$pv)
    ScoresCV[ll,] = unlist(statsAndPVals[[ll]]$TScores)
    BetasCV[ll,] = unlist(statsAndPVals[[ll]]$beta)
  }
  
  #write out the Stats and p values
  if(weightedScores == FALSE){
    if(standardizeScores == FALSE){
      write.csv(pValsCV, file = paste0(path,"/TIE_CatY_Prev",YPrev*100,"_LinHypTest_Score",score,"_GammaOnly_NoCovs_Sim",start,"to",numSims,"_PValuesCV_ForceCovFit.csv"))
      write.csv(ScoresCV, file = paste0(path,"/TIE_CatY_Prev",YPrev*100,"_LinHypTest_Score",score,"_GammaOnly_NoCovs_Sim",start,"to",numSims,"_StatsCV_ForceCovFit.csv"))
      write.csv(BetasCV, file = paste0(path,"/TIE_CatY_Prev",YPrev*100,"_LinHypTest_Score",score,"_GammaOnly_NoCovs_Sim",start,"to",numSims,"_BetasCV_ForceCovFit.csv"))
    } else { #scores are standardized
      write.csv(statsAndPValsCatPhenos, file = paste0(path,"/TIE_CatPhenos_Prev",YPrev*100,"_LinHypTest_Score",score,"_GammaOnly_NoCovs_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
    }
  } else { #scores are weighted
    write.csv(statsAndPValsCatPhenos, file = paste0(path,"/TIE_CatPhenos_Prev",YPrev*100,"_LinHypTest_Score",score,"_GammaOnly_NoCovs_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
  }
}

#continuous outcome
RunTIEPipelineLinHypTestCont_GammaTest = function(chr, gene, numPairs, standardizeScores = FALSE, weightedScores = FALSE, scoreWeights, score, start, stop){
  #function to run whole TIE pipeline, Calculates LRT test stat for Lin Hyp test method and whether p-value <0.05
  # only use when Y is continuous (Normal)
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #standardizeScores = T or F whether the scores should be standardized based on maximum score value
  #weightedScores = T or F whether the scores will be weighted
  #scoreWeights = m x 1 vector of weights, one weight for each SNP
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes TIE values to csv files
  
  library(parallel)
  
  #always the same
  numSims = stop
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_II_Sims/HapGen_Files/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/ProjectIISourceFunctions_v2.R")
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/Linear_ADMM0.r")
  
  myList = lapply(start:numSims, rep, times = 1)
  # p-value initialization
  pv <- rep(0, 6)
  Tall <- matrix(0, 1, 3)
  beta.al <- matrix(0, 1, 2)
  
  statsAndPVals = mclapply(myList, function(ii){
    #define matrix to hold all Stats and Pvalues
    statsAndPVals = matrix(NA, nrow = numSims, ncol = 6)
    
    #pull recipient and donor genotypes
    RGenos = obtainRGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    DGenos = obtainDGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    
    #calculate single snp scores
    if(score == "IBS"){
      Score.snp = calcIBSMismatch(RGenosMat = RGenos, DGenosMat = DGenos)
    } else if(score == "Incomp"){
      Score.snp = calcIncompatibilityScore(RGenosMat = RGenos, DGenosMat = DGenos)
    } else if(score == "AMS"){
      Score.snp = calcAMS(RGenosMat = RGenos, DGenosMat = DGenos)
    } else{
      Score.snp = calcBinaryMM(RGenosMat = RGenos, DGenosMat = DGenos)
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
    
    #generate covariates
    #for now, a single binary and a single continous covariate
    CovData = GenCovData(SampleSize = numPairs, BinaryValues = 1, ContinuousValues = 1)
    
    #generate phenotypes, both continuous and binary
    # ContPhenos = GenNullPhenos(SampleSize = numPairs, includeCov = FALSE, YCat = FALSE)
    ContPhenos = GenNullPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = FALSE,  Covariates = CovData)
    
    # define location of zero components
    # basically, this defines what the null hypothesis is, I think
    # so for joint null, we need the non-zero components to be for covariates
    # and the beta and gamma components to be zero
    numCov = dim(CovData)[2]
    numSNPs = dim(RGenos)[2]
    # N = c(rep(FALSE, numSNPs+1), TRUE)
    N = c(rep(FALSE, numCov+1+numSNPs), TRUE)
    
    #need to combine the covariates, r genos, and score matrices together
    #need to include column of 1s for intercept?
    intercept = matrix(1,nrow = numPairs, ncol = 1)
    # designMat = cbind(intercept, RGenos, Score.gene)
    designMat = cbind(intercept, CovData, RGenos, Score.gene)
    
    #then combine the phenos with the design matrix as a list
    Model = list(X=designMat, Y=ContPhenos)
    
    # estimate the uncontrained estimator
    beta.unre <- cv.SCAD_ADMM_unre(X=Model$X, Y=Model$Y, N=N, beta0=rep(0, dim(Model$X)[2]), err=1e-4, tune="cv", unpen = c(1,2,3))
    indice.unre <- beta.unre!=0
    
    # estimate the constrained estimator
    beta.re <- cv.SCAD_ADMM_re(X=Model$X, Y=Model$Y, N=N, beta0=rep(0, dim(Model$X)[2]), err=1e-4, tune="cv", unpen = c(1,2,3))
    indice.re <- beta.re!=0
    
    # estimate the conditional variance
    #    sig2 <- mean((Model$Y-Model$X%*%beta.unre)^2)
    n = numSims
    sig2 <- mean((Model$Y-Model$X%*%beta.unre)^2)*n/(n-sum(beta.unre!=0))
    
    # construct the likelihood ratio statistic
    TL <- sum((Model$X%*%beta.re-Model$Y)^2)-sum((Model$X%*%beta.unre-Model$Y)^2)
    if(TL < 0){
      TL = -1
    }
    
    # construct the Wald statistic
    B = crossprod(Model$X[,indice.unre|N], Model$X[,indice.unre|N])
    if(rcond(B) >= 1e-10){
      B_0 <- solve(B)
      d_0 <- sum(which(indice.unre|N)<=length(N))
      TW <- crossprod(beta.unre[N], solve(B_0[d_0,d_0], beta.unre[N]))
    } else {
      TW = -1
    }
    
    # construct the score statistic
    eps <- Model$Y-Model$X%*%beta.re
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
    
    #set pv to 1 if reject, 0 otherwise
    #LRT
    if (TL>=sig2*qchisq(0.95, df = doF)){
      pv[1] <- pv[1]+1/5000
    }
    if (TL>=qchisq(0.95, df = doF)){
      pv[2] <- pv[2]+1/5000
    }
    #Wald
    if (TW>=sig2*qchisq(0.95, df = doF)){
      pv[3] <- pv[3]+1/5000
    }
    if (TW>=qchisq(0.95, df = doF)){
      pv[4] <- pv[4]+1/5000
    }
    #Score
    if (TS>=sig2*qchisq(0.95, df = doF)){
      pv[5] <- pv[5]+1/5000
    }
    if (TS>=qchisq(0.95, df = doF)){
      pv[6] <- pv[6]+1/5000
    }
    
    print(paste0("Simulation ",ii," is complete."))
    
    Tall[1, ] <- c(TL, TW, TS)
    beta.al[1, ] <- c(sum(beta.re!=0), sum(beta.unre!=0))
    
    statsAndPVals = list(pv=pv, TScores=Tall, beta=beta.al)
    statsAndPVals
  }, mc.cores=4)
  
  n = length(statsAndPVals)
  
  pVals = matrix(nrow = n, ncol = 6)
  Scores = matrix(nrow = n, ncol = 3)
  Betas = matrix(nrow = n, ncol = 2)
  
  for(ll in 1:n){
    pVals[ll,] = unlist(statsAndPVals[[ll]]$pv)
    Scores[ll,] = unlist(statsAndPVals[[ll]]$TScores)
    Betas[ll,] = unlist(statsAndPVals[[ll]]$beta)
  }
  
  
  #write out the Stats and p values
  if(weightedScores == FALSE){
    if(standardizeScores == FALSE){
      write.csv(pVals, file = paste0(path,"/TIE_ContY_LinHypTest_Score",score,"_GammaOnly_NoCovs_Sim",start,"to",numSims,"_PValuesCV_ForceCovs.csv"))
      write.csv(Scores, file = paste0(path,"/TIE_ContY_LinHypTest_Score",score,"_GammaOnly_NoCovs_Sim",start,"to",numSims,"_StatsCV_ForceCovs.csv"))
      write.csv(Betas, file = paste0(path,"/TIE_ContY_LinHypTest_Score",score,"_GammaOnly_NoCovs_Sim",start,"to",numSims,"_BetasCV_ForceCovs.csv"))
    } else { #scores are standardized
      write.csv(statsAndPValsContPhenos, file = paste0(path,"/TIE_ContPhenos_LinHypTest_Score",score,"_GammaOnly_NoCovs_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
    }
  } else { #scores are weighted
    write.csv(statsAndPValsContPhenos, file = paste0(path,"/TIE_ContPhenos_LinHypTest_Score",score,"_GammaOnly_NoCovs_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
  }
}

#############################
## Power Pipelines
## Score true (gamma neq 0), beta = 0

#categorical outcome
RunPowerPipelineLinHypTestCat_GammaTest_Score_BetaZero = function(chr, gene, numPairs, YPrev, Gamma, TrueScore, ORSize, standardizeScores = FALSE, weightedScores = FALSE, scoreWeights, start, stop, percentageAssoc, LowLD){
  #Function to determine power of linear hyp testing method when score is associated, gamma = 0 testing only
  # Only use when Y is binary
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #Gamma = effect size for score, length 1, can be 0
  #TrueScore = IBS.gene, Incomp.gene, AMS.gene, or BinMM.gene
  #ORSize = Small, Medium, or Large for what OR was used for the associated SNP/score
  #standardizeScores = T or F whether the scores should be standardized based on maximum score value
  #weightedScores = T or F whether the scores will be weighted
  #scoreWeights = m x 1 vector of weights, one weight for each SNP
  # start - simulation number to start at
  # stop - simulation number to stop at
  #percentageAssoc = percentage of SNPs associated with outcome (either 5, 25, 50, 75, or 100) 
  #LowLD = True or FALSE whether the associated SNPs are in low LD or high LD
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes power values to csv files
  
  library(parallel)
  
  #need to define this for naming at the end
  snpOrScore = TrueScore
  
  #always the same
  numSims = stop
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_II_Sims/HapGen_Files/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/ProjectIISourceFunctions_v2.R")
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/Logistic_ADMM0.r")
  
  myList = lapply(start:numSims, rep, times = 1)
  # p-value initialization
  pv.IBS <- rep(0, 3)
  pv.Incomp <- rep(0, 3)
  pv.AMS <- rep(0, 3)
  pv.BinMM <- rep(0, 3)
  
  Tall.IBS <- matrix(0, 1, 3)
  Tall.Incomp <- matrix(0, 1, 3)
  Tall.AMS <- matrix(0, 1, 3)
  Tall.BinMM <- matrix(0, 1, 3)
  
  beta.al.IBS <- matrix(0, 1, 2)
  beta.al.Incomp <- matrix(0, 1, 2)
  beta.al.AMS <- matrix(0, 1, 2)
  beta.al.BinMM <- matrix(0, 1, 2)
  
  statsAndPVals = mclapply(myList, function(ii){
    #define matrix to hold all Stats and Pvalues
    statsAndPVals = matrix(NA, nrow = numSims, ncol = 4)
    
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
    CatPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = TRUE, YPrev = YPrev,  Covariates = CovData, RGenoData = RGenos, ScoreData = PhenoScore, Betas = Betas, Gamma = Gamma)
    
    # define location of zero components
    # basically, this defines what the null hypothesis is, I think
    # so for joint null, we need the non-zero components to be for covariates
    # and the beta and gamma components to be zero
    numCov = dim(CovData)[2]
    numSNPs = dim(RGenos)[2]
    N = c(rep(FALSE, numCov+1+numSNPs), TRUE)

    #need to combine the covariates, r genos, and score matrices together
    #need to include column of 1s for intercept?
    intercept = matrix(1,nrow = numPairs, ncol = 1)
    designMat.IBS = cbind(intercept, CovData, RGenos, IBS.gene)
    designMat.Incomp = cbind(intercept, CovData, RGenos, Incomp.gene)
    designMat.AMS = cbind(intercept, CovData, RGenos, AMS.gene)
    designMat.BinMM = cbind(intercept, CovData, RGenos, BinMM.gene)
    
    #then combine the phenos with the design matrix as a list
    Model.IBS = list(X=designMat.IBS, Y=CatPhenos)
    Model.Incomp = list(X=designMat.Incomp, Y=CatPhenos)
    Model.AMS = list(X=designMat.AMS, Y=CatPhenos)
    Model.BinMM = list(X=designMat.BinMM, Y=CatPhenos)
    
    # estimate the uncontrained estimator
    beta.unre.IBS <- cv.SCAD_ADMM_unre(X=Model.IBS$X, Y=Model.IBS$Y, N=N, beta0=rep(0, dim(Model.IBS$X)[2]), err=1e-4, tune="cv", unpen=c(1,2,3))
    indice.unre.IBS <- beta.unre.IBS!=0
    
    beta.unre.Incomp <- cv.SCAD_ADMM_unre(X=Model.Incomp$X, Y=Model.Incomp$Y, N=N, beta0=rep(0, dim(Model.Incomp$X)[2]), err=1e-4, tune="cv", unpen=c(1,2,3))
    indice.unre.Incomp <- beta.unre.Incomp!=0
    
    beta.unre.AMS <- cv.SCAD_ADMM_unre(X=Model.AMS$X, Y=Model.AMS$Y, N=N, beta0=rep(0, dim(Model.AMS$X)[2]), err=1e-4, tune="cv", unpen=c(1,2,3))
    indice.unre.AMS <- beta.unre.AMS!=0
    
    beta.unre.BinMM <- cv.SCAD_ADMM_unre(X=Model.BinMM$X, Y=Model.BinMM$Y, N=N, beta0=rep(0, dim(Model.BinMM$X)[2]), err=1e-4, tune="cv", unpen=c(1,2,3))
    indice.unre.BinMM <- beta.unre.BinMM!=0
    
    #estimate the constrained estimator (only need this for score test)
    beta.re.IBS <- cv.SCAD_ADMM_re(X=Model.IBS$X, Y=Model.IBS$Y, N=N, beta0=rep(0, dim(Model.IBS$X)[2]), err=1e-4, tune="cv", unpen=c(1,2,3))
    indice.re.IBS <- beta.re.IBS!=0
    
    beta.re.Incomp <- cv.SCAD_ADMM_re(X=Model.Incomp$X, Y=Model.Incomp$Y, N=N, beta0=rep(0, dim(Model.Incomp$X)[2]), err=1e-4, tune="cv", unpen=c(1,2,3))
    indice.re.Incomp <- beta.re.Incomp!=0
    
    beta.re.AMS <- cv.SCAD_ADMM_re(X=Model.AMS$X, Y=Model.AMS$Y, N=N, beta0=rep(0, dim(Model.AMS$X)[2]), err=1e-4, tune="cv", unpen=c(1,2,3))
    indice.re.AMS <- beta.re.AMS!=0
    
    beta.re.BinMM <- cv.SCAD_ADMM_re(X=Model.BinMM$X, Y=Model.BinMM$Y, N=N, beta0=rep(0, dim(Model.BinMM$X)[2]), err=1e-4, tune="cv", unpen=c(1,2,3))
    indice.re.BinMM <- beta.re.BinMM!=0
    
    pi.unre.AMS <- logit(Model.AMS$X%*%beta.unre.AMS)
    pi.re.AMS <- logit(Model.AMS$X%*%beta.re.AMS)
    
    pi.unre.BinMM <- logit(Model.BinMM$X%*%beta.unre.BinMM)
    pi.re.BinMM <- logit(Model.BinMM$X%*%beta.re.BinMM)
    
    pi.unre.IBS <- logit(Model.IBS$X%*%beta.unre.IBS)
    pi.re.IBS <- logit(Model.IBS$X%*%beta.re.IBS)
    
    pi.unre.Incomp <- logit(Model.Incomp$X%*%beta.unre.Incomp)
    pi.re.Incomp <- logit(Model.Incomp$X%*%beta.re.Incomp)
    
    #construct the likelihood ratio statistics
    TL.IBS <- 2*(sum(log(1+exp(Model.IBS$X%*%beta.re.IBS))-Model.IBS$Y*(Model.IBS$X%*%beta.re.IBS))-
                   sum(log(1+exp(Model.IBS$X%*%beta.unre.IBS))-Model.IBS$Y*(Model.IBS$X%*%beta.unre.IBS)))
    #if LRT stat is less than 0, there is error, so set to -1
    if(TL.IBS < 0){
      TL.IBS = -1
    }
    
    TL.Incomp <- 2*(sum(log(1+exp(Model.Incomp$X%*%beta.re.Incomp))-Model.Incomp$Y*(Model.Incomp$X%*%beta.re.Incomp))-
                      sum(log(1+exp(Model.Incomp$X%*%beta.unre.Incomp))-Model.Incomp$Y*(Model.Incomp$X%*%beta.unre.Incomp)))
    #if LRT stat is less than 0, there is error, so set to -1
    if(TL.Incomp < 0){
      TL.Incomp = -1
    }
    
    TL.AMS <- 2*(sum(log(1+exp(Model.AMS$X%*%beta.re.AMS))-Model.AMS$Y*(Model.AMS$X%*%beta.re.AMS))-
                   sum(log(1+exp(Model.AMS$X%*%beta.unre.AMS))-Model.AMS$Y*(Model.AMS$X%*%beta.unre.AMS)))
    #if LRT stat is less than 0, there is error, so set to -1
    if(TL.AMS < 0){
      TL.AMS = -1
    }
    
    TL.BinMM <- 2*(sum(log(1+exp(Model.BinMM$X%*%beta.re.BinMM))-Model.BinMM$Y*(Model.BinMM$X%*%beta.re.BinMM))-
                     sum(log(1+exp(Model.BinMM$X%*%beta.unre.BinMM))-Model.BinMM$Y*(Model.BinMM$X%*%beta.unre.BinMM)))
    #if LRT stat is less than 0, there is error, so set to -1
    if(TL.BinMM < 0){
      TL.BinMM = -1
    }
    
    # construct the Wald statistics
    #B_0 should give you Omega_a hat
    A.IBS = crossprod(Model.IBS$X[,indice.unre.IBS|N], as.vector(pi.unre.IBS*(1-pi.unre.IBS))*Model.IBS$X[,indice.unre.IBS|N])
    if(rcond(A.IBS) >= 1e-10){
      B_0.IBS <- solve(A.IBS)
      #d_0 should determine which rows to subset
      #gives the first value that is restricted
      d_0.IBS <- dim(B_0.IBS)[1]
      #in case the inverse doesn't exist, need another if else
      B.IBS = B_0.IBS[d_0.IBS,d_0.IBS]
      #so the B_0 needs to be subset to m rows and columns that are restricted based on H0
      TW.IBS <- crossprod(beta.unre.IBS[N], solve(B.IBS, beta.unre.IBS[N]))
    } else {
      TW.IBS = -1
    }
    
    A.Incomp = crossprod(Model.Incomp$X[,indice.unre.Incomp|N], as.vector(pi.unre.Incomp*(1-pi.unre.Incomp))*Model.Incomp$X[,indice.unre.Incomp|N])
    if(rcond(A.Incomp) >= 1e-10){
      B_0.Incomp <- solve(A.Incomp)
      #d_0 should determine which rows to subset
      #gives the first value that is restricted
      d_0.Incomp <- dim(B_0.Incomp)[1]
      #in case the inverse doesn't exist, need another if else
      B.Incomp = B_0.Incomp[d_0.Incomp,d_0.Incomp]
      TW.Incomp <- crossprod(beta.unre.Incomp[N], solve(B.Incomp, beta.unre.Incomp[N]))
    } else {
      TW.Incomp = -1
    }
    
    A.AMS = crossprod(Model.AMS$X[,indice.unre.AMS|N], as.vector(pi.unre.AMS*(1-pi.unre.AMS))*Model.AMS$X[,indice.unre.AMS|N])
    if(rcond(A.AMS) >= 1e-10){
      B_0.AMS <- solve(A.AMS)
      #d_0 should determine which rows to subset
      #gives the first value that is restricted
      d_0.AMS <- dim(B_0.AMS)[1]
      #in case the inverse doesn't exist, need another if else
      B.AMS = B_0.AMS[d_0.AMS,d_0.AMS]
      #so the B_0 needs to be subset to m rows and columns that are restricted based on H0
      TW.AMS <- crossprod(beta.unre.AMS[N], solve(B.AMS, beta.unre.AMS[N]))
    } else {
      TW.AMS = -1
    }
    
    A.BinMM = crossprod(Model.BinMM$X[,indice.unre.BinMM|N], as.vector(pi.unre.BinMM*(1-pi.unre.BinMM))*Model.BinMM$X[,indice.unre.BinMM|N])
    if(rcond(A.BinMM) >= 1e-10){
      B_0.BinMM <- solve(A.BinMM)
      #d_0 should determine which rows to subset
      #gives the first value that is restricted
      d_0.BinMM <- dim(B_0.BinMM)[1]
      #in case the inverse doesn't exist, need another if else
      B.BinMM = B_0.BinMM[d_0.BinMM,d_0.BinMM]
      #so the B_0 needs to be subset to m rows and columns that are restricted based on H0
      TW.BinMM <- crossprod(beta.unre.BinMM[N], solve(B.BinMM, beta.unre.BinMM[N]))
    } else {
      TW.BinMM = -1
    }
    
    # construct the score statistics
    eps.IBS <- Model.IBS$Y-pi.re.IBS #Y - e(Y)
    #this is the X^T times Y-E(Y)
    Xeps.IBS <- crossprod(Model.IBS$X[,indice.re.IBS|N], eps.IBS)
    
    C.IBS = crossprod(Model.IBS$X[,indice.re.IBS|N], as.vector(pi.re.IBS*(1-pi.re.IBS))*Model.IBS$X[,indice.re.IBS|N])
    if(rcond(C.IBS) >= 1e-10){
      TS.IBS <- crossprod(Xeps.IBS, solve(C.IBS, Xeps.IBS))
    } else {
      TS.IBS = -1
    }
    
    eps.Incomp <- Model.Incomp$Y-pi.re.Incomp #Y - e(Y)
    #this is the X^T times Y-E(Y)
    Xeps.Incomp <- crossprod(Model.Incomp$X[,indice.re.Incomp|N], eps.Incomp)
    
    C.Incomp = crossprod(Model.Incomp$X[,indice.re.Incomp|N], as.vector(pi.re.Incomp*(1-pi.re.Incomp))*Model.Incomp$X[,indice.re.Incomp|N])
    if(rcond(C.Incomp) >= 1e-10){
      TS.Incomp <- crossprod(Xeps.Incomp, solve(C.Incomp, Xeps.Incomp))
    } else {
      TS.Incomp = -1
    }
    
    eps.AMS <- Model.AMS$Y-pi.re.AMS #Y - e(Y)
    #this is the X^T times Y-E(Y)
    Xeps.AMS <- crossprod(Model.AMS$X[,indice.re.AMS|N], eps.AMS)
    
    C.AMS = crossprod(Model.AMS$X[,indice.re.AMS|N], as.vector(pi.re.AMS*(1-pi.re.AMS))*Model.AMS$X[,indice.re.AMS|N])
    if(rcond(C.AMS) >= 1e-10){
      TS.AMS <- crossprod(Xeps.AMS, solve(C.AMS, Xeps.AMS))
    } else {
      TS.AMS = -1
    }
    
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
    
    #LRT
    if (TL.IBS>=qchisq(0.95, df = doF)){
      pv.IBS[1] <- pv.IBS[1]+1/5000
    }
    if (TL.Incomp>=qchisq(0.95, df = doF)){
      pv.Incomp[1] <- pv.Incomp[1]+1/5000
    }
    if (TL.AMS>=qchisq(0.95, df = doF)){
      pv.AMS[1] <- pv.AMS[1]+1/5000
    }
    if (TL.BinMM>=qchisq(0.95, df = doF)){
      pv.BinMM[1] <- pv.BinMM[1]+1/5000
    }
    #Wald
    if (TW.IBS>=qchisq(0.95, df = doF)){
      pv.IBS[2] <- pv.IBS[2]+1/5000
    }
    if (TW.Incomp>=qchisq(0.95, df = doF)){
      pv.Incomp[2] <- pv.Incomp[2]+1/5000
    }
    if (TW.AMS>=qchisq(0.95, df = doF)){
      pv.AMS[2] <- pv.AMS[2]+1/5000
    }
    if (TW.BinMM>=qchisq(0.95, df = doF)){
      pv.BinMM[2] <- pv.BinMM[2]+1/5000
    }
    #Score
    if (TS.IBS>=qchisq(0.95, df = doF)){
      pv.IBS[3] <- pv.IBS[3]+1/5000
    }
    if (TS.Incomp>=qchisq(0.95, df = doF)){
      pv.Incomp[3] <- pv.Incomp[3]+1/5000
    }
    if (TS.AMS>=qchisq(0.95, df = doF)){
      pv.AMS[3] <- pv.AMS[3]+1/5000
    }
    if (TS.BinMM>=qchisq(0.95, df = doF)){
      pv.BinMM[3] <- pv.BinMM[3]+1/5000
    }
    
    print(paste0("Simulation ",ii," is complete."))
    
    Tall.IBS[1, ] <- c(TL.IBS, TW.IBS, TS.IBS)
    beta.al.IBS[1, ] <- c(sum(beta.re.IBS!=0), sum(beta.unre.IBS!=0))
    
    Tall.Incomp[1, ] <- c(TL.Incomp, TW.Incomp, TS.Incomp)
    beta.al.Incomp[1, ] <- c(sum(beta.re.Incomp!=0), sum(beta.unre.Incomp!=0))
    
    Tall.AMS[1, ] <- c(TL.AMS, TW.AMS, TS.AMS)
    beta.al.AMS[1, ] <- c(sum(beta.re.AMS!=0), sum(beta.unre.AMS!=0))
    
    Tall.BinMM[1, ] <- c(TL.BinMM, TW.BinMM, TS.BinMM)
    beta.al.BinMM[1, ] <- c(sum(beta.re.BinMM!=0), sum(beta.unre.BinMM!=0))
    
    statsAndPVals = list(pv.IBS=pv.IBS, TScores.IBS=Tall.IBS, beta.IBS=beta.al.IBS,
                         pv.Incomp=pv.Incomp, TScores.Incomp=Tall.Incomp, beta.Incomp=beta.al.Incomp,
                         pv.AMS=pv.AMS, TScores.AMS=Tall.AMS, beta.AMS=beta.al.AMS,
                         pv.BinMM=pv.BinMM, TScores.BinMM=Tall.BinMM, beta.BinMM=beta.al.BinMM)
    statsAndPVals
  }, mc.cores=4)
  
  n = length(statsAndPVals)
  
  pVals.IBS = pVals.Incomp = pVals.AMS = pVals.BinMM = matrix(nrow = n, ncol = 3)
  Scores.IBS = Scores.Incomp = Scores.AMS = Scores.BinMM = matrix(nrow = n, ncol = 3)
  Betas.IBS = Betas.Incomp = Betas.AMS = Betas.BinMM = matrix(nrow = n, ncol = 2)
  
  for(ll in 1:n){
    pVals.IBS[ll,] = unlist(statsAndPVals[[ll]]$pv.IBS)
    Scores.IBS[ll,] = unlist(statsAndPVals[[ll]]$TScores.IBS)
    Betas.IBS[ll,] = unlist(statsAndPVals[[ll]]$beta.IBS)
    
    pVals.Incomp[ll,] = unlist(statsAndPVals[[ll]]$pv.Incomp)
    Scores.Incomp[ll,] = unlist(statsAndPVals[[ll]]$TScores.Incomp)
    Betas.Incomp[ll,] = unlist(statsAndPVals[[ll]]$beta.Incomp)
    
    pVals.AMS[ll,] = unlist(statsAndPVals[[ll]]$pv.AMS)
    Scores.AMS[ll,] = unlist(statsAndPVals[[ll]]$TScores.AMS)
    Betas.AMS[ll,] = unlist(statsAndPVals[[ll]]$beta.AMS)
    
    pVals.BinMM[ll,] = unlist(statsAndPVals[[ll]]$pv.BinMM)
    Scores.BinMM[ll,] = unlist(statsAndPVals[[ll]]$TScores.BinMM)
    Betas.BinMM[ll,] = unlist(statsAndPVals[[ll]]$beta.BinMM)
  }
  
  pVals = cbind(pVals.IBS, pVals.Incomp, pVals.AMS, pVals.BinMM)
  Scores = cbind(Scores.IBS, Scores.Incomp, Scores.AMS, Scores.BinMM)
  Betas = cbind(Betas.IBS, Betas.Incomp, Betas.AMS, Betas.BinMM)
  
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
  write.csv(pVals, file = paste0(path,"/Power_CatY_Prev",YPrev*100,"_LinHypTest_Score",snpOrScore,"_",percentageAssoc,"SNPsAssoc_GammaOnly_OR",ORSize,"_",ld,"_Sim",start,"to",numSims,"_PValues_ForceCovFit.csv"))
  write.csv(Scores, file = paste0(path,"/Power_CatY_Prev",YPrev*100,"_LinHypTest_Score",snpOrScore,"_",percentageAssoc,"SNPsAssoc_GammaOnly_OR",ORSize,"_",ld,"_Sim",start,"to",numSims,"_Stats_ForceCovFit.csv"))
  write.csv(Betas, file = paste0(path,"/Power_CatY_Prev",YPrev*100,"_LinHypTest_Score",snpOrScore,"_",percentageAssoc,"SNPsAssoc_GammaOnly_OR",ORSize,"_",ld,"_Sim",start,"to",numSims,"_Betas_ForceCovFit.csv"))
}

#continuous outcome
RunPowerPipelineLinHypTestCont_GammaTest_Score_BetaZero = function(chr, gene, numPairs, Gamma, TrueScore, ORSize, standardizeScores = FALSE, weightedScores = FALSE, scoreWeights, start, stop, percentageAssoc, LowLD){
  # Function to determine power of linear hyp testing method when score is associated, beta = 0 only testing
  # Only use when Y is continuous
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #Gamma = effect size for score, length 1, can be 0
  #TrueScore = IBS.gene, Incomp.gene, AMS.gene, or BinMM.gene
  #ORSize = Small, Medium, or Large for what OR was used for the associated SNP/score
  #standardizeScores = T or F whether the scores should be standardized based on maximum score value
  #weightedScores = T or F whether the scores will be weighted
  #scoreWeights = m x 1 vector of weights, one weight for each SNP
  # start - simulation number to start at
  # stop - simulation number to stop at
  #percentageAssoc = percentage of SNPs associated with outcome (either 5, 25, 50, 75, or 100) 
  #LowLD = True or FALSE whether the associated SNPs are in low LD or high LD
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes power values to csv files
  
  library(parallel)
  
  #need to define this for naming at the end
  snpOrScore = TrueScore
  
  #always the same
  numSims = stop
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_II_Sims/HapGen_Files/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/ProjectIISourceFunctions_v2.R")
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/Linear_ADMM0.r")
  
  myList = lapply(start:numSims, rep, times = 1)
  # p-value initialization
  pv.IBS <- rep(0, 3)
  pv.Incomp <- rep(0, 3)
  pv.AMS <- rep(0, 3)
  pv.BinMM <- rep(0, 3)
  
  Tall.IBS <- matrix(0, 1, 3)
  Tall.Incomp <- matrix(0, 1, 3)
  Tall.AMS <- matrix(0, 1, 3)
  Tall.BinMM <- matrix(0, 1, 3)
  
  beta.al.IBS <- matrix(0, 1, 2)
  beta.al.Incomp <- matrix(0, 1, 2)
  beta.al.AMS <- matrix(0, 1, 2)
  beta.al.BinMM <- matrix(0, 1, 2)
  
  statsAndPVals = mclapply(myList, function(ii){
    #define matrix to hold all Stats and Pvalues
    statsAndPVals = matrix(NA, nrow = numSims, ncol = 6)
    
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
    ContPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = FALSE,  Covariates = CovData, RGenoData = RGenos, ScoreData = PhenoScore, Betas = Betas, Gamma = Gamma)
    
    # define location of zero components
    # basically, this defines what the null hypothesis is, I think
    # so for joint null, we need the non-zero components to be for covariates
    # and the beta and gamma components to be zero
    numCov = dim(CovData)[2]
    numSNPs = dim(RGenos)[2]
    N = c(rep(FALSE, numCov+1+numSNPs), TRUE)
    
    #need to combine the covariates, r genos, and score matrices together
    #need to include column of 1s for intercept?
    intercept = matrix(1,nrow = numPairs, ncol = 1)
    designMat.IBS = cbind(intercept, CovData, RGenos, IBS.gene)
    designMat.Incomp = cbind(intercept, CovData, RGenos, Incomp.gene)
    designMat.AMS = cbind(intercept, CovData, RGenos, AMS.gene)
    designMat.BinMM = cbind(intercept, CovData, RGenos, BinMM.gene)
    
    #then combine the phenos with the design matrix as a list
    Model.IBS = list(X=designMat.IBS, Y=ContPhenos)
    Model.Incomp = list(X=designMat.Incomp, Y=ContPhenos)
    Model.AMS = list(X=designMat.AMS, Y=ContPhenos)
    Model.BinMM = list(X=designMat.BinMM, Y=ContPhenos)
    
    # estimate the uncontrained estimator
    beta.unre.IBS <- cv.SCAD_ADMM_unre(X=Model.IBS$X, Y=Model.IBS$Y, N=N, beta0=rep(0, dim(Model.IBS$X)[2]), err=1e-4, tune="cv", unpen = c(1,2,3))
    indice.unre.IBS <- beta.unre.IBS!=0
    
    beta.unre.Incomp <- cv.SCAD_ADMM_unre(X=Model.Incomp$X, Y=Model.Incomp$Y, N=N, beta0=rep(0, dim(Model.Incomp$X)[2]), err=1e-4, tune="cv", unpen = c(1,2,3))
    indice.unre.Incomp <- beta.unre.Incomp!=0
    
    beta.unre.AMS <- cv.SCAD_ADMM_unre(X=Model.AMS$X, Y=Model.AMS$Y, N=N, beta0=rep(0, dim(Model.AMS$X)[2]), err=1e-4, tune="cv", unpen = c(1,2,3))
    indice.unre.AMS <- beta.unre.AMS!=0
    
    beta.unre.BinMM <- cv.SCAD_ADMM_unre(X=Model.BinMM$X, Y=Model.BinMM$Y, N=N, beta0=rep(0, dim(Model.BinMM$X)[2]), err=1e-4, tune="cv", unpen = c(1,2,3))
    indice.unre.BinMM <- beta.unre.BinMM!=0
    
    # estimate the constrained estimator
    beta.re.IBS <- cv.SCAD_ADMM_re(X=Model.IBS$X, Y=Model.IBS$Y, N=N, beta0=rep(0, dim(Model.IBS$X)[2]), err=1e-4, tune="cv", unpen = c(1,2,3))
    indice.re.IBS <- beta.re.IBS!=0
    
    beta.re.Incomp <- cv.SCAD_ADMM_re(X=Model.Incomp$X, Y=Model.Incomp$Y, N=N, beta0=rep(0, dim(Model.Incomp$X)[2]), err=1e-4, tune="cv", unpen = c(1,2,3))
    indice.re.Incomp <- beta.re.Incomp!=0
    
    beta.re.AMS <- cv.SCAD_ADMM_re(X=Model.AMS$X, Y=Model.AMS$Y, N=N, beta0=rep(0, dim(Model.AMS$X)[2]), err=1e-4, tune="cv", unpen = c(1,2,3))
    indice.re.AMS <- beta.re.AMS!=0
    
    beta.re.BinMM <- cv.SCAD_ADMM_re(X=Model.BinMM$X, Y=Model.BinMM$Y, N=N, beta0=rep(0, dim(Model.BinMM$X)[2]), err=1e-4, tune="cv", unpen = c(1,2,3))
    indice.re.BinMM <- beta.re.BinMM!=0
    
    # estimate the conditional variance
    #    sig2 <- mean((Model$Y-Model$X%*%beta.unre)^2)
    n = numSims
    sig2.IBS <- mean((Model.IBS$Y-Model.IBS$X%*%beta.unre.IBS)^2)*n/(n-sum(beta.unre.IBS!=0))
    sig2.Incomp <- mean((Model.Incomp$Y-Model.Incomp$X%*%beta.unre.Incomp)^2)*n/(n-sum(beta.unre.Incomp!=0))
    sig2.AMS <- mean((Model.AMS$Y-Model.AMS$X%*%beta.unre.AMS)^2)*n/(n-sum(beta.unre.AMS!=0))
    sig2.BinMM <- mean((Model.BinMM$Y-Model.BinMM$X%*%beta.unre.BinMM)^2)*n/(n-sum(beta.unre.BinMM!=0))
    
    # construct the likelihood ratio statistic
    TL.IBS <- sum((Model.IBS$X%*%beta.re.IBS-Model.IBS$Y)^2)-sum((Model.IBS$X%*%beta.unre.IBS-Model.IBS$Y)^2)
    #if LRT stat is less than 0, there is error, so set to -1
    if(TL.IBS < 0){
      TL.IBS = -1
    }
    TL.Incomp <- sum((Model.Incomp$X%*%beta.re.Incomp-Model.Incomp$Y)^2)-sum((Model.Incomp$X%*%beta.unre.Incomp-Model.Incomp$Y)^2)
    #if LRT stat is less than 0, there is error, so set to -1
    if(TL.Incomp < 0){
      TL.Incomp = -1
    }
    TL.AMS <- sum((Model.AMS$X%*%beta.re.AMS-Model.AMS$Y)^2)-sum((Model.AMS$X%*%beta.unre.AMS-Model.AMS$Y)^2)
    #if LRT stat is less than 0, there is error, so set to -1
    if(TL.AMS < 0){
      TL.AMS = -1
    }
    TL.BinMM <- sum((Model.BinMM$X%*%beta.re.BinMM-Model.BinMM$Y)^2)-sum((Model.BinMM$X%*%beta.unre.BinMM-Model.BinMM$Y)^2)
    #if LRT stat is less than 0, there is error, so set to -1
    if(TL.BinMM < 0){
      TL.BinMM = -1
    }
    
    # construct the Wald statistic
    B.IBS = crossprod(Model.IBS$X[,indice.unre.IBS|N], Model.IBS$X[,indice.unre.IBS|N])
    if(rcond(B.IBS) >= 1e-10){
      B_0.IBS <- solve(B.IBS)
      #d_0.IBS should determine which rows to subset
      #gives the first value that is restricted
      d_0.IBS <- dim(B_0.IBS)[1]
      #in case the inverse doesn't exist, need another if else
      A.IBS = B_0.IBS[d_0.IBS,d_0.IBS]
      TW.IBS <- crossprod(beta.unre.IBS[N], solve(A.IBS, beta.unre.IBS[N]))
    } else {
      TW.IBS = -1
    }
    
    B.Incomp = crossprod(Model.Incomp$X[,indice.unre.Incomp|N], Model.Incomp$X[,indice.unre.Incomp|N])
    if(rcond(B.Incomp) >= 1e-10){
      B_0.Incomp <- solve(B.Incomp)
      #d_0.Incomp should determine which rows to subset
      #gives the first value that is restricted
      d_0.Incomp <- dim(B_0.Incomp)[1]
      #in case the inverse doesn't exist, need another if else
      A.Incomp = B_0.Incomp[d_0.Incomp,d_0.Incomp]
      TW.Incomp <- crossprod(beta.unre.Incomp[N], solve(A.Incomp, beta.unre.Incomp[N]))
    } else {
      TW.Incomp = -1
    }
    
    B.AMS = crossprod(Model.AMS$X[,indice.unre.AMS|N], Model.AMS$X[,indice.unre.AMS|N])
    if(rcond(B.AMS) >= 1e-10){
      B_0.AMS <- solve(B.AMS)
      #d_0.AMS should determine which rows to subset
      #gives the first value that is restricted
      d_0.AMS <- dim(B_0.AMS)[1]
      #in case the inverse doesn't exist, need another if else
      A.AMS = B_0.AMS[d_0.AMS,d_0.AMS]
      TW.AMS <- crossprod(beta.unre.AMS[N], solve(A.AMS, beta.unre.AMS[N]))
    } else {
      TW.AMS = -1
    }
    
    B.BinMM = crossprod(Model.BinMM$X[,indice.unre.BinMM|N], Model.BinMM$X[,indice.unre.BinMM|N])
    if(rcond(B.BinMM) >= 1e-10){
      B_0.BinMM <- solve(B.BinMM)
      #d_0.BinMM should determine which rows to subset
      #gives the first value that is restricted
      d_0.BinMM <- dim(B_0.BinMM)[1]
      #in case the inverse doesn't exist, need another if else
      A.BinMM = B_0.BinMM[d_0.BinMM,d_0.BinMM]
      TW.BinMM <- crossprod(beta.unre.BinMM[N], solve(A.BinMM, beta.unre.BinMM[N]))
    } else {
      TW.BinMM = -1
    }
    
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
    
    #LRT
    if (TL.IBS>=qchisq(0.95, df = doF)){
      pv.IBS[1] <- pv.IBS[1]+1/5000
    }
    if (TL.Incomp>=qchisq(0.95, df = doF)){
      pv.Incomp[1] <- pv.Incomp[1]+1/5000
    }
    if (TL.AMS>=qchisq(0.95, df = doF)){
      pv.AMS[1] <- pv.AMS[1]+1/5000
    }
    if (TL.BinMM>=qchisq(0.95, df = doF)){
      pv.BinMM[1] <- pv.BinMM[1]+1/5000
    }
    #Wald
    if (TW.IBS>=qchisq(0.95, df = doF)){
      pv.IBS[2] <- pv.IBS[2]+1/5000
    }
    if (TW.Incomp>=qchisq(0.95, df = doF)){
      pv.Incomp[2] <- pv.Incomp[2]+1/5000
    }
    if (TW.AMS>=qchisq(0.95, df = doF)){
      pv.AMS[2] <- pv.AMS[2]+1/5000
    }
    if (TW.BinMM>=qchisq(0.95, df = doF)){
      pv.BinMM[2] <- pv.BinMM[2]+1/5000
    }
    #Score
    if (TS.IBS>=qchisq(0.95, df = doF)){
      pv.IBS[3] <- pv.IBS[3]+1/5000
    }
    if (TS.Incomp>=qchisq(0.95, df = doF)){
      pv.Incomp[3] <- pv.Incomp[3]+1/5000
    }
    if (TS.AMS>=qchisq(0.95, df = doF)){
      pv.AMS[3] <- pv.AMS[3]+1/5000
    }
    if (TS.BinMM>=qchisq(0.95, df = doF)){
      pv.BinMM[3] <- pv.BinMM[3]+1/5000
    }
    
    print(paste0("Simulation ",ii," is complete."))
    
    Tall.IBS[1, ] <- c(TL.IBS, TW.IBS, TS.IBS)
    beta.al.IBS[1, ] <- c(sum(beta.re.IBS!=0), sum(beta.unre.IBS!=0))
    
    Tall.Incomp[1, ] <- c(TL.Incomp, TW.Incomp, TS.Incomp)
    beta.al.Incomp[1, ] <- c(sum(beta.re.Incomp!=0), sum(beta.unre.Incomp!=0))
    
    Tall.AMS[1, ] <- c(TL.AMS, TW.AMS, TS.AMS)
    beta.al.AMS[1, ] <- c(sum(beta.re.AMS!=0), sum(beta.unre.AMS!=0))
    
    Tall.BinMM[1, ] <- c(TL.BinMM, TW.BinMM, TS.BinMM)
    beta.al.BinMM[1, ] <- c(sum(beta.re.BinMM!=0), sum(beta.unre.BinMM!=0))
    
    statsAndPVals = list(pv.IBS=pv.IBS, TScores.IBS=Tall.IBS, beta.IBS=beta.al.IBS,
                         pv.Incomp=pv.Incomp, TScores.Incomp=Tall.Incomp, beta.Incomp=beta.al.Incomp,
                         pv.AMS=pv.AMS, TScores.AMS=Tall.AMS, beta.AMS=beta.al.AMS,
                         pv.BinMM=pv.BinMM, TScores.BinMM=Tall.BinMM, beta.BinMM=beta.al.BinMM)
    statsAndPVals
  }, mc.cores=4)
  
  n = length(statsAndPVals)
  
  pVals.IBS = pVals.Incomp = pVals.AMS = pVals.BinMM = matrix(nrow = n, ncol = 3)
  Scores.IBS = Scores.Incomp = Scores.AMS = Scores.BinMM = matrix(nrow = n, ncol = 3)
  Betas.IBS = Betas.Incomp = Betas.AMS = Betas.BinMM = matrix(nrow = n, ncol = 2)
  
  for(ll in 1:n){
    pVals.IBS[ll,] = unlist(statsAndPVals[[ll]]$pv.IBS)
    Scores.IBS[ll,] = unlist(statsAndPVals[[ll]]$TScores.IBS)
    Betas.IBS[ll,] = unlist(statsAndPVals[[ll]]$beta.IBS)
    
    pVals.Incomp[ll,] = unlist(statsAndPVals[[ll]]$pv.Incomp)
    Scores.Incomp[ll,] = unlist(statsAndPVals[[ll]]$TScores.Incomp)
    Betas.Incomp[ll,] = unlist(statsAndPVals[[ll]]$beta.Incomp)
    
    pVals.AMS[ll,] = unlist(statsAndPVals[[ll]]$pv.AMS)
    Scores.AMS[ll,] = unlist(statsAndPVals[[ll]]$TScores.AMS)
    Betas.AMS[ll,] = unlist(statsAndPVals[[ll]]$beta.AMS)
    
    pVals.BinMM[ll,] = unlist(statsAndPVals[[ll]]$pv.BinMM)
    Scores.BinMM[ll,] = unlist(statsAndPVals[[ll]]$TScores.BinMM)
    Betas.BinMM[ll,] = unlist(statsAndPVals[[ll]]$beta.BinMM)
  }
  
  pVals = cbind(pVals.IBS, pVals.Incomp, pVals.AMS, pVals.BinMM)
  Scores = cbind(Scores.IBS, Scores.Incomp, Scores.AMS, Scores.BinMM)
  Betas = cbind(Betas.IBS, Betas.Incomp, Betas.AMS, Betas.BinMM)
  
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
  
  write.csv(pVals, file = paste0(path,"/Power_ContY_LinHypTest_Score",snpOrScore,"_GammaOnly_",percentageAssoc,"SNPsAssoc_JointTesting_OR",ORSize,"_",ld,"_Sim",start,"to",numSims,"_PValues_ForceCovs.csv"))
  write.csv(Scores, file = paste0(path,"/Power_ContY_LinHypTest_Score",snpOrScore,"_GammaOnly_",percentageAssoc,"SNPsAssoc_JointTesting_OR",ORSize,"_",ld,"_Sim",start,"to",numSims,"_Stats_ForceCovs.csv"))
  write.csv(Betas, file = paste0(path,"/Power_ContY_LinHypTest_Score",snpOrScore,"_GammaOnly_",percentageAssoc,"SNPsAssoc_JointTesting_OR",ORSize,"_",ld,"_Sim",start,"to",numSims,"_Betas_ForceCovs.csv"))
}

#################### Case 2: Testing Gamma = 0 when Beta not = 0

#categorical outcome
RunTIEPipelineLinHypTestCat_GammaTest_GLM_BetaSig = function(chr, gene, numPairs, YPrev, standardizeScores = FALSE, weightedScores = FALSE, scoreWeights, ORSize, percentageAssoc, LowLD, score, start, stop){
  #function to run whole TIE pipeline, Calculates Score test stat for Lin Hyp test method and whether p-value <0.05
  # Only use when Y is binary
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #standardizeScores = T or F whether the scores should be standardized based on maximum score value
  #weightedScores = T or F whether the scores will be weighted
  #scoreWeights = m x 1 vector of weights, one weight for each SNP
  # score - which score we are using in the main JST model
  # start - simulation number to start at
  # stop - simulation number to stop at
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes TIE values to csv files
  
  library(parallel)
  suppressMessages(library(epicalc))
  
  #always the same
  numSims = stop  
  
  #define effect based on OR size
  if(ORSize == "small"){
    effect = 0.14
  } else if(ORSize == "medium"){
    effect = 0.41
  } else {
    effect = 0.69
  }
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_II_Sims/HapGen_Files/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/ProjectIISourceFunctions_v2.R")
  
  #determine which SNPs to actually set as assoc. based on gene
  assocSNPs = DetermineAssocRSNPs(gene = gene, LowLD = LowLD, percentageAssoc = percentageAssoc)
  
  myList = lapply(start:numSims, rep, times = 1)
  nullValues =  altValues = statsAndPValsmat = finalOutput = list()
  
  statsAndPVals = mclapply(myList, function(ii){
    #pull recipient and donor genotypes
    RGenos = obtainRGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    DGenos = obtainDGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    
    #calculate single snp scores
    if(score == "IBS"){
      Score.snp = calcIBSMismatch(RGenosMat = RGenos, DGenosMat = DGenos)
    } else if(score == "Incomp"){
      Score.snp = calcIncompatibilityScore(RGenosMat = RGenos, DGenosMat = DGenos)
    } else if(score == "AMS"){
      Score.snp = calcAMS(RGenosMat = RGenos, DGenosMat = DGenos)
    } else{
      Score.snp = calcBinaryMM(RGenosMat = RGenos, DGenosMat = DGenos)
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
    
    #define null gamma
    Gamma = c(0)
    
    #generate phenotypes, both continuous and binary
    # CatPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = TRUE, YPrev = YPrev,  Covariates = CovData, RGenoData = RGenos, ScoreData = Score.gene, Betas = Betas, Gamma = Gamma)
    CatPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = FALSE, YCat = TRUE, YPrev = YPrev, RGenoData = RGenos, ScoreData = Score.gene, Betas = Betas, Gamma = Gamma)
    
    # allData = cbind(CovData, RGenos, Score.gene, CatPhenos)
    allData = cbind(RGenos, Score.gene, CatPhenos)
    
    #fit the null and alternative models
    fitNull = glm(CatPhenos~RGenos, family = binomial, epsilon = 1e-6)
    fitAlt = glm(CatPhenos~RGenos+Score.gene, family = binomial, epsilon = 1e-6)	
    
    nullValues = summary(fitNull)$coefficients
    altValues = summary(fitAlt)$coefficients
    
    if(fitNull$deviance - fitAlt$deviance >= 0){
      # construct the likelihood ratio statistics
      TL = lrtest(fitNull, fitAlt)
      
      pv = TL$p.value
      stat = TL$Chisquared
      doF = TL$df
    } else {
      pv = 1
      stat = 0
      doF = (dim(RGenos)[2]+1)
    }
    
    print(paste0("Simulation ",ii," is complete."))
    
    statsAndPValsmat = c(pv, stat)
    
    finalOutput[[ii]] = list(statsAndPValsmat, nullValues, altValues)
  }, mc.cores=4)
  
  statsAndPValsAll = list()
  nulValuesAll = list()
  altValuesAll = list()
  
  for(ii in 1:length(statsAndPVals)){
    statsAndPValsAll[[ii]] = statsAndPVals[[ii]][1]
    nulValuesAll[[ii]] =  statsAndPVals[[ii]][2]
    altValuesAll[[ii]] = statsAndPVals[[ii]][3]
  }
  
  listLength = length(statsAndPValsAll)  
  statsAndPValsAll.mat = matrix(unlist(statsAndPValsAll),nrow = listLength, ncol = 2, byrow = TRUE)
  dim1 = dim(nulValuesAll[[1]][[1]])
  dim2 = dim(altValuesAll[[1]][[1]])
  rownames1 = c()
  rownames2 = c()
  for(jj in 1:length(statsAndPVals)){
    rownames1 = c(rownames1, rownames(nulValuesAll[[jj]][[1]]))
    rownames2 = c(rownames2, rownames(altValuesAll[[jj]][[1]]))
  }
  nulValuesAll.mat = matrix(unlist(nulValuesAll), ncol = dim1[2], byrow = TRUE, dimnames = list(rownames1, colnames(nulValuesAll[[1]][[1]])))
  altValuesAll.mat = matrix(unlist(altValuesAll), ncol = dim2[2], byrow = TRUE, dimnames = list(rownames2, colnames(altValuesAll[[1]][[1]])))
  
  #rename 
  statsAndPValsCatPhenos = statsAndPValsAll.mat
  
  #write out the Stats and p values, and Summary stats
  if(weightedScores == FALSE){
    if(standardizeScores == FALSE){
      write.csv(statsAndPValsCatPhenos, file = paste0(path,"/TIE_CatPhenos_Prev",YPrev*100,"_GLM_BetaSig_Score",score,"_Sim",start,"to",numSims,"_StatsAndPValues.csv"))
      write.csv(nulValuesAll.mat, file = paste0(path,"/TIE_CatNullBetasSummary_Prev",YPrev*100,"_GLM_BetaSig_Score",score,"_Sim",start,"to",numSims,"_StatsAndPValues.csv"))
      write.csv(altValuesAll.mat, file = paste0(path,"/TIE_CatAltBetasSummary_Prev",YPrev*100,"_GLM_BetaSig_Score",score,"_Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    } else { #scores are standardized
      write.csv(statsAndPValsCatPhenos, file = paste0(path,"/TIE_CatPhenos_Prev",YPrev*100,"_GLM_BetaSig_Score",score,"_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
      write.csv(nulValuesAll.mat, file = paste0(path,"/TIE_CatBetasSummary_Prev",YPrev*100,"_GLM_BetaSig_Score",score,"_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
      write.csv(altValuesAll.mat, file = paste0(path,"/TIE_CatAltBetasSummary_Prev",YPrev*100,"_GLM_BetaSig_Score",score,"_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
    }
  } else { #scores are weighted
    write.csv(statsAndPValsCatPhenos, file = paste0(path,"/TIE_CatPhenos_Prev",YPrev*100,"_GLM_BetaSig_Score",score,"_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
    write.csv(nulValuesAll.mat, file = paste0(path,"/TIE_CatBetasSummary_Prev",YPrev*100,"_GLM_BetaSig_Score",score,"_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
    write.csv(altValuesAll.mat, file = paste0(path,"/TIE_CatAltBetasSummary_Prev",YPrev*100,"_GLM_BetaSig_Score",score,"_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
  }
}

#continuous outcome
RunTIEPipelineLinHypTestCont_GammaTest_GLM_BetaSig = function(chr, gene, numPairs, standardizeScores = FALSE, weightedScores = FALSE, scoreWeights, ORSize, percentageAssoc, LowLD, score, start, stop){
  #function to run whole TIE pipeline, Calculates LRT test stat for Lin Hyp test method and whether p-value <0.05
  # only use when Y is continuous (Normal)
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #standardizeScores = T or F whether the scores should be standardized based on maximum score value
  #weightedScores = T or F whether the scores will be weighted
  #scoreWeights = m x 1 vector of weights, one weight for each SNP
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes TIE values to csv files
  
  library(parallel)
  suppressMessages(library(epicalc))
  
  #always the same
  numSims = stop 
  
  #define effect based on OR size
  if(ORSize == "small"){
    effect = 0.14
  } else if(ORSize == "medium"){
    effect = 0.41
  } else {
    effect = 0.69
  }
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_II_Sims/HapGen_Files/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/ProjectIISourceFunctions_v2.R")
  
  #determine which SNPs to actually set as assoc. based on gene
  #check gene
  if(gene == "NAT2"){
    if(LowLD == TRUE){
      #then looking at SNPs in low LD
      if(percentageAssoc == 5){
        #then 1 SNP is associated
        assocSNPs = c(1)
      } else if(percentageAssoc == 15){
        #then 2 SNPs are associated
        assocSNPs = c(1,2)
      } else {
        #then 5 SNPs are associated
        assocSNPs = c(1,2,9,7,10)
      } 
    } else {
      #looking at SNPs in high LD 
      if(percentageAssoc == 5){
        #then 1 SNP is associated
        assocSNPs = c(13)
      } else if(percentageAssoc == 15){
        #then 2 SNPs are associated
        assocSNPs = c(13,14)
      } else {
        #then 5 SNPs are associated
        assocSNPs = c(13,14,4,5,3)
      } 
    }
  } else if(gene == "CHI3L2"){
    if(LowLD == TRUE){
      #then looking at SNPs in low LD
      if(percentageAssoc == 5){
        #then 2 SNPs are associated
        assocSNPs = c(12,2)
      } else if(percentageAssoc == 15){
        #then 5 SNPs are associated
        assocSNPs = c(12,2,3,9,8)
      } else {
        #then 8 SNPs are associated
        assocSNPs = c(12,2,3,9,8,24,28,19)
      } 
    } else {
      #looking at SNPs in high LD 
      #then looking at SNPs in low LD
      if(percentageAssoc == 5){
        #then 2 SNPs are associated
        assocSNPs = c(17,29)
      } else if(percentageAssoc == 15){
        #then 5 SNPs are associated
        assocSNPs = c(5,6,23,18,13)
      } else {
        #then 8 SNPs are associated
        assocSNPs = c(5,6,23,18,13,15,17,29)
      }
    }
  } else { 
    #otherwise we have ASAH1
    if(LowLD == TRUE){
      #then looking at SNPs in low LD
      if(percentageAssoc == 5){
        #then 2 SNPs are associated
        assocSNPs = c(31,9)
      } else if(percentageAssoc == 15){
        #then 6 SNPs are associated
        assocSNPs = c(31,9,4,8,10,5)
      } else if(percentageAssoc == 25){
        #then 10 SNPs are associated
        assocSNPs = c(31,9,4,8,10,5,3,2,7,1)
      }
    } else {
      #looking at SNPs in high LD 
      #then looking at SNPs in low LD
      if(percentageAssoc == 5){
        #then 2 SNPs are associated
        assocSNPs = c(25, 32)
      } else if(percentageAssoc == 15){
        #then 6 SNPs are associated
        assocSNPs = c(37,29,36,28,30,33)
      } else if(percentageAssoc == 25){
        #then 10 SNPs are associated
        assocSNPs = c(37,29,36,28,30,33,34,38,25,32)
      } 
    }
  }
  
  myList = lapply(start:numSims, rep, times = 1)
  nullValues =  altValues = statsAndPValsmat = finalOutput = list()
  
  statsAndPVals = mclapply(myList, function(ii){
    #pull recipient and donor genotypes
    RGenos = obtainRGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    DGenos = obtainDGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    
    #calculate single snp scores
    if(score == "IBS"){
      Score.snp = calcIBSMismatch(RGenosMat = RGenos, DGenosMat = DGenos)
    } else if(score == "Incomp"){
      Score.snp = calcIncompatibilityScore(RGenosMat = RGenos, DGenosMat = DGenos)
    } else if(score == "AMS"){
      Score.snp = calcAMS(RGenosMat = RGenos, DGenosMat = DGenos)
    } else{
      Score.snp = calcBinaryMM(RGenosMat = RGenos, DGenosMat = DGenos)
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
    
    #generate covariates
    #for now, a single binary and a single continous covariate
    # CovData = GenCovData(SampleSize = numPairs, BinaryValues = 1, ContinuousValues = 1)
    
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
    
    #define null gamma
    Gamma = c(0)
    
    #generate phenotypes, both continuous and binary
    ContPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = FALSE, YCat = FALSE,  RGenoData = RGenos, ScoreData = Score.gene, Betas = Betas, Gamma = Gamma)
    # ContPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = FALSE,  Covariates = CovData, RGenoData = RGenos, ScoreData = score, Betas = Betas, Gamma = Gamma)
    
    #fit the null and alternative models
    fitNull = glm(ContPhenos~RGenos, family = gaussian)
    fitAlt = glm(ContPhenos~RGenos+Score.gene, family = gaussian, epsilon = 1e-6)	
    
    nullValues = summary(fitNull)$coefficients
    altValues = summary(fitAlt)$coefficients
    
    if(fitNull$deviance - fitAlt$deviance >= 0){
      # construct the likelihood ratio statistics
      TL = lrtest(fitNull, fitAlt)
      
      pv = TL$p.value
      stat = TL$Chisquared
      doF = TL$df
    } else {
      pv = 1
      stat = 0
      doF = (dim(RGenos)[2]+1)
    }
    
    print(paste0("Simulation ",ii," is complete."))
    
    statsAndPValsmat = c(pv, stat)
    
    finalOutput[[ii]] = list(statsAndPValsmat, nullValues, altValues)
  }, mc.cores=1)
  
  statsAndPValsAll = list()
  nulValuesAll = list()
  altValuesAll = list()
  
  for(ii in 1:length(statsAndPVals)){
    statsAndPValsAll[[ii]] = statsAndPVals[[ii]][1]
    nulValuesAll[[ii]] =  statsAndPVals[[ii]][2]
    altValuesAll[[ii]] = statsAndPVals[[ii]][3]
  }
  
  listLength = length(statsAndPValsAll)  
  statsAndPValsAll.mat = matrix(unlist(statsAndPValsAll),nrow = listLength, ncol = 2, byrow = TRUE)
  dim1 = dim(nulValuesAll[[1]][[1]])
  dim2 = dim(altValuesAll[[1]][[1]])
  rownames1 = c()
  rownames2 = c()
  for(jj in 1:length(statsAndPVals)){
    rownames1 = c(rownames1, rownames(nulValuesAll[[jj]][[1]]))
    rownames2 = c(rownames2, rownames(altValuesAll[[jj]][[1]]))
  }
  nulValuesAll.mat = matrix(unlist(nulValuesAll), ncol = dim1[2], byrow = TRUE, dimnames = list(rownames1, colnames(nulValuesAll[[1]][[1]])))
  altValuesAll.mat = matrix(unlist(altValuesAll), ncol = dim2[2], byrow = TRUE, dimnames = list(rownames2, colnames(altValuesAll[[1]][[1]])))
  
  #separate into cat and cont phenos
  statsAndPValsContPhenos = statsAndPValsAll.mat
  
  #write out the Stats and p values
  if(weightedScores == FALSE){
    if(standardizeScores == FALSE){
      write.csv(statsAndPValsContPhenos, file = paste0(path,"/TIE_ContPhenos_LinHypTest_BetaSig_Score",score,"_Sim",start,"to",numSims,"_StatsAndPValues.csv"))
      write.csv(nulValuesAll.mat, file = paste0(path,"/TIE_ContNullBetasSummary_LinHypTest_BetaSig_Score",score,"_Sim",start,"to",numSims,"_StatsAndPValues.csv"))
      write.csv(altValuesAll.mat, file = paste0(path,"/TIE_ContAltBetasSummary_LinHypTest_BetaSig_Score",score,"_Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    } else { #scores are standardized
      write.csv(statsAndPValsContPhenos, file = paste0(path,"/TIE_ContPhenos_LinHypTest_BetaSig_Score",score,"_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
      write.csv(nulValuesAll.mat, file = paste0(path,"/TIE_ContNullBetasSummary_LinHypTest_BetaSig_Score",score,"_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
      write.csv(altValuesAll.mat, file = paste0(path,"/TIE_ContAltBetasSummary_LinHypTest_BetaSig_Score",score,"_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
    }
  } else { #scores are weighted
    write.csv(statsAndPValsContPhenos, file = paste0(path,"/TIE_ContPhenos_LinHypTest_BetaSig_Score",score,"_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
    write.csv(nulValuesAll.mat, file = paste0(path,"/TIE_ContNullBetasSummary_LinHypTest_BetaSig_Score",score,"_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
    write.csv(altValuesAll.mat, file = paste0(path,"/TIE_ContAltBetasSummary_LinHypTest_BetaSig_Score",score,"_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
  }
}

#categorical outcome
RunTIEPipelineLinHypTestCat_GammaTest_Lasso_BetaSig = function(chr, gene, numPairs, YPrev, standardizeScores = FALSE, weightedScores = FALSE, scoreWeights, ORSize, percentageAssoc, LowLD, score, start, stop){
  #function to run whole TIE pipeline, Calculates Score test stat for Lin Hyp test method and whether p-value <0.05
  # Only use when Y is binary
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #standardizeScores = T or F whether the scores should be standardized based on maximum score value
  #weightedScores = T or F whether the scores will be weighted
  #scoreWeights = m x 1 vector of weights, one weight for each SNP
  # score - which score we are using in the main JST model
  # start - simulation number to start at
  # stop - simulation number to stop at
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes TIE values to csv files
  
  library(parallel)
  suppressMessages(library(epicalc))
  suppressMessages(library(glmnet))
  
  #always the same
  numSims = stop  
  
  #define effect based on OR size
  if(ORSize == "small"){
    effect = 0.14
  } else if(ORSize == "medium"){
    effect = 0.41
  } else {
    effect = 0.69
  }
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_II_Sims/HapGen_Files/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/ProjectIISourceFunctions_v2.R")
  
  #determine which SNPs to actually set as assoc. based on gene
  #check gene
  if(gene == "NAT2"){
    if(LowLD == TRUE){
      #then looking at SNPs in low LD
      if(percentageAssoc == 5){
        #then 1 SNP is associated
        assocSNPs = c(1)
      } else if(percentageAssoc == 15){
        #then 2 SNPs are associated
        assocSNPs = c(1,2)
      } else {
        #then 5 SNPs are associated
        assocSNPs = c(1,2,9,7,10)
      } 
    } else {
      #looking at SNPs in high LD 
      if(percentageAssoc == 5){
        #then 1 SNP is associated
        assocSNPs = c(13)
      } else if(percentageAssoc == 15){
        #then 2 SNPs are associated
        assocSNPs = c(13,14)
      } else {
        #then 5 SNPs are associated
        assocSNPs = c(13,14,4,5,3)
      } 
    }
  } else if(gene == "CHI3L2"){
    if(LowLD == TRUE){
      #then looking at SNPs in low LD
      if(percentageAssoc == 5){
        #then 2 SNPs are associated
        assocSNPs = c(12,2)
      } else if(percentageAssoc == 15){
        #then 5 SNPs are associated
        assocSNPs = c(12,2,3,9,8)
      } else {
        #then 8 SNPs are associated
        assocSNPs = c(12,2,3,9,8,24,28,19)
      } 
    } else {
      #looking at SNPs in high LD 
      #then looking at SNPs in low LD
      if(percentageAssoc == 5){
        #then 2 SNPs are associated
        assocSNPs = c(17,29)
      } else if(percentageAssoc == 15){
        #then 5 SNPs are associated
        assocSNPs = c(5,6,23,18,13)
      } else {
        #then 8 SNPs are associated
        assocSNPs = c(5,6,23,18,13,15,17,29)
      }
    }
  } else { 
    #otherwise we have ASAH1
    if(LowLD == TRUE){
      #then looking at SNPs in low LD
      if(percentageAssoc == 5){
        #then 2 SNPs are associated
        assocSNPs = c(31,9)
      } else if(percentageAssoc == 15){
        #then 6 SNPs are associated
        assocSNPs = c(31,9,4,8,10,5)
      } else if(percentageAssoc == 25){
        #then 10 SNPs are associated
        assocSNPs = c(31,9,4,8,10,5,3,2,7,1)
      }
    } else {
      #looking at SNPs in high LD 
      #then looking at SNPs in low LD
      if(percentageAssoc == 5){
        #then 2 SNPs are associated
        assocSNPs = c(25, 32)
      } else if(percentageAssoc == 15){
        #then 6 SNPs are associated
        assocSNPs = c(37,29,36,28,30,33)
      } else if(percentageAssoc == 25){
        #then 10 SNPs are associated
        assocSNPs = c(37,29,36,28,30,33,34,38,25,32)
      } 
    }
  }
  
  myList = lapply(start:numSims, rep, times = 1)
  nullValues =  altValues = statsAndPValsmat = finalOutput = list()
  
  statsAndPVals = mclapply(myList, function(ii){
    #pull recipient and donor genotypes
    RGenos = obtainRGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    DGenos = obtainDGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    
    #calculate single snp scores
    if(score == "IBS"){
      Score.snp = calcIBSMismatch(RGenosMat = RGenos, DGenosMat = DGenos)
    } else if(score == "Incomp"){
      Score.snp = calcIncompatibilityScore(RGenosMat = RGenos, DGenosMat = DGenos)
    } else if(score == "AMS"){
      Score.snp = calcAMS(RGenosMat = RGenos, DGenosMat = DGenos)
    } else{
      Score.snp = calcBinaryMM(RGenosMat = RGenos, DGenosMat = DGenos)
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
    
    #generate covariates
    #for now, a single binary and a single continous covariate
    # CovData = GenCovData(SampleSize = numPairs, BinaryValues = 1, ContinuousValues = 1)
    
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
    
    #define null gamma
    Gamma = c(0)
    
    #generate phenotypes, both continuous and binary
    # CatPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = TRUE, YPrev = YPrev,  Covariates = CovData, RGenoData = RGenos, ScoreData = Score.gene, Betas = Betas, Gamma = Gamma)
    CatPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = FALSE, YCat = TRUE, YPrev = YPrev, RGenoData = RGenos, ScoreData = Score.gene, Betas = Betas, Gamma = Gamma)
    
    #combine the data for x design matrix
    allData = cbind(RGenos, Score.gene)
    
    #split into training and test data
    #train with 80% of data
    train.null.index = sample(1:nrow(RGenos), 5*(nrow(RGenos)/10))
    train.all.index = sample(1:nrow(allData), 5*(nrow(allData)/10))
    #pull the index values for the test data (20%)
    test.null.index = -train.null.index
    test.all.index = -train.all.index
    #separate the data
    x.train.null = RGenos[train.null.index,]
    x.train.all = allData[train.all.index,]
    x.test.null = RGenos[test.null.index,]
    x.test.all = allData[test.all.index,]
    
    #convert y phenos outcome to numerical variable
    CatPhenos.fac = ifelse(CatPhenos == 1,1,0)
    
    #pull the y phenos for the 
    y.train.null = CatPhenos.fac[train.null.index]
    y.train.all = CatPhenos.fac[train.all.index]
    y.test.null = CatPhenos.fac[test.null.index]
    y.test.all = CatPhenos.fac[test.all.index]
    
    #use cv to find best lambda value
    cv_output_null = cv.glmnet(x.train.null, y.train.null, family = "binomial", alpha=1, intercept=FALSE)
    cv_output_alt = cv.glmnet(x.train.all, y.train.all, family = "binomial", alpha=1, intercept=FALSE)
    cv_output_null.int = cv.glmnet(x.train.null, y.train.null, family = "binomial", alpha=1, intercept=TRUE)
    cv_output_alt.int = cv.glmnet(x.train.all, y.train.all, family = "binomial", alpha=1, intercept=TRUE)
    
    #identifying the best lambda
    best_lambda_null = cv_output_null$lambda.min
    best_lambda_alt = cv_output_alt$lambda.min
    best_lambda_null.int = cv_output_null.int$lambda.min
    best_lambda_alt.int = cv_output_alt.int$lambda.min
    
    #fit the null and alternative models with the best lambda values
    fitNull = glmnet(RGenos, CatPhenos.fac, family = "binomial", alpha=1, lambda = best_lambda_null, intercept=FALSE)
    fitAlt = glmnet(allData, CatPhenos.fac, family = "binomial", alpha=1, lambda = best_lambda_alt, intercept=FALSE)	
    fitNull.int = glmnet(RGenos, CatPhenos.fac, family = "binomial", alpha=1, lambda = best_lambda_null.int, intercept=TRUE)
    fitAlt.int = glmnet(allData, CatPhenos.fac, family = "binomial", alpha=1, lambda = best_lambda_alt.int, intercept=TRUE)	
    
    nullValues = coef(fitNull)
    altValues = coef(fitAlt)
    nullValues.int = coef(fitNull.int)
    altValues.int = coef(fitAlt.int)
    
    nullValues.lambda = fitNull$lambda
    altValues.lambda = fitAlt$lambda
    nullValues.int.lambda = fitNull.int$lambda
    altValues.int.lambda = fitAlt.int$lambda
    
    lambdas = c(nullValues.lambda, altValues.lambda, nullValues.int.lambda, altValues.int.lambda)
    
    print(paste0("Simulation ",ii," is complete."))
    
    finalOutput[[ii]] = list(nullValues, altValues, nullValues.int, altValues.int, lambdas)
  }, mc.cores=4)
  
  betasNull = list()
  betasAll = list()
  betasNull.int = list()
  betasAll.int = list()
  lambdas = list()
  
  for(ii in 1:length(statsAndPVals)){
    betasNull[[ii]] = statsAndPVals[[ii]][1]
    betasAll[[ii]] =  statsAndPVals[[ii]][2]
    betasNull.int[[ii]] = statsAndPVals[[ii]][3]
    betasAll.int[[ii]] = statsAndPVals[[ii]][4]
    lambdas[[ii]] = statsAndPVals[[ii]][5]
  }
  
  betasNull.mat =  do.call(rbind, unlist(betasNull, recursive = FALSE))
  betasAll.mat = do.call(rbind, unlist(betasAll, recursive = FALSE))
  betasNull.int.mat = do.call(rbind, unlist(betasNull.int, recursive = FALSE))
  betasAll.int.mat = do.call(rbind, unlist(betasAll.int, recursive = FALSE))
  
  betasNullmat = as.matrix(betasNull.mat)
  betasAllmat = as.matrix(betasAll.mat)
  betasNull.intmat = as.matrix(betasNull.int.mat)
  betasAll.intmat = as.matrix(betasAll.int.mat)
  lambdas.mat = matrix(unlist(lambdas), ncol = 4, byrow = TRUE)
  
  #write out the Stats and p values, and Summary stats
  if(weightedScores == FALSE){
    if(standardizeScores == FALSE){
      write.csv(betasNullmat, file = paste0(path,"/TIE_CatNullBetas_Prev",YPrev*100,"_Lasso_BetaSig_Score",score,"_Sim",start,"to",numSims,"_StatsAndPValues.csv"))
      write.csv(betasAllmat, file = paste0(path,"/TIE_CatAltBetas_Prev",YPrev*100,"_Lasso_BetaSig_Score",score,"_Sim",start,"to",numSims,"_StatsAndPValues.csv"))
      write.csv(betasNull.intmat, file = paste0(path,"/TIE_CatNullBetasIntercept_Prev",YPrev*100,"_Lasso_BetaSig_Score",score,"_Sim",start,"to",numSims,"_StatsAndPValues.csv"))
      write.csv(betasAll.intmat, file = paste0(path,"/TIE_CatAltBetasIntercept_Prev",YPrev*100,"_Lasso_BetaSig_Score",score,"_Sim",start,"to",numSims,"_StatsAndPValues.csv"))
      write.csv(lambdas.mat, file = paste0(path,"/TIE_CatLambdas_Prev",YPrev*100,"_Lasso_BetaSig_Score",score,"_Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    } else { #scores are standardized
      write.csv(betasNullmat, file = paste0(path,"/TIE_CatNullBetas_Prev",YPrev*100,"_BetaSig_Score",score,"_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
      write.csv(betasAllmat, file = paste0(path,"/TIE_CatAltBetas_Prev",YPrev*100,"_BetaSig_Score",score,"_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
      write.csv(betasNull.intmat, file = paste0(path,"/TIE_CatNullBetasIntercept_Prev",YPrev*100,"_BetaSig_Score",score,"_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
      write.csv(betasAll.intmat, file = paste0(path,"/TIE_CatAltBetasIntercept_Prev",YPrev*100,"_BetaSig_Score",score,"_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
      write.csv(lambdas.mat, file = paste0(path,"/TIE_CatLambdas_Prev",YPrev*100,"_BetaSig_Score",score,"_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
    }
  } else { #scores are weighted
    write.csv(betasNullmat, file = paste0(path,"/TIE_CatNullBetas_Prev",YPrev*100,"_BetaSig_Score",score,"_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
    write.csv(betasAllmat, file = paste0(path,"/TIE_CatAltBetas_Prev",YPrev*100,"_BetaSig_Score",score,"_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
    write.csv(betasNull.intmat, file = paste0(path,"/TIE_CatNullBetasIntercept_Prev",YPrev*100,"_BetaSig_Score",score,"_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
    write.csv(betasAll.intmat, file = paste0(path,"/TIE_CatAltBetasIntercept_Prev",YPrev*100,"_BetaSig_Score",score,"_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
    write.csv(lambdas.mat, file = paste0(path,"/TIE_CatLambdas_Prev",YPrev*100,"_BetaSig_Score",score,"_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
  }
}

#continuous outcome
RunTIEPipelineLinHypTestCont_GammaTest_Lasso_BetaSig = function(chr, gene, numPairs, standardizeScores = FALSE, weightedScores = FALSE, scoreWeights, ORSize, percentageAssoc, LowLD, score, start, stop){
  #function to run whole TIE pipeline, Calculates LRT test stat for Lin Hyp test method and whether p-value <0.05
  # only use when Y is continuous (Normal)
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #standardizeScores = T or F whether the scores should be standardized based on maximum score value
  #weightedScores = T or F whether the scores will be weighted
  #scoreWeights = m x 1 vector of weights, one weight for each SNP
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes TIE values to csv files
  
  library(parallel)
  suppressMessages(library(epicalc))
  suppressMessages(library(glmnet))

  #always the same
  numSims = stop 
  
  #define effect based on OR size
  if(ORSize == "small"){
    effect = 0.14
  } else if(ORSize == "medium"){
    effect = 0.41
  } else {
    effect = 0.69
  }
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_II_Sims/HapGen_Files/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/ProjectIISourceFunctions_v2.R")
  
  #determine which SNPs to actually set as assoc. based on gene
  #check gene
  if(gene == "NAT2"){
    if(LowLD == TRUE){
      #then looking at SNPs in low LD
      if(percentageAssoc == 5){
        #then 1 SNP is associated
        assocSNPs = c(1)
      } else if(percentageAssoc == 15){
        #then 2 SNPs are associated
        assocSNPs = c(1,2)
      } else {
        #then 5 SNPs are associated
        assocSNPs = c(1,2,9,7,10)
      } 
    } else {
      #looking at SNPs in high LD 
      if(percentageAssoc == 5){
        #then 1 SNP is associated
        assocSNPs = c(13)
      } else if(percentageAssoc == 15){
        #then 2 SNPs are associated
        assocSNPs = c(13,14)
      } else {
        #then 5 SNPs are associated
        assocSNPs = c(13,14,4,5,3)
      } 
    }
  } else if(gene == "CHI3L2"){
    if(LowLD == TRUE){
      #then looking at SNPs in low LD
      if(percentageAssoc == 5){
        #then 2 SNPs are associated
        assocSNPs = c(12,2)
      } else if(percentageAssoc == 15){
        #then 5 SNPs are associated
        assocSNPs = c(12,2,3,9,8)
      } else {
        #then 8 SNPs are associated
        assocSNPs = c(12,2,3,9,8,24,28,19)
      } 
    } else {
      #looking at SNPs in high LD 
      #then looking at SNPs in low LD
      if(percentageAssoc == 5){
        #then 2 SNPs are associated
        assocSNPs = c(17,29)
      } else if(percentageAssoc == 15){
        #then 5 SNPs are associated
        assocSNPs = c(5,6,23,18,13)
      } else {
        #then 8 SNPs are associated
        assocSNPs = c(5,6,23,18,13,15,17,29)
      }
    }
  } else { 
    #otherwise we have ASAH1
    if(LowLD == TRUE){
      #then looking at SNPs in low LD
      if(percentageAssoc == 5){
        #then 2 SNPs are associated
        assocSNPs = c(31,9)
      } else if(percentageAssoc == 15){
        #then 6 SNPs are associated
        assocSNPs = c(31,9,4,8,10,5)
      } else if(percentageAssoc == 25){
        #then 10 SNPs are associated
        assocSNPs = c(31,9,4,8,10,5,3,2,7,1)
      }
    } else {
      #looking at SNPs in high LD 
      #then looking at SNPs in low LD
      if(percentageAssoc == 5){
        #then 2 SNPs are associated
        assocSNPs = c(25, 32)
      } else if(percentageAssoc == 15){
        #then 6 SNPs are associated
        assocSNPs = c(37,29,36,28,30,33)
      } else if(percentageAssoc == 25){
        #then 10 SNPs are associated
        assocSNPs = c(37,29,36,28,30,33,34,38,25,32)
      } 
    }
  }
  
  myList = lapply(start:numSims, rep, times = 1)
  nullValues =  altValues = statsAndPValsmat = finalOutput = list()
  
  statsAndPVals = mclapply(myList, function(ii){
    #pull recipient and donor genotypes
    RGenos = obtainRGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    DGenos = obtainDGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    
    #calculate single snp scores
    if(score == "IBS"){
      Score.snp = calcIBSMismatch(RGenosMat = RGenos, DGenosMat = DGenos)
    } else if(score == "Incomp"){
      Score.snp = calcIncompatibilityScore(RGenosMat = RGenos, DGenosMat = DGenos)
    } else if(score == "AMS"){
      Score.snp = calcAMS(RGenosMat = RGenos, DGenosMat = DGenos)
    } else{
      Score.snp = calcBinaryMM(RGenosMat = RGenos, DGenosMat = DGenos)
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
    
    #generate covariates
    #for now, a single binary and a single continous covariate
    # CovData = GenCovData(SampleSize = numPairs, BinaryValues = 1, ContinuousValues = 1)
    
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
    
    #define null gamma
    Gamma = c(0)
    
    #generate phenotypes, both continuous and binary
    ContPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = FALSE, YCat = FALSE,  RGenoData = RGenos, ScoreData = Score.gene, Betas = Betas, Gamma = Gamma)
    # ContPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = FALSE,  Covariates = CovData, RGenoData = RGenos, ScoreData = score, Betas = Betas, Gamma = Gamma)
    
    #combine the data for x design matrix
    allData = cbind(RGenos, Score.gene)
    
    #split into training and test data
    #train with 80% of data
    train.null.index = sample(1:nrow(RGenos), 5*(nrow(RGenos)/10))
    train.all.index = sample(1:nrow(allData), 5*(nrow(allData)/10))
    #pull the index values for the test data (20%)
    test.null.index = -train.null.index
    test.all.index = -train.all.index
    #separate the data
    x.train.null = RGenos[train.null.index,]
    x.train.all = allData[train.all.index,]
    x.test.null = RGenos[test.null.index,]
    x.test.all = allData[test.all.index,]
    
    #pull the y phenos for the 
    y.train.null = ContPhenos[train.null.index]
    y.train.all = ContPhenos[train.all.index]
    y.test.null = ContPhenos[test.null.index]
    y.test.all = ContPhenos[test.all.index]
    
    #use cv to find best lambda value
    cv_output_null = cv.glmnet(x.train.null, y.train.null, family = "gaussian", alpha=1, intercept=FALSE)
    cv_output_alt = cv.glmnet(x.train.all, y.train.all, family = "gaussian", alpha=1, intercept=FALSE)
    cv_output_null.int = cv.glmnet(x.train.null, y.train.null, family = "gaussian", alpha=1, intercept=TRUE)
    cv_output_alt.int = cv.glmnet(x.train.all, y.train.all, family = "gaussian", alpha=1, intercept=TRUE)
    
    #identifying the best lambda
    best_lambda_null = cv_output_null$lambda.min
    best_lambda_alt = cv_output_alt$lambda.min
    best_lambda_null.int = cv_output_null.int$lambda.min
    best_lambda_alt.int = cv_output_alt.int$lambda.min
    
    #fit the null and alternative models with the best lambda values
    fitNull = glmnet(RGenos, ContPhenos, family = "gaussian", alpha=1, lambda = best_lambda_null, intercept=FALSE)
    fitAlt = glmnet(allData, ContPhenos, family = "gaussian", alpha=1, lambda = best_lambda_alt, intercept=FALSE)	
    fitNull.int = glmnet(RGenos, ContPhenos, family = "gaussian", alpha=1, lambda = best_lambda_null.int, intercept=TRUE)
    fitAlt.int = glmnet(allData, ContPhenos, family = "gaussian", alpha=1, lambda = best_lambda_alt.int, intercept=TRUE)	
    
    nullValues = coef(fitNull)
    altValues = coef(fitAlt)
    nullValues.int = coef(fitNull.int)
    altValues.int = coef(fitAlt.int)
    
    nullValues.lambda = fitNull$lambda
    altValues.lambda = fitAlt$lambda
    nullValues.int.lambda = fitNull.int$lambda
    altValues.int.lambda = fitAlt.int$lambda
    
    lambdas = c(nullValues.lambda, altValues.lambda, nullValues.int.lambda, altValues.int.lambda)
    
    print(paste0("Simulation ",ii," is complete."))
    
    finalOutput[[ii]] = list(nullValues, altValues, nullValues.int, altValues.int, lambdas)
  }, mc.cores=1)
  
  betasNull = list()
  betasAll = list()
  betasNull.int = list()
  betasAll.int = list()
  lambdas = list()
  
  for(ii in 1:length(statsAndPVals)){
    betasNull[[ii]] = statsAndPVals[[ii]][1]
    betasAll[[ii]] =  statsAndPVals[[ii]][2]
    betasNull.int[[ii]] = statsAndPVals[[ii]][3]
    betasAll.int[[ii]] = statsAndPVals[[ii]][4]
    lambdas[[ii]] = statsAndPVals[[ii]][5]
  }
  
  betasNull.mat =  do.call(rbind, unlist(betasNull, recursive = FALSE))
  betasAll.mat = do.call(rbind, unlist(betasAll, recursive = FALSE))
  betasNull.int.mat = do.call(rbind, unlist(betasNull.int, recursive = FALSE))
  betasAll.int.mat = do.call(rbind, unlist(betasAll.int, recursive = FALSE))
  lambdas.mat = matrix(unlist(lambdas), ncol = 4, byrow = TRUE)
  
  betasNullmat = as.matrix(betasNull.mat)
  betasAllmat = as.matrix(betasAll.mat)
  betasNull.intmat = as.matrix(betasNull.int.mat)
  betasAll.intmat = as.matrix(betasAll.int.mat)
  
  #write out the Stats and p values
  if(weightedScores == FALSE){
    if(standardizeScores == FALSE){
      write.csv(betasNullmat, file = paste0(path,"/TIE_ContNullBetas_LinHypTest_Lasso_BetaSig_Score",score,"_Sim",start,"to",numSims,"_StatsAndPValues.csv"))
      write.csv(betasAllmat, file = paste0(path,"/TIE_ContAltBetas_LinHypTest_Lasso_BetaSig_Score",score,"_Sim",start,"to",numSims,"_StatsAndPValues.csv"))
      write.csv(betasNull.intmat, file = paste0(path,"/TIE_ContNullBetasIntercept_LinHypTest_Lasso_BetaSig_Score",score,"_Sim",start,"to",numSims,"_StatsAndPValues.csv"))
      write.csv(betasAll.intmat, file = paste0(path,"/TIE_ContAltBetasIntercept_LinHypTest_Lasso_BetaSig_Score",score,"_Sim",start,"to",numSims,"_StatsAndPValues.csv"))
      write.csv(lambdas.mat, file = paste0(path,"/TIE_ContLambdas_LinHypTest_Lasso_BetaSig_Score",score,"_Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    } else { #scores are standardized
      write.csv(statsAndPValsContPhenos, file = paste0(path,"/TIE_ContNullBetas_LinHypTest_BetaSig_Score",score,"_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
      write.csv(nulValuesAll.mat, file = paste0(path,"/TIE_ContAltBetas_LinHypTest_BetaSig_Score",score,"_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
      write.csv(nulValuesAll.mat, file = paste0(path,"/TIE_ContNullBetasIntercept_LinHypTest_BetaSig_Score",score,"_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
      write.csv(nulValuesAll.mat, file = paste0(path,"/TIE_ContAltBetasIntercept_LinHypTest_BetaSig_Score",score,"_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
      write.csv(altValuesAll.mat, file = paste0(path,"/TIE_ContLambdas_LinHypTest_BetaSig_Score",score,"_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
    }
  } else { #scores are weighted
    write.csv(statsAndPValsContPhenos, file = paste0(path,"/TIE_ContNullBetas_LinHypTest_BetaSig_Score",score,"_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
    write.csv(nulValuesAll.mat, file = paste0(path,"/TIE_ContAltBetas_LinHypTest_BetaSig_Score",score,"_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
    write.csv(nulValuesAll.mat, file = paste0(path,"/TIE_ContNullBetasIntercept_LinHypTest_BetaSig_Score",score,"_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
    write.csv(nulValuesAll.mat, file = paste0(path,"/TIE_ContAltBetasIntercept_LinHypTest_BetaSig_Score",score,"_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
    write.csv(altValuesAll.mat, file = paste0(path,"/TIE_ContLambdas_LinHypTest_BetaSig_Score",score,"_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
  }
}

#categorical outcome
RunTIEPipelineLinHypTestCat_GammaTest_BetaSig = function(chr, gene, numPairs, YPrev, standardizeScores = FALSE, weightedScores = FALSE, scoreWeights, ORSize, percentageAssoc, LowLD, score, start, stop){
  #function to run whole TIE pipeline, Calculates Score test stat for Lin Hyp test method and whether p-value <0.05
  # Only use when Y is binary
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #standardizeScores = T or F whether the scores should be standardized based on maximum score value
  #weightedScores = T or F whether the scores will be weighted
  #scoreWeights = m x 1 vector of weights, one weight for each SNP
  # score - which score we are using in the main JST model
  # start - simulation number to start at
  # stop - simulation number to stop at
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes TIE values to csv files

  library(parallel)

  #always the same
  numSims = stop

  #define effect based on OR size
  if(ORSize == "small"){
    effect = 0.14
  } else if(ORSize == "medium"){
    effect = 0.41
  } else {
    effect = 0.69
  }

  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_II_Sims/HapGen_Files/",gene,"_Results_",numPairs,"Pairs")

  #source the needed functions
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/ProjectIISourceFunctions_v2.R")
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/Logistic_ADMM0.r")

  #determine which SNPs to actually set as assoc. based on gene
  assocSNPs = DetermineAssocRSNPs(gene = gene, LowLD = LowLD, percentageAssoc = percentageAssoc)

  myList = lapply(start:numSims, rep, times = 1)
  
  # p-value initialization
  pv <- rep(0, 3)
  Tall <- matrix(0, 1, 3)
  beta.al <- matrix(0, 1, 2)
  
  statsAndPVals = mclapply(myList, function(ii){
    #define matrix to hold all Stats and Pvalues
    statsAndPVals = matrix(NA, nrow = numSims, ncol = 4)

    #pull recipient and donor genotypes
    RGenos = obtainRGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    DGenos = obtainDGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)

    #calculate single snp scores
    if(score == "IBS"){
      Score.snp = calcIBSMismatch(RGenosMat = RGenos, DGenosMat = DGenos)
    } else if(score == "Incomp"){
      Score.snp = calcIncompatibilityScore(RGenosMat = RGenos, DGenosMat = DGenos)
    } else if(score == "AMS"){
      Score.snp = calcAMS(RGenosMat = RGenos, DGenosMat = DGenos)
    } else{
      Score.snp = calcBinaryMM(RGenosMat = RGenos, DGenosMat = DGenos)
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

    #define null gamma
    Gamma = c(0)

    #generate phenotypes, both continuous and binary
    CatPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = TRUE, YPrev = YPrev,  Covariates = CovData, RGenoData = RGenos, ScoreData = Score.gene, Betas = Betas, Gamma = Gamma)

    # define location of zero components
    # basically, this defines what the null hypothesis is, I think
    # so for joint null, we need the non-zero components to be for covariates
    # and the beta and gamma components to be zero
    numCov = dim(CovData)[2]
    numSNPs = dim(RGenos)[2]
    N = c(rep(FALSE, numCov+1+numSNPs), TRUE)

    #need to combine the covariates, r genos, and score matrices together
    #need to include column of 1s for intercept?
    intercept = matrix(1,nrow = numPairs, ncol = 1)
    designMat = cbind(intercept, CovData, RGenos, Score.gene)

    #then combine the phenos with the design matrix as a list
    Model = list(X=designMat, Y=CatPhenos)

    # estimate the uncontrained estimator
    beta.unre <- cv.SCAD_ADMM_unre(X=Model$X, Y=Model$Y, N=N, beta0=rep(0, dim(Model$X)[2]), err=1e-4, tune="cv", unpen=c(1,2,3))
    indice.unre <- beta.unre!=0

    # estimate the constrained estimator (only need this for score test)
    beta.re <- cv.SCAD_ADMM_re(X=Model$X, Y=Model$Y, N=N, beta0=rep(0, dim(Model$X)[2]), err=1e-4, tune="cv", unpen=c(1,2,3))
    indice.re <- beta.re!=0

    pi.unre <- logit(Model$X%*%beta.unre)
    pi.re <- logit(Model$X%*%beta.re)

    # construct the likelihood ratio statistics
    TL <- 2*(sum(log(1+exp(Model$X%*%beta.re))-Model$Y*(Model$X%*%beta.re))-
               sum(log(1+exp(Model$X%*%beta.unre))-Model$Y*(Model$X%*%beta.unre)))
    #if LRT stat is less than 0, there is error, so set to -1
    if(TL < 0){
      TL = -1
    }
    
    # construct the Wald statistics
    #B_0 should give you Omega_a hat
    A = crossprod(Model$X[,indice.unre|N], as.vector(pi.unre*(1-pi.unre))*Model$X[,indice.unre|N])
    if(rcond(A) >= 1e-10){
      B_01 <- solve(A)
      #d_0 should determine which rows to subset
      #gives the first value that is restricted
      d_01 <- sum(which(indice.unre|N)<=length(N))
      #so the B_0 needs to be subset to m rows and columns that are restricted based on H0
      TW <- crossprod(beta.unre[N], solve(B_01[d_01,d_01], beta.unre[N]))
    } else {
      TW = -1
    }
    
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

    #set pv to 1 if reject, 0 otherwise
    #LRT
    if (TL>=qchisq(0.95, df = doF)){
      pv[1] <- pv[1]+1/5000
    }
    #Wald
    if (TW>=qchisq(0.95, df = doF)){
      pv[2] <- pv[2]+1/5000
    }
    #Score
    if (TS>=qchisq(0.95, df = doF)){
      pv[3] <- pv[3]+1/5000
    }

    print(paste0("Simulation ",ii," is complete."))

    Tall[1, ] <- c(TL, TW, TS)
    beta.al[1, ] <- c(sum(beta.re!=0), sum(beta.unre!=0))
    
    statsAndPVals = list(pv=pv, TScores=Tall, beta=beta.al)
    statsAndPVals
  }, mc.cores=4)

  n = length(statsAndPVals)
  
  pValsCV = matrix(nrow = n, ncol = 3)
  ScoresCV = matrix(nrow = n, ncol = 3)
  BetasCV = matrix(nrow = n, ncol = 2)
  
  for(ll in 1:n){
    pValsCV[ll,] = unlist(statsAndPVals[[ll]]$pv)
    ScoresCV[ll,] = unlist(statsAndPVals[[ll]]$TScores)
    BetasCV[ll,] = unlist(statsAndPVals[[ll]]$beta)
  }
  
  #write out the Stats and p values
  if(weightedScores == FALSE){
    if(standardizeScores == FALSE){
      write.csv(pValsCV, file = paste0(path,"/TIE_CatY_Prev",YPrev*100,"_LinHypTest_Score",score,"_GammaOnly_BetaSig_",ORSize,"OR_",percentageAssoc,"AssocSNPs_",LowLD,"_Sim",start,"to",numSims,"_PValuesCV_ForceCovFit.csv"))
      write.csv(ScoresCV, file = paste0(path,"/TIE_CatY_Prev",YPrev*100,"_LinHypTest_Score",score,"_GammaOnly_BetaSig_",ORSize,"OR_",percentageAssoc,"AssocSNPs_",LowLD,"_Sim",start,"to",numSims,"_StatsCV_ForceCovFit.csv"))
      write.csv(BetasCV, file = paste0(path,"/TIE_CatY_Prev",YPrev*100,"_LinHypTest_Score",score,"_GammaOnly_BetaSig_",ORSize,"OR_",percentageAssoc,"AssocSNPs_",LowLD,"_Sim",start,"to",numSims,"_BetasCV_ForceCovFit.csv"))
    } else { #scores are standardized
      write.csv(statsAndPValsCatPhenos, file = paste0(path,"/TIE_CatPhenos_Prev",YPrev*100,"_LinHypTest_Score",score,"_GammaOnly_NoCovs_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
    }
  } else { #scores are weighted
    write.csv(statsAndPValsCatPhenos, file = paste0(path,"/TIE_CatPhenos_Prev",YPrev*100,"_LinHypTest_Score",score,"_GammaOnly_NoCovs_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
  }
}

#continuous outcome
RunTIEPipelineLinHypTestCont_GammaTest_BetaSig = function(chr, gene, numPairs, standardizeScores = FALSE, weightedScores = FALSE, scoreWeights, ORSize, percentageAssoc, LowLD, score, start, stop){
  #function to run whole TIE pipeline, Calculates LRT test stat for Lin Hyp test method and whether p-value <0.05
  # only use when Y is continuous (Normal)
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #standardizeScores = T or F whether the scores should be standardized based on maximum score value
  #weightedScores = T or F whether the scores will be weighted
  #scoreWeights = m x 1 vector of weights, one weight for each SNP
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes TIE values to csv files

  library(parallel)

  #always the same
  numSims = stop

  #define effect based on OR size
  if(ORSize == "small"){
    effect = 0.14
  } else if(ORSize == "medium"){
    effect = 0.41
  } else {
    effect = 0.69
  }

  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_II_Sims/HapGen_Files/",gene,"_Results_",numPairs,"Pairs")

  #source the needed functions
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/ProjectIISourceFunctions_v2.R")
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/Linear_ADMM0.r")

  #determine which SNPs to actually set as assoc. based on gene
  assocSNPs = DetermineAssocRSNPs(gene = gene, LowLD = LowLD, percentageAssoc = percentageAssoc)

  myList = lapply(start:numSims, rep, times = 1)
  # p-value initialization
  pv <- rep(0, 6)
  Tall <- matrix(0, 1, 3)
  beta.al <- matrix(0, 1, 2)

  statsAndPVals = mclapply(myList, function(ii){
    #define matrix to hold all Stats and Pvalues
    statsAndPVals = matrix(NA, nrow = numSims, ncol = 6)

    #pull recipient and donor genotypes
    RGenos = obtainRGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)
    DGenos = obtainDGenotypes(chr = chr, numSamples = numPairs, simNum = ii, gene = gene, path = path)

    #calculate single snp scores
    if(score == "IBS"){
      Score.snp = calcIBSMismatch(RGenosMat = RGenos, DGenosMat = DGenos)
    } else if(score == "Incomp"){
      Score.snp = calcIncompatibilityScore(RGenosMat = RGenos, DGenosMat = DGenos)
    } else if(score == "AMS"){
      Score.snp = calcAMS(RGenosMat = RGenos, DGenosMat = DGenos)
    } else{
      Score.snp = calcBinaryMM(RGenosMat = RGenos, DGenosMat = DGenos)
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

    #define null gamma
    Gamma = c(0)

    #generate phenotypes, both continuous and binary
    ContPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = FALSE,  Covariates = CovData, RGenoData = RGenos, ScoreData = Score.gene, Betas = Betas, Gamma = Gamma)

    # define location of zero components
    # basically, this defines what the null hypothesis is, I think
    # so for joint null, we need the non-zero components to be for covariates
    # and the beta and gamma components to be zero
    numCov = dim(CovData)[2]
    numSNPs = dim(RGenos)[2]
    N = c(rep(FALSE, numCov+1+numSNPs), TRUE)

    #need to combine the covariates, r genos, and score matrices together
    #need to include column of 1s for intercept?
    intercept = matrix(1,nrow = numPairs, ncol = 1)
    designMat = cbind(intercept, CovData, RGenos, Score.gene)

    #then combine the phenos with the design matrix as a list
    Model = list(X=designMat, Y=ContPhenos)

    # estimate the uncontrained estimator
    beta.unre <- cv.SCAD_ADMM_unre(X=Model$X, Y=Model$Y, N=N, beta0=rep(0, dim(Model$X)[2]), err=1e-4, tune="cv", unpen=c(1,2,3))
    indice.unre <- beta.unre!=0

    # estimate the constrained estimator
    beta.re <- cv.SCAD_ADMM_re(X=Model$X, Y=Model$Y, N=N, beta0=rep(0, dim(Model$X)[2]), err=1e-4, tune="cv", unpen=c(1,2,3))
    indice.re <- beta.re!=0

    # estimate the conditional variance
    #    sig2 <- mean((Model$Y-Model$X%*%beta.unre)^2)
    n = numSims
    sig2 <- mean((Model$Y-Model$X%*%beta.unre)^2)*n/(n-sum(beta.unre!=0))

    # construct the likelihood ratio statistic
    TL <- sum((Model$X%*%beta.re-Model$Y)^2)-sum((Model$X%*%beta.unre-Model$Y)^2)
    if(TL < 0){
      TL = -1
    }
    
    # construct the Wald statistic
    B = crossprod(Model$X[,indice.unre|N], Model$X[,indice.unre|N])
    if(rcond(B) >= 1e-10){
      B_0 <- solve(B)
      d_0 <- sum(which(indice.unre|N)<=length(N))
      TW <- crossprod(beta.unre[N], solve(B_0[d_0,d_0], beta.unre[N]))
    } else {
      TW = -1
    }
    
    # construct the score statistic
    eps <- Model$Y-Model$X%*%beta.re
    Xeps <- crossprod(Model$X[,indice.re|N], eps)
    D = crossprod(Model$X[,indice.re|N], Model$X[,indice.re|N])
    if(rcond(D) >= 1e-10){
      TS <- crossprod(Xeps, solve(D, Xeps))
    } else {
      TS = -1
    }
    
    #degrees of freedom is equal to the number of restricted values
    doF = sum(N == TRUE)

    #set pv to 1 if reject, 0 otherwise
    #LRT
    if (TL>=sig2*qchisq(0.95, df = doF)){
      pv[1] <- pv[1]+1/5000
    }
    if (TL>=qchisq(0.95, df = doF)){
      pv[2] <- pv[2]+1/5000
    }
    #Wald
    if (TW>=sig2*qchisq(0.95, df = doF)){
      pv[3] <- pv[3]+1/5000
    }
    if (TW>=qchisq(0.95, df = doF)){
      pv[4] <- pv[4]+1/5000
    }
    #Score
    if (TS>=sig2*qchisq(0.95, df = doF)){
      pv[5] <- pv[5]+1/5000
    }
    if (TS>=qchisq(0.95, df = doF)){
      pv[6] <- pv[6]+1/5000
    }

    print(paste0("Simulation ",ii," is complete."))

    Tall[1, ] <- c(TL, TW, TS)
    beta.al[1, ] <- c(sum(beta.re!=0), sum(beta.unre!=0))
    
    statsAndPVals = list(pv=pv, TScores=Tall, beta=beta.al)
    statsAndPVals
  }, mc.cores=4)

  n = length(statsAndPVals)
  
  pVals = matrix(nrow = n, ncol = 6)
  Scores = matrix(nrow = n, ncol = 3)
  Betas = matrix(nrow = n, ncol = 2)
  
  for(ll in 1:n){
    pVals[ll,] = unlist(statsAndPVals[[ll]]$pv)
    Scores[ll,] = unlist(statsAndPVals[[ll]]$TScores)
    Betas[ll,] = unlist(statsAndPVals[[ll]]$beta)
  }
  
  #write out the Stats and p values
  if(weightedScores == FALSE){
    if(standardizeScores == FALSE){
      write.csv(pVals, file = paste0(path,"/TIE_ContY_LinHypTest_Score",score,"_GammaOnly_BetaSig_",ORSize,"OR_",percentageAssoc,"AssocSNPs_",LowLD,"_Sim",start,"to",numSims,"_PValuesCV_ForceCovs.csv"))
      write.csv(Scores, file = paste0(path,"/TIE_ContY_LinHypTest_Score",score,"_GammaOnly_BetaSig_",ORSize,"OR_",percentageAssoc,"AssocSNPs_",LowLD,"_Sim",start,"to",numSims,"_StatsCV_ForceCovs.csv"))
      write.csv(Betas, file = paste0(path,"/TIE_ContY_LinHypTest_Score",score,"_GammaOnly_BetaSig_",ORSize,"OR_",percentageAssoc,"AssocSNPs_",LowLD,"_Sim",start,"to",numSims,"_BetasCV_ForceCovs.csv"))
    } else { #scores are standardized
      write.csv(statsAndPValsContPhenos, file = paste0(path,"/TIE_ContPhenos_LinHypTest_Score",score,"_GammaOnly_NoCovs_Sim",start,"to",numSims,"_StandardizedScores_StatsAndPValues.csv"))
    }
  } else { #scores are weighted
    write.csv(statsAndPValsContPhenos, file = paste0(path,"/TIE_ContPhenos_LinHypTest_Score",score,"_GammaOnly_NoCovs_Sim",start,"to",numSims,"_WeightedScores_StatsAndPValues.csv"))
  }
}

#############################
## Power Pipelines

## Score is assoc (gamma neq 0), beta associated (neq 0)
#categorical outcome
RunPowerPipelineLinHypTestCat_GammaTest_Score_BetaNZero = function(chr, gene, numPairs, YPrev, standardizeScores = FALSE, weightedScores = FALSE, scoreWeights, ORSize, percentageAssocGene, percentageAssocSNP, LowLD, TrueScore, start, stop){
  #Function to determine power of linear hyp testing method when score is associated, gamma = 0 testing only
  # Only use when Y is binary
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #TrueScore = IBS.gene, Incomp.gene, AMS.gene, or BinMM.gene
  #ORSize = Small, Medium, or Large for what OR was used for the associated SNP/score
  #standardizeScores = T or F whether the scores should be standardized based on maximum score value
  #weightedScores = T or F whether the scores will be weighted
  #scoreWeights = m x 1 vector of weights, one weight for each SNP
  # start - simulation number to start at
  # stop - simulation number to stop at
  #percentageAssocGene = percentage of SNPs associated with outcome (either 5, 25, 50, 75, or 100) 
  #percentageAssocSNP = percentage of SNPs associated with outcome (5, 15, or 25)
  #LowLD = True or FALSE whether the associated SNPs are in low LD or high LD
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes power values to csv files
  
  library(parallel)
  
  #need to define this for naming at the end
  snpOrScore = TrueScore
  
  #always the same
  numSims = stop
  
  #define effect based on OR size
  if(ORSize == "small"){
    effect = 0.14
  } else if(ORSize == "medium"){
    effect = 0.41
  } else {
    effect = 0.69
  }
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_II_Sims/HapGen_Files/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/ProjectIISourceFunctions_v2.R")
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/Logistic_ADMM0.r")
  
  #determine which SNPs to actually set as assoc. based on gene
  assocSNPs = DetermineAssocRSNPs(gene = gene, LowLD = LowLD, percentageAssoc = percentageAssocSNP)
  
  myList = lapply(start:numSims, rep, times = 1)
  # p-value initialization
  pv.IBS <- rep(0, 3)
  pv.Incomp <- rep(0, 3)
  pv.AMS <- rep(0, 3)
  pv.BinMM <- rep(0, 3)
  
  Tall.IBS <- matrix(0, 1, 3)
  Tall.Incomp <- matrix(0, 1, 3)
  Tall.AMS <- matrix(0, 1, 3)
  Tall.BinMM <- matrix(0, 1, 3)
  
  beta.al.IBS <- matrix(0, 1, 2)
  beta.al.Incomp <- matrix(0, 1, 2)
  beta.al.AMS <- matrix(0, 1, 2)
  beta.al.BinMM <- matrix(0, 1, 2)
  
  statsAndPVals = mclapply(myList, function(ii){
    #define matrix to hold all Stats and Pvalues
    statsAndPVals = matrix(NA, nrow = numSims, ncol = 4)
    
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
    IBS.gene.PercentOfSNPs = calcGeneScorePercentOfSNPs(SingleSNPKernel = IBS.snp, gene = gene, percentageAssoc = percentageAssocGene, LowLD = LowLD, standardize = FALSE, useWeights = FALSE)
    Incomp.gene.PercentOfSNPs = calcGeneScorePercentOfSNPs(SingleSNPKernel = Incomp.snp, gene =  gene, percentageAssoc = percentageAssocGene, LowLD = LowLD, standardize = FALSE, useWeights = FALSE)
    AMS.gene.PercentOfSNPs = calcGeneScorePercentOfSNPs(SingleSNPKernel = AMS.snp, gene =  gene, percentageAssoc = percentageAssocGene, LowLD = LowLD, standardize = FALSE, useWeights = FALSE)
    BinMM.gene.PercentOfSNPs = calcGeneScorePercentOfSNPs(SingleSNPKernel = BinMM.snp, gene =  gene, percentageAssoc = percentageAssocGene, LowLD = LowLD, standardize = FALSE, useWeights = FALSE)
    
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
    nullBetas = rep(0,nSNP) #generate null beta values
    Betas = nullBetas
    
    #set assoc Betas
    #all betas have same effect for now
    for(jj in assocSNPs){
      Betas[jj] = effect
      Betas = as.matrix(Betas, ncol = 1)
    }
    
    Gamma = effect
    
    #generate phenotypes, both continuous and binary
    CatPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = TRUE, YPrev = YPrev,  Covariates = CovData, RGenoData = RGenos, ScoreData = PhenoScore, Betas = Betas, Gamma = Gamma)
    
    # define location of zero components
    # basically, this defines what the null hypothesis is, I think
    # so for joint null, we need the non-zero components to be for covariates
    # and the beta and gamma components to be zero
    numCov = dim(CovData)[2]
    numSNPs = dim(RGenos)[2]
    N = c(rep(FALSE, numCov+1+numSNPs), TRUE)

    #need to combine the covariates, r genos, and score matrices together
    #need to include column of 1s for intercept?
    intercept = matrix(1,nrow = numPairs, ncol = 1)
    designMat.IBS = cbind(intercept, CovData, RGenos, IBS.gene)
    designMat.Incomp = cbind(intercept, CovData, RGenos, Incomp.gene)
    designMat.AMS = cbind(intercept, CovData, RGenos, AMS.gene)
    designMat.BinMM = cbind(intercept, CovData, RGenos, BinMM.gene)
    
    #then combine the phenos with the design matrix as a list
    Model.IBS = list(X=designMat.IBS, Y=CatPhenos)
    Model.Incomp = list(X=designMat.Incomp, Y=CatPhenos)
    Model.AMS = list(X=designMat.AMS, Y=CatPhenos)
    Model.BinMM = list(X=designMat.BinMM, Y=CatPhenos)
    
    # estimate the uncontrained estimator
    beta.unre.IBS <- cv.SCAD_ADMM_unre(X=Model.IBS$X, Y=Model.IBS$Y, N=N, beta0=rep(0, dim(Model.IBS$X)[2]), err=1e-4, tune="cv", unpen=c(1,2,3))
    indice.unre.IBS <- beta.unre.IBS!=0
    
    beta.unre.Incomp <- cv.SCAD_ADMM_unre(X=Model.Incomp$X, Y=Model.Incomp$Y, N=N, beta0=rep(0, dim(Model.Incomp$X)[2]), err=1e-4, tune="cv", unpen=c(1,2,3))
    indice.unre.Incomp <- beta.unre.Incomp!=0
    
    beta.unre.AMS <- cv.SCAD_ADMM_unre(X=Model.AMS$X, Y=Model.AMS$Y, N=N, beta0=rep(0, dim(Model.AMS$X)[2]), err=1e-4, tune="cv", unpen=c(1,2,3))
    indice.unre.AMS <- beta.unre.AMS!=0
    
    beta.unre.BinMM <- cv.SCAD_ADMM_unre(X=Model.BinMM$X, Y=Model.BinMM$Y, N=N, beta0=rep(0, dim(Model.BinMM$X)[2]), err=1e-4, tune="cv", unpen=c(1,2,3))
    indice.unre.BinMM <- beta.unre.BinMM!=0
    
    #estimate the constrained estimator (only need this for score test)
    beta.re.IBS <- cv.SCAD_ADMM_re(X=Model.IBS$X, Y=Model.IBS$Y, N=N, beta0=rep(0, dim(Model.IBS$X)[2]), err=1e-4, tune="cv", unpen=c(1,2,3))
    indice.re.IBS <- beta.re.IBS!=0
    
    beta.re.Incomp <- cv.SCAD_ADMM_re(X=Model.Incomp$X, Y=Model.Incomp$Y, N=N, beta0=rep(0, dim(Model.Incomp$X)[2]), err=1e-4, tune="cv", unpen=c(1,2,3))
    indice.re.Incomp <- beta.re.Incomp!=0
    
    beta.re.AMS <- cv.SCAD_ADMM_re(X=Model.AMS$X, Y=Model.AMS$Y, N=N, beta0=rep(0, dim(Model.AMS$X)[2]), err=1e-4, tune="cv", unpen=c(1,2,3))
    indice.re.AMS <- beta.re.AMS!=0
    
    beta.re.BinMM <- cv.SCAD_ADMM_re(X=Model.BinMM$X, Y=Model.BinMM$Y, N=N, beta0=rep(0, dim(Model.BinMM$X)[2]), err=1e-4, tune="cv", unpen=c(1,2,3))
    indice.re.BinMM <- beta.re.BinMM!=0
    
    pi.unre.AMS <- logit(Model.AMS$X%*%beta.unre.AMS)
    pi.re.AMS <- logit(Model.AMS$X%*%beta.re.AMS)
    
    pi.unre.BinMM <- logit(Model.BinMM$X%*%beta.unre.BinMM)
    pi.re.BinMM <- logit(Model.BinMM$X%*%beta.re.BinMM)
    
    pi.unre.IBS <- logit(Model.IBS$X%*%beta.unre.IBS)
    pi.re.IBS <- logit(Model.IBS$X%*%beta.re.IBS)
    
    pi.unre.Incomp <- logit(Model.Incomp$X%*%beta.unre.Incomp)
    pi.re.Incomp <- logit(Model.Incomp$X%*%beta.re.Incomp)
    
    #construct the likelihood ratio statistics
    TL.IBS <- 2*(sum(log(1+exp(Model.IBS$X%*%beta.re.IBS))-Model.IBS$Y*(Model.IBS$X%*%beta.re.IBS))-
                   sum(log(1+exp(Model.IBS$X%*%beta.unre.IBS))-Model.IBS$Y*(Model.IBS$X%*%beta.unre.IBS)))
    #if LRT stat is less than 0, there is error, so set to -1
    if(TL.IBS < 0){
      TL.IBS = -1
    }
    
    TL.Incomp <- 2*(sum(log(1+exp(Model.Incomp$X%*%beta.re.Incomp))-Model.Incomp$Y*(Model.Incomp$X%*%beta.re.Incomp))-
                      sum(log(1+exp(Model.Incomp$X%*%beta.unre.Incomp))-Model.Incomp$Y*(Model.Incomp$X%*%beta.unre.Incomp)))
    #if LRT stat is less than 0, there is error, so set to -1
    if(TL.Incomp < 0){
      TL.Incomp = -1
    }
    
    TL.AMS <- 2*(sum(log(1+exp(Model.AMS$X%*%beta.re.AMS))-Model.AMS$Y*(Model.AMS$X%*%beta.re.AMS))-
                   sum(log(1+exp(Model.AMS$X%*%beta.unre.AMS))-Model.AMS$Y*(Model.AMS$X%*%beta.unre.AMS)))
    #if LRT stat is less than 0, there is error, so set to -1
    if(TL.AMS < 0){
      TL.AMS = -1
    }
    
    TL.BinMM <- 2*(sum(log(1+exp(Model.BinMM$X%*%beta.re.BinMM))-Model.BinMM$Y*(Model.BinMM$X%*%beta.re.BinMM))-
                     sum(log(1+exp(Model.BinMM$X%*%beta.unre.BinMM))-Model.BinMM$Y*(Model.BinMM$X%*%beta.unre.BinMM)))
    #if LRT stat is less than 0, there is error, so set to -1
    if(TL.BinMM < 0){
      TL.BinMM = -1
    }
    
    # construct the Wald statistics
    #B_0 should give you Omega_a hat
    A.IBS = crossprod(Model.IBS$X[,indice.unre.IBS|N], as.vector(pi.unre.IBS*(1-pi.unre.IBS))*Model.IBS$X[,indice.unre.IBS|N])
    if(rcond(A.IBS) >= 1e-10){
      B_0.IBS <- solve(A.IBS)
      #d_0 should determine which rows to subset
      #gives the first value that is restricted
      d_0.IBS <- dim(B_0.IBS)[1]
      #in case the inverse doesn't exist, need another if else
      B.IBS = B_0.IBS[d_0.IBS,d_0.IBS]
      #so the B_0 needs to be subset to m rows and columns that are restricted based on H0
      TW.IBS <- crossprod(beta.unre.IBS[N], solve(B.IBS, beta.unre.IBS[N]))
    } else {
      TW.IBS = -1
    }
    
    A.Incomp = crossprod(Model.Incomp$X[,indice.unre.Incomp|N], as.vector(pi.unre.Incomp*(1-pi.unre.Incomp))*Model.Incomp$X[,indice.unre.Incomp|N])
    if(rcond(A.Incomp) >= 1e-10){
      B_0.Incomp <- solve(A.Incomp)
      #d_0 should determine which rows to subset
      #gives the first value that is restricted
      d_0.Incomp <- dim(B_0.Incomp)[1]
      #in case the inverse doesn't exist, need another if else
      B.Incomp = B_0.Incomp[d_0.Incomp,d_0.Incomp]
      TW.Incomp <- crossprod(beta.unre.Incomp[N], solve(B.Incomp, beta.unre.Incomp[N]))
    } else {
      TW.Incomp = -1
    }
    
    A.AMS = crossprod(Model.AMS$X[,indice.unre.AMS|N], as.vector(pi.unre.AMS*(1-pi.unre.AMS))*Model.AMS$X[,indice.unre.AMS|N])
    if(rcond(A.AMS) >= 1e-10){
      B_0.AMS <- solve(A.AMS)
      #d_0 should determine which rows to subset
      #gives the first value that is restricted
      d_0.AMS <- dim(B_0.AMS)[1]
      #in case the inverse doesn't exist, need another if else
      B.AMS = B_0.AMS[d_0.AMS,d_0.AMS]
      #so the B_0 needs to be subset to m rows and columns that are restricted based on H0
      TW.AMS <- crossprod(beta.unre.AMS[N], solve(B.AMS, beta.unre.AMS[N]))
    } else {
      TW.AMS = -1
    }
    
    A.BinMM = crossprod(Model.BinMM$X[,indice.unre.BinMM|N], as.vector(pi.unre.BinMM*(1-pi.unre.BinMM))*Model.BinMM$X[,indice.unre.BinMM|N])
    if(rcond(A.BinMM) >= 1e-10){
      B_0.BinMM <- solve(A.BinMM)
      #d_0 should determine which rows to subset
      #gives the first value that is restricted
      d_0.BinMM <- dim(B_0.BinMM)[1]
      #in case the inverse doesn't exist, need another if else
      B.BinMM = B_0.BinMM[d_0.BinMM,d_0.BinMM]
      #so the B_0 needs to be subset to m rows and columns that are restricted based on H0
      TW.BinMM <- crossprod(beta.unre.BinMM[N], solve(B.BinMM, beta.unre.BinMM[N]))
    } else {
      TW.BinMM = -1
    }
    
    # construct the score statistics
    eps.IBS <- Model.IBS$Y-pi.re.IBS #Y - e(Y)
    #this is the X^T times Y-E(Y)
    Xeps.IBS <- crossprod(Model.IBS$X[,indice.re.IBS|N], eps.IBS)
    
    C.IBS = crossprod(Model.IBS$X[,indice.re.IBS|N], as.vector(pi.re.IBS*(1-pi.re.IBS))*Model.IBS$X[,indice.re.IBS|N])
    if(rcond(C.IBS) >= 1e-10){
      TS.IBS <- crossprod(Xeps.IBS, solve(C.IBS, Xeps.IBS))
    } else {
      TS.IBS = -1
    }
    
    eps.Incomp <- Model.Incomp$Y-pi.re.Incomp #Y - e(Y)
    #this is the X^T times Y-E(Y)
    Xeps.Incomp <- crossprod(Model.Incomp$X[,indice.re.Incomp|N], eps.Incomp)
    
    C.Incomp = crossprod(Model.Incomp$X[,indice.re.Incomp|N], as.vector(pi.re.Incomp*(1-pi.re.Incomp))*Model.Incomp$X[,indice.re.Incomp|N])
    if(rcond(C.Incomp) >= 1e-10){
      TS.Incomp <- crossprod(Xeps.Incomp, solve(C.Incomp, Xeps.Incomp))
    } else {
      TS.Incomp = -1
    }
    
    eps.AMS <- Model.AMS$Y-pi.re.AMS #Y - e(Y)
    #this is the X^T times Y-E(Y)
    Xeps.AMS <- crossprod(Model.AMS$X[,indice.re.AMS|N], eps.AMS)
    
    C.AMS = crossprod(Model.AMS$X[,indice.re.AMS|N], as.vector(pi.re.AMS*(1-pi.re.AMS))*Model.AMS$X[,indice.re.AMS|N])
    if(rcond(C.AMS) >= 1e-10){
      TS.AMS <- crossprod(Xeps.AMS, solve(C.AMS, Xeps.AMS))
    } else {
      TS.AMS = -1
    }
    
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
    
    #LRT
    if (TL.IBS>=qchisq(0.95, df = doF)){
      pv.IBS[1] <- pv.IBS[1]+1/5000
    }
    if (TL.Incomp>=qchisq(0.95, df = doF)){
      pv.Incomp[1] <- pv.Incomp[1]+1/5000
    }
    if (TL.AMS>=qchisq(0.95, df = doF)){
      pv.AMS[1] <- pv.AMS[1]+1/5000
    }
    if (TL.BinMM>=qchisq(0.95, df = doF)){
      pv.BinMM[1] <- pv.BinMM[1]+1/5000
    }
    #Wald
    if (TW.IBS>=qchisq(0.95, df = doF)){
      pv.IBS[2] <- pv.IBS[2]+1/5000
    }
    if (TW.Incomp>=qchisq(0.95, df = doF)){
      pv.Incomp[2] <- pv.Incomp[2]+1/5000
    }
    if (TW.AMS>=qchisq(0.95, df = doF)){
      pv.AMS[2] <- pv.AMS[2]+1/5000
    }
    if (TW.BinMM>=qchisq(0.95, df = doF)){
      pv.BinMM[2] <- pv.BinMM[2]+1/5000
    }
    #Score
    if (TS.IBS>=qchisq(0.95, df = doF)){
      pv.IBS[3] <- pv.IBS[3]+1/5000
    }
    if (TS.Incomp>=qchisq(0.95, df = doF)){
      pv.Incomp[3] <- pv.Incomp[3]+1/5000
    }
    if (TS.AMS>=qchisq(0.95, df = doF)){
      pv.AMS[3] <- pv.AMS[3]+1/5000
    }
    if (TS.BinMM>=qchisq(0.95, df = doF)){
      pv.BinMM[3] <- pv.BinMM[3]+1/5000
    }
    
    print(paste0("Simulation ",ii," is complete."))
    
    Tall.IBS[1, ] <- c(TL.IBS, TW.IBS, TS.IBS)
    beta.al.IBS[1, ] <- c(sum(beta.re.IBS!=0), sum(beta.unre.IBS!=0))
    
    Tall.Incomp[1, ] <- c(TL.Incomp, TW.Incomp, TS.Incomp)
    beta.al.Incomp[1, ] <- c(sum(beta.re.Incomp!=0), sum(beta.unre.Incomp!=0))
    
    Tall.AMS[1, ] <- c(TL.AMS, TW.AMS, TS.AMS)
    beta.al.AMS[1, ] <- c(sum(beta.re.AMS!=0), sum(beta.unre.AMS!=0))
    
    Tall.BinMM[1, ] <- c(TL.BinMM, TW.BinMM, TS.BinMM)
    beta.al.BinMM[1, ] <- c(sum(beta.re.BinMM!=0), sum(beta.unre.BinMM!=0))
    
    statsAndPVals = list(pv.IBS=pv.IBS, TScores.IBS=Tall.IBS, beta.IBS=beta.al.IBS,
                         pv.Incomp=pv.Incomp, TScores.Incomp=Tall.Incomp, beta.Incomp=beta.al.Incomp,
                         pv.AMS=pv.AMS, TScores.AMS=Tall.AMS, beta.AMS=beta.al.AMS,
                         pv.BinMM=pv.BinMM, TScores.BinMM=Tall.BinMM, beta.BinMM=beta.al.BinMM)
    statsAndPVals
  }, mc.cores=4)
  
  n = length(statsAndPVals)
  
  pVals.IBS = pVals.Incomp = pVals.AMS = pVals.BinMM = matrix(nrow = n, ncol = 3)
  Scores.IBS = Scores.Incomp = Scores.AMS = Scores.BinMM = matrix(nrow = n, ncol = 3)
  Betas.IBS = Betas.Incomp = Betas.AMS = Betas.BinMM = matrix(nrow = n, ncol = 2)
  
  for(ll in 1:n){
    pVals.IBS[ll,] = unlist(statsAndPVals[[ll]]$pv.IBS)
    Scores.IBS[ll,] = unlist(statsAndPVals[[ll]]$TScores.IBS)
    Betas.IBS[ll,] = unlist(statsAndPVals[[ll]]$beta.IBS)
    
    pVals.Incomp[ll,] = unlist(statsAndPVals[[ll]]$pv.Incomp)
    Scores.Incomp[ll,] = unlist(statsAndPVals[[ll]]$TScores.Incomp)
    Betas.Incomp[ll,] = unlist(statsAndPVals[[ll]]$beta.Incomp)
    
    pVals.AMS[ll,] = unlist(statsAndPVals[[ll]]$pv.AMS)
    Scores.AMS[ll,] = unlist(statsAndPVals[[ll]]$TScores.AMS)
    Betas.AMS[ll,] = unlist(statsAndPVals[[ll]]$beta.AMS)
    
    pVals.BinMM[ll,] = unlist(statsAndPVals[[ll]]$pv.BinMM)
    Scores.BinMM[ll,] = unlist(statsAndPVals[[ll]]$TScores.BinMM)
    Betas.BinMM[ll,] = unlist(statsAndPVals[[ll]]$beta.BinMM)
  }
  
  pVals = cbind(pVals.IBS, pVals.Incomp, pVals.AMS, pVals.BinMM)
  Scores = cbind(Scores.IBS, Scores.Incomp, Scores.AMS, Scores.BinMM)
  Betas = cbind(Betas.IBS, Betas.Incomp, Betas.AMS, Betas.BinMM)
  
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
  write.csv(pVals, file = paste0(path,"/Power_CatY_Prev",YPrev*100,"_LinHypTest_ScoreAndRSNPs_Score",snpOrScore,"_",percentageAssocGene,"Gene_",percentageAssocSNP,"SNPs_OR",ORSize,"_",ld,"_Sim",start,"to",numSims,"_PValues_ForceCovFit.csv"))
  write.csv(Scores, file = paste0(path,"/Power_CatY_Prev",YPrev*100,"_LinHypTest_ScoreAndRSNPs_Score",snpOrScore,"_",percentageAssocGene,"Gene_",percentageAssocSNP,"SNPs_OR",ORSize,"_",ld,"_Sim",start,"to",numSims,"_Stats_ForceCovFit.csv"))
  write.csv(Betas, file = paste0(path,"/Power_CatY_Prev",YPrev*100,"_LinHypTest_ScoreAndRSNPs_Score",snpOrScore,"_",percentageAssocGene,"Gene_",percentageAssocSNP,"SNPs_OR",ORSize,"_",ld,"_Sim",start,"to",numSims,"_Betas_ForceCovFit.csv"))
}

#continuous outcome
RunPowerPipelineLinHypTestCont_GammaTest_Score_BetaNZero = function(chr, gene, numPairs, standardizeScores = FALSE, weightedScores = FALSE, scoreWeights, ORSize, percentageAssocGene, percentageAssocSNP, LowLD, TrueScore, start, stop){
  #function to run whole TIE pipeline, Calculates LRT test stat for Lin Hyp test method and whether p-value <0.05
  # only use when Y is continuous (Normal)
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #standardizeScores = T or F whether the scores should be standardized based on maximum score value
  #weightedScores = T or F whether the scores will be weighted
  #scoreWeights = m x 1 vector of weights, one weight for each SNP
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes TIE values to csv files
  
  library(parallel)
  
  #need to define this for naming at the end
  snpOrScore = TrueScore
  
  #always the same
  numSims = stop
  
  #define effect based on OR size
  if(ORSize == "small"){
    effect = 0.14
  } else if(ORSize == "medium"){
    effect = 0.41
  } else {
    effect = 0.69
  }
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_II_Sims/HapGen_Files/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/ProjectIISourceFunctions_v2.R")
  source("/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/Linear_ADMM0.r")
  
  #determine which SNPs to actually set as assoc. based on gene
  assocSNPs = DetermineAssocRSNPs(gene = gene, LowLD = LowLD, percentageAssoc = percentageAssocSNP)
  
  myList = lapply(start:numSims, rep, times = 1)
  # p-value initialization
  pv.IBS <- rep(0, 3)
  pv.Incomp <- rep(0, 3)
  pv.AMS <- rep(0, 3)
  pv.BinMM <- rep(0, 3)
  
  Tall.IBS <- matrix(0, 1, 3)
  Tall.Incomp <- matrix(0, 1, 3)
  Tall.AMS <- matrix(0, 1, 3)
  Tall.BinMM <- matrix(0, 1, 3)
  
  beta.al.IBS <- matrix(0, 1, 2)
  beta.al.Incomp <- matrix(0, 1, 2)
  beta.al.AMS <- matrix(0, 1, 2)
  beta.al.BinMM <- matrix(0, 1, 2)
  
  statsAndPVals = mclapply(myList, function(ii){
    #define matrix to hold all Stats and Pvalues
    statsAndPVals = matrix(NA, nrow = numSims, ncol = 6)
    
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
    IBS.gene.PercentOfSNPs = calcGeneScorePercentOfSNPs(SingleSNPKernel = IBS.snp, gene = gene, percentageAssoc = percentageAssocGene, LowLD = LowLD, standardize = FALSE, useWeights = FALSE)
    Incomp.gene.PercentOfSNPs = calcGeneScorePercentOfSNPs(SingleSNPKernel = Incomp.snp, gene =  gene, percentageAssoc = percentageAssocGene, LowLD = LowLD, standardize = FALSE, useWeights = FALSE)
    AMS.gene.PercentOfSNPs = calcGeneScorePercentOfSNPs(SingleSNPKernel = AMS.snp, gene =  gene, percentageAssoc = percentageAssocGene, LowLD = LowLD, standardize = FALSE, useWeights = FALSE)
    BinMM.gene.PercentOfSNPs = calcGeneScorePercentOfSNPs(SingleSNPKernel = BinMM.snp, gene =  gene, percentageAssoc = percentageAssocGene, LowLD = LowLD, standardize = FALSE, useWeights = FALSE)
    
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
    
    #set assoc Betas
    #all betas have same effect for now
    for(jj in assocSNPs){
      Betas[jj] = effect
      Betas = as.matrix(Betas, ncol = 1)
    }
    
    Gamma = effect
    
    #generate phenotypes, both continuous and binary
    ContPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = FALSE,  Covariates = CovData, RGenoData = RGenos, ScoreData = PhenoScore, Betas = Betas, Gamma = Gamma)
    
    # define location of zero components
    # basically, this defines what the null hypothesis is, I think
    # so for joint null, we need the non-zero components to be for covariates
    # and the beta and gamma components to be zero
    numCov = dim(CovData)[2]
    numSNPs = dim(RGenos)[2]
    N = c(rep(FALSE, numCov+1+numSNPs), TRUE)
    
    #need to combine the covariates, r genos, and score matrices together
    #need to include column of 1s for intercept?
    intercept = matrix(1,nrow = numPairs, ncol = 1)
    designMat.IBS = cbind(intercept, CovData, RGenos, IBS.gene)
    designMat.Incomp = cbind(intercept, CovData, RGenos, Incomp.gene)
    designMat.AMS = cbind(intercept, CovData, RGenos, AMS.gene)
    designMat.BinMM = cbind(intercept, CovData, RGenos, BinMM.gene)
    
    #then combine the phenos with the design matrix as a list
    Model.IBS = list(X=designMat.IBS, Y=ContPhenos)
    Model.Incomp = list(X=designMat.Incomp, Y=ContPhenos)
    Model.AMS = list(X=designMat.AMS, Y=ContPhenos)
    Model.BinMM = list(X=designMat.BinMM, Y=ContPhenos)
    
    # estimate the uncontrained estimator
    beta.unre.IBS <- cv.SCAD_ADMM_unre(X=Model.IBS$X, Y=Model.IBS$Y, N=N, beta0=rep(0, dim(Model.IBS$X)[2]), err=1e-4, tune="cv", unpen = c(1,2,3))
    indice.unre.IBS <- beta.unre.IBS!=0
    
    beta.unre.Incomp <- cv.SCAD_ADMM_unre(X=Model.Incomp$X, Y=Model.Incomp$Y, N=N, beta0=rep(0, dim(Model.Incomp$X)[2]), err=1e-4, tune="cv", unpen = c(1,2,3))
    indice.unre.Incomp <- beta.unre.Incomp!=0
    
    beta.unre.AMS <- cv.SCAD_ADMM_unre(X=Model.AMS$X, Y=Model.AMS$Y, N=N, beta0=rep(0, dim(Model.AMS$X)[2]), err=1e-4, tune="cv", unpen = c(1,2,3))
    indice.unre.AMS <- beta.unre.AMS!=0
    
    beta.unre.BinMM <- cv.SCAD_ADMM_unre(X=Model.BinMM$X, Y=Model.BinMM$Y, N=N, beta0=rep(0, dim(Model.BinMM$X)[2]), err=1e-4, tune="cv", unpen = c(1,2,3))
    indice.unre.BinMM <- beta.unre.BinMM!=0
    
    # estimate the constrained estimator
    beta.re.IBS <- cv.SCAD_ADMM_re(X=Model.IBS$X, Y=Model.IBS$Y, N=N, beta0=rep(0, dim(Model.IBS$X)[2]), err=1e-4, tune="cv", unpen = c(1,2,3))
    indice.re.IBS <- beta.re.IBS!=0
    
    beta.re.Incomp <- cv.SCAD_ADMM_re(X=Model.Incomp$X, Y=Model.Incomp$Y, N=N, beta0=rep(0, dim(Model.Incomp$X)[2]), err=1e-4, tune="cv", unpen = c(1,2,3))
    indice.re.Incomp <- beta.re.Incomp!=0
    
    beta.re.AMS <- cv.SCAD_ADMM_re(X=Model.AMS$X, Y=Model.AMS$Y, N=N, beta0=rep(0, dim(Model.AMS$X)[2]), err=1e-4, tune="cv", unpen = c(1,2,3))
    indice.re.AMS <- beta.re.AMS!=0
    
    beta.re.BinMM <- cv.SCAD_ADMM_re(X=Model.BinMM$X, Y=Model.BinMM$Y, N=N, beta0=rep(0, dim(Model.BinMM$X)[2]), err=1e-4, tune="cv", unpen = c(1,2,3))
    indice.re.BinMM <- beta.re.BinMM!=0
    
    # estimate the conditional variance
    #    sig2 <- mean((Model$Y-Model$X%*%beta.unre)^2)
    n = numSims
    sig2.IBS <- mean((Model.IBS$Y-Model.IBS$X%*%beta.unre.IBS)^2)*n/(n-sum(beta.unre.IBS!=0))
    sig2.Incomp <- mean((Model.Incomp$Y-Model.Incomp$X%*%beta.unre.Incomp)^2)*n/(n-sum(beta.unre.Incomp!=0))
    sig2.AMS <- mean((Model.AMS$Y-Model.AMS$X%*%beta.unre.AMS)^2)*n/(n-sum(beta.unre.AMS!=0))
    sig2.BinMM <- mean((Model.BinMM$Y-Model.BinMM$X%*%beta.unre.BinMM)^2)*n/(n-sum(beta.unre.BinMM!=0))
    
    # construct the likelihood ratio statistic
    TL.IBS <- sum((Model.IBS$X%*%beta.re.IBS-Model.IBS$Y)^2)-sum((Model.IBS$X%*%beta.unre.IBS-Model.IBS$Y)^2)
    #if LRT stat is less than 0, there is error, so set to -1
    if(TL.IBS < 0){
      TL.IBS = -1
    }
    TL.Incomp <- sum((Model.Incomp$X%*%beta.re.Incomp-Model.Incomp$Y)^2)-sum((Model.Incomp$X%*%beta.unre.Incomp-Model.Incomp$Y)^2)
    #if LRT stat is less than 0, there is error, so set to -1
    if(TL.Incomp < 0){
      TL.Incomp = -1
    }
    TL.AMS <- sum((Model.AMS$X%*%beta.re.AMS-Model.AMS$Y)^2)-sum((Model.AMS$X%*%beta.unre.AMS-Model.AMS$Y)^2)
    #if LRT stat is less than 0, there is error, so set to -1
    if(TL.AMS < 0){
      TL.AMS = -1
    }
    TL.BinMM <- sum((Model.BinMM$X%*%beta.re.BinMM-Model.BinMM$Y)^2)-sum((Model.BinMM$X%*%beta.unre.BinMM-Model.BinMM$Y)^2)
    #if LRT stat is less than 0, there is error, so set to -1
    if(TL.BinMM < 0){
      TL.BinMM = -1
    }
    
    # construct the Wald statistic
    B.IBS = crossprod(Model.IBS$X[,indice.unre.IBS|N], Model.IBS$X[,indice.unre.IBS|N])
    if(rcond(B.IBS) >= 1e-10){
      B_0.IBS <- solve(B.IBS)
      #d_0.IBS should determine which rows to subset
      #gives the first value that is restricted
      d_0.IBS <- dim(B_0.IBS)[1]
      #in case the inverse doesn't exist, need another if else
      A.IBS = B_0.IBS[d_0.IBS,d_0.IBS]
      TW.IBS <- crossprod(beta.unre.IBS[N], solve(A.IBS, beta.unre.IBS[N]))
    } else {
      TW.IBS = -1
    }
    
    B.Incomp = crossprod(Model.Incomp$X[,indice.unre.Incomp|N], Model.Incomp$X[,indice.unre.Incomp|N])
    if(rcond(B.Incomp) >= 1e-10){
      B_0.Incomp <- solve(B.Incomp)
      #d_0.Incomp should determine which rows to subset
      #gives the first value that is restricted
      d_0.Incomp <- dim(B_0.Incomp)[1]
      #in case the inverse doesn't exist, need another if else
      A.Incomp = B_0.Incomp[d_0.Incomp,d_0.Incomp]
      TW.Incomp <- crossprod(beta.unre.Incomp[N], solve(A.Incomp, beta.unre.Incomp[N]))
    } else {
      TW.Incomp = -1
    }
    
    B.AMS = crossprod(Model.AMS$X[,indice.unre.AMS|N], Model.AMS$X[,indice.unre.AMS|N])
    if(rcond(B.AMS) >= 1e-10){
      B_0.AMS <- solve(B.AMS)
      #d_0.AMS should determine which rows to subset
      #gives the first value that is restricted
      d_0.AMS <- dim(B_0.AMS)[1]
      #in case the inverse doesn't exist, need another if else
      A.AMS = B_0.AMS[d_0.AMS,d_0.AMS]
      TW.AMS <- crossprod(beta.unre.AMS[N], solve(A.AMS, beta.unre.AMS[N]))
    } else {
      TW.AMS = -1
    }
    
    B.BinMM = crossprod(Model.BinMM$X[,indice.unre.BinMM|N], Model.BinMM$X[,indice.unre.BinMM|N])
    if(rcond(B.BinMM) >= 1e-10){
      B_0.BinMM <- solve(B.BinMM)
      #d_0.BinMM should determine which rows to subset
      #gives the first value that is restricted
      d_0.BinMM <- dim(B_0.BinMM)[1]
      #in case the inverse doesn't exist, need another if else
      A.BinMM = B_0.BinMM[d_0.BinMM,d_0.BinMM]
      TW.BinMM <- crossprod(beta.unre.BinMM[N], solve(A.BinMM, beta.unre.BinMM[N]))
    } else {
      TW.BinMM = -1
    }
    
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
    
    #LRT
    if (TL.IBS>=qchisq(0.95, df = doF)){
      pv.IBS[1] <- pv.IBS[1]+1/5000
    }
    if (TL.Incomp>=qchisq(0.95, df = doF)){
      pv.Incomp[1] <- pv.Incomp[1]+1/5000
    }
    if (TL.AMS>=qchisq(0.95, df = doF)){
      pv.AMS[1] <- pv.AMS[1]+1/5000
    }
    if (TL.BinMM>=qchisq(0.95, df = doF)){
      pv.BinMM[1] <- pv.BinMM[1]+1/5000
    }
    #Wald
    if (TW.IBS>=qchisq(0.95, df = doF)){
      pv.IBS[2] <- pv.IBS[2]+1/5000
    }
    if (TW.Incomp>=qchisq(0.95, df = doF)){
      pv.Incomp[2] <- pv.Incomp[2]+1/5000
    }
    if (TW.AMS>=qchisq(0.95, df = doF)){
      pv.AMS[2] <- pv.AMS[2]+1/5000
    }
    if (TW.BinMM>=qchisq(0.95, df = doF)){
      pv.BinMM[2] <- pv.BinMM[2]+1/5000
    }
    #Score
    if (TS.IBS>=qchisq(0.95, df = doF)){
      pv.IBS[3] <- pv.IBS[3]+1/5000
    }
    if (TS.Incomp>=qchisq(0.95, df = doF)){
      pv.Incomp[3] <- pv.Incomp[3]+1/5000
    }
    if (TS.AMS>=qchisq(0.95, df = doF)){
      pv.AMS[3] <- pv.AMS[3]+1/5000
    }
    if (TS.BinMM>=qchisq(0.95, df = doF)){
      pv.BinMM[3] <- pv.BinMM[3]+1/5000
    }
    
    print(paste0("Simulation ",ii," is complete."))
    
    Tall.IBS[1, ] <- c(TL.IBS, TW.IBS, TS.IBS)
    beta.al.IBS[1, ] <- c(sum(beta.re.IBS!=0), sum(beta.unre.IBS!=0))
    
    Tall.Incomp[1, ] <- c(TL.Incomp, TW.Incomp, TS.Incomp)
    beta.al.Incomp[1, ] <- c(sum(beta.re.Incomp!=0), sum(beta.unre.Incomp!=0))
    
    Tall.AMS[1, ] <- c(TL.AMS, TW.AMS, TS.AMS)
    beta.al.AMS[1, ] <- c(sum(beta.re.AMS!=0), sum(beta.unre.AMS!=0))
    
    Tall.BinMM[1, ] <- c(TL.BinMM, TW.BinMM, TS.BinMM)
    beta.al.BinMM[1, ] <- c(sum(beta.re.BinMM!=0), sum(beta.unre.BinMM!=0))
    
    statsAndPVals = list(pv.IBS=pv.IBS, TScores.IBS=Tall.IBS, beta.IBS=beta.al.IBS,
                         pv.Incomp=pv.Incomp, TScores.Incomp=Tall.Incomp, beta.Incomp=beta.al.Incomp,
                         pv.AMS=pv.AMS, TScores.AMS=Tall.AMS, beta.AMS=beta.al.AMS,
                         pv.BinMM=pv.BinMM, TScores.BinMM=Tall.BinMM, beta.BinMM=beta.al.BinMM)
    statsAndPVals
  }, mc.cores=4)
  
  n = length(statsAndPVals)
  
  pVals.IBS = pVals.Incomp = pVals.AMS = pVals.BinMM = matrix(nrow = n, ncol = 3)
  Scores.IBS = Scores.Incomp = Scores.AMS = Scores.BinMM = matrix(nrow = n, ncol = 3)
  Betas.IBS = Betas.Incomp = Betas.AMS = Betas.BinMM = matrix(nrow = n, ncol = 2)
  
  for(ll in 1:n){
    pVals.IBS[ll,] = unlist(statsAndPVals[[ll]]$pv.IBS)
    Scores.IBS[ll,] = unlist(statsAndPVals[[ll]]$TScores.IBS)
    Betas.IBS[ll,] = unlist(statsAndPVals[[ll]]$beta.IBS)
    
    pVals.Incomp[ll,] = unlist(statsAndPVals[[ll]]$pv.Incomp)
    Scores.Incomp[ll,] = unlist(statsAndPVals[[ll]]$TScores.Incomp)
    Betas.Incomp[ll,] = unlist(statsAndPVals[[ll]]$beta.Incomp)
    
    pVals.AMS[ll,] = unlist(statsAndPVals[[ll]]$pv.AMS)
    Scores.AMS[ll,] = unlist(statsAndPVals[[ll]]$TScores.AMS)
    Betas.AMS[ll,] = unlist(statsAndPVals[[ll]]$beta.AMS)
    
    pVals.BinMM[ll,] = unlist(statsAndPVals[[ll]]$pv.BinMM)
    Scores.BinMM[ll,] = unlist(statsAndPVals[[ll]]$TScores.BinMM)
    Betas.BinMM[ll,] = unlist(statsAndPVals[[ll]]$beta.BinMM)
  }
  
  pVals = cbind(pVals.IBS, pVals.Incomp, pVals.AMS, pVals.BinMM)
  Scores = cbind(Scores.IBS, Scores.Incomp, Scores.AMS, Scores.BinMM)
  Betas = cbind(Betas.IBS, Betas.Incomp, Betas.AMS, Betas.BinMM)
  
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
  write.csv(pVals, file = paste0(path,"/Power_ContY_LinHypTest_ScoreAndRSNP_Score",snpOrScore,"_",percentageAssocGene,"Gene_",percentageAssocSNP,"SNPs_OR",ORSize,"_",ld,"_Sim",start,"to",numSims,"_PValues_ForceCovs.csv"))
  write.csv(Scores, file = paste0(path,"/Power_ContY_LinHypTest_ScoreAndRSNP_Score",snpOrScore,"_",percentageAssocGene,"Gene_",percentageAssocSNP,"SNPs_OR",ORSize,"_",ld,"_Sim",start,"to",numSims,"_Stats_ForceCovs.csv"))
  write.csv(Betas, file = paste0(path,"/Power_ContY_LinHypTest_ScoreAndRSNP_Score",snpOrScore,"_",percentageAssocGene,"Gene_",percentageAssocSNP,"SNPs_OR",ORSize,"_",ld,"_Sim",start,"to",numSims,"_Betas_ForceCovs.csv"))
}
