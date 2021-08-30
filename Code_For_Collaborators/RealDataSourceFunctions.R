################################################
## Functions to source for real data analyses
## 
## Written by : Victoria Arthur
## Edited: 01/28/2021
################################################

##Functions to calculate gene-based scores
##First calculate individual scores for each pair,
##then calculate score for gene region
 
#IBS mismatch score - SNP based
calcIBSMismatch = function(RGenosMat = matrix(), DGenosMat = matrix()){
  #Function to calculate the IBS mismatch score for all D/R pairs
  #Inputs: A matrix of Recipient genotypes for all individuals (N x m)
  #        A matrix of Donor genotypes for all individuals (N x m)
  #Returns: A matrix of IBS mismatch scores  (N x m)
  
  #number of SNPs is the number of columns in D Genotypes matrix
  ncolPlink = ncol(DGenosMat)
  
  #calculate the difference between the two subjects
  diffsPlink = abs(DGenosMat - RGenosMat)
  
  #######################
  ## IBS Score
  ####################### 
  #if diff = 0, score = 0
  #if diff = 1, score is unchanged
  #if diff = 2, score = 2
  IBSMismatch = diffsPlink
  
  #save IBS scores
  RIBSScoresMat = matrix(unlist(IBSMismatch), ncol = ncolPlink, byrow = F)
  return(RIBSScoresMat)
}

#Incompatibility score - SNP based
calcIncompatibilityScore = function(RGenosMat = matrix(), DGenosMat = matrix()){
  #Function to calculate the incompatibility score for all D/R pairs
  #Inputs: A matrix of Recipient genotypes for all individuals (N x m)
  #        A matrix of Donor genotypes for all individuals (N x m)
  #Returns: A matrix of incompatibility scores  (N x m)
  
  #number of SNPs is the number of columns in D Genotypes matrix
  ncolPlink = ncol(DGenosMat)
  
  #calculate the difference between the two subjects
  diffsPlink = abs(DGenosMat - RGenosMat)
  
  #######################
  ## Incomp Score
  #######################
  #initialize a list of empty dfs with same number of columns as original
  incomp = diffsPlink
  for(ii in 1:ncol(diffsPlink)){
    incomp[diffsPlink[,ii] == 0,ii] = 0
  }
  for(ii in 1:ncol(diffsPlink)){
    incomp[diffsPlink[,ii] != 0,ii] = 1
  }
  
  #save incomp scores
  RIncompScoresMat = matrix(unlist(incomp), ncol = ncolPlink, byrow = F)
  return(RIncompScoresMat)
}

#AMS score - SNP based
calcAMS = function(RGenosMat = matrix(), DGenosMat = matrix()){
  #function to calculate the AMS score for all D/R pairs 
  #Inputs: A matrix of Recipient genotypes for all individuals (N x m)
  #        A matrix of Donor genotypes for all individuals (N x m)
  #Returns: A matrix of AMS scores (N x m)
  
  #number of SNPs is the number of columns in D Genotypes matrix
  ncolPlink = ncol(DGenosMat)
  
  #calculate the difference between the two subjects
  diffsPlink = abs(DGenosMat - RGenosMat)
  
  #######################
  ## AMS
  #######################
  #mismatch if D has allele not in R
  #sum across both alleles in genotype
  #Score is either 0, 1, or 2
  alloMismatch = matrix(0, nrow = nrow(diffsPlink), ncol = ncol(diffsPlink)) #make default value 0
  alloMismatch[(DGenosMat == 0) & (RGenosMat == 2)] = 2 #Donor AA, Recip aa
  alloMismatch[(DGenosMat == 2) & (RGenosMat == 0)] = 2 #Donor aa, Recip AA
  alloMismatch[(DGenosMat == 1) & (RGenosMat == 2)] = 1 #Donor Aa, Recip aa
  alloMismatch[(DGenosMat == 1) & (RGenosMat == 0)] = 1 #Donor Aa, Recip AA
  alloMismatch[is.na(DGenosMat) | is.na(RGenosMat)] = NA #make sure NAs are preserved
  #match row and column names from the original data set
  #row names should be R Ids
  rownames(alloMismatch) = rownames(DGenosMat)
  colnames(alloMismatch) = colnames(DGenosMat)
  
  #save AMS scores
  RAlloMismatchMat = matrix(unlist(alloMismatch), ncol = ncolPlink, byrow = F)
  return(RAlloMismatchMat)
}

#Binary mismatch score - SNP based
calcBinaryMM = function(RGenosMat = matrix(), DGenosMat = matrix()){
  #Function to calculate the binary MM score for all D/R pairs
  #Inputs: A matrix of Recipient genotypes for all individuals (N x m)
  #        A matrix of Donor genotypes for all individuals (N x m)
  #Returns: A matrix of binary MM scores  (N x m)
  
  #number of SNPs is the number of columns in D Genotypes matrix
  ncolPlink = ncol(DGenosMat)
  
  #calculate the difference between the two subjects
  diffsPlink = abs(DGenosMat - RGenosMat)
  
  #######################
  ## Binary  Mismatch
  #######################
  #mismatch if D has allele not in R
  #Score is either 0 or 1
  binMismatch = matrix(0, nrow = nrow(diffsPlink), ncol = ncol(diffsPlink)) #make default value 0
  binMismatch[(DGenosMat == 1) & (RGenosMat == 2)] = 1 #Donor Aa, Recip aa
  binMismatch[(DGenosMat == 0) & (RGenosMat == 2)] = 1 #Donor AA, Recip aa
  binMismatch[(DGenosMat == 2) & (RGenosMat == 0)] = 1 #Donor aa, Recip AA
  binMismatch[(DGenosMat == 1) & (RGenosMat == 0)] = 1 #Donor Aa, Recip AA
  binMismatch[is.na(DGenosMat) | is.na(RGenosMat)] = NA #make sure NAs are preserved
  #match row and column names from the original data set
  #row names should be R Ids
  rownames(binMismatch) = rownames(DGenosMat)
  colnames(binMismatch) = colnames(DGenosMat)
  
  #save Binary mismatch scores
  RBinMismatchMat = matrix(unlist(binMismatch), ncol = ncolPlink, byrow = F)
  return(RBinMismatchMat)
}

#Then use any of the single SNP kernel functions to calculate gene-based score
calcGeneScore = function(SingleSNPKernel = matrix(), standardize = FALSE, useWeights = FALSE, weights){
  #Function to calculate the gene based score based on the single SNP kernels and optional weights
  #Inputs: SingleSNPKernel: A m x N vector of single SNP scores (either IBS, Incompatitibility, AMS or Bin MM)
  #         standardize: True or False for if data should be standardized by maximum score value
  #                       if TRUE, then the summed gene score is divided by the number of SNPs (for 0,1 scores)
  #                       or the number of SNPs * 2 (for 0,1,2 scores) (AKA divided by maximum summed score value possible)
  #         useWeights: True or False for if the individual scores are weighted or unweighted, if TRUE, weights must be specified 
  #         weights: A m x 1 vector of weights, allows individual SNPs to be weighted differently if wanted
  #Returns: N x 1 vector, with single gene score for each individual
  
  if(useWeights == FALSE){
    #unweighted sum of kernel
    if(standardize == FALSE){
      #then no division by anything, just a simple sum of kernels
      #sum over all the SNPs
      geneScore = rowSums(SingleSNPKernel, na.rm = TRUE)
    } 
    else {
      #unweighted sum but standardized by maximum score value
      #sum over all the SNPs
      geneScore.raw = rowSums(SingleSNPKernel, na.rm = TRUE)
      #figure out if we need to multiply number of SNPs by 2
      #if the max score is larger than the total number of SNPs,
      #we know the score can go up to 2, so we multiply denominator by 2
      if(max(geneScore.raw) > ncol(SingleSNPKernel)){
        geneScore = geneScore.raw/(2*ncol(SingleSNPKernel))
      } else {
        #otherwise just divide by total number of SNPs
        geneScore = geneScore.raw/(ncol(SingleSNPKernel))
      }
    }
  }
  #otherwise we are using weights
  else {
    #need the sum of all weights
    weightTotal = sum(weights, na.rm = TRUE)
    #need to multiply kernel by weight and sum
    #make weights a vector (m x 1)
    weights.vec = as.matrix(weights)
    sum.w.Kernel = SingleSNPKernel %*% weights.vec
    geneScore = sum.w.Kernel/weightTotal
  }
  geneScore.mat = as.matrix(geneScore)
  return(geneScore.mat)
}

#Calculate Scores for either R genotypes OR genetic matching score
#Calculate individual score values (for either R geno or mismatch score)
CalcUScore = function(SampleSize, includeCov = FALSE, CovData, CalcUR = TRUE, RGenoData, ScoreData, Phenos, BinPhenos = TRUE){
  
  #define variables
  #SampleSize is number of D/R pairs (N)
  #includeCov: T or F whether or not covariates are included in the modeling
  #CovData: Covariate matrix (N x K) if covariates are included
  #CalcUR: T if we are calculating U score for R geno, F if we are calculating U for gene-based score
  #RGenoData: Matrix of R genotype data (N x m)
  #ScoreData: Vector of Gene-based score data (N x 1)
  #Phenos: Vector of generated phenotypes from either GenNullPhenos or GenAltPhenos
  #BinPhenos: T if phenotypes are binary, F is they are continuous
  
  #define variables
  ss = SampleSize
  if(includeCov == TRUE){
    W = CovData
    K = ncol(W)
  } else {
    K = 0
    W = matrix(1, nrow = ss, ncol = 1)
  }
  if(CalcUR == TRUE){
    #then we have R geno data
    XZ = RGenoData
  } else {
    #otherwise we have score data
    XZ = ScoreData
  }
  Y = Phenos
  
  if(BinPhenos == TRUE){
    nulllogitreg=glm(Y~W,family=binomial)
    # p: expected value of Y from the logistic regression
    p1=fitted(nulllogitreg)
    p0 = 1 - p1
    
    # X or Z hat: expected genotype from weighted linear regression
    XZhat=XZ
    for(j in 1:ncol(XZ)){
      xz=XZ[,j]
      linearreg=lm(xz~W,weights=p1*p0) # weighted linear regression
      XZhat[,j]=fitted(linearreg)
    }
    
    U_score_XZhat = as.matrix(t(XZ-XZhat)%*%(Y-p1))
  }  else {
    #phenotypes are continuous, so slightly different approach
    # X or Z hat: expected genotype from unweighted linear regression
    XZhat=XZ
    for(j in 1:ncol(XZ)){
      xz=XZ[,j]
      linearreg=lm(xz~W) # unweighted linear regression
      XZhat[,j]=fitted(linearreg)
    }
    
    #need expected value of Ys
    nulllinearreg = glm(Y~W,family=gaussian)
    # expected value of Y from the linear regression
    Yhat = fitted(nulllinearreg)
    
    U_score_XZhat = as.matrix(t(XZ-XZhat)%*%(Y-Yhat))
  }
  return(U_score_XZhat)
}

CalcQValues = function(SampleSize, includeCov = FALSE, CovData, CalcUR = TRUE, RGenoData, ScoreData, Phenos, BinPhenos = TRUE){
  
  #define variables
  #SampleSize is number of D/R pairs (N)
  #includeCov: T or F whether or not covariates are included in the modeling
  #CovData: Covariate matrix (N x K) if covariates are included
  #CalcUR: T if we are calculating U for R geno, F if we are calculating U for Score
  #RGenoData: Matrix of R genotype data (N x m)
  #ScoreData: Vector of Gene-based score data (Nx1)
  #Phenos: Vector of generated phenotypes from either GenNullPhenos or GenAltPhenos
  #BinPhenos: T if phenotypes are binary, F is they are continuous
  
  #output should be matrix of dim n x m for QR, and n x 1 for QS
  
  #define variables
  ss = SampleSize
  if(includeCov == TRUE){
    W = CovData
    K = ncol(W)
  } else {
    K = 0
    W = matrix(1, nrow = ss, ncol = 1)
  }
  if(CalcUR == TRUE){
    #then we have R geno data
    XZ = RGenoData
    m = ncol(XZ)
  } else {
    #otherwise we have score data
    XZ = ScoreData
    m = 1
  }
  Y = Phenos
  
  if(BinPhenos == TRUE){
    nulllogitreg=glm(Y~W,family=binomial)
    # p: expected value of Y from the logistic regression
    p1=fitted(nulllogitreg)
    p0 = 1 - p1
    
    # X or Z hat: expected genotype from weighted linear regression
    XZhat=XZ
    for(j in 1:ncol(XZ)){
      xz=XZ[,j]
      linearreg=lm(xz~W,weights=p1*p0) # weighted linear regression
      XZhat[,j]=fitted(linearreg)
    }
    
    #define Q matrix
    Q=matrix(ncol=m,nrow=ss)
    
    #populate Q differently for R geno vs Score
    for (i in 1:nrow(Q)){
      for (j in 1:ncol(Q)){
        Q[i,j]=(XZ[i,j]-XZhat[i,j])*(Y[i]-p1[i]) 
      }
    }
  }
  else{
    #phenotypes are continuous, so slightly different approach
    # X or Z hat: expected genotype from unweighted linear regression
    XZhat=XZ
    for(j in 1:ncol(XZ)){
      xz=XZ[,j]
      linearreg=lm(xz~W) # unweighted linear regression
      XZhat[,j]=fitted(linearreg)
    }
    
    #need expected value of Ys
    nulllinearreg = glm(Y~W,family=gaussian)
    # expected value of Y from the linear regression
    Yhat = fitted(nulllinearreg)
    
    #define Q matrix
    Q=matrix(ncol=m,nrow=ss)
    
    for (i in 1:nrow(Q)){
      for (j in 1:ncol(Q)){
        Q[i,j]=(XZ[i,j]-XZhat[i,j])*(Y[i]-Yhat[i]) 
      }
    }
  }
  return(Q)
}

#Calculate original variance
CalcVariance = function(SampleSize, QValues){
  #Calculates original variance/cov matrix for score
  #define variables
  #SampleSize is number of D/R pairs (N)
  #QValues is a matrix of combined Q values for R geno and Score
  ### dim should be n x (m+1)
  Q = QValues
  ss = SampleSize
  #variance is nQ'Q
  V = ss*t(Q)%*%Q
  #return variance matrix
  return(V)
}

#Calc final statistic and p-value
CalcStatisticPVal = function(SampleSize, Variance, UscoresR, UscoreS, s){
  #define variables
  #SampleSize is number of D/R pairs (N)
  #Variance is the original Var/Cov matrix for the scores calculated  using CalcVariance()
  #### dim should be (m+1) x (m+1)
  #UscoresR is the m x 1 vector of scores for the r genos
  #UscoreS should be the 1x1 score value for the gene-based score
  #s is the number of PCs we want to keep
  
  ss = SampleSize
  VFull = Variance
  UR = UscoresR
  US = UscoreS
  s = s
  
  #Decompose the VFull matrix
  m = dim(VFull)[1] - 1
  VR = VFull[1:m, 1:m]
  CRS = matrix(VFull[1:m, m+1], nrow = m, ncol = 1)
  CSR = matrix(VFull[m+1, 1:m], nrow = 1, ncol = m)
  VS = VFull[m+1, m+1]
  
  #Eigen Decomp of VR
  A = eigen(VR)$vectors
  lambdas = eigen(VR)$values
  
  # keep PCs that explain s% of the variance
  # eigen_percent = cumsum(lambdas/sum(lambdas))
  #determine if s is too small
  #if s is smaller than smallest PVE, use 1 PC
  # if(eigen_percent[1] >= s){
    # num_pc = 1
  # } else {
    # num_pc = which(eigen_percent >= s)[1]
  # }
  
  num_pc = s #number of PCs is directly equal to s
  
  A_s = A[,1:num_pc]
  lambda_s = lambdas[1:num_pc]
  
  #define UPR vector
  UPR = matrix(nrow = num_pc, ncol = 1)
  #if num_pc = 1, then this errors, need to split into cases
  if(num_pc == 1){
    UPR = t(UR) %*% A_s %*% (1/sqrt(lambda_s))
  } else {
    for(i in 1:num_pc){
      UPR[i,] = t(UR) %*% A_s[,i] %*% (1/sqrt(lambda_s[i]))
    }
  }
  
  #construct UP vector
  UP = rbind(UPR, US)
  UP.t = t(UP)
  
  #construct new variance/cov matrix
  Ident = diag(1, nrow = num_pc, ncol = num_pc)
  
  Cov_UPR_US = matrix(nrow = num_pc, ncol = 1)
  #if num_pc = 1, then this errors, need to split into cases
  if(num_pc == 1){
    Cov_UPR_US = (1/sqrt(lambda_s)) * t(A_s) %*% CRS
  } else {
    for(i in 1:num_pc){
      Cov_UPR_US[i,] = (1/sqrt(lambda_s[i])) * t(A_s[,i]) %*% CRS 
    }
  }
  Cov_US_UPR = t(Cov_UPR_US)
  
  LHS = rbind(Ident, Cov_US_UPR)
  RHS = rbind(Cov_UPR_US, VS)
  
  NewVar = cbind(LHS, RHS)
  #have to increase tolerance or else it says matrix singular for binary phenos data
  NewVarInv = solve(NewVar) #, tol = 1.0e-30
  
  finalStat = ss * (UP.t %*% NewVarInv %*% UP)
  
  #calculate p value
  #stat should be dist as chi-sq s+1
  pval = 1 - pchisq(finalStat, df = num_pc + 1)
  
  finalOutput = cbind(finalStat, pval)
  return(finalOutput)
}
