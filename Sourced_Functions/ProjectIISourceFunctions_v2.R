############################################
## Source functions for Project II       ###
############################################
###################################   
DetermineAssocRSNPs = function(gene, LowLD, percentageAssoc){
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
  return(assocSNPs)
}  
##Step 1) Generate paired genotype data, get this data into correct format
#pull R genotype information
obtainRGenotypes = function(chr = c(), numSamples = c(), simNum = c(), gene = "", path=paste0("/path/to/data/")){
  #function of obtain Recipient genotypes for D/R transplant pairs
  #output is a matrix of R genotypes for all individuals (N x m)
  #load needed packages
  suppressMessages(require(ARTP2,lib.loc='/home/vlynn/R/library'))
  suppressMessages(require(dplyr))
  
  #make sure number of samples is integer value
  numSamples = as.integer(numSamples) 
  
  #setwd to location of plink data
  setwd(path)
  
  #make lists of names of plink data
  #for HapGen generated data
  bedFile = paste0(gene,"_",numSamples,"Pairs_SimNum",simNum,"_Subset.bed")
  bimFile = paste0(gene,"_",numSamples,"Pairs_SimNum",simNum,"_Subset.bim")
  famFile = paste0(gene,"_",numSamples,"Pairs_SimNum",simNum,"_Subset.fam")
  
  #read in plink binary files
  plinkFile = read.bed(bed = bedFile, bim = bimFile, fam = famFile)
  
  #matched pairs are even and odd columns of each file
  #define variables for number of columns and rows of the df
  nrowPlink = nrow(plinkFile) #this is number of generated subjects (2*SS)
  ncolPlink = ncol(plinkFile) #this is number of SNPs
  
  recipientGenotypes = plinkFile[seq(2,nrowPlink,2),]
  RGenosMat = matrix(unlist(recipientGenotypes), ncol = ncolPlink, byrow = F)
  
  return(RGenosMat)
}
#pull D genotype information
obtainDGenotypes = function(chr = c(), numSamples = c(), simNum = c(), gene = "", path=paste0("/path/to/data/")){
  #function of obtain Donor genotypes for D/R transplant pairs
  #output is a matrix of D genotypes for all individuals (m x N)
  #load needed packages
  suppressMessages(require(ARTP2,lib.loc='/home/vlynn/R/library'))
  suppressMessages(require(dplyr))
  
  #make sure number of samples is integer value
  numSamples = as.integer(numSamples) 
  
  #setwd to location of plink data
  setwd(path)
  
  #make lists of names of plink data
  #for HapGen generated data
  bedFile = paste0(gene,"_",numSamples,"Pairs_SimNum",simNum,"_Subset.bed")
  bimFile = paste0(gene,"_",numSamples,"Pairs_SimNum",simNum,"_Subset.bim")
  famFile = paste0(gene,"_",numSamples,"Pairs_SimNum",simNum,"_Subset.fam")
    
  #read in plink binary files
  plinkFile = read.bed(bed = bedFile, bim = bimFile, fam = famFile)
  
  #matched pairs are even and odd columns of each file
  #define variables for number of columns and rows of the df
  nrowPlink = nrow(plinkFile) #this is number of generated subjects (2*SS)
  ncolPlink = ncol(plinkFile) #this is number of SNPs
  
  #Ds are odd numbered rows, Rs are even numbered rows
  donorGenotypes = plinkFile[seq(1,nrowPlink,2),]
  DGenosMat = matrix(unlist(donorGenotypes), ncol = ncolPlink, byrow = F)
 
  return(DGenosMat)
  
}

###################################       
##Step 2) Calculate individual scores for each pair,
##then calc score for gene region
#use R and D genotypes to calc IBS mismatch score
calcIBSMismatch = function(RGenosMat = matrix(), DGenosMat = matrix()){
  #function to calculate the IBS mismatch score for all D/R pairs
  #returns a matrix of IBS mismatch scores  (m x N)
  
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
#use R and D genotypes to calc Incompatibility score
calcIncompatibilityScore = function(RGenosMat = matrix(), DGenosMat = matrix()){
  #function to calculate the incompatibility score for all D/R pairs
  #returns a matrix of incompatibility scores  (m x N)
  
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
#use R and D genotypes to calc AMS score
calcAMS = function(RGenosMat = matrix(), DGenosMat = matrix()){
  #function to calculate the AMS score for all D/R pairs 
  #returns a matrix of AMS scores (m x N)
  
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
#use R and D genotypes to calc binary mismatch score
calcBinaryMM = function(RGenosMat = matrix(), DGenosMat = matrix()){
  #function to calculate the binary MM score for all D/R pairs
  #returns a matrix of binary MM scores  (m x N) 
  
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

#use single SNP kernel functions to calc gene-based score
calcGeneScore = function(SingleSNPKernel = matrix(), standardize = FALSE, useWeights = FALSE, weights){
  #function to calculate the gene based score based on the single SNP
  #kernels and optional weights
  #Returns 1 x N vector, with single gene score for each individual
  if(useWeights == FALSE){
  #unweighted sum of kernel
    if(standardize == FALSE){
    #then no division by anything, just a simple sum of kernels
    #sum over all the SNPs
      geneScore = rowSums(SingleSNPKernel)
    } 
    else {
      #unweighted sum but standardized by maximum score value
      #sum over all the SNPs
      geneScore.raw = rowSums(SingleSNPKernel)
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
    #need to determine number of weights to pull
    nSNP = dim(SingleSNPKernel)[2]
    
    #pull nSNP weights from the main weight vector
    weights_subset = weights[1:nSNP]
    
    #need the sum of all weights
    weightTotal = sum(weights_subset)
    #need to multiply kernel by weight and sum
    #make weights a vector (m x 1)
    weights.vec = as.matrix(weights_subset)
    sum.w.Kernel = SingleSNPKernel %*% weights.vec
    geneScore = sum.w.Kernel/weightTotal
  }
  geneScore.mat = as.matrix(geneScore)
  return(geneScore.mat)
}

#use single SNP kernel functions to calc gene-based score, using percent of SNPs
calcGeneScorePercentOfSNPs = function(SingleSNPKernel = matrix(), gene =  "", percentageAssoc = 100, LowLD = TRUE, standardize = FALSE, useWeights = FALSE, weights){
  #function to calculate the gene based score based on the single SNP
  #Only uses a percentage of SNPs in the gene region, assuming not all SNPs are associated with outcome
  #kernels and optional weights
  #Returns 1 x N vector, with single gene score for each individual
  #check gene
  if(gene == "NAT2"){
    if(LowLD == TRUE){
      #then looking at SNPs in low LD
      if(percentageAssoc == 5){
        #then 1 SNP is associated
        SingleSNPKernel = as.matrix(SingleSNPKernel[,1],ncol=1)
      } else if(percentageAssoc == 15){
		#then 2 SNPs are associated
		SingleSNPKernel = SingleSNPKernel[,c(1,2)]
	  } else if(percentageAssoc == 25){
        #then 5 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(1,2,9,7,10)]
      } else if(percentageAssoc == 50){
        #then 7 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(1,2,9,7,10,6,8)]
      } else if(percentageAssoc == 75){
        #then 11 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(1,2,9,7,10,6,8,11,12,3,5)]
      } else {
        #otherwise use all SNPs
        SingleSNPKernel = SingleSNPKernel
      }
    } else {
      #looking at SNPs in high LD 
      if(percentageAssoc == 5){
        #then 1 SNP is associated
        SingleSNPKernel = as.matrix(SingleSNPKernel[,13],ncol=1)
      } else if(percentageAssoc == 15){
		#then 2 SNPs are associated
		SingleSNPKernel = SingleSNPKernel[,c(13,14)]
	  } else if(percentageAssoc == 25){
        #then 5 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(13,14,4,5,3)]
      } else if(percentageAssoc == 50){
        #then 7 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(11,12,3,5,4,14,13)]
      } else if(percentageAssoc == 75){
        #then 11 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(7,10,6,8,11,12,3,5,4,14,13)]
      } else {
        #otherwise use all SNPs
        SingleSNPKernel = SingleSNPKernel
      }
    }
  } else if(gene == "CHI3L2"){
    if(LowLD == TRUE){
      #then looking at SNPs in low LD
      if(percentageAssoc == 5){
        #then 2 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(12,2)]
      } else if(percentageAssoc == 15){
		#then 5 SNPs are associated
		SingleSNPKernel = SingleSNPKernel[,c(12,2,3,9,8)]
	  } else if(percentageAssoc == 25){
        #then 8 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(12,2,3,9,8,24,28,19)]
      } else if(percentageAssoc == 50){
        #then 17 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(12,2,3,9,8,24,28,19,21,22,20,11,1,10,16,4,31)]
      } else if(percentageAssoc == 75){
        #then 25 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(12,2,3,9,8,24,28,19,21,22,20,11,1,10,16,4,31,26,25,32,14,27,33,30,7)]
      } else {
        #otherwise use all SNPs
        SingleSNPKernel = SingleSNPKernel
      }
    } else {
      #looking at SNPs in high LD 
      if(percentageAssoc == 5){
        #then 2 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(17,29)]
      } else if(percentageAssoc == 15){
		#then 5 SNPs are associated
		SingleSNPKernel = SingleSNPKernel[,c(5,6,23,18,13)]
	  } else if(percentageAssoc == 25){
        #then 8 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(5,6,23,18,13,15,17,29)]
      } else if(percentageAssoc == 50){
        #then 17 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(31,26,25,32,14,27,33,30,7,5,6,23,18,13,15,17,29)]
      } else if(percentageAssoc == 75){
        #then 25 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(21,22,20,11,1,10,16,4,31,26,25,32,14,27,33,30,7,5,6,23,18,13,15,17,29)]
      } else {
        #otherwise use all SNPs
        SingleSNPKernel = SingleSNPKernel
      }
    }
  } else { 
    #otherwise we have ASAH1
    if(LowLD == TRUE){
      #then looking at SNPs in low LD
      if(percentageAssoc == 5){
        #then 2 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(31, 9)]
      } else if(percentageAssoc == 15){
		#then 6 SNPs are associated
		SingleSNPKernel = SingleSNPKernel[,c(31,9,4,8,10,40)]
	  } else if(percentageAssoc == 25){
        #then 10 SNPs are associated
		if(dim(SingleSNPKernel)[2] < 40){
			SingleSNPKernel = SingleSNPKernel[,c(31,9,4,8,10,5,3,2,7,1)]
		} else {
			SingleSNPKernel = SingleSNPKernel[,c(31,9,4,8,10,40,5,3,2,7)]
		}
      } else if(percentageAssoc == 50){
        #then 20 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(31,9,4,8,10,40,5,3,2,7,1,35,21,22,11,6,15,12,13,14)]
      } else if(percentageAssoc == 75){
        #then 30 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(31,9,4,8,10,40,5,3,2,7,1,35,21,22,11,6,15,12,13,14,18,39,16,19,20,24,17,26,23,27)]
      } else {
        #otherwise use all SNPs
        SingleSNPKernel = SingleSNPKernel
      }
    } else {
      #looking at SNPs in high LD 
      if(percentageAssoc == 5){
        #then 2 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(25, 32)]
      } else if(percentageAssoc == 15){
		#then 6 SNPs are associated
		SingleSNPKernel = SingleSNPKernel[,c(30,33,34,38,25,32)]
	  } else if(percentageAssoc == 25){
        #then 10 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(27, 37, 29, 36, 28, 30, 33, 34, 38, 25, 32)]
      } else if(percentageAssoc == 50){
        #then 20 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(18,39,16,19,20,24,17,26,23,27,37,29,36,28,30,33,34,38,25,32)]
      } else if(percentageAssoc == 75){
        #then 30 SNPs are associated
        SingleSNPKernel = SingleSNPKernel[,c(1,35,21,22,11,6,15,12,13,14,18,39,16,19,20,24,17,26,23,27,37,29,36,28,30,33,34,38,25,32)]
      } else {
        #otherwise use all SNPs
        SingleSNPKernel = SingleSNPKernel
      }
    }
  }
  
  if(useWeights == FALSE){
    #unweighted sum of kernel
    if(standardize == FALSE){
      #then no division by anything, just a simple sum of kernels
      geneScore = rowSums(SingleSNPKernel)
    } else {
      #unweighted sum but standardized by maximum score value
      #sum over all the SNPs
      geneScore.raw = rowSums(SingleSNPKernel)
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
  } else {
    #otherwise we are using weights
    #need to determine number of weights to pull
    nSNP = dim(SingleSNPKernel)[2]
    
    #pull nSNP weights from the main weight vector
    weights_subset = weights[1:nSNP]
    
    #need the sum of all weights
    weightTotal = sum(weights_subset)
    #need to multiply kernel by weight and sum
    #make weights a vector (m x 1)
    weights.vec = as.matrix(weights_subset)
    sum.w.Kernel = SingleSNPKernel %*% weights.vec
    geneScore = sum.w.Kernel/weightTotal
  }
  geneScore.mat = as.matrix(geneScore)
  return(geneScore.mat)
}

###################################       
##Step 3) Generate covariates
#generate covariate data
GenCovData = function(SampleSize, BinaryValues, ContinuousValues){
  #####################################
  # Simulate covariates
  # returns a matrix of size N x K (N = sample size)
  # where K is number of covariates generated
  # allows input of number of binary and
  # continuous covariate values
  #####################################
  # Binary covariates are drawn from a random
  # binimial distribution with p=0.5
  # Continuous covariate sare drawn from a standard
  # Normal distribution
  
  K = BinaryValues + ContinuousValues
  ss = SampleSize
  
  #define W matrix
  W = matrix(NA,nrow = ss, ncol = K)
  
  for(ii in 1:BinaryValues){
    W[,ii] = rbinom(ss,1,0.5)
  }
  
  for(ii in (BinaryValues+1):K){
    W[,ii] = rnorm(ss, mean = 0, sd = 1)
  }
  
  return(W)
}

###################################       
##Step 4) Generate true phenotype data using covariates only (TIE)
## or using covariates, R genotype, and score (power)
#generate phenotype data - null model, with covariates
GenNullPhenos = function(SampleSize, includeCov = FALSE, YCat = TRUE, YPrev,  Covariates){
  # generates null phenotypes to be used as "true" Y values
  # returns N x 1 vector of phenotypes
  # SampleSize (N) is number of D/R pairs
  # includeCov is FALSE if no covariates are being used, TRUE if they are being used
  # YCat = TRUE if outcome is categorical, otherwise Y is continuous
  # YPrev is the probability that Y=1, only needed if YCat=TRUE
  # Covariates is N x K matrix of covariate values (if they are being used)
  
  ss = SampleSize
  #only define W and eff sizes if we have covariates
  if(includeCov == TRUE){
    W = Covariates
    K = ncol(W)
    #define all effect sizes as 0.5, since this is what most people do in 
    #similar papers for TIE simulations
    eff_sizes = matrix(0.5, nrow = K, ncol = 1)
  }
  if(YCat == TRUE){
    prev = YPrev
  }
  
  #if not including covariates
  if(includeCov == FALSE){
    #if Y is catergorical, and no covariates
    # then Y is pulled from a binom(prev) for each individual
    #logit (pi) = alpha_0
    if(YCat == TRUE){
      nullPhenos = rbinom(ss, 1, prev)
    } 
    #otherwise Y is continuous, and no covs
    #then Y is pulled from a standard normal for each individual
    # basically, Y = epsilon (random error)
    else {
      nullPhenos = rnorm(ss, 0, 1)
    }
  }
  # otherwise we are including covariates
  else{
    #if Y is categorical, and covariates
    #then logit Pr(Y=1) = a0 + a*Covariates
    if(YCat == TRUE){
      #alpha 1 
      alphas = W %*% eff_sizes
      #set alpha0 value to start
      alpha_0 = rep(1,ss)
      #calc alpha_0 + alphas
      lin_predict = alpha_0 + alphas
      #calculate p values for binomials, want colMeans of this to be around prevalence
      prob_case = exp(lin_predict)/(1 + exp(lin_predict))
      #change alpha_0 such that prob of case is close to prevalence
      while((mean(prob_case) > prev) == TRUE) {
        #set alpha_0 value such that prob of case is approx prevalence
        alpha_0 = alpha_0 - 0.01
        #calc alpha_0 + alphas
        lin_predict = alpha_0 + alphas
        #calculate p values for binomials, want colMeans of this to be around prevalence
        prob_case = exp(lin_predict)/(1 + exp(lin_predict))
        #use binomial dist to get phenos
        nullPhenos = rbinom(ss,1,prob_case)
      } 
    }
    #otherwise Y is continuous 
    else {
        #alphas 
        alphas = W %*% eff_sizes
        #epsilon
        epsilon = rnorm(ss,0,1)
        #linear model, added in normally dist error
        lin_predict = alphas + epsilon
        #add to null phenos
        nullPhenos = lin_predict
      }
  }
  nullPhenos = as.matrix(nullPhenos, nrow = ss, ncol = 1)
  return(nullPhenos)
}

#generate phenotype data - alt model, with covariates
GenAltPhenos = function(SampleSize, includeCov = FALSE, YCat = TRUE, YPrev, Covariates, RGenoData = matrix(), ScoreData = matrix(), Betas = matrix(), Gamma = c()){
  # generates "true" phenotypes for power analysis using the alternative
  # hypothesis model
  # returns N x 1 vector of phenotypes
  
  # SampleSize: number of D/R pairs
  # includeCov: TRUE or FALSE for if covariates are included in the modeling
  # YCat: TRUE of FALSE for whether Y is categorical or continuous
  # YPrev: if Y is categorical, the p(y=1)
  # Covariates is K x N matrix of covariate values
  # RGenoData is N x m matrix of R genotype values
  # ScoreData is N x 1 vector of gene-based score values
  # Betas: cov effect sizes for R geno, could all be 0 (m x 1)
  # Gamma: cov effect size for score, could all be 0 (1 x 1)
  
  ss = SampleSize
  #only define W and eff sizes if we have covariates
  if(includeCov == TRUE){
    W = Covariates
    K = ncol(W)
    #define all effect sizes as 0.5, since this is what most people do in 
    #similar papers for TIE simulations
    eff_sizes = matrix(0.5, nrow = K, ncol = 1)
  }
  if(YCat == TRUE){
    prev = YPrev
  }
  
  X = RGenoData
  B = as.matrix(Betas, ncol=1)
  Z = ScoreData
  G = matrix(Gamma, nrow = 1, ncol = 1)
  
  #if not including covariates
  if(includeCov == FALSE){
    #if Y is catergorical, and no covariates
    # then Y is calculated using only XBeta and ZGamma values
    #logit (pi) = alpha_0 + XBeta + ZGamma
    if(YCat == TRUE){
      #beta terms
      betas = X %*% B
      #set alpha0 value to start
      alpha_0 = rep(1,ss)
      #gamma term
      gamma = ScoreData %*% Gamma
      #calc alpha_0 + Betas*SNPGeno + Gamma*Score
      lin_predict = alpha_0 + betas + gamma
      #calculate p values for binomials, want colMeans of this to be around prevalence
      prob_case = exp(lin_predict)/(1 + exp(lin_predict))
      #change beta_0 such that prob of case is close to prevalence
      while((mean(prob_case) > prev) == TRUE) {
        #set beta0 value such that prob of case is approx prevalence
        alpha_0 = alpha_0 - 0.01
        #calc alpha_0 + Betas*SNPGeno + Gamma*Score
        lin_predict = alpha_0 + betas + gamma
        #calculate p values for binomials, want colMeans of this to be around prevalence
        prob_case = exp(lin_predict)/(1 + exp(lin_predict))
        #use binomial dist to get phenos
        altPhenos = rbinom(ss,1,prob_case)
      } 
    } else {
	#otherwise Y is continuous, and no covs
    #then Y is pulled from a standard normal for each individual
    # basically, Y = epsilon (random error)
      #beta terms 
      betas = X %*% B
      #gamma term
      gamma = ScoreData %*% Gamma
      #epsilon
      epsilon = rnorm(ss,0,1)
      #linear model, added in normally dist error
      lin_predict = betas + gamma + epsilon
      #add to null phenos
      altPhenos = lin_predict
    }
  } else { 
  # otherwise we are including covariates
    #if Y is categorical, and covariates
    #then logit Pr(Y=1) = a0 + a*Covariates + XBeta + ZGamma
    if(YCat == TRUE){
      #alpha terms
      alphas = W %*% eff_sizes
      #beta terms
      betas = X %*% B
      #set alpha0 value to start
      alpha_0 = rep(1,ss)
      #gamma term
      gamma = ScoreData %*% G
      #calc alpha_0 + Betas*SNPGeno + Gamma*Score
      lin_predict = alpha_0 + alphas + betas + gamma
      #calculate p values for binomials, want colMeans of this to be around prevalence
      prob_case = exp(lin_predict)/(1 + exp(lin_predict))
      #change beta_0 such that prob of case is close to prevalence
      while((mean(prob_case) > prev) == TRUE) {
        #set beta0 value such that prob of case is approx prevalence
        alpha_0 = alpha_0 - 0.01
        #calc alpha_0 + Betas*SNPGeno + Gamma*Score
        lin_predict = alpha_0 + alphas + betas + gamma
        #calculate p values for binomials, want colMeans of this to be around prevalence
        prob_case = exp(lin_predict)/(1 + exp(lin_predict))
        #use binomial dist to get phenos
        altPhenos = rbinom(ss,1,prob_case)
      } 
    } else {
	  #otherwise Y is continuous 
      #alpha terms
      alphas = W %*% eff_sizes
      #beta terms 
      betas = X %*% B
      #gamma term
      gamma = ScoreData %*% G
      #epsilon
      epsilon = rnorm(ss,0,1)
      #linear model, added in normally dist error
      lin_predict = alphas + betas + gamma + epsilon
      #add to null phenos
      altPhenos = lin_predict
    }
  }
  altPhenos = as.matrix(altPhenos, nrow = ss, ncol = 1)
  return(altPhenos)
}

###################################       
##Step 5) Calculate Scores for either R genotypes OR genetic matching score
#calculate individual score values (for either R geno or mismatch score)
CalcUScore = function(SampleSize, includeCov = FALSE, CovData, CalcUR = TRUE, RGenoData, ScoreData, Phenos, BinPhenos = TRUE){
  
  #SampleSize is number of D/R pairs (N)
  #includeCov: T or F whether or not covariates are included in the modeling
  #CovData: Covariate matrix (N x K) if covariates are included
  #CalcUR: T if we are calculating U for R geno, F if we are calculating U for Score
  #RGenoData: Matrix of R genotype data (N x m)
  #ScoreData: Vector of Gene-based score data (Nx1)
  #Phenos: Vector of generated phenotypes from either GenNullPhenos or GenAltPhenos
  #BinPhenos: T if phenotypes are binary, F is they are continuous
  
  #define variables
  ss = SampleSize
  if(includeCov == TRUE){
    W = CovData
    K = ncol(W)
  } else {
    K = 0
  }
  if(CalcUR == TRUE){
    #then we have R geno data
    XZ = RGenoData
  } else {
    #otherwise we have score data
    XZ = ScoreData
  }
  Y = Phenos
  
  #######################################################
  #For binary phenotypes
  if(BinPhenos == TRUE){
    #need to calculate p_hat (predicted prob of Y=1 under H0)
    nulllogitreg = glm(Y~W,family=binomial)
    # p: expected value of Y from the logistic regression
    p1 = fitted(nulllogitreg)
    p0 = 1 - p1
    
    #calculate U stat using Taylor expanded equation
    #we are summing over N (rows)
    #define (1,W) as N x (k+1) vector
    if(includeCov == FALSE){
      #if no covariates, then (1,W) will be a N x 1 vector of 1s
      OneW = matrix(1, nrow = nrow(Y), ncol = 1)
    } else {
      #otherwise we have covariates, so (1,W) will be a N x (K + 1) vector 
      OneW = cbind(matrix(1, nrow = nrow(Y), ncol = 1), W)
    }
    OneW = as.matrix(OneW)
    
    #calculate the inverse term
    tosum = list()
    for(i in 1:ss){
      tosum[[i]] = t(t(OneW[i,])) %*% t(OneW[i,]) * p1[i] * p0[i]
    }
    SummedMats = Reduce('+', tosum)
    invTerm = solve(SummedMats)
    
    #calculate the summed term
    sumterm.raw = list()
    for(i in 1:ss){
      sumterm.raw[[i]] = t(t(XZ[i,])) %*% t(OneW[i,]) * p1[i] * p0[i]
    }
    summedTerm = Reduce('+', sumterm.raw)
    
    #combine with other terms and sum to get score(s)
    finalTermsToSum = list()
    for(i in  1:ss){
      finalTermsToSum[[i]] = (t(t(XZ[i,])) - summedTerm %*% invTerm %*% t(t(OneW[i,])))*(Y[i]-p1[i])
    }
    UScores = Reduce('+', finalTermsToSum)
  } 
  else {
    #otherwise we have continuous phenotypes so the equation changes slightly
    #calculate U stat using Taylor expanded equation
    #we are summing over N (rows)
    nulllinearreg = glm(Y~W,family=gaussian)
    # expected value of Y from the linear regression
    Yhat = fitted(nulllinearreg)
    
    #define (1,W) as N x (k+1) vector
    if(includeCov == FALSE){
      #if no covariates, then (1,W) will be a N x 1 vector of 1s
      OneW = matrix(1, nrow = nrow(Y), ncol = 1)
    } else {
      #otherwise we have covariates, so (1,W) will be a N x (K + 1) vector 
      OneW = cbind(matrix(1, nrow = nrow(Y), ncol = 1), W)
    }
    OneW = as.matrix(OneW)
    
    #calculate the inverse term
    tosum = list()
    for(i in 1:ss){
      tosum[[i]] = t(t(OneW[i,])) %*% t(OneW[i,])
    }
    SummedMats = Reduce('+', tosum)
    invTerm = solve(SummedMats)
    
    #calculate the summed term
    sumterm.raw = list()
    for(i in 1:ss){
      sumterm.raw[[i]] = t(t(XZ[i,])) %*% t(OneW[i,])
    }
    summedTerm = Reduce('+', sumterm.raw)
    
    #combine with other terms and sum to get score(s)
    finalTermsToSum = list()
    for(i in  1:ss){
      finalTermsToSum[[i]] = (t(t(XZ[i,])) - summedTerm %*% invTerm %*% t(t(OneW[i,])))*(Y[i]-Yhat[i])
    }
    UScores = Reduce('+', finalTermsToSum)
  } 
  return (UScores)
}

CalcUScoreGhat = function(SampleSize, includeCov = FALSE, CovData, CalcUR = TRUE, RGenoData, ScoreData, Phenos, BinPhenos = TRUE){
  
  #define variables
  #SampleSize is number of D/R pairs (N)
  #includeCov: T or F whether or not covariates are included in the modeling
  #CovData: Covariate matrix (N x K) if covariates are included
  #CalcUR: T if we are calculating U for R geno, F if we are calculating U for Score
  #RGenoData: Matrix of R genotype data (N x m)
  #ScoreData: Vector of Gene-based score data (Nx1)
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
  }
  else {
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

#also need to calc Q values
CalcQValues = function(SampleSize, includeCov = FALSE, CovData, CalcUR = TRUE, RGenoData, ScoreData, Phenos, BinPhenos = TRUE){
  
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
  }
  if(CalcUR == TRUE){
    #then we have R geno data
    XZ = RGenoData
    m = ncol(XZ)
  } else {
    #otherwise we have score data
    XZ = ScoreData
    m=1
  }
  Y = Phenos
  
  #######################################################
  #For binary phenotypes
  if(BinPhenos == TRUE){
    #need to calculate p_hat (predicted prob of Y=1 under H0)
    nulllogitreg = glm(Y~W,family=binomial)
    # p: expected value of Y from the logistic regression
    p1 = fitted(nulllogitreg)
    p0 = 1 - p1
    
    #calculate U stat using Taylor expanded equation
    #we are summing over N (rows)
    #define (1,W) as N x (k+1) vector
    if(includeCov == FALSE){
      #if no covariates, then (1,W) will be a N x 1 vector of 1s
      OneW = matrix(1, nrow = nrow(Y), ncol = 1)
    } else {
      #otherwise we have covariates, so (1,W) will be a N x (K + 1) vector 
      OneW = cbind(matrix(1, nrow = nrow(Y), ncol = 1), W)
    }
    OneW = as.matrix(OneW)
    
    #calculate the inverse term
    tosum = list()
    for(i in 1:ss){
      tosum[[i]] = t(t(OneW[i,])) %*% t(OneW[i,]) * p1[i] * p0[i]
    }
    SummedMats = Reduce('+', tosum)
    invTerm = solve(SummedMats)
    
    #calculate the summed term
    sumterm.raw = list()
    for(i in 1:ss){
      sumterm.raw[[i]] = t(t(XZ[i,])) %*% t(OneW[i,]) * p1[i] * p0[i]
    }
    summedTerm = Reduce('+', sumterm.raw)
    
    #calculate Q values
    finalTermsToSum = list()
    for(i in  1:ss){
      finalTermsToSum[[i]] = (t(t(XZ[i,])) - summedTerm %*% invTerm %*% t(t(OneW[i,])))*(Y[i]-p1[i])
    }
    
    Qmat = matrix(unlist(finalTermsToSum), nrow = ss, ncol = m, byrow = TRUE)
    
  } 
  else {
    #otherwise we have continuous phenotypes so the equation changes slightly
    #calculate U stat using Taylor expanded equation
    #we are summing over N (rows)
    nulllinearreg = glm(Y~W,family=gaussian)
    # expected value of Y from the linear regression
    Yhat = fitted(nulllinearreg)
    
    #define (1,W) as N x (k+1) vector
    if(includeCov == FALSE){
      #if no covariates, then (1,W) will be a N x 1 vector of 1s
      OneW = matrix(1, nrow = nrow(Y), ncol = 1)
    } else {
      #otherwise we have covariates, so (1,W) will be a N x (K + 1) vector 
      OneW = cbind(matrix(1, nrow = nrow(Y), ncol = 1), W)
    }
    OneW = as.matrix(OneW)
    
    #calculate the inverse term
    tosum = list()
    for(i in 1:ss){
      tosum[[i]] = t(t(OneW[i,])) %*% t(OneW[i,])
    }
    SummedMats = Reduce('+', tosum)
    invTerm = solve(SummedMats)
    
    #calculate the summed term
    sumterm.raw = list()
    for(i in 1:ss){
      sumterm.raw[[i]] = t(t(XZ[i,])) %*% t(OneW[i,])
    }
    summedTerm = Reduce('+', sumterm.raw)
    
    #combine with other terms and sum to get score(s)
    finalTermsToSum = list()
    for(i in  1:ss){
      finalTermsToSum[[i]] = (t(t(XZ[i,])) - summedTerm %*% invTerm %*% t(t(OneW[i,])))*(Y[i]-Yhat[i])
    }
    Qmat = matrix(unlist(finalTermsToSum), nrow = ss, ncol = m, byrow = TRUE)
  } 
  return(Qmat)
}

CalcQValuesGhat = function(SampleSize, includeCov = FALSE, CovData, CalcUR = TRUE, RGenoData, ScoreData, Phenos, BinPhenos = TRUE){
  
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
###################################       
##Step 6) Calculate original variance
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

###################################       
##Step 7) Calc final statistic
CalcStatisticPVal = function(SampleSize, Variance, UscoresR, UscoreS, s){
  #define variables
  #SampleSize is number of D/R pairs (N)
  #Variance is the original Var/Cov matrix for the scores calculated  using CalcVariance()
  #### dim should be (m+1) x (m+1)
  #UscoresR is the m x 1 vector of scores for the r genos
  #UscoreS should be the 1x1 score value for the gene-based score
  #s is the percent of variance we want explained by the PCs we are keeping
  
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
  eigen_percent = cumsum(lambdas/sum(lambdas))
  #determine if s is too small
  #if s is smaller than smallest PVE, use 1 PC
  s_pct = s/100
  if(eigen_percent[1] >= s_pct){
	num_pc = 1
  } else {
    num_pc = which(eigen_percent >= s_pct)[1]
  }
  
  # num_pc = s #number of PCs is directly equal to s
  
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
  dof = num_pc + 1
  
  finalOutput = cbind(finalStat, pval, dof)
  return(finalOutput)
}

################################################
## Putting some functions together
CalcNullPhenotypeData = function(chr, numPairs, simNum, YPrev, gene, path, weightedScores = FALSE, standardizeScores = FALSE){
  #pull recipient and donor genotypes
  RGenos = obtainRGenotypes(chr = chr, numSamples = numPairs, simNum = simNum, gene = gene, path = path)
  DGenos = obtainDGenotypes(chr = chr, numSamples = numPairs, simNum = simNum, gene = gene, path = path)
  
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
  
  Phenos = list(CatPhenos, ContPhenos)
  
  return(Phenos)
}

CalcAltPhenotypeData_Scores = function(chr, numPairs, simNum, YPrev, gene, path, percentageAssoc, LowLD, TrueScore, Gamma){
  #pull recipient and donor genotypes
  RGenos = obtainRGenotypes(chr = chr, numSamples = numPairs, simNum = simNum, gene = gene, path = path)
  DGenos = obtainDGenotypes(chr = chr, numSamples = numPairs, simNum = simNum, gene = gene, path = path)
  
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
  
  AltScorePhenos = list(RGenos, IBS.gene, Incomp.gene, AMS.gene, BinMM.gene, CovData, CatPhenos, ContPhenos)
  return(AltScorePhenos)
}

CalcAltPhenotypeData_RSNPs = function(chr, numPairs, simNum, YPrev, gene, path, ORSize, LowLD, percentAssoc, TrueScore){
    #define effect based on OR size
  if(ORSize == "small"){
    effect = 0.14
  } else if(ORSize == "medium"){
    effect = 0.41
  } else {
    effect = 0.69
  }
  
  assocSNPs = DetermineAssocRSNPs(gene = gene, LowLD = LowLD, percentageAssoc = percentageAssoc)
  
  #pull recipient and donor genotypes
  RGenos = obtainRGenotypes(chr = chr, numSamples = numPairs, simNum = simNum, gene = gene, path = path)
  DGenos = obtainDGenotypes(chr = chr, numSamples = numPairs, simNum = simNum, gene = gene, path = path)
  
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
  
  nullGamma = c(0)
  
  #generate phenotypes, both continuous and binary
  #Based on single true score
  CatPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = TRUE, YPrev = YPrev,  Covariates = CovData, RGenoData = RGenos, ScoreData = PhenoScore, Betas = Betas, Gamma = nullGamma)
  ContPhenos = GenAltPhenos(SampleSize = numPairs, includeCov = TRUE, YCat = FALSE,  Covariates = CovData, RGenoData = RGenos, ScoreData = PhenoScore, Betas = Betas, Gamma = nullGamma)
  
  AltPhenos_RSNPs = list(RGenos, IBS.gene, Incomp.gene, AMS.gene, BinMM.gene, CovData, CatPhenos, ContPhenos)
}

RunSKATAnalysis = function(PhenoList, kernel, kernelWeights = c()){
  RGenos = PhenoList[[1]]
  IBS.gene = PhenoList[[2]]
  Incomp.gene = PhenoList[[3]]
  AMS.gene = PhenoList[[4]]
  BinMM.gene = PhenoList[[5]]
  CovData = PhenoList[[6]]
  CatPhenos = PhenoList[[7]]
  ContPhenos = PhenoList[[8]]
  
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
  
  return(statsAndPValsMat) 
}

RunSKATBinaryAnalysis = function(PhenoList, kernel, kernelWeights = c()){
  RGenos = PhenoList[[1]]
  IBS.gene = PhenoList[[2]]
  Incomp.gene = PhenoList[[3]]
  AMS.gene = PhenoList[[4]]
  BinMM.gene = PhenoList[[5]]
  CovData = PhenoList[[6]]
  CatPhenos = PhenoList[[7]]
  ContPhenos = PhenoList[[8]]
  
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
    Stat_IBS_CatPhenos = SKATBinary(RGeno.IBS.snp, obj_dich, kernel = kernel, is_check_genotype = FALSE)
    Stat_Incomp_CatPhenos = SKATBinary(RGeno.Incomp.snp, obj_dich, kernel = kernel, is_check_genotype = FALSE)
    Stat_AMS_CatPhenos = SKATBinary(RGeno.AMS.snp, obj_dich, kernel = kernel, is_check_genotype = FALSE)
    Stat_BinMM_CatPhenos = SKATBinary(RGeno.BinMM.snp, obj_dich, kernel = kernel, is_check_genotype = FALSE)
    
    Stat_IBS_ContPhenos = SKAT(RGeno.IBS.snp, obj_cont, kernel = kernel, is_check_genotype = FALSE)
    Stat_Incomp_ContPhenos = SKAT(RGeno.Incomp.snp, obj_cont, kernel = kernel, is_check_genotype = FALSE)
    Stat_AMS_ContPhenos = SKAT(RGeno.AMS.snp, obj_cont, kernel = kernel, is_check_genotype = FALSE)
    Stat_BinMM_ContPhenos = SKAT(RGeno.BinMM.snp, obj_cont, kernel = kernel, is_check_genotype = FALSE)
  } else {
    Stat_IBS_CatPhenos = SKATBinary(RGeno.IBS.snp, obj_dich, kernel = kernel, weights = kernelWeights, is_check_genotype = FALSE)
    Stat_Incomp_CatPhenos = SKATBinary(RGeno.Incomp.snp, obj_dich, kernel = kernel, weights = kernelWeights, is_check_genotype = FALSE)
    Stat_AMS_CatPhenos = SKATBinary(RGeno.AMS.snp, obj_dich, kernel = kernel, weights = kernelWeights, is_check_genotype = FALSE)
    Stat_BinMM_CatPhenos = SKATBinary(RGeno.BinMM.snp, obj_dich, kernel = kernel, weights = kernelWeights, is_check_genotype = FALSE)
    
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
  
  return(statsAndPValsMat) 
}
