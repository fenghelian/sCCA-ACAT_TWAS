####################################################
#Example code to perform sCCA+ACAT TWAS
#Author: Helian Feng
#Date: 2019-06-07
###################################################
library(dplyr)
library(tidyr)

#ACAT function
CCT.pval<-function(Pvals,Weights=NULL){
  #### check if there is NA
  if (sum(is.na(Pvals))>0){
    stop("Cannot have NAs in the p-values!")
  }
  #### check if Pvals are between 0 and 1
  if ((sum(Pvals<0)+sum(Pvals>1))>0){
    warning("P-values must be between 0 and 1!")
  }
  #### check if there are pvals that are either exactly 0 or 1.
  is.zero<-(sum(Pvals==0)>=1)
  is.one<-(sum(Pvals==1)>=1)
  if (is.zero && is.one){
    warning("Cannot have both 0 and 1 p-values!")
  }
  if (is.zero){
    return(0)
  }
  if (is.one){
    Pvals[Pvals==1]<-1-1/length(Pvals)
    warning("There are p-values that are exactly 1!")
  }

  #### Default: equal weights. If not, check the validity of the user supplied weights and standadize them.
  if (is.null(Weights)){
    Weights<-rep(1/length(Pvals),length(Pvals))
  }else if (length(Weights)!=length(Pvals)){
    stop("The length of weights should be the same as that of the p-values")
  }else if (sum(Weights<0)>0){
    stop("All the weights must be positive!")
  }else{
    Weights<-Weights/sum(Weights)
  }


  #### check if there are very small non-zero p values
  is.small<-(Pvals<1e-16)
  if (sum(is.small)==0){
    cct.stat<-sum(Weights*tan((0.5-Pvals)*pi))
  }else{
    cct.stat<-sum((Weights[is.small]/Pvals[is.small])/pi)
    cct.stat<-cct.stat+sum(Weights[!is.small]*tan((0.5-Pvals[!is.small])*pi))
  }
  #### check if the test statistic is very large.
  if (cct.stat>1e+15){
    pval<-(1/cct.stat)/pi
  }else{
    pval<-1-pcauchy(cct.stat)
  }
  return(pval)
}
##################################################

#unzip sCCA weights and single-tissue weights from TWAS HUB

mkdir WEIGHTS
cd WEIGHTS
tar xjf GTEx.sCCA.tar.bz2

wget https://data.broadinstitute.org/alkesgroup/FUSION/WGT/GTEx.Whole_Blood.tar.bz2
tar xjf GTEx.Whole_Blood.tar.bz2


# run TWAS one chromosome at a time for sCCA and single tissue weights

Rscript FUSION.assoc_test.R \
--sumstats PGC2.SCZ.sumstats \
--weights ./WEIGHTS/GTEx.sCCA.pos \
--weights_dir ./WEIGHTS/ \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 22 \
--out PGC2.SCZ.22.dat

Rscript FUSION.assoc_test.R \
--sumstats PGC2.SCZ.sumstats \
--weights ./WEIGHTS/GTEx.Whole_Blood.pos \
--weights_dir ./WEIGHTS/ \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 22 \
--out PGC2.SCZ.22.dat

# Combine test results with ACAT
#all_p is the TWAS test results matrix with each gene in rows and test results 
#for each tissue in column 3 to 50
acat<-apply(all_p[,3:50],FUN=CCT.pval,1)
#acat is the test statistics vector for the genes
