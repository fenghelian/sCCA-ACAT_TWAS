####################################################
#Example code to create sCCA cross-tissue weight
#Author: Helian Feng
#Date: 2019-06-07
###################################################
require(dplyr)
require(tidyr)
require(PMA)
require(mice)
suppressMessages(require('plink2R'))

# read in genotype file
geno.file="../data/sample/2.keep"
ref.geno.file="../data/sample/2.ref"
genos = read_plink(geno.file,impute="avg")
ref.genos = read_plink(ref.geno.file,impute="avg")
m = match(  genos$bim[,2] , ref.genos$bim[,2] )
m.keep = !is.na(m)
genos$bed = scale(genos$bed[,m[m.keep]])
genos$bim = genos$bim[m[m.keep],]
genos$bed = scale(genos$bed)
geno<-as.data.frame(genos$bed)
ref.genos$bed = scale(ref.genos$bed[,m[m.keep]])
ref.genos$bim = ref.genos$bim[m[m.keep],]
ref.genos$bed = scale(ref.genos$bed)
ref.geno<-as.data.frame(ref.genos$bed)
geno$ID<-matrix(unlist(strsplit(rownames(geno),split=":")), ncol=2, byrow=TRUE)[,1]

# read in expression matrix(each row for an individual, column for tissue)
p_tab <- read.csv("/data/p_tab.csv", stringsAsFactors=FALSE)
exp<-p_tab[,3:24]
exp_imp <-mice(exp,m=1)
exp_imp<-as.matrix(complete(exp_imp))
geno_i<-as.matrix(geno%>%select(-ID))

#sCCA step
perm.out <- PMA::CCA.permute(exp_imp,geno_i,typex="standard",typez="ordered",nperms=15,trace=F)
out <- PMA::CCA(exp_imp,geno_i,typex="standard",typez="ordered",K=3,penaltyx=perm.out$bestpenaltyx,penaltyz=perm.out$bestpenaltyz, v=perm.out$v.init,trace=F)
for( j in 1:3){
  expr_l<-exp_imp%*%out$u[,j]
  pheno<-cbind(p_tab[,1:2],0,0,0,expr_l)
  write.table(pheno,paste0("/n/pmage1/hlfeng/helian/BSLMM/",OUT,NUM,".keep.fam"),row.names = F,col.names = F,quote = F)

  cmd<-paste0("Rscript /Tools/fusion/FUSION.compute_weights.R --bfile /n/",
              OUT,NUM,".keep --tmp /n/",OUT,NUM,
              ".tmp --out /gtex_weight/",
              NUM,".",as.character(i),"_scca",as.character(j),
              " --verbose 1 --hsq_p 0.01 --save_hsq --PATH_plink /Tools/Plink2/plink --PATH_gcta Tools/gcta_nr_robust --PATH_gemma /Tools/gemma.linux --models lasso,top1,enet")
  system(cmd)
  cmd=paste0("cat /CROSS_TIS/gtex_weight/",NUM,".",as.character(i),"_scca",as.character(j),".hsq >> /CROSS_TIS/hsq/",NUM,".hsq")
  system(cmd)
  cmd=paste0("rm -f /CROSS_TIS/gtex_weight/",NUM,".",as.character(i),"_scca",as.character(j),".hsq")
  system(cmd)
  }
  cmd=paste0("rm /n/",OUT,NUM,"*")
  system(cmd)
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
