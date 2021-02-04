####################################################
#Script to create sCCA cross-tissue weight
#Author: Helian Feng
#Date: 2019-06-07
###################################################
require(dplyr)
require(tidyr)
require(PMA)
require(mice)
suppressMessages(require('plink2R'))

option_list = list(
  make_option("--geno_file", action="store", default=NA, type='character',
              help="Path to genotype file with PLINK format"),
  make_option("--ref_geno_file", action="store", default=NA, type='character',
              help="Path to reference genotype file with PLINK format"),
  make_option("--expression_file", action="store", default=NA, type='character',
              help="Path to gene expression in csv"),
  make_option("--fusion_location", action="store", default=NA, type='character',
              help="Path to FUSION package"),
  make_option("--plink_location", action="store", default=NA, type='character',
              help="Path to plink"),
  make_option("--gcta_location", action="store", default=NA, type='character',
              help="Path to gcta"),
  make_option("--gemma_location", action="store", default=NA, type='character',
              help="Path to gemma"),
  make_option("--tmp_location", action="store", default=NA, type='character',
              help="Path to save tmp file"),
  make_option("--out", action="store", default=NA, type='character',
              help="Path to output files [required]"),
)

opt = parse_args(OptionParser(option_list=option_list))


# read in genotype file

genos = read_plink(opt$geno_file,impute="avg")
ref.genos = read_plink(opt$ref_geno_file,impute="avg")
m = match(  genos$bim[,2] , opt$ref_genos$bim[,2] )
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
p_tab <- read.csv(opt$expression_file, stringsAsFactors=FALSE)
exp<-p_tab[,2:]
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

  cmd<-paste0("Rscript",opt$fusion_location, "FUSION.compute_weights.R --bfile ",
              opt$geno_file," --tmp ",opt$tmp_location,
              " --out ",opt$out,
              " --verbose 1 --hsq_p 0.01 --save_hsq --PATH_plink ", opt$plink_location, " --PATH_gcta ",opt$gcta_location," --PATH_gemma ",opt$gemma_location , " --models lasso,top1,enet")
system(cmd)
