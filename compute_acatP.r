suppressMessages(library('dplyr'))
suppressMessages(library("tidyr"))
suppressMessages(library("readr"))
suppressMessages(library("optparse"))

option_list = list(
  make_option("--file_location", action="store", default=NA, type='character',
              help="Path to TWAS result files from FUSION [required]"),
  make_option("--out", action="store", default=NA, type='character',
              help="Path to output files [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

# ACAT function
CCT.pval<-function(Pvals,Weights=NULL){
  #### check if there is NA
  if (sum(is.na(Pvals))>0){
    if(all(is.na(Pvals))){
      return(NA)
    }else{
      Pvals=Pvals[!is.na(Pvals)]
      if(length(Pvals)==1){
        return(Pvals)
      }
    }
    #stop("Cannot have NAs in the p-values!")
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
    Pvals[Pvals==0]<-as.numeric(1.507e-143)
    warning("There are p-values that are exactly 1!")
    #return(0)
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

#Read in data
file_list <- list.files(path=opt$file_location)
#file_list <- list.files(path="/Users/Helian/OneDrive/biosta/TWAS/cancer_twas/data/V8/sCCA/brca/")
dat<-NULL
for (i in 1:length(file_list)){
  filename<-paste0(opt$file_location,file_list[i])
  #filename<-paste0("/Users/Helian/OneDrive/biosta/TWAS/cancer_twas/data/V8/sCCA/brca/",file_list[i])
  tmp <- read_delim(filename,"\t", escape_double = FALSE, trim_ws = TRUE)
  dat<-rbind(dat,tmp)
}
# Get ACAT Pvalue
file_tab <- dat%>%select(c(FILE))
file_tab <-file_tab%>%separate(col=FILE,into=paste0("x",0:12),sep="/")
panel_list <-file_tab$x11
dat$PANEL <- ifelse(is.na(dat$PANEL),panel_list,dat$PANEL)
g_list<-dat%>%group_by(ID,PANEL)%>%summarize(n=n())
g_list<-g_list$ID[g_list$n>1]
dat<-dat%>%filter(!ID %in% g_list)
dat<-dat%>%select(c(PANEL,ID,TWAS.P))
dat[,3]<-as.data.frame(sapply(dat[,3], as.numeric))
dat <- dat%>%filter(!is.na(TWAS.P))
dat_spread <- dat%>% spread(key = PANEL,value = TWAS.P, drop = T)
n_col <-dim(dat_spread)[2]
acat<-apply(dat_spread[,2:n_col],FUN=CCT.pval,1)
dat_spread$acat<-acat

write.table(dat_spread,opt$out,sep = "\t",quote = F,row.names = F,col.names = T)
