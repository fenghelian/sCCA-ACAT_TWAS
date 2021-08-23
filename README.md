# sCCA-ACAT_TWAS
Repository for conducting sCCA+ACAT TWAS with Fusion (1,2) and ACAT (3). The precomputed sCCA TWAS weights for GTEx gene expression version 6 and 8 (5) can be found at TWAS HUB (http://gusevlab.org/projects/fusion/#reference-functional-data) as well as FUSION weight page on for Alkes Group Page (https://alkesgroup.broadinstitute.org/FUSION/WGT/) (4). We also provided code to generate sCCA TWAS weights with raw gene expression and genotype data (Note that this is not required for conducting sCCA-TWAS). 

# Installation 
* Download and unpack the  FUSION software package from github:

wget https://github.com/gusevlab/fusion_twas/archive/master.zip

unzip master.zip

* Download and unpack the sCCA+ACAT package from github:

wget https://github.com/fenghelian/sCCA-ACAT_TWAS/archive/master.zip

unzip master.zip

cd sCCA-ACAT_TWAS-master

* Download and unpack the (1000 Genomes)  LD reference data:

wget https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2

tar xjvf LDREF.tar.bz2

* Launch R and install required libraries:

install.packages(c('optparse','RColorBrewer'))

install.packages('plink2R-master/plink2R/',repos=NULL)

install.packages(c('dplyr','tidyr'))

# Steps to conduct sCCA+ACAT TWAS
To conduct sCCA+ACAT TWAS, we first need to conduct single tissue and sCCA TWAS with single-tissue and sCCA cross-tissue weights with FUSION first. The typical TWAS analysis takes pre-computed gene expression weights (below) together with disease GWAS summary statistics to estimate the association of each gene to disease. For example, we could use the PGC Schizophrenia summary statistics to perform a TWAS with interested single GTEx tissue weights data as well as sCCA cross-tissue weights (for simplicity, we are showing the example with 2 tissues, but the workflow can be easily generalized to more than 2 tissues). The example assumes you have setup FUSION and LD reference data as above and are in the FUSION directory with an LDREF subdirectory.

1. First, download and prepare the GWAS and GTEx whole blood data:

wget https://data.broadinstitute.org/alkesgroup/FUSION/SUM/PGC2.SCZ.sumstats

mkdir WEIGHTS

cd WEIGHTS

wget https://data.broadinstitute.org/alkesgroup/FUSION/WGT/GTEx.Whole_Blood.tar.bz2

tar xjf GTEx.Whole_Blood.tar.bz2

wget https://data.broadinstitute.org/alkesgroup/FUSION/WGT/GTEx.Brain_Caudate_basal_ganglia.tar.bz2

tar xjf GTEx.Brain_Caudate_basal_ganglia.tar.bz2

wget http://gusevlab.org/projects/fusion/weights/sCCA_weights_v8_2.zip

unzip sCCA_weights_v8_2.zip

The WEIGHTS directory should contain a subdirectory of expression weights (which you can inspect in R), as well as several report files that describe the data (see below for details). The following sections describe the inputs in detail.

2. Perform TWAS on each tissue

More detailed instruction are provided at FUSION webpage (http://gusevlab.org/projects/fusion/)

* Input: GWAS summary statistics

The primary input is genome-wide summary statistics in LD-score format. At minimum, this is a flat file with a header row containing the following fields:

SNP – SNP identifier (rsID)
A1 – first allele (effect allele)
A2 – second allele (other allele)
Z – Z-scores, sign with respect to A1.
and subsequent data rows for each SNP (all white-space separated). Additional columns are allowed and will be ignored. We recommend using the LDSC munge_stats.py utility for converting GWAS summary data into this format, which detects and reports many common pitfalls.


* Input: Expression weights
The functional reference data are loaded using the ./WEIGHTS/#TISSUE_NAME#.pos which points to the individual *.RDat weight files, their gene identifiers, physical positions (and an optional PANEL column for the panel/study name). We have pre-computed *.pos files for all reference data downloadable below. Only weights in the file will be evaluated. The physical positions should correspond to the feature (e.g. TSS and TES) and will be used in the final output and for plotting.

Finally, we run FUSION.test.R using this data on chromosome 22:
for tissue in GTEx.Whole_Blood GTEx.Brain_Caudate_basal_ganglia sCCA_weights_v8_2
do 
Rscript fusion_twas-master/FUSION.assoc_test.R \
--sumstats PGC2.SCZ.sumstats \
--weights ./WEIGHTS/$tissue.pos \
--weights_dir ./WEIGHTS/ \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 22 \
--out ./out/PGC2.$tissue.22.dat
done

3. Combine TWAS test results cross tissue with ACAT

Rscript ./compute_acatP.r \
--file_location ./out/ \
--out ACAT_pvalue.tsv

The output file would contain the Pvalue for each gene in each gene expression panel and the last column with column name "acat" contains the ACAT Pvalue for all the tissues.

# To compute your own sCCA cross tissue weights 

Note that this step is not required for sCCA+ACAT TWAS, as the pre-computed weights are available at TWAS HUB.

The script for computing expression weights,it works one gene at a time, taking as input a standard binary PLINK format file (bed/bim/fam) for genotype and reference geno file, which contains only the desired SNPs in the cis-locus of the gene (or any other desired SNPs) and expression matrix for the gene (individual  in rows and tissue in columns, with first column as individual ID). Note that, to compute the TWAS weights you also need to have the FUSION package (1,2) and the required software and package by FUSION installed installed. 

A typical run looks like this:

Rscript FUSION.compute_weights.R \
--geno_file $INP \
--ref_geno_file $REF \ 
--expression_file $EXP \
--fusion_location $FUSION \
--tmp $TMP \
--out $OUT \

# To cite this work

Feng, H., Mancuso, N., Gusev, A., Majumdar, A., Major, M., Pasaniuc, B. and Kraft, P., 2021. Leveraging expression from multiple tissues using sparse canonical correlation analysis and aggregate tests improves the power of transcriptome-wide association studies. PLoS genetics, 17(4), p.e1008973.

# Reference
1. FUSION Tool (http://gusevlab.org/projects/fusion/)
2. Gusev, A., Ko, A., Shi, H., Bhatia, G., Chung, W., Penninx, B.W., Jansen, R., De Geus, E.J., Boomsma, D.I., Wright, F.A. and Sullivan, P.F., 2016. Integrative approaches for large-scale transcriptome-wide association studies. Nature genetics, 48(3), pp.245-252.
3. Liu, Y., Chen, S., Li, Z., Morrison, A.C., Boerwinkle, E. and Lin, X., 2019. Acat: A fast and powerful p value combination method for rare-variant analysis in sequencing studies. The American Journal of Human Genetics, 104(3), pp.410-421.
4. TWAS HUB (www.http://twas-hub.org/)
5. GTEx Consortium, 2017. Genetic effects on gene expression across human tissues. Nature, 550(7675), p.204.
