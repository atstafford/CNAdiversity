Packages <- c("readxl","data.table",
              "tidyverse",
              "ggrepel","ggplot2","ggpubr","cowplot",
              "pvclust",
              "BSgenome","rtracklayer","chromoMap",
              "fastDummies","betareg","frmselection","car",
              "stats",
              "foreach","doParallel","ranger","palmerpenguins","kableExtra","Metrics",
              "corrplot","survival","survminer",
              "rrvgo",
              "limma")
lapply(Packages, library, character.only = TRUE)

# Datasets
# Training
ad.raw <- readRDS("~/Documents/CNA/Data/ad.raw.rds")
ad.info <- readRDS("~/Documents/CNA/Data/ad.info.rds")
ad.clonality <- readRDS("~/Documents/CNA/Data/ad.clonality.rds")
ad.diversity <- readRDS("~/Documents/CNA/Data/ad.diversity.rds")
car.raw <- readRDS("~/Documents/CNA/Data/car.raw.rds")
car.info <- readRDS("~/Documents/CNA/Data/car.info.rds")
car.clonality <- readRDS("~/Documents/CNA/Data/car.clonality.rds")
car.diversity <- readRDS("~/Documents/CNA/Data/car.diversity.rds")
actualIth_train <- readRDS("~/Documents/CNA/Data/actualIth_train.rds")
actualIth_ad <- readRDS("~/Documents/CNA/Data/actualIth_ad.rds")
ad.matrices <- readRDS("~/Documents/CNA/Data/ad.matrices.rds")
car.matrices <- readRDS("~/Documents/CNA/Data/car.matrices.rds")
predictedIth_train <- readRDS("~/Documents/CNA/Data/predictedIth_train.rds")

# Validation
val.info <- readRDS("~/Documents/CNA/Data/val.info.rds")
val.clonality <- readRDS("~/Documents/CNA/Data/val.clonality.rds")
val.diversity <- readRDS("~/Documents/CNA/Data/val.diversity.rds")
testcnBinned <- readRDS("~/Documents/CNA/Data/testcnBinned.rds")
testcnBinned.list <- readRDS("~/Documents/CNA/Data/testcnBinned.list.rds")
actualIth_test <- readRDS("~/Documents/CNA/Data/actualIth_test.rds")
predictedIth_test <- readRDS("~/Documents/CNA/Data/predictedIth_test.rds")

# epicc?
epicc.info <- readRDS("~/Documents/CNA/Data/epicc.info.rds")
epicc.clonality <- readRDS("~/Documents/CNA/Data/epicc.clonality.rds")
epicc.diversity <- readRDS("~/Documents/CNA/Data/epicc.diversity.rds")
epiccbinned <- readRDS("~/Documents/CNA/Data/epiccbinned.rds")
epiccbinned.list <- readRDS("~/Documents/CNA/Data/epiccbinned.list.rds")

TCGA.info.list <- readRDS("~/Documents/CNA/Data/TCGA.info.list.rds")
TCGA.byCan.list <- readRDS("~/Documents/CNA/Data/TCGA.byCan.list.rds")
TCGA_clinical <- readRDS("~/Documents/CNA/Data/TCGA_clinical.rds")
TCGA.cn <- readRDS("~/Documents/CNA/Data/TCGA.cn.rds")
TCGA_raw <- readRDS("~/Documents/CNA/Data/TCGA_raw.rds")


# COAD
TCGA.info <- readRDS("~/Documents/CNA/Data/TCGA.info.rds")
COADbinned <- readRDS("~/Documents/CNA/Data/COADbinned.rds")
COAD_clinical <- readRDS("~/Documents/CNA/Data/COAD_clinical.rds")

# LAUD
LAUD.info <- readRDS("~/Documents/CNA/Data/LAUD.info.rds")
LAUDbinned <- readRDS("~/Documents/CNA/Data/LAUDbinned.rds")
LAUD_clin <- readRDS("~/Documents/CNA/Data/LAUD_clin.rds")

# TracerX
tracerx.info <- readRDS("~/Documents/CNA/Data/tracerx.info.rds")
tracerx.clonality <- readRDS("~/Documents/CNA/Data/tracerx.clonality.rds")
tracerx.diversity <- readRDS("~/Documents/CNA/Data/tracerx.diversity.rds")
lung.matrices <- readRDS("~/Documents/CNA/Data/lung.matrices.rds")
tracerxBinned <- readRDS("~/Documents/CNA/Data/tracerxBinned.rds")
actualIth_tracerx <- readRDS("~/Documents/CNA/Data/actualIth_tracerx.rds")
uniReg.out.list.lung <- readRDS("~/Documents/CNA/Data/uniReg.out.list.lung.rds")
sig.hclust.lung <- readRDS("~/Documents/CNA/Data/sig.hclust.lung.rds")
candidate.bins.lung <- readRDS("~/Documents/CNA/Data/candidate.bins.lung.rds")
representitive.bins.lung <- readRDS("~/Documents/CNA/Data/representitive.bins.lung.rds")
selection_for.lung <- readRDS("~/Documents/CNA/Data/selection_for.lung.rds")
selection_back.lung <- readRDS("~/Documents/CNA/Data/selection_back.lung.rds")
lung.inbackrem <- readRDS("~/Documents/CNA/Data/lung.inbackrem.rds")
beta_lung.backrem <- readRDS("~/Documents/CNA/Data/beta_lung.backrem.rds")
hg19predictors.lung <- readRDS("~/Documents/CNA/Data/hg19predictors.lung.rds")


# Clusters
load("~/Documents/CNA/HPC/hiclust_diploid.aneu.rda")
hiclust <- hiclust_diploid.aneu
sig.hclust <- readRDS("~/Documents/CNA/Data/sig.hclust.rds")

# Load singlebin models
uniReg.out.list <- readRDS("~/Documents/CNA/Data/uniReg.out.list.rds")

# betareg
candidate.bins <- readRDS("~/Documents/CNA/Data/candidate.bins.rds")
representitive.bins<- readRDS("~/Documents/CNA/Data/representitive.bins.rds")
beta_52backrem <- readRDS("~/Documents/CNA/Data/beta_52backrem.rds")
hg19predictors <- readRDS("~/Documents/CNA/Data/hg19predictors.rds")
predictedIth_train <- readRDS("~/Documents/CNA/Data/predictedIth_train.rds")
predictedIth_test <- readRDS("~/Documents/CNA/Data/predictedIth_test.rds")

# Load gene annotation of bins that are individually predictive
gene.anno <- readRDS("~/Documents/CNA/Data/gene.anno.rds")
GOanno <- readRDS("~/Documents/CNA/Data/GOanno.rds")
reducedTerms_all <- readRDS("~/Documents/CNA/Data/reducedTerms_all.rds")





new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages) > 0){
  install.packages(new.packages, dep=TRUE)
}
for(package.i in list.of.packages){
  suppressPackageStartupMessages(
    library(
      package.i, 
      character.only = TRUE
    )
  )
}