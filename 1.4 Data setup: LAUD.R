# TCGA PROCESSING ####

#load
load("~/Dropbox/CNA diversity/Data/TCGA/chrom_start_stop_copynumbermean.rData")
LAUD_clin <- as.data.frame(fread("~/Dropbox/CNA diversity/Data/TCGA/clinical.project-TCGA-LUAD.2022-02-23/clinical.tsv"))
colnames(LAUD_clin)[c(1,2)] <- c("UUID_copyN","TCGA_barcode")

# Remove duplicates
LAUD_clin <- LAUD_clin[!duplicated(LAUD_clin$TCGA_barcode), ]

# Remove columns with no information
LAUD_clin <- LAUD_clin[,colSums(LAUD_clin=="'--") < nrow(LAUD_clin)]
LAUD_clin <- LAUD_clin[,-c(3,14,15,19,18,20,23,28,32)]

# Change stage names to remove A/B
colnames(LAUD_clin)[13] <- 'stage'
for ( i in 1:nrow(LAUD_clin) ) {
  LAUD_clin$stage[i] <- ifelse(LAUD_clin$stage[i]=="Stage IA", "Stage I", 
                               ifelse(LAUD_clin$stage[i]=="Stage IB", "Stage I",
                               ifelse(LAUD_clin$stage[i]=="Stage IIA", "Stage II",
                                      ifelse(LAUD_clin$stage[i]=="Stage IIB", "Stage II",
                                             ifelse(LAUD_clin$stage[i]=="Stage IIIA", "Stage III",
                                                    ifelse(LAUD_clin$stage[i]=="Stage IIIB", "Stage III",
                                                           ifelse(LAUD_clin$stage[i]=="Stage IIIC", "Stage III",
                                                                  ifelse(LAUD_clin$stage[i]=="Stage IVA", "Stage IV",
                                                                         ifelse(LAUD_clin$stage[i]=="Stage IVB", "Stage IV",LAUD_clin$stage[i])))))))))
}

# Keep patient which has cn and clinical data
LAUD_cn <- chrom_start_stop

LAUD_cn <- cbind(LAUD_cn[,c(1:3)], LAUD_cn[ ,which(str_sub(colnames(LAUD_cn),1,str_length(colnames(LAUD_cn))-3) %in% LAUD_clin$TCGA_barcode)])
LAUD_cn <- LAUD_cn[,colSums(is.na(LAUD_cn))<nrow(LAUD_cn)]

LAUD_clin <- LAUD_clin[which(LAUD_clin$TCGA_barcode %in% str_sub(colnames(LAUD_cn),1,str_length(colnames(LAUD_cn))-3)),] #2 patients have 2 samples

# Pull cn data into element corresponding to sample ID
names <- (colnames(LAUD_cn)[-c(1:3)])
LAUD.list <- list()
for ( i in 1:length(names) ) {
  LAUD.list[[i]] <- cbind(LAUD_cn[,c(1:3)], LAUD_cn[,which(colnames(LAUD_cn) == names[i])])
  colnames(LAUD.list[[i]])[c(1:3)] <- c("chr","start","stop")
}

# Assess ploidy and recentre if average CN is 2.8 or more
LAUD.list <- ploidyRecentre(LAUD.list,skipcol = 3, multi_cn = FALSE)

# Bin CN data to match training dataset
# LAUD doesnt have data on 898, 1074, 1090, 1374, 1449, 1552
LAUDbinned.list <- newAlignBins(bins = car.info$start.stop, cn.list = LAUD.list)

# If absolute copy number >=3, its a gain, if <=2 its a loss
LAUDbinned.list <- lapply(LAUDbinned.list, function(x) {
  x[,-c(1:4)] <- ifelse(x[,-c(1:4)] >= 3, 3, ifelse(x[,-c(1:4)] < 2, 1, 2))
  x
})

# Rename columns with sampleID and combine all patients into one dataframe
for ( i in 1:length(LAUDbinned.list) ) {
  for ( k in 5:ncol(LAUDbinned.list[[i]]) ) {
    colnames(LAUDbinned.list[[i]])[k] <- names[i]
  }
}
LAUDbinned <- LAUDbinned.list %>% purrr::reduce(full_join, by = c("chr","start","stop","bin"))
LAUDbinned <- LAUDbinned[,-1] # remove bin columns

# Use PullData functions to create necessary dataframes for analysis
LAUD.info <- PullDataInfo(rawdata = LAUDbinned)

saveRDS(LAUD.info, "~/Documents/CNA/Data/LAUD.info.rds")
saveRDS(LAUDbinned, "~/Documents/CNA/Data/LAUDbinned.rds")
saveRDS(LAUD_clin, "~/Documents/CNA/Data/LAUD_clin.rds")



