# TCGA PROCESSING ####

# Load clinical data for CRC (COAD)
COAD_clinical <- read_excel("~/Documents/CNA/Data/TCGA/TCGAcn/Clinical COAD.xlsx")
colnames(COAD_clinical)[c(1,2)] <- c("UUID_copyN","TCGA_barcode")

# Remove duplicates
COAD_clinical <- COAD_clinical[!duplicated(COAD_clinical$TCGA_barcode), ]

# Remove columns with no information
COAD_clinical <- COAD_clinical[,colSums(COAD_clinical=="'--") < nrow(COAD_clinical)]
COAD_clinical <- COAD_clinical[,-c(3,14,15,19,18,20,23,28,32,33)]

# Load additional data and merge
COAD_MSI <- read_excel("~/Documents/CNA/Data/TCGA/TCGAcn/TCGA_MSIstatus.xlsx", col_names = TRUE, skip =1)
COAD_MSI <- COAD_MSI[COAD_MSI$Organ=="COAD",]
COAD_clinical$MSI <- unlist(COAD_MSI[match(COAD_clinical$TCGA_barcode, COAD_MSI$`TCGA Participant Barcode`), 14])
COAD_clinical$CMS <- unlist(COAD_MSI[match(COAD_clinical$TCGA_barcode, COAD_MSI$`TCGA Participant Barcode`), 27])
COAD_clinical$MSI[is.na(COAD_clinical$MSI)] <- 'NA'

# Change stage names to remove A/B
colnames(COAD_clinical)[13] <- 'stage'
for ( i in 1:nrow(COAD_clinical) ) {
  COAD_clinical$stage[i] <- ifelse(COAD_clinical$stage[i]=="Stage IA", "Stage I", 
                                   ifelse(COAD_clinical$stage[i]=="Stage IIA", "Stage II",
                                          ifelse(COAD_clinical$stage[i]=="Stage IIB", "Stage II",
                                                 ifelse(COAD_clinical$stage[i]=="Stage IIC", "Stage II",
                                                 ifelse(COAD_clinical$stage[i]=="Stage IIIA", "Stage III",
                                                        ifelse(COAD_clinical$stage[i]=="Stage IIIB", "Stage III",
                                                               ifelse(COAD_clinical$stage[i]=="Stage IIIC", "Stage III",
                                                                      ifelse(COAD_clinical$stage[i]=="Stage IVA", "Stage IV",
                                                                             ifelse(COAD_clinical$stage[i]=="Stage IVB", "Stage IV",COAD_clinical$stage[i])))))))))
}

# Load TCGA copy number data, and keep only COAD
COAD_copyN <- read.table("~/Documents/CNA/Data/TCGA/TCGAcn/TCGA_mastercalls.abs_segtabs.fixed.txt", header=T, skip=0, sep="\t")
COAD_copyN <- COAD_copyN[which((str_sub(COAD_copyN$Sample,1,str_length(COAD_copyN$Sample)-3)) %in% COAD_clinical$TCGA_barcode),c(1:4,9)]

# Remove patients missing CN data
COAD_copyN <- COAD_copyN[!is.na(COAD_copyN$Modal_Total_CN),]

# Rename columns
colnames(COAD_copyN) <- c('sample','chr','start','stop','cn')

# Pull cn data into element corresponding to sample ID
names <- unique(COAD_copyN$sample)
COAD.list <- list()
for ( i in 1:length(names) ) {
  COAD.list[[i]] <- COAD_copyN[which(COAD_copyN$sample == names[i]),]
  COAD.list[[i]] <- COAD.list[[i]][-1] #remove col holding sample as its sample per list element
}

# Assess ploidy and recentre if average CN is 2.8 or more
COAD.list <- ploidyRecentre(COAD.list,skipcol = 3, multi_cn = FALSE)

# Bin CN data to match training dataset
COADbinned.list <- newAlignBins(bins = car.info$start.stop, cn.list = COAD.list)

# AlignBins allows for MR sampling, so remove col 5
COADbinned.list <- lapply(COADbinned.list, function(x) {
  x[] <- apply(x,2,as.numeric)
  x
})

# If absolute copy number >=3, its a gain, if <=2 its a loss
COADbinned.list <- lapply(COADbinned.list, function(x) {
  x[,-c(1:4)] <- ifelse(x[,-c(1:4)] >= 3, 3, ifelse(x[,-c(1:4)] < 2, 1, 2))
  x
})

# Rename columns with sampleID and combine all patients into one dataframe
for ( i in 1:length(COADbinned.list) ) {
  for ( k in 5:ncol(COADbinned.list[[i]]) ) {
    colnames(COADbinned.list[[i]])[k] <- names[i]
  }
}
COADbinned <- COADbinned.list %>% purrr::reduce(full_join, by = c("chr","start","stop","bin"))
COADbinned <- COADbinned[,-1] # remove bin columns

# Use PullData functions to create necessary dataframes for analysis
TCGA.info <- PullDataInfo(rawdata = COADbinned)

saveRDS(TCGA.info, "~/Documents/CNA/Data/TCGA.info.rds")
saveRDS(COADbinned, "~/Documents/CNA/Data/COADbinned.rds")
saveRDS(COAD_clinical, "~/Documents/CNA/Data/COAD_clinical.rds")
