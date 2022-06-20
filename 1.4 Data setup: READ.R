# TCGA PROCESSING ####

# Load clinical data for CRC (COAD)
READ_clinical <- read_excel("~/Documents/CNA/Data/TCGA/TCGAcn/Clinical COAD.xlsx")
colnames(READ_clinical)[c(1,2)] <- c("UUID_copyN","TCGA_barcode")

# Remove duplicates
READ_clinical <- READ_clinical[!duplicated(READ_clinical$TCGA_barcode), ]

# Remove columns with no information
READ_clinical <- READ_clinical[,colSums(READ_clinical=="'--") < nrow(READ_clinical)]
READ_clinical <- READ_clinical[,-c(3,14,15,19,18,20,23,28,32,33)]

# Load additional data and merge
READ_MSI <- read_excel("~/Documents/CNA/Data/TCGA/TCGAcn/TCGA_MSIstatus.xlsx", col_names = TRUE, skip =1)
READ_MSI <- READ_MSI[READ_MSI$Organ=="READ",]
READ_clinical$MSI <- unlist(READ_MSI[match(READ_clinical$TCGA_barcode, READ_MSI$`TCGA Participant Barcode`), 14])
READ_clinical$CMS <- unlist(READ_MSI[match(READ_clinical$TCGA_barcode, READ_MSI$`TCGA Participant Barcode`), 27])
READ_clinical$MSI[is.na(READ_clinical$MSI)] <- 'NA'

# Change stage names to remove A/B
colnames(READ_clinical)[13] <- 'stage'
for ( i in 1:nrow(READ_clinical) ) {
  READ_clinical$stage[i] <- ifelse(READ_clinical$stage[i]=="Stage IA", "Stage I", 
                                   ifelse(READ_clinical$stage[i]=="Stage IIA", "Stage II",
                                          ifelse(READ_clinical$stage[i]=="Stage IIB", "Stage II",
                                                 ifelse(READ_clinical$stage[i]=="Stage IIC", "Stage II",
                                                 ifelse(READ_clinical$stage[i]=="Stage IIIA", "Stage III",
                                                        ifelse(READ_clinical$stage[i]=="Stage IIIB", "Stage III",
                                                               ifelse(READ_clinical$stage[i]=="Stage IIIC", "Stage III",
                                                                      ifelse(READ_clinical$stage[i]=="Stage IVA", "Stage IV",
                                                                             ifelse(READ_clinical$stage[i]=="Stage IVB", "Stage IV",READ_clinical$stage[i])))))))))
}

# Load TCGA copy number data, and keep only COAD
READ_copyN <- read.table("~/Documents/CNA/Data/TCGA/TCGAcn/TCGA_mastercalls.abs_segtabs.fixed.txt", header=T, skip=0, sep="\t")
READ_copyN <- READ_copyN[which((str_sub(READ_copyN$Sample,1,str_length(READ_copyN$Sample)-3)) %in% READ_clinical$TCGA_barcode),c(1:4,9)]

# Remove patients missing CN data
READ_copyN <- READ_copyN[!is.na(READ_copyN$Modal_Total_CN),]

# Rename columns
colnames(READ_copyN) <- c('sample','chr','start','stop','cn')

# Pull cn data into element corresponding to sample ID
names <- unique(READ_copyN$sample)
READ.list <- list()
for ( i in 1:length(names) ) {
  READ.list[[i]] <- READ_copyN[which(READ_copyN$sample == names[i]),]
  READ.list[[i]] <- READ.list[[i]][-1] #remove col holding sample as its sample per list element
}

# Assess ploidy and recentre if average CN is 2.8 or more
#COAD.list <- ploidyRecentre(COAD.list,skipcol = 3, multi_cn = FALSE)

# Assess ploidy and recentre 
for ( i in 1:length(READ.list) ) {
  weights <- READ.list[[i]][,3] - READ.list[[i]][,2]
  sumweight <- sum(weights)
  
  for ( j in 4:ncol(READ.list[[i]]) ) {
    ploidy <- round(sum(READ.list[[i]][ ,j]*weights, na.rm = T)/sumweight, 1)
    
    if ( ploidy > 2.5 ) {
      distance <- abs(2.5 - ploidy)
      print(paste(i, ploidy))
      READ.list[[i]][ ,j] <- READ.list[[i]][ ,j] - distance                
    }
    else {next}
  }
}
x <- READ.list[[1]]

# Bin CN data to match training dataset
READbinned.list <- newAlignBins(bins = car.info$start.stop, cn.list = READ.list)
donttouch <- READbinned.list

READbinned.list <- lapply(READbinned.list, function(x) {
  x[] <- apply(x,2,as.numeric)
  x
})

# Rename columns with sampleID and combine all patients into one dataframe
for ( i in 1:length(READbinned.list) ) {
  for ( k in 5:ncol(READbinned.list[[i]]) ) {
    colnames(READbinned.list[[i]])[k] <- names[i]
  }
}

# make dataframe
READbinned.list <- lapply(READbinned.list, function(x) {
  x <- as.data.frame(x)
  x[,5][is.nan(x[,5])] <- NA
  x
})

# Fill in missing
for ( i in 1:length(READbinned.list) ) {
  if ( sum(is.na(READbinned.list[[i]][hg19predictors$bin,-c(1:4)])) == 0 ) { # no nas in predictors
    next
  } else {
    for ( j in 1:length(hg19predictors$bin) ) {
      bin <- hg19predictors$bin[j]
      if ( is.na(READbinned.list[[i]][bin,5]) ) {
        #print(bin)
        READbinned.list[[i]][bin,5] <- mean( c(READbinned.list[[i]][(bin-1),5], READbinned.list[[i]][(bin+1),5]) , na.rm=T)
      } else {next}
    }
  }
}

x <- READbinned.list %>% purrr::reduce(full_join, by = c("chr","start","stop","bin"))
x <- x[c(hg19predictors$bin), -c(1:4)]
x <- x[ , colSums(is.na(x)) > 0]


# If absolute copy number >=3, its a gain, if <=2 its a loss
READbinned.list <- lapply(READbinned.list, function(x) {
  x[,-c(1:4)] <- ifelse(x[,-c(1:4)] >= 3, 3, ifelse(x[,-c(1:4)] < 2, 1, 2))
  x <- as.data.frame(x)
  x
})

READbinned <- READbinned.list %>% purrr::reduce(full_join, by = c("chr","start","stop","bin"))
READbinned <- READbinned[,-1] # remove bin columns

# Use PullData functions to create necessary dataframes for analysis
READ.info <- PullDataInfo(rawdata = READbinned)

saveRDS(READ.info, "~/Documents/CNA/Data/READ.info.rds")
saveRDS(READbinned, "~/Documents/CNA/Data/READbinned.rds")
saveRDS(READ_clinical, "~/Documents/CNA/Data/READ_clinical.rds")
