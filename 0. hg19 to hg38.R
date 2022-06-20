# Pull hg19 IDs (bin names)
#qtlids <- predictors$bin

# Create input data: chromosomes must be chr1, chr2 etc
iddf <- car.info$start.stop
iddf$chr <- as.character(paste('chr',iddf$chr,sep = "", collapse = NULL))
idrange <- makeGRangesFromDataFrame(iddf)

# Inport chain for hg19 to hg38
chain <- import.chain("~/Documents/CNA/Data/hg19Tohg38.over.chain")

# Liftover
hg38ids <- liftOver(idrange, chain)
numfound <- unlist(lapply(hg38ids,length))
hg38ids <- as.data.frame(hg38ids)
hg38ids$group_name <- rep(car.info$start.stop$bin, numfound)

# Create table showing hg38 loci
hg38predictors <- data.frame(bin=car.info$start.stop$bin, chr=NA, start=NA, stop=NA)

for ( i in 1:length(car.info$start.stop$bin) ) {
  wd <- hg38ids[which(hg38ids$group_name==car.info$start.stop$bin[i]),]
  
  hg38predictors$chr[i] <- sub('...','',wd$seqnames[1])
  hg38predictors$start[i] <- wd$start[1]
  hg38predictors$stop[i] <- wd$end[nrow(wd)]
}

saveRDS(hg38predictors, "~/Documents/CNA/Data/hg38predictors.rds")