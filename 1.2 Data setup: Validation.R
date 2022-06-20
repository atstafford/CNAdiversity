# VALIDATION COHORT ####

# In order to use upcoming functions, data must be in format of a list with a dataframe per patient, with columns: 
# chr | start | stop | sample1 | sample2...

# Read in the datasets
test_car1 <- read.table("~/Documents/CNA/Data/Validation/testSet_car/Set.01.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_car2 <- read.table("~/Documents/CNA/Data/Validation/testSet_car/Set.02.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_car3 <- read.table("~/Documents/CNA/Data/Validation/testSet_car/Set.03.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_car4 <- read.table("~/Documents/CNA/Data/Validation/testSet_car/Set.04.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_car5 <- read.table("~/Documents/CNA/Data/Validation/testSet_car/Set.05.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_car6 <- read.table("~/Documents/CNA/Data/Validation/testSet_car/Set.06.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_car7 <- read.table("~/Documents/CNA/Data/Validation/testSet_car/Set.07.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_car8 <- read.table("~/Documents/CNA/Data/Validation/testSet_car/Set.08.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_car9p <- read.table("~/Documents/CNA/Data/Validation/testSet_car/Set.09.Proximal.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_car9d <- read.table("~/Documents/CNA/Data/Validation/testSet_car/Set.09.Distal.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_car10 <- read.table("~/Documents/CNA/Data/Validation/testSet_car/Set.10.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")

# Load patient datasets into a list
testraw.list <- list(test_car1,test_car2,test_car3,test_car4,test_car5,test_car6,test_car7,test_car8,test_car9p,test_car9d,test_car10)
rm(test_car1,test_car2,test_car3,test_car4,test_car5,test_car6,test_car7,test_car8,test_car9p,test_car9d,test_car10)

# Convert data to numeric
testraw.list <- lapply(testraw.list, function(x) {
  x <- apply(x, 2, function(x) as.numeric(as.character(x)))
  x
})

# For each patient, create dataframe to hold absolute copy number per chr/start/stop
testabs.list <- lapply(testraw.list, function(x) {
  x <- x[,c(1,2,4)]
  colnames(x) <- c('chr','start','stop')
  x <- as.data.frame(x)
  x
})

# Generate absolute copy number (minor + major)
for ( j in 1:length(testraw.list)) {
  i <- 5
  l <- 1
  abscn <- list()
  
  while ( i < ncol(testraw.list[[j]]) ) {
    
    for ( k in 1:nrow(testraw.list[[j]]) ) {
      abscn[[l]] <- sum(testraw.list[[j]][k,i], testraw.list[[j]][k,i+1])
      l <- l + 1
    }
    i <- i + 2
  }
  
  testabs.list[[j]] <- cbind(testabs.list[[j]], matrix(unlist(abscn), nrow = nrow(testraw.list[[j]])) )
}

# Assess ploidy and recentre 
for ( i in 1:length(testabs.list) ) {
  weights <- testabs.list[[i]][,3] - testabs.list[[i]][,2]
  sumweight <- sum(weights)
  
  for ( j in 4:ncol(testabs.list[[i]]) ) {
    ploidy <- round(sum(testabs.list[[i]][ ,j]*weights, na.rm = T)/sumweight, 1)
    
    if ( ploidy > 2.5 ) {
      distance <- abs(2.5 - ploidy)
      print(paste(i, ploidy))
      testabs.list[[i]][ ,j] <- testabs.list[[i]][ ,j] - distance                
    }
    else {next}
  }
}
x <- testabs.list[[1]]

# Bin CN data to match training dataset ONLY IF HG19
testcnBinned.list <- newAlignBins(bins = car.info$start.stop, cn.list = testabs.list)
testcnBinned.list <- lapply(testcnBinned.list, function(x) {
  x <- apply(x, 2, function(x) as.numeric(as.character(x)))
  x
})
x <- testcnBinned.list[[1]]

# Rename columns with sampleID and combine all patients into one dataframe
for ( i in 1:length(testcnBinned.list) ) {
  for ( k in 5:ncol(testcnBinned.list[[i]]) ) {
    colnames(testcnBinned.list[[i]])[k] <- paste(i, k-4, sep = '.')
  }
}

# make dataframe
testcnBinned.list <- lapply(testcnBinned.list, function(x) {
  x <- as.data.frame(x)
  x
})

# Fill in missing
x <- testcnBinned.list %>% purrr::reduce(full_join, by = c("chr","start","stop","bin"))
x <- x[c(hg19predictors$bin,2375,2370,2373), ]
x <- x[ , colSums(is.na(x)) > 0]
testcnBinned.list[[1]][2374, -c(1:4)] <- testcnBinned.list[[1]][2375, -c(1:4)]
testcnBinned.list[[3]][2374, -c(1:4)] <- testcnBinned.list[[3]][2375, -c(1:4)]

# If absolute copy number >=3, its a gain, if <=2 its a loss
testcnBinned.list <- lapply(testcnBinned.list, function(x) {
  x[,-c(1:4)] <- ifelse(x[,-c(1:4)] >= 3, 3, ifelse(x[,-c(1:4)] < 2, 1, 2))
  x <- as.data.frame(x)
  x
})

testcnBinned <- testcnBinned.list %>% purrr::reduce(full_join, by = c("chr","start","stop","bin"))

# Use PullData functions to create necessary dataframes for analysis
val.info <- PullDataInfo(rawdata = testcnBinned[,-1])
val.clonality <- PullDataClonality(rawdata = testcnBinned[,-1], dataInfo = val.info)
val.diversity <- PullDataDiversity(rawdata = testcnBinned[,-1], dataInfo = val.info)

actualIth_test <- data.frame(patient = lapply(data.frame(patient=val.info$patientIDs), rep, val.info$sampPerPatient),
                              sample = val.info$sampleIDs,
                              actual = lapply(data.frame(actual=val.diversity$pic.frac$pic.frac), rep, val.info$sampPerPatient))

# Save
saveRDS(val.info, "~/Documents/CNA/Data/val.info.rds")
saveRDS(val.clonality, "~/Documents/CNA/Data/val.clonality.rds")
saveRDS(val.diversity, "~/Documents/CNA/Data/val.diversity.rds")
saveRDS(testcnBinned, "~/Documents/CNA/Data/testcnBinned.rds")
saveRDS(testcnBinned.list, "~/Documents/CNA/Data/testcnBinned.list.rds")
saveRDS(testabs.list, "~/Documents/CNA/Data/testabs.list.rds")
saveRDS(actualIth_test, "~/Documents/CNA/Data/actualIth_test.rds")

