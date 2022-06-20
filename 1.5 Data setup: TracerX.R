# VALIDATION COHORT ####

# In order to use upcoming functions, data must be in format of a list with a dataframe per patient, with columns: 
# chr | start | stop | sample1 | sample2...

# Read in the datasets
tracerx <- read_excel("~/Documents/CNA/Data/TracerX/cn_data.xlsx", sheet = 5)

# Pull cn data into element corresponding to patient ID
names <- unique(tracerx$sample)
patient <- unique(sub("\\-.*", "", names))

tracerx.list <- list()
for ( i in 1:length(names) ) {
  tracerx.list[[i]] <- tracerx[which(tracerx$sample == names[i]),]
  tracerx.list[[i]] <- tracerx.list[[i]][ ,c(2,3,4,6)]
  colnames(tracerx.list[[i]]) <- c("chr", "start", "stop", names[i])
}

x <- tracerx.list[[1]]

# Convert data to numeric
tracerx.list <- lapply(tracerx.list, function(x) {
  x <- data.frame(apply(x, 2, function(x) as.numeric(as.character(x))))
  x
})
x <- tracerx.list[[1]]

# Assess ploidy and recentre 
for ( i in 1:length(tracerx.list) ) {
  weights <- tracerx.list[[i]][,3] - tracerx.list[[i]][,2]
  sumweight <- sum(weights)
  
  for ( j in 4:ncol(tracerx.list[[i]]) ) {
    ploidy <- round(sum(tracerx.list[[i]][ ,j]*weights, na.rm = T)/sumweight, 1)
    
    if ( ploidy > 2.5 ) {
      distance <- abs(2.5 - ploidy)
      tracerx.list[[i]][ ,j] <- tracerx.list[[i]][ ,j] - distance                
    }
    else {next}
  }
}
x <- tracerx.list[[1]]

# Bin CN data to match training dataset ONLY IF HG19
tracerxBinned.list <- newAlignBins(bins = car.info$start.stop, cn.list = tracerx.list)
tracerxBinned.list <- lapply(tracerxBinned.list, function(x) {
  x <- apply(x, 2, function(x) as.numeric(as.character(x)))
  x <- as.data.frame(x)
  x
})
x <- tracerxBinned.list[[1]]
donttouch <- tracerxBinned.list


# Pull samples from same patient into dataframe per patient
tracerxcn.list <- list()
for ( i in 1:length(patient) ) {
  x <- which(sub("\\-.*","",names) == patient[i])
  tracerxcn.list[[i]] <- tracerxBinned.list[x] %>% purrr::reduce(full_join, by = c("chr","start","stop","bin"))
  tracerxcn.list[[i]] <- tracerxcn.list[[i]][-1]
}

# Convert data to numeric
tracerxcn.list <- lapply(tracerxcn.list, function(x) {
  x <- data.frame(apply(x, 2, function(x) as.numeric(as.character(x))))
  x
})

# Drop patient without MR
drop <- which(sapply(tracerxcn.list, function(x) ncol(x) < 5) == TRUE)
tracerxcn.list <- tracerxcn.list[-drop]
x <- tracerxcn.list[[1]]

# Fill in missing
x <- tracerxcn.list %>% purrr::reduce(full_join, by = c("chr","start","stop"))
x <- x[ ,-c(1:3)]
for ( i in 1:ncol(x) ) {
  for ( j in 1:nrow(x) ) {
    if ( !is.na(x[j,i]) ) {
      next
    }
    else {
      x[j,i] <- mean(x[c(j-1,j+1),i], na.rm = T)
    }
  }
}

# If absolute copy number >=3, its a gain, if <=2 its a loss
x <- ifelse(x >= 3, 3, ifelse(x < 2, 1, 2))
tracerxBinned <- cbind(car.info$start.stop, x)
tracerxBinned <- tracerxBinned[-1]

# Use PullData functions to create necessary dataframes for analysis
tracerx.info <- PullDataInfo(rawdata = tracerxBinned)
tracerx.clonality <- PullDataClonality(rawdata = tracerxBinned, dataInfo = tracerx.info)
tracerx.diversity <- PullDataDiversity(rawdata = tracerxBinned, dataInfo = tracerx.info)

actualIth_tracerx <- data.frame(patient = lapply(data.frame(patient=tracerx.info$patientIDs), rep, tracerx.info$sampPerPatient),
                             sample = tracerx.info$sampleIDs,
                             actual = lapply(data.frame(actual=tracerx.diversity$pic.frac$pic.frac), rep, tracerx.info$sampPerPatient))

lung.matrices <- genMatrices(tracerxBinned)

# Save
saveRDS(tracerx.info, "~/Documents/CNA/Data/tracerx.info.rds")
saveRDS(tracerx.clonality, "~/Documents/CNA/Data/tracerx.clonality.rds")
saveRDS(tracerx.diversity, "~/Documents/CNA/Data/tracerx.diversity.rds")
saveRDS(lung.matrices, "~/Documents/CNA/Data/lung.matrices.rds")
saveRDS(tracerxBinned, "~/Documents/CNA/Data/tracerxBinned.rds")
saveRDS(actualIth_tracerx, "~/Documents/CNA/Data/actualIth_tracerx.rds")

