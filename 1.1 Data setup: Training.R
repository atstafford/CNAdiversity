# TRAINING COHORT ####

# The original dataset from Cross et al has CN data across 2694 bins for 19 adenomas and 81 carcinomas
# Two sample from each lesion we taken
# Columns are: chr | start | stop | sample1 | sample2...
# Rows will contain the raw copy number assignment per bin (1=loss, 2=diploid, 3=gain). 
# From columns 4 onwards, sample identifiers must be constructed as: patientName.sampleNum

# Load data
rawdata <- read.table("~/Documents/CNA/Data/Training/sWGS Cross rename with type.txt", header=T, skip=2, sep="\t")

# Split into 2 separate raw data files (ad.raw and car.raw) 
ad.raw<- rawdata[,c(1:3,which(sub("^([[:alpha:]]*).*", "\\1", colnames(rawdata))=="A"))]
car.raw <- rawdata[,c(1:3,which(sub("^([[:alpha:]]*).*", "\\1", colnames(rawdata))=="C"))]

# Use PullData functions to create necessary dataframes for analysis
ad.info <- PullDataInfo(rawdata = ad.raw)
ad.clonality <- PullDataClonality(rawdata = ad.raw, dataInfo = ad.info)
ad.diversity <- PullDataDiversity(rawdata = ad.raw, dataInfo = ad.info)

car.info <- PullDataInfo(rawdata = car.raw)
car.clonality <- PullDataClonality(rawdata = car.raw, dataInfo = car.info)
car.diversity <- PullDataDiversity(rawdata = car.raw, dataInfo = car.info)

actualIth_train <- data.frame(patient = lapply(data.frame(patient=car.info$patientIDs), rep, car.info$sampPerPatient),
                              sample = car.info$sampleIDs,
                              actual = lapply(data.frame(actual=car.diversity$pic.frac$pic.frac), rep, car.info$sampPerPatient))

actualIth_ad <- data.frame(patient = lapply(data.frame(patient=ad.info$patientIDs), rep, ad.info$sampPerPatient),
                              sample = ad.info$sampleIDs,
                              actual = lapply(data.frame(actual=ad.diversity$pic.frac$pic.frac), rep, ad.info$sampPerPatient))

# Use genMatrices to generate four matrices from the raw data (where 1=loss, 2=noCNA, 3=gain).
ad.matrices <- genMatrices(ad.raw)
car.matrices <- genMatrices(car.raw)

# Save
saveRDS(ad.raw, "~/Documents/CNA/Data/ad.raw.rds")
saveRDS(ad.info, "~/Documents/CNA/Data/ad.info.rds")
saveRDS(ad.clonality, "~/Documents/CNA/Data/ad.clonality.rds")
saveRDS(ad.diversity, "~/Documents/CNA/Data/ad.diversity.rds")

saveRDS(car.raw, "~/Documents/CNA/Data/car.raw.rds")
saveRDS(car.info, "~/Documents/CNA/Data/car.info.rds")
saveRDS(car.clonality, "~/Documents/CNA/Data/car.clonality.rds")
saveRDS(car.diversity, "~/Documents/CNA/Data/car.diversity.rds")
saveRDS(actualIth_train, "~/Documents/CNA/Data/actualIth_train.rds")
saveRDS(actualIth_ad, "~/Documents/CNA/Data/actualIth_ad.rds")

saveRDS(ad.matrices, "~/Documents/CNA/Data/ad.matrices.rds")
saveRDS(car.matrices, "~/Documents/CNA/Data/car.matrices.rds")

