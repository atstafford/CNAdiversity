# APPLY STEP.MODEL TO ADENOMAS
# We have 19 adenomas from the orinigal cohort and an additional 5 from the validation cohort

# Prepare input for predicting adenomas in original cohort

# Requires bin per col
multiReg.in_ad <- t(ad.raw[,-c(1:3)])

# Change to factor
multiReg.in_ad <- data.frame(apply(multiReg.in_ad, 2, as.character), check.names = FALSE)
multiReg.in_ad <- data.frame(lapply(multiReg.in_ad, factor, levels=c(2,1,3), labels=c('diploid','loss','gain')), check.names = FALSE)
names(multiReg.in_ad) = paste("bin_", names(multiReg.in_ad), sep="")
multiReg.in_ad <- dummy_cols(multiReg.in_ad, remove_selected_columns = TRUE, remove_first_dummy = TRUE)

# Predict ITH using step.model
predictedIth_ad <- predict(beta, multiReg.in_ad) 

# Enter actual and predicted ITH into a comparative dataframe
predictedIth_ad <- data.frame(cbind(actualIth_ad, predicted=predictedIth_ad))
predictedIth_ad$type <- 'ad'

# Load the 10 test validation set adenoma patients
test_ad1 <- read.table("~/Documents/CNA/Data/Validation/testSet_ad/Polyp.08.WGS.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_ad2 <- read.table("~/Documents/CNA/Data/Validation/testSet_ad/Polyp.02.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_ad3 <- read.table("~/Documents/CNA/Data/Validation/testSet_ad/Polyp.05.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_ad4 <- read.table("~/Documents/CNA/Data/Validation/testSet_ad/Polyp.09.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_ad5 <- read.table("~/Documents/CNA/Data/Validation/testSet_ad/Polyp.03.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")

testraw.list <- list(test_ad1,test_ad2,test_ad3,test_ad4,test_ad5)
rm(test_ad1,test_ad2,test_ad3,test_ad4,test_ad5)

# Convert rawdata to numeric
testraw.list <- lapply(testraw.list, function(x) {
  x <- apply(x, 2, function(x) as.numeric(as.character(x)))
  x
})

x <- testraw.list[[1]]

# Generate absolute copy number (minor + major)
testcn.list_ad <- lapply(testraw.list, function(x) {
  x <- x[,c(1,2,4)]
  colnames(x) <- c('chr','start','stop')
  x <- as.data.frame(x)
  x
})

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
  
  testcn.list_ad[[j]] <- cbind(testcn.list_ad[[j]], matrix(unlist(abscn), nrow = nrow(testraw.list[[j]])) )
}

# Assess ploidy and recentre if average CN is 2.8 or more
testcn.list_ad <- ploidyRecentre(testcn.list_ad, skipcol = 3, multi_cn = TRUE)

# Bin cn data to match training dataset
testcnBinned.list_ad <- newAlignBins(bins = ad.info$start.stop, cn.list = testcn.list_ad)

# If absolute copy number >=3, its a gain, if <=2 its a loss
testcnBinned.list_ad <- lapply(testcnBinned.list_ad, function(x) {
  x[,-c(1:4)] <- ifelse(x[,-c(1:4)] >= 3, 3, ifelse(x[,-c(1:4)] < 2, 1, 2))
  x
})

# Rename columns with sampleID and combine all patients into one dataframe
for ( i in 1:length(testcnBinned.list_ad) ) {
  for ( k in 5:ncol(testcnBinned.list_ad[[i]]) ) {
    colnames(testcnBinned.list_ad[[i]])[k] <- paste(i, k-4, sep = '.')
  }
}

testcnBinned_ad <- testcnBinned.list_ad %>% purrr::reduce(full_join, by = c("chr","start","stop","bin"))
testcnBinned_ad <- testcnBinned_ad[-1]

valad.info <- PullDataInfo(testcnBinned_ad)
valad.diversity <- PullDataDiversity(testcnBinned_ad, dataInfo = valad.info)
valad.clonality <- PullDataClonality(testcnBinned_ad, dataInfo = valad.info)

# Prepare input for prediction

# Transpose as input table requires cols=bins and convert to character
test.input <- data.frame(t(testcnBinned[,-c(1:3)]), check.names = FALSE)
test.input <- ifelse(test.input == 3, "gain", ifelse(test.input == 1, "loss", "diploid"))
names(test.input) <- 1:ncol(test.input)
test.input <- dummy_cols(test.input, remove_selected_columns = TRUE, remove_first_dummy = TRUE)
names(test.input) = paste("bin_", names(test.input), sep="")
test.input <- data.frame(test.input, check.names = FALSE)
#test.input$bin_2385_gain <- 0
#test.input$bin_2371_gain <- 0
#test.input$bin_2374_gain <- 0

# List actual patient ITH, repeat each value twice as there are two samples
actualIth_valad <- data.frame(patient = lapply(data.frame(patient=valad.info$patientIDs), rep, valad.info$sampPerPatient),
                           sample = valad.info$sampleIDs,
                           actual = lapply(data.frame(actual=valad.diversity$pic.frac$pic.frac), rep, valad.info$sampPerPatient))


# Predict ITH using step.model
new_model = beta$model[,names(beta$model) != 'bin_2385_gain']
fit.new = betareg(bin_ITH ~., data=new_model)
new_model = fit.new$model[,names(fit.new$model) != 'bin_2371_gain']
fit.new = betareg(bin_ITH ~., data=new_model)
new_model = fit.new$model[,names(fit.new$model) != 'bin_2374_gain']
fit.new = betareg(bin_ITH ~., data=new_model)
predictedIth_valad <- predict(fit.new, test.input)

predictedIth_valad <- data.frame(cbind(actualIth_valad, predicted=predictedIth_valad))
predictedIth_valad$type <- 'valad'

x <- rbind(predictedIth_ad, predictedIth_valad)

# PLOT
plot <- ggplot(data = x, aes(x = actual, y = predicted, color=type)) +
  geom_line(aes(group = patient), size=0.2, colour = "#003366") +
  geom_point(size= 5) +
  #geom_point(fill = "#003366", colour = "#003366", size = 5) +
  geom_line(aes(x = actual, y = actual), linetype = "dashed") +
  #xlim(0,0.34) +
  #ylim(0,0.34) +
  xlab("Patient CNA diversity") +
  ylab("Patient CNA diversity (predicted)") +
  theme(plot.margin = unit(c(t=0,r=0,b=0,l=0), "cm"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.title = element_text(size=24, colour='black'),
        axis.text = element_text(size=24, colour='black'),
        legend.position = "top")

jpeg('tempfig.jpeg', width = 1000, height = 1000)
plot
dev.off()

