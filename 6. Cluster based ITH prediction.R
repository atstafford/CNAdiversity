# ASSESS IF ANEUPLOIDY PER CLUSTER CAN PREDICT ITH

# For each cluster count how many chromosomes are within, and identify chromosome
list1 <- list()
list2 <- list()
for ( i in 1:max(sig.hclust$cluster) ) {
  data <- sig.hclust[which(sig.hclust$cluster==i),]
  list1[[i]] <- paste(unique(data$chr), collapse = ',')
  list2[[i]] <- length(unique(data$chr))
}
cluster.info <- data.frame(cluster = 1:max(sig.hclust$cluster),
                           no.ofChr = do.call('rbind',list2),
                           chr = do.call('rbind',list1) )

# Add cluster group for facet
cluster.info$chrGroup <- as.numeric(as.character(cluster.info$chr))
cluster.info$chrGroup[is.na(cluster.info$chrGroup)] <- 23


# Add % of chromosome represented by each cluster

# For each cluster holding bins from one chromosome, claculate % of chromosome represented
list <- list()
for ( i in 1:length(sig.hclust$cluster) ) {
  chr <- as.numeric(cluster.info[which(cluster.info$cluster==i),]$chr)
  chrlength <- car.info$chr.end$end[chr] - car.info$chr.end$start[chr]
  number <- nrow(sig.hclust[which(sig.hclust$cluster==i),])
  list[[i]] <- number / chrlength
}

# Add % of chromosome into cluster.info dataframe
cluster.info <- cbind(cluster.info, data.frame(pcChrPerClus = unlist(list)))
cluster.info$pcChrPerClus <- cluster.info$pcChrPerClus*100


# Calculate the % aneuploidy per cluster per patient
list <- list()
l <- 1

# For a given sample
for ( k in 1:car.info$noSamples ) { 
  
  # For a given cluster, calculate the pc aneuploidy across the bins in the cluster
  for ( i in 1:max(sig.hclust$cluster) ) { 
    bins <- sig.hclust$bin[which(sig.hclust$cluster == i)]
    data <- car.matrices$diploid.aneu[bins, k]
    list[[l]] <- (sum(data[data==1])) / (length(data))
    l <- l + 1
  }
}

# Unlist into a dataframe
pcClusAneu <- data.frame(matrix(unlist(list), ncol = car.info$noSamples), stringsAsFactors = FALSE)
colnames(pcClusAneu) <- car.info$sampleIDs

# Add to cluster.info
cluster.info <- cbind(cluster.info, 
                      pcClusAneu = (apply(pcClusAneu, 1, mean))*100)


# Use % aneuploidy per cluster to predict diversity per patient

# Transpose to make sample per row, and a cluster per col
clustReg.in <- t(pcClusAneu)

# Add in ITH per patient
clustReg.in <- cbind(actualIth_train$actual, clustReg.in)
colnames(clustReg.in) <- c("actual",1:129)

# To store output
clustReg.out <- data.frame(cluster=1:max(cluster.info$cluster), coeff=NA, pval=NA, adjR2=NA, Modelpval=NA)

# For a each cluster, consider if % aneuploidy can predict ITH
for ( k in 2:ncol(clustReg.in) ) { 
  data <- data.frame(clustReg.in[,c(1,k)])
  data <- na.omit(data)
  colnames(data) <- c('ITH','predictor')
  
  # Regression can only be run when there are >1 factors present
  if ( length(unique(as.character(data$predictor))) == 1 ) { 
    next
  }
  reg <- lm(ITH ~ ., data)
  
  coeffs <- data.frame(t(summary(reg)$coefficients[,1])) 
  pvals <- data.frame(t(summary(reg)$coefficients[,4]))
  adjR2 <- round(summary(reg)$adj.r.squared,3)
  pval <- signif(lmp(reg),3)
  
  clustReg.out[(clustReg.out$cluster==(k-1)),c(2:5)] <- c(coeffs$predictor,pvals$predictor,adjR2,pval)
}

# Add col to indicate if aneuploidy in cluster is predictive (p<0.05)
clustReg.out$sig <- 'insig'
clustReg.out$sig[which(clustReg.out$Modelpval <= 0.05)] <- 'sig'

# Add FDR
clustReg.out$FDR = p.adjust(clustReg.out$Modelpval, method = "fdr")

# Designate significance of each cluster
clustReg.out$sigFDR <- 'insig'
clustReg.out$sigFDR[which(clustReg.out$BH <= 0.05)] <- 'sig'

# add coeffs and pvalx
cluster.info <- cbind(cluster.info, clustReg.out[,c(2:8)])
