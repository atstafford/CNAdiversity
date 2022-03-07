# BUILD BETA REGRESSION MODELS: A MULTI-BIN MODEL TO PREDICT ITH FROM SINGLE SAMPLE----

# Perform stepwise selection to maximise Akaike information criteria, using one bin per cluster

# Bind the three dataframes holding the univariate regression data
candidate.bins <- do.call('rbind',uniReg.out.list[c(2:3)])

# Extract bins that are significant (p<0.05) before the FDR was applied as so few bins was sig after FDR
candidate.bins <- candidate.bins[which(candidate.bins$pval <= 0.05),]

# pull clonality
candidate.bins$pcSubclonal <- NA
for ( i in 1:nrow(candidate.bins) ) {
  bin <- candidate.bins$bin[i]
  CNA <- candidate.bins$CNA[i]
  
  if ( CNA=='gain') {
    candidate.bins$pcSubclonal[i] <- car.clonality$pcSubclonal$gain[which(car.clonality$pcSubclonal$bin==bin)]
  }
  if ( CNA=='loss') {
    candidate.bins$pcSubclonal[i] <- car.clonality$pcSubclonal$loss[which(car.clonality$pcSubclonal$bin==bin)]
  }
}

candidate.bins$cluster <- sig.hclust[match(candidate.bins$bin, sig.hclust$bin), 2]
rownames(candidate.bins) <- 1:nrow(candidate.bins)

# Assign unique clusters to bins with no cluster
new.clust <- max(candidate.bins$cluster[is.finite(candidate.bins$cluster)]) + 1
for (i in 1:nrow(candidate.bins) ) {
  if ( is.na(candidate.bins$cluster[i]) ) {
    candidate.bins$cluster[i] <- new.clust
    new.clust <- new.clust + 1
  }
}

# Select one representive bin per cluster
list <- list()
for ( i in 1:max(candidate.bins$cluster) ) {
  data <- candidate.bins[which(candidate.bins$cluster == i),]
  list[[i]] <- data[which.max(abs(data$coeff)),]
}
representitive.bins <- do.call('rbind',list)


# BUILD FOR 616 BINS------------------------------------------------------------
# Build input
multiReg.in_616 <- t(car.raw[,-c(1:3)]) # The input matrix requires data on loss/gain/diploid, with a bin per column
multiReg.in_616 <- data.frame(apply(multiReg.in_616, 2, as.character), check.names = FALSE)
multiReg.in_616 <- data.frame(lapply(multiReg.in_616, factor, levels=c(2,1,3), labels=c('diploid','loss','gain')), check.names = FALSE)
rownames(multiReg.in_616) <- car.info$sampleIDs
multiReg.in_616 <- multiReg.in_616[,c(candidate.bins$bin)] # Keep only the candidate bins (columns)

names(multiReg.in_616) = paste("bin_", names(multiReg.in_616), sep="")
multiReg.in_616 = multiReg.in_616 %>%
  mutate(across(everything(), as.character))

# Make dummy for beta
'%!in%' <- function(x,y)!('%in%'(x,y))
multiReg.in_616 <- dummy_cols(multiReg.in_616, remove_selected_columns = TRUE, remove_first_dummy = TRUE)
multiReg.in_616$bin_ITH <- actualIth_train$actual
remove <- attributes(alias(lm(bin_ITH ~ ., data = multiReg.in_616))$Complete)$dimnames[[1]]
multiReg.in_616 <- multiReg.in_616[ ,colnames(multiReg.in_616) %!in% remove]

# Betareg
beta_616 <- betareg(bin_ITH ~., data = multiReg.in_616)
summary(beta_616, type = "pearson")
x <- data.frame(vif(beta_616))


# BULD FOR 52-------------------------------------------------------------------
# Build input
multiReg.in_52 <- t(car.raw[,-c(1:3)]) # The input matrix requires data on loss/gain/diploid, with a bin per column
multiReg.in_52 <- data.frame(apply(multiReg.in_52, 2, as.character), check.names = FALSE)
multiReg.in_52 <- data.frame(lapply(multiReg.in_52, factor, levels=c(2,1,3), labels=c('diploid','loss','gain')), check.names = FALSE)
rownames(multiReg.in_52) <- car.info$sampleIDs
multiReg.in_52 <- multiReg.in_52[,c(representitive.bins$bin)] # Keep only the candidate bins (columns)

names(multiReg.in_52) = paste("bin_", names(multiReg.in_52), sep="")
multiReg.in_52 = multiReg.in_52 %>%
  mutate(across(everything(), as.character))

# Make dummy for beta
multiReg.in_52 <- dummy_cols(multiReg.in_52, remove_selected_columns = TRUE, remove_first_dummy = TRUE)
multiReg.in_52$bin_ITH <- actualIth_train$actual
remove <- attributes(alias(lm(bin_ITH ~ ., data = multiReg.in_52))$Complete)$dimnames[[1]]
multiReg.in_52 <- multiReg.in_52[ ,colnames(multiReg.in_52) %!in% remove]

# Betareg
beta_52 <- betareg(bin_ITH ~., data = multiReg.in_52)
summary(beta_52, type = "pearson")
x <- data.frame(vif(beta_52))


# PERFORM SELCTION ON 52--------------------------------------------------------
x <- multiReg.in_52[-85]
y <- actualIth_train$actual
selection_for <- betaselect(x, y, criterion = "AIC",link = "logit", method = "forward", plotit = FALSE)
keep <- selection_for$variable
multiReg.in_52for <- multiReg.in_52[ ,colnames(multiReg.in_52) %in% keep]
multiReg.in_52for$bin_ITH <- actualIth_train$actual
beta_52for <- betareg(bin_ITH ~., data = multiReg.in_52for)
summary(beta_52for, type = "pearson")
x <- data.frame(vif(beta_52for))

x <- multiReg.in_52[-85]
y <- actualIth_train$actual
selection_back <- betaselect(x, y, criterion = "AIC",link = "logit", method = "backward", plotit = FALSE)
keep <- selection_back$variable
multiReg.in_52back <- multiReg.in_52[ ,colnames(multiReg.in_52) %in% keep]
multiReg.in_52back$bin_ITH <- actualIth_train$actual
beta_52back <- betareg(bin_ITH ~., data = multiReg.in_52back)
summary(beta_52back, type = "pearson")
x <- data.frame(vif(beta_52back))

remove <- c("bin_1413_gain","bin_1418_gain","bin_1112_gain","bin_1173_gain","bin_1083_gain","bin_1085_gain")
multiReg.in_52backrem <- multiReg.in_52back[ ,colnames(multiReg.in_52back) %!in% remove]
beta_52backrem <- betareg(bin_ITH ~., data = multiReg.in_52backrem)
summary(beta_52backrem, type = "pearson")
x <- data.frame(vif(beta_52backrem))
round(AIC(beta_52backrem),2)

# ASSESS------------------------------------------------------------------------
# pull predictors
beta <- beta_52backrem
multiReg.in <- multiReg.in_52backrem

x <- rownames(data.frame(summary(beta)$coefficients))[-1]
hg19predictors <- data.frame(bin=as.numeric(str_extract_all(x, "[0-9]+")), cna=str_sub(x,-4,-1))
hg19predictors <- merge(unique(hg19predictors), car.info$start.stop)
hg19predictors$cluster <- candidate.bins[match(hg19predictors$bin, candidate.bins$bin), 11] 
hg19predictors$pcSubclonal <- candidate.bins[match(hg19predictors$bin, candidate.bins$bin), 10] 
both <- hg19predictors$bin[which(duplicated(hg19predictors$bin)==T)]
hg19predictors$cna[which(hg19predictors$bin %in% hg19predictors$bin[which(duplicated(hg19predictors$bin)==T)] )] <- 'both'

hg19predictors <- hg19predictors[c(1,3:7,2)]
hg19predictors <- hg19predictors[!duplicated(hg19predictors),]

psuedoR2 <- round(beta$pseudo.r.squared,3)
loglik <- round(beta$loglik,2)
aic <- round(AIC(beta),2)
predictedIth_train <- data.frame(predicted = predict(beta, multiReg.in))

# Enter actual and predicted ITH into a comparative dataframe
predictedIth_train <- data.frame(cbind(actualIth_train, predictedIth_train))
predictedIth_train$type <- 'train'

# Add two labels for plot
ggplot(data = predictedIth_train, aes(x = actual, y = predicted)) +
  geom_line(aes(group = patient), size=0.2, colour = "#003366") +
  geom_point(fill = "#003366", colour = "#003366", size = 6) +
  geom_line(aes(x = actual, y = actual), linetype = "dashed") +
  xlab("Patient CNA diversity") +
  ylab("Patient CNA diversity (predicted)") +
  #xlim(0,0.5) +
  #ylim(0,0.5) +
  theme(plot.margin = unit(c(t=0,r=0,b=0,l=0), "cm"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.title = element_text(size=24, colour='black'),
        axis.text = element_text(size=24, colour='black'),
        legend.position = "none") +
  annotate('text', x = c(0), y = c(0.500,0.450), label = c(paste('pseudoR[adj]^2 ==', psuedoR2), paste('AIC ==', aic)), parse=TRUE, size = 8, hjust = 0, vjust = 1)

jpeg('tempfig.jpeg', width = 1000, height = 1000)
plot
dev.off()

# save
saveRDS(candidate.bins, "~/Documents/CNA/Data/candidate.bins.rds")
saveRDS(representitive.bins, "~/Documents/CNA/Data/representitive.bins.rds")
saveRDS(beta_52backrem, "~/Documents/CNA/Data/beta_52backrem.rds")
saveRDS(hg19predictors, "~/Documents/CNA/Data/hg19predictors.rds")
saveRDS(predictedIth_train, "~/Documents/CNA/Data/predictedIth_train.rds")

