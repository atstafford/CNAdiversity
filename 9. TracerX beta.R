# CNA cooccurance ---------------------------------------------------------------
x <- lung.matrices$diploid.aneu
names <- colnames(x)
patient <- unique(sub("\\..*","",names))

list <- list()
for ( i in 1:length(patient) ) {
  wd <- x[ ,which(sub("\\..*","",names) == patient[i])]
  list[[i]] <- data.frame(rowMeans(wd, na.rm = T))
  colnames(list[[i]]) <- patient[i]
}

x <- do.call("cbind",list)
x <- data.matrix(x, rownames.force = NA)
x <- t(x)
lung.hiclust.input <- cor(x, method = 'pearson', use = 'pairwise.complete.obs')
colnames(x) <- 1:ncol(x)
rownames(x) <- 1:nrow(x)
#saveRDS(lung.hiclust.input, "~/Documents/CNA/Data/lung.hiclust.input.rds")
lung.hiclust.input <- readRDS("~/Documents/CNA/Data/lung.hiclust.input.rds")

# Melt into long form for geom_tile
meltedcormat.lung <- melt(lung.hiclust.input, na.rm = FALSE)
colnames(meltedcormat.lung) <- c('first_bin','second_bin','correlation')
meltedcormat.lung$chr1 <- as.numeric(car.info$start.stop[match(meltedcormat.lung$first_bin, car.info$start.stop$bin), 2])
meltedcormat.lung$chr2 <- as.factor(car.info$start.stop[match(meltedcormat.lung$second_bin, car.info$start.stop$bin), 2])
meltedcormat.lung$chr2 <- factor(meltedcormat.lung$chr2, levels = rev(c(seq(1,22,1))) )

# Prepare heatmaps 
hm <- ggplot(meltedcormat.lung, aes(first_bin, second_bin, fill = correlation)) +
      geom_tile() +
      scale_fill_distiller(palette ="Spectral", na.value="grey85", limits = c(-1,1), breaks=c(-1,0,1), guide=guide_colorbar()) +
      facet_grid(chr2 ~ chr1, scales='free', space = 'free') +
      scale_x_continuous(name='Chromosome', breaks = c(round(car.info$chr.mid, digits = 0)),  
                         labels = c(1:18,'\n19','20','\n21','22')) + 
      scale_y_continuous(name='Chromosome', breaks = c(round(car.info$chr.mid, digits = 0)),
                         labels = c(1:18,'19    ','20','21    ','22')) +
      guides(fill = guide_colourbar(title = 'Correlation', title.position = 'left', title.hjust = 0, title.vjust = 0.9,
                                    label.position = "bottom", label.vjust = 0)) +
      theme(
        legend.position = "none",
        plot.background=element_blank(),
        plot.margin=margin(t=0.5,r=0.5,b=0,l=0,"cm"),
        panel.spacing = unit(0.2, "lines"),
        panel.background = element_blank(),
        panel.border=element_blank(),
        strip.text = element_blank(),
        axis.text = element_text(size=24, colour='black'),
        axis.title = element_text(size=24, colour='black'),
        axis.ticks = element_line(size=0.4),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.ticks.length=unit(0.2, "cm"),) 

hm.legend <- as_ggplot(cowplot::get_legend(hm + 
                                             guides(fill = guide_colorbar(title="Correlation", label.position = "bottom",
                                                                          title.position = "left", title.vjust = 0.9)) +
                                             theme(plot.margin = unit(c(t=-100,r=0,b=-100,l=0), "cm"),
                                                   legend.position = "bottom",legend.direction="horizontal",
                                                   legend.margin = margin(grid::unit(c(t=-100,r=0,b=-100,l=0),"cm")),
                                                   legend.key.height = grid::unit(0.8,"cm"),
                                                   legend.key.width = grid::unit(1.4,"cm")) ))   


hmPlot <- cowplot::plot_grid(hm.legend, hm, ncol = 1, align = 'v', rel_heights = c(0.09, 1))

jpeg('tempfig.jpeg', width = 1000, height = 1000)
hmPlot
dev.off()

# split data -------------------------------------------------------------------
n <- round(tracerx.info$noPatients*0.6)
sample <- sample(tracerx.info$patientIDs, n)

tracerxBinned_train <- cbind(tracerxBinned[, c(1:3)], tracerxBinned[, sub("\\..*", "", colnames(tracerxBinned)) %in% sample])
tracerx_train.info <- PullDataInfo(rawdata = tracerxBinned_train)
tracerx_train.clonality <- PullDataClonality(rawdata = tracerxBinned_train, dataInfo = tracerx_train.info)
tracerx_train.diversity <- PullDataDiversity(rawdata = tracerxBinned_train, dataInfo = tracerx_train.info)

actualIth_tracerx_train <- data.frame(patient = lapply(data.frame(patient=tracerx_train.info$patientIDs), rep, tracerx_train.info$sampPerPatient),
                                sample = tracerx_train.info$sampleIDs,
                                actual = lapply(data.frame(actual=tracerx_train.diversity$pic.frac$pic.frac), rep, tracerx_train.info$sampPerPatient))

lung.matrices_train <- genMatrices(tracerxBinned_train)

# single bin models -------------------------------------------------------------
uniReg.in.list.lung <- lung.matrices_train[c(1,2,3)]

uniReg.in.list.lung <- lapply(uniReg.in.list.lung, function(x) {
  x <- t(x)
  x <- data.frame(apply(x, 2, as.character), row.names = rownames(x), check.names = FALSE)
  x <- data.frame(lapply(x, factor, levels=c(0,1), labels=c('diploid','CNA')), row.names = rownames(x), check.names = FALSE)
  x <- cbind(PIC.frac = actualIth_tracerx_train$actual, x)
  x
})

uniReg.out.list.lung <- list(data.frame(CNA='aneuploid', chr=tracerx_train.info$start.stop$chr, bin=tracerx_train.info$start.stop$bin, coeff=NA, pval=NA),
                             data.frame(CNA='gain', chr=tracerx_train.info$start.stop$chr, bin=tracerx_train.info$start.stop$bin, coeff=NA, pval=NA),
                             data.frame(CNA='loss', chr=tracerx_train.info$start.stop$chr, bin=tracerx_train.info$start.stop$bin, coeff=NA, pval=NA))
names(uniReg.out.list.lung) = c('diploid.aneu', 'diploid.gain', 'diploid.loss')

for ( i in 1:length(uniReg.in.list.lung) ) { 
  for ( k in 2:ncol(uniReg.in.list.lung[[i]]) ) {
    
    data <- uniReg.in.list.lung[[i]][,c(1,k)]
    data <- na.omit(data)
    colnames(data) <- c('ITH','predictor')
    
    data <- data %>%
      mutate(predictor = relevel(predictor, ref = 'diploid')) 
    
    if ( length(unique(as.character(data$predictor))) == 1 ) { 
      next
    }
    
    reg <- lm(ITH ~ ., data)
    
    coeffs <- data.frame(t(summary(reg)$coefficients[,1])) 
    names(coeffs) <- unlist(reg$xlevels)
    pvals <- data.frame(t(summary(reg)$coefficients[,4]))
    names(pvals) <- unlist(reg$xlevels)
    
    if ('CNA' %in% names(coeffs)) {
      uniReg.out.list.lung[[i]][(uniReg.out.list.lung[[i]]$bin==(k-1)),c(4:5)] <- c(coeffs$CNA,pvals$CNA)
    }
  }
}

# cluster -----------------------------------------------------------------------
load("~/Documents/CNA/Data/lung.hiclust.rda")
sig.hclust.lung <- pvpick(lung.hiclust, alpha=0.95, pv="au", max.only = T)

# For each cluster extract the corresponding bins
list <- list()
for ( i in 1:length(sig.hclust.lung$clusters) ) { 
  list[[i]] <- data.frame(bin = unlist(sig.hclust.lung$clusters[[i]]), cluster = i)
}

# Create a dataframe with columns: bin | cluster
sig.hclust.lung <- do.call('rbind', list)

# Convert dataframe to numeric
sig.hclust.lung <- data.frame(apply(sig.hclust.lung, 2, function(x) as.numeric(as.character(x))))

# Add chromosome data based on bin number
sig.hclust.lung$chr <- car.info$start.stop[match(sig.hclust.lung$bin, car.info$start.stop$bin), 2]


# Candidates and reps -------------------------------------------------------------
candidate.bins.lung <- do.call('rbind',uniReg.out.list.lung[c(2:3)])
candidate.bins.lung <- candidate.bins.lung[which(candidate.bins.lung$pval <= 0.05),]

# pull clonality
candidate.bins.lung$pcSubclonal <- NA
for ( i in 1:nrow(candidate.bins.lung) ) {
  bin <- candidate.bins.lung$bin[i]
  CNA <- candidate.bins.lung$CNA[i]
  
  if ( CNA=='gain') {
    candidate.bins.lung$pcSubclonal[i] <- tracerx_train.clonality$pcSubclonal$gain[which(tracerx_train.clonality$pcSubclonal$bin==bin)]
  }
  if ( CNA=='loss') {
    candidate.bins.lung$pcSubclonal[i] <- tracerx_train.clonality$pcSubclonal$loss[which(tracerx_train.clonality$pcSubclonal$bin==bin)]
  }
}

candidate.bins.lung$cluster <- sig.hclust.lung[match(candidate.bins.lung$bin, sig.hclust.lung$bin), 2]
rownames(candidate.bins.lung) <- 1:nrow(candidate.bins.lung)

# Assign unique clusters to bins with no cluster
new.clust <- max(candidate.bins.lung$cluster[is.finite(candidate.bins.lung$cluster)]) + 1
for (i in 1:nrow(candidate.bins.lung) ) {
  if ( is.na(candidate.bins.lung$cluster[i]) ) {
    candidate.bins.lung$cluster[i] <- new.clust
    new.clust <- new.clust + 1
  }
}

# Select one representive bin per cluster
list <- list()
for ( i in 1:max(candidate.bins.lung$cluster) ) {
  data <- candidate.bins.lung[which(candidate.bins.lung$cluster == i),]
  list[[i]] <- data[which.max(abs(data$coeff)),]
}
representitive.bins.lung <- do.call('rbind',list)

# Build full betareg -----------------------------------------------------------------
lung.in <- t(tracerxBinned_train[,-c(1:3)]) # The input matrix requires data on loss/gain/diploid, with a bin per column
lung.in <- data.frame(apply(lung.in, 2, as.character), check.names = FALSE)
lung.in <- data.frame(lapply(lung.in, factor, levels=c(2,1,3), labels=c('diploid','loss','gain')), check.names = FALSE)
rownames(lung.in) <- tracerx_train.info$sampleIDs
lung.in <- lung.in[,c(representitive.bins.lung$bin)] # Keep only the candidate bins (columns)

names(lung.in) = paste("bin_", names(lung.in), sep="")
lung.in = lung.in %>%
  mutate(across(everything(), as.character))

# Make dummy for beta
lung.in <- dummy_cols(lung.in, remove_selected_columns = TRUE, remove_first_dummy = TRUE)
lung.in$bin_ITH <- actualIth_tracerx_train$actual
remove <- attributes(alias(lm(bin_ITH ~ ., data = lung.in))$Complete)$dimnames[[1]]
lung.in <- lung.in[ ,colnames(lung.in) %!in% remove]

# Betareg
beta.lung <- betareg(bin_ITH ~., data = lung.in)
summary(beta.lung, type = "pearson")
x <- data.frame(vif(beta.lung))

# selection on 96 --------------------------------------------------------
x <- lung.in[-185]
y <- actualIth_tracerx_train$actual
selection_for.lung <- betaselect(x, y, criterion = "AIC",link = "logit", method = "forward", plotit = FALSE)
keep <- selection_for.lung$variable
lung.infor <- lung.in[ ,colnames(lung.in) %in% keep]
lung.infor$bin_ITH <- actualIth_tracerx_train$actual
beta_lung.for <- betareg(bin_ITH ~., data = lung.infor)
summary(beta_lung.for, type = "pearson")
x <- data.frame(vif(beta_lung.for))

x <- lung.in[-136]
x$pga <- tracerx_train.diversity$pga$pic.frac
y <- actualIth_tracerx_train$actual
selection_back.lungpga <- betaselect(x, y, criterion = "AIC",link = "logit", method = "backward", plotit = FALSE)
keep <- selection_back.lungpga$variable
lung.inback <- lung.in[ ,colnames(lung.in) %in% keep]
lung.inback$bin_ITH <- actualIth_tracerx_train$actual
beta_lung.back <- betareg(bin_ITH ~., data = lung.inback)
summary(beta_lung.back, type = "pearson")
y <- data.frame(vif(beta_lung.back))

x <- lung.in[-136]
y <- actualIth_tracerx_train$actual
selection_back.lung <- betaselect(x, y, criterion = "AIC",link = "logit", method = "backward", plotit = FALSE)
keep <- selection_back.lung$variable
lung.inback <- lung.in[ ,colnames(lung.in) %in% keep]
lung.inback$bin_ITH <- actualIth_tracerx_train$actual
beta_lung.back <- betareg(bin_ITH ~., data = lung.inback)
summary(beta_lung.back, type = "pearson")
y <- data.frame(vif=vif(beta_lung.back))

#remove <- c("bin_2331_loss","bin_142_gain","bin_2327_gain","bin_1366_gain","bin_2291_gain",
#            "bin_2393_loss","bin_1366_loss","bin_156_gain","bin_145_gain","bin_2335_loss",
#            "bin_1390_loss","bin_2392_gain","bin_2322_gain","bin_2391_loss","bin_2294_loss",
#            "bin_1440_gain","bin_2325_loss","bin_139_gain","bin_2402_gain","bin_1388_loss")
#lung.inbackrem <- lung.inback[ ,colnames(lung.inback) %!in% remove]
lung.inbackrem <- lung.inback
for ( i in 1:50 ) {
  max <- max(y$vif)
  remove <- rownames(y)[which(y$vif== max)]
  
  if (max <= 10) {
    break
  }
  
  lung.inbackrem <- lung.inbackrem[ ,colnames(lung.inbackrem) %!in% remove]
  beta_lung.backrem <- betareg(bin_ITH ~., data = lung.inbackrem)
  summary(beta_lung.backrem, type = "pearson")
  y <- data.frame(vif=vif(beta_lung.backrem))
}
round(AIC(beta_lung.backrem),2)

saveRDS(uniReg.out.list.lung, "~/Documents/CNA/Data/uniReg.out.list.lung.rds")
saveRDS(sig.hclust.lung, "~/Documents/CNA/Data/sig.hclust.lung.rds")
saveRDS(candidate.bins.lung, "~/Documents/CNA/Data/candidate.bins.lung.rds")
saveRDS(representitive.bins.lung, "~/Documents/CNA/Data/representitive.bins.lung.rds")
saveRDS(selection_for.lung, "~/Documents/CNA/Data/selection_for.lung.rds")
saveRDS(selection_back.lung, "~/Documents/CNA/Data/selection_back.lung.rds")
saveRDS(lung.inbackrem, "~/Documents/CNA/Data/lung.inbackrem.rds")
saveRDS(beta_lung.backrem, "~/Documents/CNA/Data/beta_lung.backrem.rds")
saveRDS(hg19predictors.lung, "~/Documents/CNA/Data/hg19predictors.lung.rds")


# assess on tracerX------------------------------------------------------------------------
# pull predictors
beta <- beta_lung.backrem
multiReg.in <- lung.inbackrem

x <- rownames(data.frame(summary(beta)$coefficients))[-1]
hg19predictors.lung <- data.frame(bin=as.numeric(str_extract_all(x, "[0-9]+")), cna=str_sub(x,-4,-1))
hg19predictors.lung <- merge(unique(hg19predictors.lung), car.info$start.stop)
hg19predictors.lung$cluster <- candidate.bins.lung[match(hg19predictors.lung$bin, candidate.bins.lung$bin), 7] 
hg19predictors.lung$pcSubclonal <- candidate.bins.lung[match(hg19predictors.lung$bin, candidate.bins.lung$bin), 6] 
both <- hg19predictors.lung$bin[which(duplicated(hg19predictors.lung$bin)==T)]
hg19predictors.lung$cna[which(hg19predictors.lung$bin %in% hg19predictors.lung$bin[which(duplicated(hg19predictors.lung$bin)==T)] )] <- 'both'
hg19predictors.lung <- hg19predictors.lung[!duplicated(hg19predictors.lung),]

psuedoR2 <- round(beta$pseudo.r.squared,3)
loglik <- round(beta$loglik,2)
aic <- round(AIC(beta),2)
predictedIth_tracerx_train <- data.frame(predicted = predict(beta_lung.backrem, lung.in))

# Enter actual and predicted ITH into a comparative dataframe
predictedIth_tracerx_train <- data.frame(cbind(actualIth_tracerx_train, predictedIth_tracerx_train))
predictedIth_tracerx_train$type <- 'tracerx'

# Add two labels for plot
plot <- ggplot(data = predictedIth_tracerx_train, aes(x = actual, y = predicted)) +
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

png('tempfig.jpeg', width = 1000, height = 1000)
plot
dev.off()

# bootstrapon random 27 bins----------------------------------------------------------------------
# Are these 33 bins actually special? try random 33 from candidate and assess model

test33 <- t(tracerxBinned_train[,-c(1:3)])
test33 <- data.frame(apply(test33, 2, as.character), check.names = FALSE)
test33 <- data.frame(lapply(test33, factor, levels=c(2,1,3), labels=c('diploid','loss','gain')), check.names = FALSE)
rownames(test33) <- tracerx_train.info$sampleIDs

test33 <- cbind(ITH=actualIth_tracerx_train$actual, test33)
names(test33) = paste("bin_", names(test33), sep="")
test33 = test33 %>%
  mutate(across(everything(), as.character))

test33 <- dummy_cols(test33[-1], remove_selected_columns = TRUE, remove_first_dummy = TRUE) 
test33$bin_ITH <- actualIth_tracerx_train$actual

bins <- unique(as.numeric(gsub("[^\\d]+", "", colnames(test33), perl=TRUE))) #pull bins to select from
bins <- bins[!is.na(bins)]

set.seed(12345)
pR2 <- list()
ll <- list()
tb <- list()
aic2 <- list()
foreach(i = 1:100) %do% {
  testbins <- sample(bins, 53)
  tb[[i]] <- testbins
  input <- test33[ , which(as.numeric(gsub("[^\\d]+", "", colnames(test33), perl=TRUE)) %in% testbins)]
  input$bin_ITH <- test33$bin_ITH
  try({
    beta <- betareg(bin_ITH ~., data = input)
    pR2[[i]] <- round(beta$pseudo.r.squared,3)
    ll[[i]] <- round(beta$loglik,2)
    aic2[[i]] <- round(AIC(beta),2)
  }, silent = TRUE)
}

length(unlist(aic2))
z <- unlist(aic2)
mean(z)
min(z)

min <- which(z==min(z))
tb[[min]][order(tb[[min]])]
z[order(z)]

jpeg('tempfig.jpeg', width = 1000, height = 1000)
hist(z, xlim = c(-500,-150))
abline(v = aic, col="blue")
dev.off()


# LOOCV ---------------------------------------------------------------
loo <- t(tracerxBinned[,-c(1:3)])
loo <- data.frame(apply(loo, 2, as.character), check.names = FALSE)
loo <- data.frame(lapply(loo, factor, levels=c(2,1,3), labels=c('diploid','loss','gain')), check.names = FALSE)
rownames(loo) <- tracerx.info$sampleIDs
loo <- loo[, which(colnames(loo) %in% hg19predictors.lung$bin)]

loo <- cbind(actualIth_tracerx$actual, loo)
names(loo) = paste("bin_", names(loo), sep="")
loo = loo %>%
  mutate(across(everything(), as.character))

loo <- dummy_cols(loo[-1], remove_selected_columns = TRUE, remove_first_dummy = TRUE) 
loo$bin_ITH <- actualIth_tracerx$actual
rownames(loo) <- tracerx.info$sampleIDs

set.seed(123456)
pR2 <- list()
RMSE <- list()
MAE <- list()
foreach(i = 1:295) %do% {
  input <- sample_n(loo, 294)
  data <- data.frame(predicted=predict(beta_lung.backrem, input), actual=input$bin_ITH)
  data <- na.omit(data)
  
  r <- lm(actual ~ predicted, data=data)
  pR2[[i]] <- round(summary(r)$r.squared ,3)
  RMSE[[i]] <- round(sqrt(mean((data$actual - data$predicted)^2)), 2)
  MAE[[i]] <- round(mae(data$actual, predict(r)) ,2)
}

mean(unlist(pR2))
mean(unlist(RMSE))
mean(unlist(MAE))

# CHROMOMAP just cna-------------------------------------------------------------
hg19 <- makeHg19()
hg19$start <- 1
hg19 <- hg19[c(1,6,2,3)]
hg19$chrom <- as.character(as.numeric(hg19$chrom))
hg19 <- hg19[which(hg19$chrom %in% hg19predictors.lung$chr),]

# make annotation file
annotation <- hg19predictors.lung[,c(1,3,4,5,2)]
annotation$chr <- as.character((annotation$chr))
annotation$cna <- as.character((annotation$cna))
chromoMap(list(hg19), list(annotation),
          #cna
          data_based_color_map = T, 
          data_type = "categorical",
          data_colors = list(c("#CC0033","#330099","#009999")),
          #labels
          labels = T, label_font = 9, label_angle = -70,
          #features
          canvas_width = 1200, canvas_height = 1000,
          chr_color = c("#CCCCCC"),
          chr_width = 12, chr_length = 9, ch_gap = -15, 
          y_chr_scale = -30, top_margin = -20,
          legend = T, lg_x = 10, lg_y = 100)
