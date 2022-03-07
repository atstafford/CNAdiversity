# PREDICT TCGA DIVERSITY

# Prepare input for prediction

# Transpose as input table requires cols=bins
TCGA.input <- t(COADbinned[,-c(1:3)])

# Convert matrix data to character then factor
TCGA.input <- data.frame(apply(TCGA.input, 2, as.character), check.names = FALSE)
TCGA.input <- data.frame(lapply(TCGA.input, factor, levels=c(2,1,3), labels=c('diploid','loss','gain')), check.names = FALSE)
rownames(TCGA.input) <- TCGA.info$patientIDs

names(TCGA.input) = paste("bin_", names(TCGA.input), sep="")
TCGA.input <- dummy_cols(TCGA.input, remove_selected_columns = TRUE, remove_first_dummy = TRUE)
rownames(TCGA.input) <- TCGA.info$sampleIDs

# Predict using step.model
TCGA.predict <- data.frame(TCGA_barcode = substr(rownames(TCGA.input),1,12), predictedITH = predict(beta_52backrem, TCGA.input))

# Add predicted ITH back into COAD_clinical
COAD_clinical$ITH <- TCGA.predict[match(COAD_clinical$TCGA_barcode, TCGA.predict$TCGA_barcode), 2]

# Keep only those samples whose ITH could be predicted
COAD_clinical <- COAD_clinical[!is.na(COAD_clinical$ITH),]

# Mark those patients that didnt have a CNA in any of the predictive bins
COAD_clinical$intercept <- "no"
COAD_clinical$intercept[which(COAD_clinical$ITH==coef(beta_52backrem)["(Intercept)"])] <- "yes"


# SURVIVAL

# ASSESS CORRELATION BETWEEN PREDICTED ITH AND SURVIVAL IN TCGA

# Remove samples missing stage data
COAD_survival <- COAD_clinical[which(COAD_clinical$stage!="'--"),]

# Keep only MSS
COAD_survival <- COAD_survival[which(COAD_survival$MSI=="MSS"),]

# Split patients into 2 groups based on ITH
COAD_survival$divGroup2 <- 'low' 
COAD_survival$divGroup2[COAD_survival$ITH>= median(COAD_survival$ITH) ] <- "high"

# Split patients into 4 groups based on ITH
COAD_survival$divGroup4 <- ntile(COAD_survival$ITH, 4)
#COAD_survival$divGroup4[which(COAD_survival$divGroup4!=1)] <- 4

# Perform chi squared test
chisq2 <- chisq.test(table(COAD_survival$divGroup2, COAD_survival$stage))
chisq4 <- chisq.test(table(COAD_survival$divGroup4, COAD_survival$stage))
corrplot(chisq2$residuals, is.cor = FALSE)
corrplot(chisq4$residuals, is.cor = FALSE)
chisq2$p.value
chisq4$p.value

# Pull diversity by stage into summary table
COAD_survival_summary <- as.data.frame(table(COAD_survival$stage, COAD_survival$divGroup2))

# Plot
divByStagePlot <- ggplot(COAD_survival_summary, aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_bar(position="fill", stat="identity", width=0.8) +
  #scale_fill_manual(values = c('#00CC66','#FF9933'), labels = c('High','Low')) + 
  ylab('Fraction') +
  labs(fill = 'Predicted CNA diversity') +
  guides(fill = guide_legend(direction = 'horizontal', label.position = "left", label.hjust = 1)) +
  scale_x_discrete(labels=c('Stage I\n(n=63)','Stage II\n(n=138)','Stage III\n(n=96)','Stage IV\n(n=45)')) +
  theme(plot.background=element_blank(),
        plot.margin=margin(t=0.4,r=0,b=1,l=0,"cm"),
        panel.background = element_blank(),
        legend.position = "top",
        legend.text = element_text(size=24, colour='black'),
        legend.title = element_text(size=24, colour='black',face='bold'),
        plot.title = element_blank(),
        axis.title.y = element_text(size=24, colour='black'),
        axis.title.x = element_blank(),
        axis.text = element_text(size=24, colour='black'),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.ticks.length=unit(0.2, "cm"),)

jpeg('tempfig.jpeg', width = 1000, height = 1000)
divByStagePlot
dev.off()

# Survival of high/low diversity patients (KM) 

# split by stage
COAD_survival$divSGroup2 <- "low"
for (i in 1:length(unique(COAD_survival$stage))) {
  stage <- unique(COAD_survival$stage)[i]
  median <- median(COAD_survival$ITH[which(COAD_survival$stage==stage)])
  patients <- COAD_survival$TCGA_barcode[which(COAD_survival$stage==stage & COAD_survival$ITH>=median)]
  COAD_survival$divSGroup2[which(COAD_survival$TCGA_barcode %in% patients)] <- "high"
}

# Calculate serial time between death and censor
#add censore
COAD_survival$OS_cens <- 0
COAD_survival$OS_cens[which(COAD_survival$vital_status=="Dead")] <- 1

#manual edits for those who dies same year as diag
COAD_survival$days_to_death <- as.numeric(as.character(COAD_survival$days_to_death))
COAD_survival$days_to_death[which(COAD_survival$TCGA_barcode=="TCGA-AA-3852")] <- 182
COAD_survival$days_to_death[which(COAD_survival$TCGA_barcode=="TCGA-AA-3845")] <- 182
COAD_survival$days_to_death[which(COAD_survival$TCGA_barcode=="TCGA-AA-3850")] <- 182
COAD_survival <- COAD_survival[-which(COAD_survival$days_to_death==1 | COAD_survival$days_to_death==0),]

COAD_survival$days_to_birth <- as.numeric(as.character(COAD_survival$days_to_birth))
COAD_survival <- COAD_survival[!is.na(COAD_survival$days_to_birth),]
COAD_survival$OS_time <- NA

for (i in 1:nrow(COAD_survival)) {
  if (COAD_survival$OS_cens[i] == 0) {
    COAD_survival$OS_time[i] <- -1 * COAD_survival$days_to_birth[i]
  }
  else if (COAD_survival$OS_cens[i] == 1) {
    COAD_survival$OS_time[i] <- COAD_survival$days_to_death[i]
  }
}

# Fit survival data using the Kaplan-Meier method
x <- data.frame(table(COAD_survival$stage, COAD_survival$divSGroup2))
stage <- as.character(unique(x$Var1))
n <- list()
for (i in 1:length(stage)) {
  n[[i]] <- paste(stage[i],
                  paste("(n=", sum(x$Freq[which(x$Var1==stage[i])]), ")",sep = ""),
                  paste("\nhigh=", x$Freq[which(x$Var1==stage[i] & x$Var2=="high")]," low=", x$Freq[which(x$Var1==stage[i] & x$Var2=="low")], sep = ""))
}

stage = c("Stage I", "Stage II", "Stage III", "Stage IV")
plot.list <- list()
for (i in 1:length(stage)) {
  surv_object <- Surv(time = COAD_survival$OS_time[which(COAD_survival$stage==stage[i])], 
                      event = COAD_survival$OS_cens[which(COAD_survival$stage==stage[i])])
  fit <- survfit(surv_object ~ divSGroup2, data = COAD_survival[which(COAD_survival$stage==stage[i]),])
  plot.list[[i]] <- ggsurvplot(fit, data = COAD_survival[which(COAD_survival$stage==stage[i]),], pval = TRUE, 
                               title=n[i], legend="none",
                               font.main=24, font.x=NA, font.y=NA, font.tickslab=24, font.legend=24, pval.size=9, 
                               #palette = c('#00CC66','#FF9933'), legend.title='Predicted CNA diversity', legend.labs=c('high','low')
                               )
}

# Plot
plot <- plot_grid(plot.list[[1]]$plot, plot.list[[2]]$plot, plot.list[[3]]$plot, plot.list[[4]]$plot, nrow = 2)

plot <- plot_grid(get_legend(divByStagePlot), plot, nrow = 2, rel_heights = c(1,10))

plot <- annotate_figure(plot, 
                        left = text_grob('Survival probability', rot = 90, size = 24),
                        bottom = text_grob('Time (years)', size = 24))
plot

jpeg('tempfig.jpeg', width = 1000, height = 1000)
plot
dev.off()

