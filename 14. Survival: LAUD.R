# PREDICT TCGA DIVERSITY

# Prepare input for prediction

# Transpose as input table requires cols=bins
names <- colnames(LAUDbinned)[-c(1:3)]
TCGA.input <- t(LAUDbinned[,-c(1:3)])

# Convert matrix data to character then factor
TCGA.input <- data.frame(apply(TCGA.input, 2, as.character), check.names = FALSE)
TCGA.input <- data.frame(lapply(TCGA.input, factor, levels=c(2,1,3), labels=c('diploid','loss','gain')), check.names = FALSE)
rownames(TCGA.input) <- names
x <- c(hg19predictors$bin,2370,2375)
TCGA.input <- TCGA.input[,x]

#bin 2371 missing, so use 2370
TCGA.input$`2371` <- TCGA.input$`2370`
TCGA.input <- TCGA.input[,-which(colnames(TCGA.input)=="2370")]

#bin 2374 missing, so use 2375
TCGA.input$`2374` <- TCGA.input$`2375`
TCGA.input <- TCGA.input[,-which(colnames(TCGA.input)=="2375")]

# make dummy
names(TCGA.input) = paste("bin_", names(TCGA.input), sep="")
TCGA.input <- dummy_cols(TCGA.input, remove_selected_columns = TRUE, remove_first_dummy = TRUE)
rownames(TCGA.input) <- names

# Predict using step.model
TCGA.predict <- data.frame(TCGA_barcode = substr(rownames(TCGA.input),1,12), predictedITH = predict(beta_52backrem, TCGA.input))

# Add predicted ITH back into LAUD_clinical
LAUD_clin$ITH <- TCGA.predict[match(LAUD_clin$TCGA_barcode, TCGA.predict$TCGA_barcode), 2]

# Keep only those samples whose ITH could be predicted
LAUD_clin <- LAUD_clin[!is.na(LAUD_clin$ITH),]

# Mark those patients that didnt have a CNA in any of the predictive bins
LAUD_clin$intercept <- "no"
LAUD_clin$intercept[which(LAUD_clin$ITH==coef(beta_52backrem)["(Intercept)"])] <- "yes"

# SURVIVAL

# ASSESS CORRELATION BETWEEN PREDICTED ITH AND SURVIVAL IN TCGA

# Remove samples missing stage data
LAUD_survival <- LAUD_clin[which(LAUD_clin$stage!="'--"),]
#LAUD_survival <- LAUD_survival[-c(1,6,7,8,12,14,16:23,25:28)]

# Split patients into 2 groups based on ITH
LAUD_survival$divGroup2 <- 'low' 
LAUD_survival$divGroup2[LAUD_survival$ITH>= median(LAUD_survival$ITH) ] <- "high"

# Split patients into 4 groups based on ITH
LAUD_survival$divGroup4 <- ntile(LAUD_survival$ITH, 4)
#LAUD_survival$divGroup4[which(LAUD_survival$divGroup4!=1)] <- 4

# Perform chi squared test
chisq2 <- chisq.test(table(LAUD_survival$divGroup2, LAUD_survival$stage))
chisq4 <- chisq.test(table(LAUD_survival$divGroup4, LAUD_survival$stage))
corrplot(chisq2$residuals, is.cor = FALSE)
corrplot(chisq4$residuals, is.cor = FALSE)
chisq2$p.value
chisq4$p.value

# Pull diversity by stage into summary table
LAUD_survival_summary <- as.data.frame(table(LAUD_survival$stage, LAUD_survival$divGroup2))

# Plot
divByStagePlot <- ggplot(LAUD_survival_summary, aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_bar(position="fill", stat="identity", width=0.8) +
  #scale_fill_manual(values = c('#00CC66','#FF9933'), labels = c('High','Low')) + 
  ylab('Fraction') +
  labs(fill = 'Predicted CNA diversity') +
  guides(fill = guide_legend(direction = 'horizontal', label.position = "left", label.hjust = 1)) +
  scale_x_discrete(labels=c('Stage I\n(n=152)','Stage II\n(n=146)','Stage III\n(n=103)','Stage IV\n(n=12)')) +
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
LAUD_survival$divSGroup2 <- "low"
for (i in 1:length(unique(LAUD_survival$stage))) {
  stage <- unique(LAUD_survival$stage)[i]
  median <- median(LAUD_survival$ITH[which(LAUD_survival$stage==stage)])
  patients <- LAUD_survival$TCGA_barcode[which(LAUD_survival$stage==stage & LAUD_survival$ITH>=median)]
  LAUD_survival$divSGroup2[which(LAUD_survival$TCGA_barcode %in% patients)] <- "high"
}


# Calculate serial time between death and censor

# Convert year of birth/death to numeric
#LAUD_survival <- LAUD_survival[c(1, 3,6,2, 9,10, 4,7,5, 8,11:14)]

#add censore
LAUD_survival$OS_cens <- 0
LAUD_survival$OS_cens[which(LAUD_survival$vital_status=="Dead")] <- 1

#manual edits
LAUD_survival$days_to_death <- as.numeric(as.character(LAUD_survival$days_to_death))
LAUD_survival <- LAUD_survival[-which(LAUD_survival$days_to_death==1 | LAUD_survival$days_to_death==0),]

LAUD_survival$days_to_birth <- as.numeric(as.character(LAUD_survival$days_to_birth))
LAUD_survival <- LAUD_survival[!is.na(LAUD_survival$days_to_birth),]

for (i in 1:nrow(LAUD_survival)) {
  if (LAUD_survival$OS_cens[i] == 0) {
    LAUD_survival$OS_time[i] <- -1 * LAUD_survival$days_to_birth[i]
  }
  else if (LAUD_survival$OS_cens[i] == 1) {
    LAUD_survival$OS_time[i] <- LAUD_survival$days_to_death[i]
  }
}

# Fit survival data using the Kaplan-Meier method
x <- data.frame(table(LAUD_survival$stage, LAUD_survival$divGroup2))
stage <- as.character(unique(x$Var1))
n <- list()
for (i in 1:length(stage)) {
  n[[i]] <- paste(stage[i],
             paste("(n=", sum(x$Freq[which(x$Var1==stage[i])]), ")",sep = ""),
             paste("\nhigh=", x$Freq[which(x$Var1==stage[i] & x$Var2=="high")]," low=", x$Freq[which(x$Var1==stage[i] & x$Var2=="low")], sep = ""))
}

plot.list <- list()
for (i in 1:length(stage)) {
  surv_object <- Surv(time = LAUD_survival$OS_time[which(LAUD_survival$stage==stage[i])], 
                      event = LAUD_survival$OS_cens[which(LAUD_survival$stage==stage[i])])
  fit <- survfit(surv_object ~ divGroup2, data = LAUD_survival[which(LAUD_survival$stage==stage[i]),])
  plot.list[[i]] <- ggsurvplot(fit, data = LAUD_survival[which(LAUD_survival$stage==stage[i]),], pval = TRUE, 
                               title=n[i], legend="none",
                               font.main=24, font.x=NA, font.y=NA, font.tickslab=24, font.legend=24, pval.size=9, 
                               legend.title='Predicted CNA diversity', 
                               #legend.labs=c('high','low'), palette = c('#00CC66','#FF9933')
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

