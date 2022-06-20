# Predict --------------------------------------------------------------

# Prepare input for prediction

# Transpose as input table requires cols=bins
TCGA.input <- t(READbinned[,-c(1:3)])

# Convert matrix data to character then factor
TCGA.input <- data.frame(apply(TCGA.input, 2, as.character), check.names = FALSE)
TCGA.input <- data.frame(lapply(TCGA.input, factor, levels=c(2,1,3), labels=c('diploid','loss','gain')), check.names = FALSE)
rownames(TCGA.input) <- READ.info$patientIDs

names(TCGA.input) = paste("bin_", names(TCGA.input), sep="")
TCGA.input <- dummy_cols(TCGA.input, remove_selected_columns = TRUE, remove_first_dummy = TRUE)
rownames(TCGA.input) <- READ.info$sampleIDs

# Predict using step.model
TCGA.predict <- data.frame(TCGA_barcode = substr(rownames(TCGA.input),1,12), predictedITH = predict(beta_52backrem, TCGA.input))

# Add predicted ITH back into COAD_clinical
READ_clinical$ITH <- TCGA.predict[match(READ_clinical$TCGA_barcode, TCGA.predict$TCGA_barcode), 2]

# Keep only those samples whose ITH could be predicted
READ_clinical <- READ_clinical[!is.na(READ_clinical$ITH),]

# Mark those patients that didnt have a CNA in any of the predictive bins
READ_clinical$intercept <- "no"
READ_clinical$intercept[which(READ_clinical$ITH==coef(beta_52backrem)["(Intercept)"])] <- "yes"

# Remove samples missing stage data
READ_survival <- READ_clinical[which(READ_clinical$stage!="'--"),]

# Keep only MSS
#READ_survival <- READ_clinical[which(READ_clinical$MSI=="MSS"),]

# ITH per stage ----------------------------------------------------------------
# Perform chi squared test
chisq2 <- chisq.test(table(READ_survival$divGroup2, COAD_survival$stage))
chisq4 <- chisq.test(table(READ_survival$divGroup4, COAD_survival$stage))
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

# Set up KM --------------------------------------------------------------

# Calculate serial time between death and censor
#add censore
READ_survival$OS_cens <- 0
READ_survival$OS_cens[which(READ_survival$vital_status=="Dead")] <- 1

#manual edits for those who dies same year as diag
READ_survival$days_to_death <- as.numeric(as.character(READ_survival$days_to_death))
#READ_survival$days_to_death[which(READ_survival$TCGA_barcode=="TCGA-AA-3852")] <- 182
#READ_survival$days_to_death[which(READ_survival$TCGA_barcode=="TCGA-AA-3845")] <- 182
#READ_survival$days_to_death[which(READ_survival$TCGA_barcode=="TCGA-AA-3850")] <- 182
READ_survival <- READ_survival[-which(READ_survival$days_to_death==1 | READ_survival$days_to_death==0),]

READ_survival$days_to_birth <- as.numeric(as.character(READ_survival$days_to_birth))
#COAD_survival <- COAD_survival[!is.na(COAD_survival$days_to_birth),]
READ_survival$OS_time <- NA

for (i in 1:nrow(READ_survival)) {
  if (READ_survival$OS_cens[i] == 0) {
    READ_survival$OS_time[i] <- READ_survival$age_at_index - (READ_survival$year_of_diagnosis - READ_survival$year_of_birth)
  }
  else if (READ_survival$OS_cens[i] == 1) {
    READ_survival$OS_time[i] <- READ_survival$days_to_death[i]/365
  }
}


# Split into groups --------------------------------------------------------------

# Split patients into 2 groups based on ITH
READ_survival$divGroup4 <- ntile(READ_survival$ITH, 4)
READ_survival$divGroup3 <- "med"
READ_survival$divGroup3[which(READ_survival$divGroup4==4)] <- "top25"
READ_survival$divGroup3[which(READ_survival$divGroup4==1)] <- "bottom25"

cut <- surv_cutpoint(READ_survival, time = "OS_time", event = "OS_cens", variables = "ITH")
res.cat <- surv_categorize(cut)
READ_survival$optim <- res.cat$ITH

length(READ_survival$ITH[which(READ_survival$optim=="low")])/nrow(READ_survival)
length(READ_survival$ITH[which(READ_survival$optim=="high")])/nrow(READ_survival)

# optim  ---------------------------------
chooseCut <- COAD_survival$optim

fit <- survfit(Surv(OS_time, OS_cens) ~ chooseCut, data = COAD_survival)
plot <- ggsurvplot(fit, data = COAD_survival, pval = TRUE, legend="none",
                   font.main=24, font.x=NA, font.y=NA, font.tickslab=24, font.legend=24, pval.size=9)

hist <- ggplot(COAD_survival, aes(ITH, fill = chooseCut)) + 
  geom_histogram(binwidth = 0.001) +
  #scale_fill_manual(values = c("red", "blue"), name = "ITH", breaks = c("high","low")) +
  theme(panel.background = element_blank(), legend.position = "top",axis.text = element_text(size=24, colour='black'),
        axis.line = element_line(colour = "black", size = 0.5),axis.title = element_text(size=24, colour='black'))

plot_grid(hist , plot$plot)

x <- data.frame(table(COAD_survival$stage, chooseCut))
stage <- as.character(unique(x$Var1))
n <- list()
for (i in 1:length(stage)) {
  n[[i]] <- paste(stage[i],
                  paste("(n=", sum(x$Freq[which(x$Var1==stage[i])]), ")",sep = ""),
                  paste("\nhigh=", x$Freq[which(x$Var1==stage[i] & x$chooseCut=="high")]," low=", x$Freq[which(x$Var1==stage[i] & x$chooseCut=="low")], sep = ""))
}

plot.list <- list()
for (i in 1:length(stage)) {
  surv_object <- Surv(time = COAD_survival$OS_time[which(COAD_survival$stage==stage[i])], 
                      event = COAD_survival$OS_cens[which(COAD_survival$stage==stage[i])])
  fit <- survfit(surv_object ~ optim, data = COAD_survival[which(COAD_survival$stage==stage[i]),])
  plot.list[[i]] <- ggsurvplot(fit, data = COAD_survival[which(COAD_survival$stage==stage[i]),], pval = TRUE, 
                               title=n[i], legend="none",
                               font.main=24, font.x=NA, font.y=NA, font.tickslab=24, font.legend=24, pval.size=9)
}
legend <- get_legend(ggsurvplot(fit, data = COAD_survival[which(COAD_survival$stage==stage[i]),])$plot)
plot <- plot_grid(plot.list[[1]]$plot, plot.list[[2]]$plot, plot.list[[3]]$plot, plot.list[[4]]$plot, nrow = 2)
plot <- plot_grid(legend, plot, nrow = 2, rel_heights = c(1,10))
plot <- annotate_figure(plot, 
                        left = text_grob('Survival probability', rot = 90, size = 24),
                        bottom = text_grob('Time (years)', size = 24))
plot
# top and bottom 25% -----
z <- READ_survival
z$COADcut <- "med"
z$COADcut[which(z$ITH<0.086)] <-  "bot25"
z$COADcut[which(z$ITH>0.317)] <- "top25"
chooseCut <- z$COADcut

legend <- as_ggplot(cowplot::get_legend(ggplot(data.frame(color=c("Top 25%","Bottom 25%"),x=c(1,2), y=c(1,2))) +
                                          geom_bar(aes(x=x, y=y, fill=(color)), stat="identity") + 
                                          scale_fill_manual(values = c("#330099","#CC0033")) +
                                          guides(color = guide_legend(direction = 'horizontal', label.position = "left", label.hjust = 1)) +
                                          theme(plot.margin = unit(c(t=-100,r=0,b=-100,l=0), "cm"),
                                                legend.position = "bottom",legend.direction="horizontal",
                                                legend.title = element_blank(),
                                                legend.margin = margin(grid::unit(c(t=-100,r=0,b=-100,l=0),"cm")),
                                                legend.text = element_text(size=24, colour='black'),
                                                legend.key.height = grid::unit(0.8,"cm"),
                                                legend.key.width = grid::unit(1.4,"cm"))  ))

hist <- ggplot(z, aes(round(ITH,2), fill = COADcut)) + 
  geom_histogram(binwidth = 0.005) +
  xlab("Predicted CNA diversity") +
  scale_fill_manual(values = c("#330099","#CCCCCC","#CC0033")) +
  #ggtitle(paste("","\n"))+
  theme(panel.background = element_blank(), legend.position = "none",axis.text = element_text(size=24, colour='black'),
        axis.line = element_line(colour = "black", size = 0.5),axis.title = element_text(size=24, colour='black'))
hist <- plot_grid(legend, hist, ncol=1, rel_heights = c(1,5))

z <-z[which(z$COADcut!="med"),]

title <- paste(paste("Pan-stage ","(n=", nrow(z),")" ,sep = ""),
               paste("\nhigh=", nrow(z[which(z$divGroup3=="top25"),])," low=", nrow(z[which(z$divGroup3=="bottom25"),]), sep = ""))

fit <- survfit(Surv(OS_time, OS_cens) ~ COADcut, data = z)
panOS <- ggsurvplot(fit, data = z, pval = TRUE, legend="none", palette = c("#330099","#CC0033"),
                    font.main=24, font.x=24, font.y=24, font.tickslab=24, font.legend=24, pval.size=9, title=title,
                    ylab="Survival Probability", xlab="Time (years)")

plot1 <- plot_grid(hist, panOS$plot, nrow = 1)


x <- data.frame(table(z$stage, z$COADcut))
stage <- as.character(unique(x$Var1))
n <- list()
for (i in 1:length(stage)) {
  n[[i]] <- paste(stage[i],
                  paste("(n=", sum(x$Freq[which(x$Var1==stage[i])]), ")",sep = ""),
                  paste("\nhigh=", x$Freq[which(x$Var1==stage[i] & x$chooseCut=="top25")]," low=", x$Freq[which(x$Var1==stage[i] & x$chooseCut=="bottom25")], sep = ""))
}
plot.list <- list()
for (i in 1:length(stage)) {
  surv_object <- Surv(time = z$OS_time[which(z$stage==stage[i])], 
                      event = z$OS_cens[which(z$stage==stage[i])])
  fit <- survfit(surv_object ~ COADcut, data = z[which(z$stage==stage[i]),])
  plot.list[[i]] <- ggsurvplot(fit, data = z[which(z$stage==stage[i]),], pval = TRUE, 
                               title=n[i], legend="none", palette = c("#330099","#CC0033"),
                               font.main=24, font.x=NA, font.y=NA, font.tickslab=24, font.legend=24, pval.size=9)
}

plot2 <- plot_grid(plot.list[[1]]$plot, plot.list[[2]]$plot, plot.list[[3]]$plot, plot.list[[4]]$plot, nrow = 2)
plot2 <- annotate_figure(plot2, 
                         left = text_grob('Survival probability', rot = 90, size = 24),
                         bottom = text_grob('Time (years)', size = 24))

plot <- plot_grid(plot1, plot2, ncol = 1, rel_widths = c(1,1), rel_heights = c(1,2))

jpeg('tempfig.jpeg', width = 1000, height = 1300)
plot
dev.off()

# Progression -------------------------
TCGA_clinical <- read_excel("~/Documents/CNA/Data/TCGA/survival.dijk.xlsx", sheet = 66, skip = 1)
z$PFI <- TCGA_clinical[match(z$TCGA_barcode, str_sub(TCGA_clinical$Samplename,1,str_length(TCGA_clinical$Samplename)-16)),c(20)]
z$PFI_event <- TCGA_clinical[match(z$TCGA_barcode, str_sub(TCGA_clinical$Samplename,1,str_length(TCGA_clinical$Samplename)-16)),c(21)]

chooseCut <- z$divGroup3

x <- data.frame(table(z$stage, chooseCut))
stage <- as.character(unique(x$Var1))
n <- list()
for (i in 1:length(stage)) {
  n[[i]] <- paste(stage[i],
                  paste("(n=", sum(x$Freq[which(x$Var1==stage[i])]), ")",sep = ""),
                  paste("\nhigh=", x$Freq[which(x$Var1==stage[i] & x$chooseCut=="top25")]," low=", x$Freq[which(x$Var1==stage[i] & x$chooseCut=="bottom25")], sep = ""))
}
plot.list <- list()
for (i in 1:length(stage)) {
  surv_object <- Surv(time = z$prog$PFI[which(z$stage==stage[i])], 
                      event = z$prog$PFI_event[which(z$stage==stage[i])])
  fit <- survfit(surv_object ~ divGroup3, data = z[which(z$stage==stage[i]),])
  plot.list[[i]] <- ggsurvplot(fit, data = z[which(z$stage==stage[i]),], pval = TRUE, 
                               title=n[i], legend="none",
                               font.main=24, font.x=NA, font.y=NA, font.tickslab=24, font.legend=24, pval.size=9)
}
legend <- get_legend(ggsurvplot(fit, data = z[which(z$stage==stage[i]),])$plot)
plot <- plot_grid(plot.list[[1]]$plot, plot.list[[2]]$plot, plot.list[[3]]$plot, plot.list[[4]]$plot, nrow = 2)
plot <- plot_grid(legend, plot, nrow = 2, rel_heights = c(1,10))
plot <- annotate_figure(plot, 
                        left = text_grob('Survival probability', rot = 90, size = 24),
                        bottom = text_grob('Time (years)', size = 24))
plot

fit <- survfit(Surv(time = z$prog$PFI, event = z$prog$PFI_event) ~ chooseCut, data = z)
plot <- ggsurvplot(fit, data = z, pval = TRUE, legend="none",
                   font.main=24, font.x=NA, font.y=NA, font.tickslab=24, font.legend=24, pval.size=9)
hist <- ggplot(z, aes(ITH, fill = chooseCut)) + 
  geom_histogram(binwidth = 0.001) +
  #scale_fill_manual(values = c("red", "blue"), name = "ITH", breaks = c("high","low")) +
  theme(panel.background = element_blank(), legend.position = "top",axis.text = element_text(size=24, colour='black'),
        axis.line = element_line(colour = "black", size = 0.5),axis.title = element_text(size=24, colour='black'))

plot <- plot_grid(hist , plot$plot)

jpeg('tempfig.jpeg', width = 1000, height = 1000)
plot
dev.off()


# location in bowel-------------------------
y <- COAD_survival[c(2,23,29)]
y$tissue_or_organ_of_origin <- as.factor(y$tissue_or_organ_of_origin)
levels(y$tissue_or_organ_of_origin) <- c("Cecum","Ascending colon","Hepatic flexure of colon",
                                         "Transverse colon","Splenic flexure of colon","Descending colon",
                                         "Sigmoid colon","Rectosigmoid junction","Colon, NOS")



ggplot(data = y, aes(x = tissue_or_organ_of_origin, y = ITH)) +
  geom_violin() +
  geom_jitter(fill = "#003366", colour = "#003366", size = 6, width = 0.2, height = 0) +
  xlab("Patient CNA diversity") +
  ylab("Patient CNA diversity (predicted)") +
  theme(plot.margin = unit(c(t=0,r=0,b=0,l=0), "cm"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.title = element_text(size=24, colour='black'),
        axis.text.x = element_text(size=24, colour='black', angle = 90),
        axis.text.y = element_text(size=24, colour='black'),
        legend.position = "none") 


