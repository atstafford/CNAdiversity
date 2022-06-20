# SURVIVAL----------------------------------------------------------------------
tracerx_clin <- read_excel("~/Documents/CNA/Data/TracerX/cn_data.xlsx", sheet = 1, skip = 1)

#use mean of samples per patient
x <- aggregate(x = predictedIth_tracerx$predicted,                # Specify data column
          by = list(predictedIth_tracerx$patient),              # Specify group indicator
          FUN = median, na.rm =TRUE)

tracerx_clin$ITH <- x[match(tracerx_clin$TRACERxID, x$Group.1),2]
tracerx_clin$actualITH <- predictedIth_tracerx[match(tracerx_clin$TRACERxID, predictedIth_tracerx$patient),3]

# surviuval per samples
x <- predictedIth_tracerx[c(1,4)]
x$`Time_to_recurrence_or_death (months)` <- unlist(tracerx_clin[match(x$patient, tracerx_clin$TRACERxID), 17])
x$`Recurrence or death` <- unlist(tracerx_clin[match(x$patient, tracerx_clin$TRACERxID), 16])
x$Stage <- unlist(tracerx_clin[match(x$patient, tracerx_clin$TRACERxID), 2])
colnames(x)[2] <- "ITH"

tracerx_clin <- x

# Keep only those samples whose ITH could be predicted
tracerx_clin <- tracerx_clin[!is.na(tracerx_clin$ITH),]

#remove MSI

#change stage names
for ( i in 1:nrow(tracerx_clin) ) {
  tracerx_clin$Stage[i] <- ifelse(tracerx_clin$Stage[i]=="2a", "Stage II", 
                                   ifelse(tracerx_clin$Stage[i]=="1b", "Stage I",
                                          ifelse(tracerx_clin$Stage[i]=="3a", "Stage III",
                                                 ifelse(tracerx_clin$Stage[i]=="2b", "Stage II",
                                                        ifelse(tracerx_clin$Stage[i]=="1a", "Stage I",
                                                               ifelse(tracerx_clin$Stage[i]=="3b", "Stage III",tracerx_clin$Stage[i]))))))
}

# ITH per stage ----------------------------------------------------------------
# Perform chi squared test
chisq2 <- chisq.test(table(tracerx_clin$preddivGroup2, tracerx_clin$Stage))
corrplot(chisq2$residuals, is.cor = FALSE)
chisq2$p.value

# Pull diversity by stage into summary table
train_survival_summary <- as.data.frame(table(tracerx_clin$Stage, tracerx_clin$preddivGroup2))

# Plot
divByStagePlot <- ggplot(train_survival_summary, aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_bar(position="fill", stat="identity", width=0.8) +
  scale_fill_manual(values = c('#00CC66','#FF9933'), labels = c('High','Low')) + 
  ylab('Fraction') +
  labs(fill = 'Predicted CNA diversity') +
  guides(fill = guide_legend(direction = 'horizontal', label.position = "left", label.hjust = 1)) +
  #scale_x_discrete(labels=c('Stage I\n(n=61)','Stage II\n(n=133)','Stage III\n(n=99)','Stage IV\n(n=47)')) +
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

# Survival of high/low diversity patients (KM) 

# Split into groups --------------------------------------------------------------
tracerx_clin$preddivGroup2 <- "low"
tracerx_clin$preddivGroup2[tracerx_clin$ITH>= median(tracerx_clin$ITH) ] <- "high"
tracerx_clin$actdivGroup2 <- "low"
tracerx_clin$actdivGroup2[tracerx_clin$actualITH>= median(tracerx_clin$actualITH) ] <- "high"


# split by stage
tracerx_clin$preddivSGroup2 <- "low"
for (i in 1:length(unique(tracerx_clin$Stage))) {
  stage <- unique(tracerx_clin$Stage)[i]
  median <- median(tracerx_clin$ITH[which(tracerx_clin$Stage==stage)])
  patients <- tracerx_clin$TRACERxID[which(tracerx_clin$Stage==stage & tracerx_clin$ITH>=median)]
  tracerx_clin$preddivSGroup2[which(tracerx_clin$TRACERxID %in% patients)] <- "high"
}

tracerx_clin$divGroup4 <- ntile(tracerx_clin$ITH, 4)
tracerx_clin$divGroup3 <- "med"
tracerx_clin$divGroup3[which(tracerx_clin$divGroup4==4)] <- "top25"
tracerx_clin$divGroup3[which(tracerx_clin$divGroup4==1)] <- "bottom25"

cut <- surv_cutpoint(tracerx_clin, time = "Time_to_recurrence_or_death (months)", event = "Recurrence or death", variables = "ITH")
res.cat <- surv_categorize(cut)
tracerx_clin$optim <- res.cat$ITH

# top and bottom 25% -----
z <- tracerx_clin[which(tracerx_clin$divGroup3!="med"),]

x <- data.frame(table(z$Stage, z$divGroup3))
stage <- as.character(unique(x$Var1))
n <- list()
for (i in 1:length(stage)) {
  n[[i]] <- paste(stage[i],
                  paste("(n=", sum(x$Freq[which(x$Var1==stage[i])]), ")",sep = ""),
                  paste("\nhigh=", x$Freq[which(x$Var1==stage[i] & x$Var2=="top25")]," low=", x$Freq[which(x$Var1==stage[i] & x$Var2=="bottom25")], sep = ""))
}
plot.list <- list()
for (i in 1:length(stage)) {
  surv_object <- Surv(time = z$`Time_to_recurrence_or_death (months)`[which(z$Stage==stage[i])], 
                      event = z$`Recurrence or death`[which(z$Stage==stage[i])])
  fit <- survfit(surv_object ~ divGroup3, data = z[which(z$Stage==stage[i]),])
  plot.list[[i]] <- ggsurvplot(fit, data = z[which(z$Stage==stage[i]),], pval = TRUE, 
                               title=n[i], legend="none",
                               font.main=24, font.x=NA, font.y=NA, font.tickslab=24, font.legend=24, pval.size=9)
}
legend <- get_legend(ggsurvplot(fit, data = z[which(z$Stage==stage[i]),])$plot)
plot <- plot_grid(plot.list[[1]]$plot, plot.list[[2]]$plot, plot.list[[3]]$plot, nrow = 2)
plot <- plot_grid(legend, plot, nrow = 2, rel_heights = c(1,10))
plot <- annotate_figure(plot, 
                        left = text_grob('Survival probability', rot = 90, size = 24),
                        bottom = text_grob('Time (years)', size = 24))
plot

fit <- survfit(Surv(time = z$`Time_to_recurrence_or_death (months)`, event = z$`Recurrence or death`) ~ divGroup3, data = z)
plot <- ggsurvplot(fit, data = z, pval = TRUE, legend="none",
                   font.main=24, font.x=NA, font.y=NA, font.tickslab=24, font.legend=24, pval.size=9)
hist <- ggplot(z, aes(ITH, fill = divGroup3)) + 
  geom_histogram(binwidth = 0.001) +
  #scale_fill_manual(values = c("red", "blue"), name = "ITH", breaks = c("high","low")) +
  theme(panel.background = element_blank(), legend.position = "top",axis.text = element_text(size=24, colour='black'),
        axis.line = element_line(colour = "black", size = 0.5),axis.title = element_text(size=24, colour='black'))

plot <- plot_grid(hist , plot$plot)
plot

# optim  ---------------------------------
fit <- survfit(Surv(tracerx_clin$`Time_to_recurrence_or_death (months)`, tracerx_clin$`Recurrence or death`) ~ optim, data = tracerx_clin)
plot <- ggsurvplot(fit, data = tracerx_clin, pval = TRUE, legend="none",
                   font.main=24, font.x=NA, font.y=NA, font.tickslab=24, font.legend=24, pval.size=9)

hist <- ggplot(tracerx_clin, aes(ITH, fill = optim)) + 
  geom_histogram(binwidth = 0.001) +
  #scale_fill_manual(values = c("red", "blue"), name = "ITH", breaks = c("high","low")) +
  theme(panel.background = element_blank(), legend.position = "top",axis.text = element_text(size=24, colour='black'),
        axis.line = element_line(colour = "black", size = 0.5),axis.title = element_text(size=24, colour='black'))

plot_grid(hist , plot$plot)

x <- data.frame(table(tracerx_clin$Stage, tracerx_clin$optim))
stage <- as.character(unique(x$Var1))
n <- list()
for (i in 1:length(stage)) {
  n[[i]] <- paste(stage[i],
                  paste("(n=", sum(x$Freq[which(x$Var1==stage[i])]), ")",sep = ""),
                  paste("\nhigh=", x$Freq[which(x$Var1==stage[i] & x$Var2=="high")]," low=", x$Freq[which(x$Var1==stage[i] & x$Var2=="low")], sep = ""))
}

plot.list <- list()
for (i in 1:length(stage)) {
  surv_object <- Surv(time = tracerx_clin$`Time_to_recurrence_or_death (months)`[which(tracerx_clin$Stage==stage[i])], 
                      event = tracerx_clin$`Recurrence or death`[which(tracerx_clin$Stage==stage[i])])
  fit <- survfit(surv_object ~ optim, data = tracerx_clin[which(tracerx_clin$Stage==stage[i]),])
  plot.list[[i]] <- ggsurvplot(fit, data = tracerx_clin[which(tracerx_clin$Stage==stage[i]),], pval = TRUE, 
                               title=n[i], legend="none",
                               font.main=24, font.x=NA, font.y=NA, font.tickslab=24, font.legend=24, pval.size=9)
}
legend <- get_legend(ggsurvplot(fit, data = tracerx_clin[which(tracerx_clin$Stage==stage[i]),])$plot)
plot <- plot_grid(plot.list[[1]]$plot, plot.list[[2]]$plot, plot.list[[3]]$plot, nrow = 2)
plot <- plot_grid(legend, plot, nrow = 2, rel_heights = c(1,10))
plot <- annotate_figure(plot, 
                        left = text_grob('Survival probability', rot = 90, size = 24),
                        bottom = text_grob('Time (years)', size = 24))
plot
