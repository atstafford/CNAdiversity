# SURVIVAL----------------------------------------------------------------------
train_clin <- read_excel("~/Documents/CNA/Data/Table_S3.Fig1.Clinical_Pathology.xlsx")
train_clin$ITH <- predictedIth_train[match(train_clin$`sampleID 2`, predictedIth_train$patient),4]
train_clin$actualITH <- predictedIth_train[match(train_clin$`sampleID 2`, predictedIth_train$patient),3]

#remove MSI
train_clin <- train_clin[which(train_clin$MSI!="MSI"),]

train_clin$preddivGroup2 <- "low"
train_clin$preddivGroup2[train_clin$ITH>= median(train_clin$ITH) ] <- "high"
train_clin$actdivGroup2 <- "low"
train_clin$actdivGroup2[train_clin$actualITH>= median(train_clin$actualITH) ] <- "high"

# Perform chi squared test
chisq2 <- chisq.test(table(train_clin$actdivGroup2, train_clin$T))
corrplot(chisq2$residuals, is.cor = FALSE)
chisq2$p.value

# Pull diversity by stage into summary table
train_survival_summary <- as.data.frame(table(train_clin$Stage, train_clin$preddivGroup2))

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

# split by stage
train_clin$preddivSGroup2 <- "low"
for (i in 1:length(unique(train_clin$N))) {
  stage <- unique(train_clin$N)[i]
  median <- median(train_clin$ITH[which(train_clin$N==stage)])
  patients <- train_clin$`sampleID 2`[which(train_clin$N==stage & train_clin$ITH>=median)]
  train_clin$preddivSGroup2[which(train_clin$`sampleID 2` %in% patients)] <- "high"
}

train_clin$actdivSGroup2 <- "low"
for (i in 1:length(unique(train_clin$N))) {
  stage <- unique(train_clin$N)[i]
  median <- median(train_clin$actualITH[which(train_clin$N==stage)])
  patients <- train_clin$`sampleID 2`[which(train_clin$N==stage & train_clin$actualITH>=median)]
  train_clin$actdivSGroup2[which(train_clin$`sampleID 2` %in% patients)] <- "high"
}

# Fit survival data using the Kaplan-Meier method
x <- data.frame(table(train_clin$N, train_clin$preddivSGroup2))
stage <- as.character(unique(x$Var1))
n <- list()
for (i in 1:length(stage)) {
  n[[i]] <- paste(stage[i],
                  paste("(n=", sum(x$Freq[which(x$Var1==stage[i])]), ")",sep = ""),
                  paste("\nhigh=", x$Freq[which(x$Var1==stage[i] & x$Var2=="high")]," low=", x$Freq[which(x$Var1==stage[i] & x$Var2=="low")], sep = ""))
}


plot.list <- list()
for (i in 1:length(stage)) {
  surv_object <- Surv(time = train_clin$OS_time[which(train_clin$N==stage[i])], 
                      event = train_clin$OS_cens[which(train_clin$N==stage[i])])
  fit <- survfit(surv_object ~ preddivSGroup2, data = train_clin[which(train_clin$N==stage[i]),])
  plot.list[[i]] <- ggsurvplot(fit, data = train_clin[which(train_clin$N==stage[i]),], pval = TRUE, 
                               title=n[i], legend="none",
                               font.main=24, font.x=NA, font.y=NA, font.tickslab=24, font.legend=24, pval.size=9, 
                               #palette = c('#00CC66','#FF9933'), legend.title='Predicted CNA diversity', legend.labs=c('high','low')
                               )
}

# Plot
plot <- plot_grid(plot.list[[1]]$plot, plot.list[[2]]$plot, plot.list[[3]]$plot, plot.list[[4]]$plot)
plot <- plot_grid(plot.list[[1]]$plot, plot.list[[2]]$plot, plot.list[[3]]$plot)

plot <- plot_grid(get_legend(divByStagePlot), plot, nrow = 2, rel_heights = c(1,10))

plot <- annotate_figure(plot, 
                        left = text_grob('Survival probability', rot = 90, size = 24),
                        bottom = text_grob('Time (years)', size = 24))
plot

jpeg('tempfig.jpeg', width = 1000, height = 1000)
plot
dev.off()



