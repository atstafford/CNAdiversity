# SPREAD ACROSS SAMPLES---------------------------------------------------------
# Are the beta bins spread evenly across both samples? What is the difference when using sample 1 vs 2
x <- rbind(predictedIth_train, predictedIth_test)
y <- gather(x[x$type=="train",c(1:5)], 'CNA diversity', value, -sample, -type, -patient)


plot <- ggplot(y, aes(x=fct_reorder(patient, value, mean) , y=value, shape =`CNA diversity`, colour=`CNA diversity`)) +
  #stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.05, colour='black', size = 0.5) +
  geom_line(aes(group = patient), size=1.2, colour = "#003366") +
  geom_point(size = 4) +
  scale_shape_manual(values=c(15, 16)) +
  scale_color_manual(values=c("#003366", alpha("#FF6666",1))) +
  ylab("CNA diversity") +
  xlab("Patient") +
  coord_flip() +
  theme(plot.margin = unit(c(t=0,r=0,b=0,l=0), "cm"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.title = element_text(size=24, colour='black'),
        axis.text.y = element_text(size=12, colour='black'),
        axis.text.x = element_text(size=20, colour='black', angle = 0),
        legend.position = "top", legend.text = element_text(size=24), legend.title = element_text(size=24)) +
  guides(color = guide_legend(direction = 'horizontal', label.position = "left", label.hjust = 1))
  
jpeg('tempfig.jpeg', width = 500, height = 1000)
plot
dev.off()

# difference between two samples
diff <- list()
error <- list()
for ( i in 1:length(unique(x$patient[which(x$type=="train")])) ) {
  patient <- unique(x$patient[which(x$type=="train")])[i]
  wd <- x[which(x$patient ==  patient),]
  diff[[i]] <- max(wd$predicted)-min(wd$predicted)
}
sampleDiff <- unlist(diff)
mean(unlist(diff))
median(unlist(diff))
max(unlist(diff))
min(unlist(diff))

# difference between actual and predicted
diff <- list()
wd <- x[which(x$type == "train"),]
for ( i in 1:nrow(wd) ) {
  if ( wd$actual[i] >= wd$predicted[i] ) {
    diff[[i]] <- wd$actual[i] - wd$predicted[i]
  }
  else if ( wd$actual[i] < wd$predicted[i] ) {
    diff[[i]] <- wd$predicted[i] - wd$actual[i]
  }
}
error <- unlist(diff)
mean(unlist(diff))
median(unlist(diff))

# difference between actual and predicted ave per patient
diff <- list()
meandiff <- list()
wd <- x[which(x$type == "train"),]
for ( k in 1:length(unique(wd$patient)) ) {
  wd2 <- wd[which(wd$patient == unique(wd$patient)[k]),]
  for ( i in 1:nrow(wd2) ) {
    if ( wd2$actual[i] >= wd2$predicted[i] ) {
      diff[[i]] <- wd2$actual[i] - wd2$predicted[i]
    }
    else if ( wd2$actual[i] < wd2$predicted[i] ) {
      diff[[i]] <- wd2$predicted[i] - wd2$actual[i]
    }
  }
  meandiff[[k]] <- mean(unlist(diff))
}

meanerror <- unlist(meandiff)

# whats corr with error
x <- predictedIth_train
y <- car.clonality$patientClo
y <- (car.clonality$patientClo$subclonal/car.clonality$patientClo$CNA)
y <- data.frame(ITH=lapply(data.frame(ITH=y), rep, car.info$sampPerPatient))
x$pcsubclonal <- y$ITH

y <- car.diversity$pga$prop.aneu
x$pga <- y

x$sampleerror <- unlist(error)
x$patienterror <- unlist(data.frame(meanerror=lapply(data.frame(meanerror=meanerror), rep, car.info$sampPerPatient)))
x$sampleDiff <- unlist(data.frame(sampleDiff=lapply(data.frame(sampleDiff=sampleDiff), rep, car.info$sampPerPatient)))
  

ggplot(data = x, aes(x = pcsubclonal, y = sampleerror)) +
  geom_line(aes(group = patient), size=0.2, colour = "#003366") +
  geom_point(fill = "#003366", colour = "#003366", size = 6) +
  theme(plot.margin = unit(c(t=0,r=0,b=0,l=0), "cm"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.title = element_text(size=24, colour='black'),
        axis.text = element_text(size=24, colour='black'),
        legend.position = "none") +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'right', label.y = 'top', size=8) 


# cut off should be 0.3, above which it is inaccurate








