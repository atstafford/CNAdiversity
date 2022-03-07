# SPREAD ACROSS SAMPLES---------------------------------------------------------
# Are the beta bins spread evenly across both samples? What is the difference when using sample 1 vs 2
x <- rbind(predictedIth_train, predictedIth_test)
y <- gather(x[x$type=="train",c(1:5)], color, value, -sample, -type, -patient)
y$patient = with(y[which(y$color=="actual"),], reorder(patient, value, median))


plot <- ggplot(y, aes(x=patient, y=value, colour=color)) +
  #geom_violin(trim=FALSE, size=1) +
  stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.05, colour='black', size = 0.5)+
  geom_point(shape=16, size = 4) +
  scale_color_manual(values=c("#999999", alpha("#E69F00",0.7))) +
  theme(plot.margin = unit(c(t=0,r=0,b=0,l=0), "cm"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.title = element_text(size=24, colour='black'),
        axis.text.y = element_text(size=24, colour='black'),
        axis.text.x = element_text(size=20, colour='black', angle = 90),
        legend.position = "none") 

jpeg('tempfig.jpeg', width = 1000, height = 1000)
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
mean(unlist(diff))
median(unlist(diff))
