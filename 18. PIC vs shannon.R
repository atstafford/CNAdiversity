library(vegan)
library(SciViews)

# Define the maximum shan depending on the number of samples
max.pics <- list()
for ( d in 1:50 ) {
  if ( d %% 3 == 0 ) {
    n1 <- n2 <- n3 <- d/3
  }
  else if ( d %% 3 == 1 ) {
    n1 <- ((d-1)/3) + 1
    n2 <- ((d-1)/3)
    n3 <- ((d-1)/3)
  }
  else if ( d %% 3 == 2 ) {
    n1 <- ((d-2)/3) + 1
    n2 <- ((d-2)/3) + 1
    n3 <- ((d-2)/3)
  }
  max.pics[[d]] <- pic.score <- 1 - ( (n1/d)^2 + (n2/d)^2 + (n3/d)^2 )
}

# Define the maximum shan score depending on the number of samples (up to 13)
max.shan <- list()
for ( d in 1:50 ) {
  if ( d %% 3 == 0 ) {
    n1 <- n2 <- n3 <- d/3
  }
  else if ( d %% 3 == 1 ) {
    n1 <- ((d-1)/3) + 1
    n2 <- ((d-1)/3)
    n3 <- ((d-1)/3)
  }
  else if ( d %% 3 == 2 ) {
    n1 <- ((d-2)/3) + 1
    n2 <- ((d-2)/3) + 1
    n3 <- ((d-2)/3)
  }
  max.shan[[d]] <- -1 * ( ((n1/d)*ln(n1/d)) + ((n2/d)*ln(n2/d)) + ((n3/d)*ln(n3/d)) )
}

# Calculate shannon using all samples for each patient
shan <- list()
for ( i in 1:length(testcnBinned.list)) {
  # Set working data as the cols holding the samples
  wd <- testcnBinned.list[[i]][,-c(1:4)]
  
  # Calculate freq of loss/dip/gain per bin
  x <- wd
  x$loss <- rowSums(wd == 1)
  x$diploid <- rowSums(wd == 2)
  x$gain <- rowSums(wd == 3)
  
  # Calculate Shannon div
  shan[[i]] <- (diversity(x[,c((ncol(x)-2):ncol(x))], index="shannon"))
  shan[[i]] <- na.omit(shan[[i]])
  
  # Store as a dataframe
  shan[[i]] <- data.frame(patient = i, shan = sum(shan[[i]], na.rm = TRUE), sampleNo = ncol(wd))
}

testShan <- do.call('rbind',shan)

# Calculate shan.frac using all samples for each patient
shan <- list()
for ( i in 1:length(testcnBinned.list)) {
  # Set working data as the cols holding the samples
  wd <- testcnBinned.list[[i]][,-c(1:4)]
  
  # Calculate freq of loss/dip/gain per bin
  x <- wd
  x$loss <- rowSums(wd == 1)
  x$diploid <- rowSums(wd == 2)
  x$gain <- rowSums(wd == 3)
  
  # Calculate Shannon div
  shan[[i]] <- (diversity(x[,c((ncol(x)-2):ncol(x))], index="shannon"))
  shan[[i]] <- na.omit(shan[[i]])
  
  # Record the number of samples
  upto <- ncol(wd)
  
  # Define the max possible diversity given the number of sample, as maxPIC*number of bins
  max.ITH <- max.shan[[upto]] * length(shan[[i]])
  
  # Store as a dataframe
  shan[[i]] <- data.frame(patient = i, shan.frac = sum(shan[[i]], na.rm = TRUE)/max.ITH, sampleNo = ncol(wd))
}

testShan.frac <- do.call('rbind',shan)

# Calculate PIC using all samples for each patient
pic <- list()
for ( i in 1:length(testcnBinned.list)) {
  # Set working data as the cols holding the samples
  wd <- testcnBinned.list[[i]][,-c(1:4)]
  
  # Record the number of samples
  upto <- ncol(wd)
  
  # Use PIC function on wd
  pic[[i]] <- PIC(wd, upto, c(1:upto))
  pic[[i]] <- na.omit(pic[[i]])
  
  # Store as a dataframe
  pic[[i]] <- data.frame(patient = i, PIC = sum(pic[[i]], na.rm = TRUE), sampleNo = ncol(wd))
}

testPIC <- do.call('rbind',pic)

# Calculate PIC.frac using all samples for each patient
pic <- list()
for ( i in 1:length(testcnBinned.list)) {
  # Set working data as the cols holding the samples
  wd <- testcnBinned.list[[i]][,-c(1:4)]
  
  # Record the number of samples
  upto <- ncol(wd)
  
  # Use PIC function on wd
  pic[[i]] <- PIC(wd, upto, c(1:upto))
  pic[[i]] <- na.omit(pic[[i]])
  
  # Define the max possible diversity given the number of sample, as maxPIC*number of bins
  max.ITH <- max.pics[[upto]] * length(pic[[i]])
  
  # Store as a dataframe
  pic[[i]] <- data.frame(patient = i, PIC.frac = sum(pic[[i]], na.rm = TRUE)/max.ITH, sampleNo = ncol(wd))
}

testPIC.frac <- do.call('rbind',pic)

test <- merge(testPIC, testPIC.frac, by=c('patient','sampleNo'))
test <- merge(test, testShan, by=c('patient','sampleNo'))
test <- merge(test, testShan.frac, by=c('patient','sampleNo'))
test$PIC.frac<- test$PIC.frac*800
test$shan.frac <- test$shan.frac*800
test$type <- "test"

x <- rbind(test)

comp <- gather(x, Metric, count, -patient, -type, -sampleNo)
comp$Value <- "Raw"
comp$Value[which(comp$Metric=="shan.frac" | comp$Metric=="PIC.frac")] <- "As fraction" 
comp$Value <- as.factor(comp$Value)
comp$Value <- factor(comp$Value, levels = c("Raw","As fraction"))
comp$Method <- "PIC"
comp$Method[which(comp$Metric=="shan.frac" | comp$Metric=="shan")] <- "Shannon" 

labels <- list()
for (i in 1:length(unique(x$patient))) {
  s <- x$sampleNo[which(x$patient==unique(x$patient)[i])]
  labels[[i]] <- paste(i,"(n=",s,")", sep = "", collapse = "")
}
labels <- unlist(labels)


plot <- ggplot(data = comp, aes(x=reorder(patient, count) , y=count, group=Metric)) +
  geom_point(size=3, aes(color=Method)) +
  geom_line(aes(linetype=Value, color=Method)) +
  scale_y_continuous(name="Raw diversity values", 
                     sec.axis = sec_axis(trans = ~./800, name="Values as fraction of maximum diversity")) +
  scale_x_discrete(name="Patient (sample number)", labels=labels)+
  scale_linetype_manual(values=c("solid","dotted")) +
  theme(plot.margin=margin(t=0,r=0.5,b=0.5,l=0.5,"cm"),
        panel.background = element_blank(),
        plot.title = element_text(size=24, colour='black',face='bold'),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.title = element_text(size=24, colour='black'),
        axis.text.y = element_text(size=24, colour='black', hjust=1),
        axis.text.x = element_text(size=15, colour='black', angle = 90),
        axis.ticks.length=unit(0.2, "cm"),
        legend.position = "top",
        legend.text = element_text(size=24, colour='black'),
        legend.title = element_blank())  

jpeg('tempfig.jpeg', width = 600, height = 600)
plot
dev.off()



