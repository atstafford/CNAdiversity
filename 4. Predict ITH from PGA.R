#PGA

# Predict ITH from PGA per sample -----------------------------------------------
pgaPredict.list <- list(ad.diversity$pga, car.diversity$pga)
names(pgaPredict.list) <- c('adenoma','carcinoma')
pgaPredict.list <- lapply(pgaPredict.list, function(x) {
  x$sample <- rownames(x)
  x$patient <- sub('\\..*', '',rownames(x))
  x <- gather(x,CNA,proportion,-sample, -patient, -pic.frac)
  x$mark <- 'no'
  x$mark[which(x$sample=="C01.1" | x$sample=="C20.1" | x$sample=="C19.2")] <- "yes"
  x$mark2 <- 'no'
  x$mark2[which(x$patient=="C01" | x$patient=="C20" | x$patient=="C19")] <- "yes"
  x
})

plot.list <- list()
for ( i in 1:length(pgaPredict.list) ) {
  
  if (i == 1) {tit <- 'Adenoma'} #adenoma
  if (i == 2) {tit <- 'Carcinoma'} #carcinoma
  
  plot.list[[i]] <- ggplot(data=pgaPredict.list[[i]][which(pgaPredict.list[[i]]$CNA=="prop.aneu"),], aes(y=pic.frac, x=proportion)) + 
    geom_smooth(method = "lm", formula = y ~ x, color="#336666") +
    geom_point(size=3, color="#003333") +
    geom_point(data=pgaPredict.list[[i]][which(pgaPredict.list[[i]]$mark=="yes" & pgaPredict.list[[i]]$CNA=="prop.aneu"),],size=4, shape=1, color="black",fill="#003333") +
    geom_text(data=pgaPredict.list[[i]][which(pgaPredict.list[[i]]$mark=="yes" & pgaPredict.list[[i]]$CNA=="prop.aneu"),], aes(label = sample), size=5, vjust=-1, check_overlap = T) +
    ggpmisc::stat_poly_eq(formula = y ~ x, 
                          aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                          parse = TRUE, label.x = 'right', label.y = 'top', size=7) +
    ggtitle(tit) +
    scale_y_continuous(expand = c(0.01,0), limits = c(0,0.6)) +
    scale_x_continuous(expand = c(0.01,0), limits = c(0,0.6), breaks = c(seq(0,0.6,0.1)), labels = c(seq(0,0.6,0.1))) +
    theme(plot.margin=margin(t=0,r=0.5,b=0.5,l=0.5,"cm"),
          panel.background = element_blank(),
          plot.title = element_text(size=24, colour='black',face='bold'),
          axis.line = element_line(colour = "black", size = 0.5),
          axis.title = element_blank(),
          axis.text = element_text(size=24, colour='black'),
          axis.ticks.length=unit(0.2, "cm"),
          legend.position = "none")  
}



pgaPredict.plot <- cowplot::plot_grid(plot.list[[1]], plot.list[[2]], ncol = 1, align = "v")

pgaPredict.plot <- annotate_figure(pgaPredict.plot, 
                                   left = text_grob('Patient CNA diversity', rot = 90, size = 24),
                                   bottom = text_grob('Proportion of genome altered', size = 24, vjust = 0))

jpeg('tempfig.jpeg', width = 600, height = 1000)
pgaPredict.plot
dev.off()





# Predict ITH from PGA in tracerX -----------------------------------------------
pgaPredict.list <- tracerx.diversity$pga
pgaPredict.list$sample <- rownames(pgaPredict.list)
pgaPredict.list$patient <- sub('\\..*', '',rownames(pgaPredict.list))
pgaPredict.list <- gather(pgaPredict.list,CNA,proportion,-sample, -patient, -pic.frac)

ggplot(data=pgaPredict.list, aes(y=pic.frac, x=proportion)) + 
    geom_smooth(method = "lm", formula = y ~ x, color="#336666") +
    geom_point(size=3, color="#003333") +
    #geom_point(data=pgaPredict.list,size=4, shape=1, color="black",fill="#003333") +
    #geom_text(data=pgaPredict.list, aes(label = sample), size=5, vjust=-1, check_overlap = T) +
    ggpmisc::stat_poly_eq(formula = y ~ x, 
                          aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                          parse = TRUE, label.x = 'right', label.y = 'top', size=7) +
    ggtitle("TRACERx") +
    scale_y_continuous(expand = c(0.01,0), limits = c(0,0.6)) +
    scale_x_continuous(expand = c(0.01,0), limits = c(0,0.6), breaks = c(seq(0,0.6,0.1)), labels = c(seq(0,0.6,0.1))) +
    theme(plot.margin=margin(t=0,r=0.5,b=0.5,l=0.5,"cm"),
          panel.background = element_blank(),
          plot.title = element_text(size=24, colour='black',face='bold'),
          axis.line = element_line(colour = "black", size = 0.5),
          axis.title = element_blank(),
          axis.text = element_text(size=24, colour='black'),
          axis.ticks.length=unit(0.2, "cm"),
          legend.position = "none")  

pgaPredict.plot <- cowplot::plot_grid(plot.list[[1]], ncol = 1, align = "v")

pgaPredict.plot <- annotate_figure(pgaPredict.plot, 
                                   left = text_grob('Patient CNA diversity', rot = 90, size = 24),
                                   bottom = text_grob('Proportion of genome altered', size = 24, vjust = 0))

jpeg('tempfig.jpeg', width = 600, height = 1000)
pgaPredict.plot
dev.off()


# Predicting ITH from average PGA -----------------------------------------------
pgaPredict.list <- list(ad.diversity$pga, car.diversity$pga)
pgaPredict.list <- lapply(pgaPredict.list, function(x) {
  x$patient <- sub('\\..*', '',rownames(x))
  x
})

for ( j in 1:length(pgaPredict.list) ) {
  for ( i in 1:nrow(pgaPredict.list[[j]]) ) {
    patient <- pgaPredict.list[[j]]$patient[i]
    wd <- pgaPredict.list[[j]]$prop.aneu[which(pgaPredict.list[[j]]$patient==patient)]
    pgaPredict.list[[j]]$avepga[i] <- mean(wd)
  }
}

pgaPredict.list <- lapply(pgaPredict.list, function(x) {
  rownames(x) <- NULL
  x <- unique(x[,c(1,5,6)])
  #x$mark <- "no"
  #x$mark[which(x$avepga < quantile(x$avepga, 0.5) & x$pic.frac > quantile(x$pic.frac, 0.5))] <- "red"
  #x$mark[which(x$avepga > quantile(x$avepga, 0.5) & x$pic.frac < quantile(x$pic.frac, 0.5))] <- "blue"
  x$mark <- 'no'
  x$mark[which(x$patient=="C01" | x$patient=="C20" | x$patient=="C19")] <- "yes"
  x
})

plot.list <- list()
for ( i in 1:length(pgaPredict.list) ) {
  
  if (i == 1) {tit <- 'Adenoma'} #adenoma
  if (i == 2) {tit <- 'Carcinoma'} #carcinoma
  
  plot.list[[i]] <- ggplot(data=pgaPredict.list[[i]], aes(y=pic.frac, x=avepga)) + 
    geom_smooth(method = "lm", formula = y ~ x, color="#336666") +
    geom_point(size=3, color="#003333") +
    geom_point(data=pgaPredict.list[[i]][which(pgaPredict.list[[i]]$mark=="yes"),],size=4, shape=1, color="black",fill="#003333") +
    geom_text(data=pgaPredict.list[[i]][which(pgaPredict.list[[i]]$mark=="yes"),], aes(label = patient), size=5, vjust=-1, check_overlap = T) +
    ggpmisc::stat_poly_eq(formula = y ~ x, 
                          aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                          parse = TRUE, label.x = 'right', label.y = 'top', size=7) +
    ggtitle(tit) +
    scale_y_continuous(expand = c(0.01,0), limits = c(0,0.6)) +
    scale_x_continuous(expand = c(0.01,0), limits = c(0,0.6), breaks = c(seq(0,0.6,0.1)), labels = c(seq(0,0.6,0.1))) +
    theme(plot.margin=margin(t=0,r=0.5,b=0.5,l=0.5,"cm"),
          panel.background = element_blank(),
          plot.title = element_text(size=24, colour='black',face='bold'),
          axis.line = element_line(colour = "black", size = 0.5),
          axis.title = element_blank(),
          axis.text = element_text(size=24, colour='black'),
          axis.ticks.length=unit(0.2, "cm"),
          legend.position = "none")  
}

pgaPredict.plot <- cowplot::plot_grid(plot.list[[1]], plot.list[[2]], ncol = 1, align = "v")

pgaPredict.plot <- annotate_figure(pgaPredict.plot, 
                                   left = text_grob('Patient CNA diversity', rot = 90, size = 24),
                                   bottom = text_grob('Average proportion of genome altered', size = 24, vjust = 0))

jpeg('tempfig.jpeg', width = 600, height = 1000)
pgaPredict.plot
dev.off()



# Mini clonality plots for C01 and C19 -----------------------------------------------
pat <- c("C01","C20","C19")
tit <- c("C01: ave. PGA=0.20, ITH=0.37", "C20: ave.PGA=0.21, ITH=0.17", "C19: ave.PGA=0.45, ITH=0.02")
plot.list <- list()
for (j in 1:length(pat) ) {
  x <- data.frame(bin=1:2694, car.raw[,which(sub('\\..*', '', colnames(car.raw))==pat[j])])
  x$val <- 0
  x$col <- "none"
  
  for ( i in 1:nrow(x) ) {
    if (x[i,3]==1 & x[i,2]==1) {
      x$val[i] <- -1
      x$col[i] <- "clonal"
    }
    if (x[i,3]==3 & x[i,2]==3) {
      x$val[i] <- 1
      x$col[i] <- "clonal"
    }
    if ((x[i,3]==1 & x[i,2]==2) | (x[i,3]==2 & x[i,2]==1)) {
      x$val[i] <- -1
      x$col[i] <- "subclonal"
    }
    if ((x[i,3]==3 & x[i,2]==2) | (x[i,3]==2 & x[i,2]==3)) {
      x$val[i] <- 1
      x$col[i] <- "subclonal"
    }
  }
  
  x <- x[which(x$col!="none"),]
  
  plot.list[[j]] <- ggplot() + 
    geom_bar(data=x, aes(fill=col, y=val, x=bin),
             position="stack", stat="identity", width = 1) +
    scale_fill_manual(values = c("#1B7837",'#762A83'), labels = c('clonal','subclonal')) + 
    ggtitle(tit[j]) +
    scale_x_continuous(expand = c(0,0), name="chromosome", breaks=car.info$chr.mid, labels = c(1:18,'\n19','20','\n21','22')) +
    scale_y_continuous(limits = c(-1,1), breaks = c(seq(-1,1,1)), labels = c("Loss","","Gain")) +
    geom_hline(yintercept = 0, size=0.3, color="black") +
    geom_rect(data = car.info$chr.end[which(car.info$chr.end$col=='W'),], 
              aes(NULL,NULL,xmin=start, xmax=end),
              fill = alpha("#CC9966", 0.1),
              ymin = -1,
              ymax = 1) +
    theme(plot.margin=margin(t=0,r=0.5,b=0.5,l=0.5,"cm"),
          panel.background = element_blank(),
          plot.title = element_text(size=24, colour='black',face='bold'),
          axis.line = element_line(colour = "black", size = 0.5),
          axis.title = element_blank(),
          axis.text = element_text(size=, colour='black'),
          axis.ticks.length=unit(0.2, "cm"),
          legend.position = "none") 
}

pgaPredict.plot <- cowplot::plot_grid(plot.list[[1]], plot.list[[2]], plot.list[[3]], ncol = 1, align = "v")

jpeg('tempfig.jpeg', width = 600, height = 600)
pgaPredict.plot
dev.off()







# PGA is similar between samples in patients with low ITH -----------------------------------------------
pgaPredict.list <- list(ad.diversity$pga, car.diversity$pga)
pgaPredict.list <- lapply(pgaPredict.list, function(x) {
  x$patient <- sub('\\..*', '',rownames(x))
  x
})

for ( j in 1:length(pgaPredict.list) ) {
  for ( i in 1:nrow(pgaPredict.list[[j]]) ) {
    patient <- pgaPredict.list[[j]]$patient[i]
    wd <- pgaPredict.list[[j]]$prop.aneu[which(pgaPredict.list[[j]]$patient==patient)]
    pgaPredict.list[[j]]$pgadiff[i] <- max(wd) - min(wd)
  }
}

pgaPredict.list <- lapply(pgaPredict.list, function(x) {
  rownames(x) <- NULL
  x <- unique(x[,c(1,5,6)])
  x
})

plot.list <- list()
for ( i in 1:length(pgaPredict.list) ) {
  
  if (i == 1) {tit <- 'Adenoma'} #adenoma
  if (i == 2) {tit <- 'Carcinoma'} #carcinoma
  
  plot.list[[i]] <- ggplot(data=pgaPredict.list[[i]], aes(y=pic.frac, x=pgadiff)) + 
    geom_smooth(method = "lm", formula = y ~ x, color="#336666") +
    geom_point(size=3, color="#003333") +
    ggpmisc::stat_poly_eq(formula = y ~ x, 
                          aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                          parse = TRUE, label.x = 'right', label.y = 'top', size=7) +
    ggtitle(tit) +
    
    scale_y_continuous(expand = c(0.01,0), limits = c(0,0.6)) +
    scale_x_continuous(expand = c(0.01,0), limits = c(0,0.3), breaks = c(seq(0,0.3,0.1)), labels = c(seq(0,0.3,0.1))) +
    theme(plot.margin=margin(t=0,r=0.5,b=0.5,l=0.5,"cm"),
          panel.background = element_blank(),
          plot.title = element_text(size=24, colour='black',face='bold'),
          axis.line = element_line(colour = "black", size = 0.5),
          axis.title = element_blank(),
          axis.text = element_text(size=24, colour='black'),
          axis.ticks.length=unit(0.2, "cm"),
          legend.position = "none")  
}

pgaPredict.plot <- cowplot::plot_grid(plot.list[[1]], plot.list[[2]], nrow = 2, align = "v")

pgaPredict.plot <- annotate_figure(pgaPredict.plot, 
                                   left = text_grob('Patient CNA diversity', rot = 90, size = 24),
                                   bottom = text_grob('Difference in pga between samples', size = 24, vjust = 0))

jpeg('tempfig.jpeg', width = 600, height = 1000)
pgaPredict.plot
dev.off()

