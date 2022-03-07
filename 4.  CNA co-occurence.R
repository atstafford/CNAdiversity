# HEATMAP OF BIN-BIN PAIRWISE CORRELATION WITH HISTOGRAMS

# Prepare two histograms to display frequency of clonal and subclonal gains and losses
# Prepare long form of GaLo.Clo data for histograms
data <- gather(car.clonality$CNA.clo.counts,variable,freq,-bin,-chr1,-chr2)


# Prepare data for heatmaps (x4) showing genome-wide correlations between:
# diploid/aneuploid, diploid/gain, diploid/loss, loss/gain, and loss/diploid/gain

# Having two multiregion samples per patient will drive up the correlation. 
# Therefore, we will consider the bin-bin correlations using the average copy number per patient
cor.inputs <- list()

x <- c(1,5)
# For each matrix 
for ( d in 1:length(car.matrices) ) {
  i <- 1
  l <- 1
  list <- list()
  
  # For a given bin
  for ( k in 1:nrow(car.matrices[[d]]) ) { 
    i <- 1
    
    # For a each patient, average the copy number between samples
    while ( i < ncol(car.matrices[[d]]) ) { 
      list[[l]] <- ((car.matrices[[d]][k,i]) + (car.matrices[[d]][k,i+1])) / 2 
      i <- i + 2
      l <- l + 1
    }
  }
  cor.inputs[[d]] <- data.frame(t(matrix(unlist(list), ncol=car.info$noBins)), check.names = FALSE)
}

# Add list titles
names(cor.inputs) = c('diploid.aneu', 'diploid.gain', 'diploid.loss', 'loss.gain', 'loss.dip.gain')

# Convert to numeric matrices
cor.inputs <- lapply(cor.inputs, function(x) data.matrix(x,rownames.force = NA))


# Generate bin-bin corrleation matrices

# Calculate pairwise correlations
cormat.list <- lapply(cor.inputs, function(x) t(x))
cormat.list <- lapply(cormat.list, function(x) cor(x, method = 'pearson', use = 'pairwise.complete.obs'))
cormat.list <- lapply(cormat.list, function(x) {
  colnames(x) <- 1:ncol(x)
  rownames(x) <- 1:nrow(x)
  x
})

# Melt into long form for geom_tile
meltedcormat.list <- cormat.list
meltedcormat.list <- lapply(meltedcormat.list, function(x) melt(x, na.rm = FALSE))
meltedcormat.list <- lapply(meltedcormat.list, "colnames<-", c('first_bin','second_bin','correlation'))

# add 2x chr data for faceting
meltedcormat.list <- lapply(meltedcormat.list, function(x) {
  x <- cbind(x, chr1 = as.numeric(car.info$start.stop[match(x$first_bin, car.info$start.stop$bin), 2])) 
  x <- cbind(x, chr2 = as.factor(car.info$start.stop[match(x$second_bin, car.info$start.stop$bin), 2]))
  x$chr2 <- factor(x$chr2, levels = rev(c(seq(1,22,1))) ) # to order xaxis facets
  x
} )


# Prepare heatmaps 
hm.list <- list()
for ( i in 1:length(meltedcormat.list) ) {
  hm.list[[i]] <- ggplot(meltedcormat.list[[i]], aes(first_bin, second_bin, fill = correlation)) +
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
}

hm.legend <- as_ggplot(cowplot::get_legend(hm.list[[1]] + 
                                             guides(fill = guide_colorbar(title="Correlation", label.position = "bottom",
                                                                          title.position = "left", title.vjust = 0.9)) +
                                             theme(plot.margin = unit(c(t=-100,r=0,b=-100,l=0), "cm"),
                                                   legend.position = "bottom",legend.direction="horizontal",
                                                   legend.title = text.bold,
                                                   legend.margin = margin(grid::unit(c(t=-100,r=0,b=-100,l=0),"cm")),
                                                   legend.text = text,
                                                   legend.key.height = grid::unit(0.8,"cm"),
                                                   legend.key.width = grid::unit(1.4,"cm")) ))   


hmPlot <- cowplot::plot_grid(hm.legend, hm.list[[i]], ncol = 1, align = 'v', rel_heights = c(0.09, 1))

jpeg('tempfig.jpeg', width = 1000, height = 1000)
hmPlot
dev.off()
