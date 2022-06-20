hg19 <- makeHg19()
hg19$start <- 1
hg19 <- hg19[c(1,6,2,3)]
hg19$chrom <- as.character(as.numeric(hg19$chrom))
hg19 <- hg19[which(hg19$chrom!=23),]
hg19 <- hg19[which(hg19$chrom!=24),]

annotation2 <- car.info$start.stop
annotation2$cluster <- sig.hclust[match(annotation2$bin, sig.hclust$bin), 2]

# define new cluster
new.clust <- max(annotation2$cluster[is.finite(annotation2$cluster)]) + 1
for (i in 1:nrow(annotation2) ) {
  if ( is.na(annotation2$cluster[i]) ) {
    annotation2$cluster[i] <- new.clust
    new.clust <- new.clust + 1
  }
}

# define breakpoints
i <- 1
l <- 1
breakpoint<- list()

for ( k in 1:nrow(annotation2) ) {
  if ( k == 1 ) {
    breakpoint[[l]] <- i
    l <- l + 1
  }
  else {
    if ( annotation2$cluster[k] == annotation2$cluster[k-1] ) {
      breakpoint[[l]] <- i
      l <- l + 1
    }
    else {
      i <- i + 1
      breakpoint[[l]] <- i
      l <- l + 1
    }
  }
}
annotation2$breakpoint <- unlist(breakpoint)

# link
list <- list()
l <- 1
for (i in 1:length(unique(annotation2$cluster))) {
  cluster <- unique(annotation2$cluster)[i]
  wd <- annotation2[which(annotation2$cluster==cluster),]
  if ( length(unique(wd$breakpoint)) == 2 ) {
    list[[l]] <- unique(wd$breakpoint)
    l <- l + 1
  }
  else if ( length(unique(wd$breakpoint)) == 3 ) {
    list[[l]] <- c(unique(wd$breakpoint)[1], unique(wd$breakpoint)[2])
    l <- l + 1
    list[[l]] <- c(unique(wd$breakpoint)[1], unique(wd$breakpoint)[3])
    l <- l + 1
    list[[l]] <- c(unique(wd$breakpoint)[2], unique(wd$breakpoint)[3])
    l <- l + 1
  }
  else {next}
}

# outline significant clusters that can predict ITH
x <- cluster.info$cluster[which(cluster.info$sig=="sig")]
annotation2$x <- NA
annotation2$x[which((annotation2$breakpoint %% 2) != 0 & annotation2$cluster %!in% x)] <- "1insig"
annotation2$x[which((annotation2$breakpoint %% 2) == 0 & annotation2$cluster %!in% x)] <- "2insig"
annotation2$x[which((annotation2$breakpoint %% 2) != 0 & annotation2$cluster %in% x)] <- "1sig"
annotation2$x[which((annotation2$breakpoint %% 2) == 0 & annotation2$cluster %in% x)] <- "2sig"

link <- data.frame(matrix(unlist(list),byrow = T,ncol = 2), check.names = F)
link <- data.frame(ann1=link$`1`, 1 , ann2=link$`2`, 1) 
link <- link[-c(2,12), ] #small=19,132
link$ann1 <- as.character(link$ann1)
link$ann2 <- as.character(link$ann2)
link$ann1[which(link$ann1==15)] <- 231
link$ann2[which(link$ann2==17)] <- 255
link$ann1[which(link$ann1==19)] <- 262
link$ann2[which(link$ann2==21)] <- 275
link$ann1[which(link$ann1==66)] <- 1201
link$ann2[which(link$ann2==150)] <- 2589
link$ann1[which(link$ann1==72)] <- 1404
link$ann2[which(link$ann2==74)] <- 1420
link$ann1[which(link$ann1==94)] <- 1554
link$ann2[which(link$ann2==96)] <- 1579
link$ann1[which(link$ann1==104)] <- 1754
link$ann2[which(link$ann2==106)] <- 1787
link$ann1[which(link$ann1==110)] <- 1807
link$ann1[which(link$ann1==131)] <- 2274
link$ann2[which(link$ann2==147)] <- 2517
link$ann1[which(link$ann1==116)] <- 1902
link$ann2[which(link$ann2==118)] <- 1949
link$ann1[which(link$ann1==122)] <- 2106
link$ann2[which(link$ann2==124)] <- 2141
link$ann1[which(link$ann1==132)] <- 2277
link$ann2[which(link$ann2==134)] <- 2303
link$ann1[which(link$ann1==135)] <- 2313
link$ann2[which(link$ann2==137)] <- 2331

annotation2 <- annotation2[c(1,2,3,4,7,5)]
x <- c("#cccccc", "#666666", "#ff9933", "#993366")


#lung
annotation2$x <- NA
annotation2$x[which((annotation2$breakpoint %% 2) != 0)] <- "1insig"
annotation2$x[which((annotation2$breakpoint %% 2) == 0)] <- "2insig"

link <- data.frame(matrix(unlist(list),byrow = T,ncol = 2), check.names = F)
link <- data.frame(ann1=link$`1`, 1 , ann2=link$`2`, 1) 
#link <- link[-c(2,12), ] #small=19,132

link$ann1[which(link$ann1==4)] <- 125
link$ann2[which(link$ann2==6)] <- 131
link$ann1[which(link$ann1==15)] <- 282
link$ann2[which(link$ann1==16)] <- 208
link$ann1[which(link$ann2==17)] <- 210
link$ann2[which(link$ann2==18)] <- 216
link$ann2[which(link$ann1==53)] <- 1429
link$ann2[which(link$ann2==55)] <- 1435
link$ann2[which(link$ann1==62)] <- 1874
link$ann2[which(link$ann2==64)] <- 1900
link$ann2[which(link$ann1==69)] <- 2105
link$ann2[which(link$ann2==71)] <- 2108
link$ann2[which(link$ann1==103)] <- 2395
link$ann2[which(link$ann1==104)] <- 2397
link$ann2[which(link$ann2==105)] <- 2398
link$ann2[which(link$ann2==106)] <- 2401
link$ann2[which(link$ann1==112)] <- 2566
link$ann2[which(link$ann2==114)] <- 2588

annotation2 <- annotation2[c(1,2,3,4,7,5)]
x <- c("#cccccc", "#666666")

chromoMap(list(hg19), list(annotation2),
          #cna
          data_based_color_map = T, 
          data_type = "categorical",
          data_colors = list(c(x)),
          #link
          show.links = T,
          loci_links = link,
          links.colors = "black", 
          #features
          canvas_width = 800, 
          canvas_height = 500,
          chr_color = c("#CCCCCC"),
          chr_width = 10, chr_length = 4, ch_gap = 4, 
          #y_chr_scale = -30, top_margin = -20,
          legend = F)

jpeg('tempfig.jpeg', width = 1000, height = 1000)
plot
dev.off()