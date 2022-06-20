hg19 <- makeHg19()
hg19$start <- 1
hg19 <- hg19[c(1,6,2,3)]
hg19$chrom <- as.character(as.numeric(hg19$chrom))
hg19 <- hg19[which(hg19$chrom %in% hg19predictors$chr),]

# CHROMOMAP just cna-------------------------------------------------------------
# make annotation file
annotation <- hg19predictors[,c(1:4,7)]
annotation$bin[annotation$bin==1157] <- "1157 g-/l+"
annotation$bin[annotation$bin==1468] <- "1468 -"
annotation$bin[annotation$bin==1455] <- "1455 -"
annotation$bin[annotation$bin==1848] <- "1848 -"
annotation$bin[annotation$bin==1470] <- "1470 +"
annotation$bin[annotation$bin==1131] <- "1131 +"
annotation$bin[annotation$bin==406] <- "406 +"
annotation$bin[annotation$bin==1449] <- "1449 +"
annotation$bin[annotation$bin==1090] <- "1090 +"
annotation$bin[annotation$bin==1491] <- "1491 g-/l+"
annotation$bin[annotation$bin==1139] <- "1139 g+/l-"
annotation$bin[annotation$bin==418] <- "418 +"
annotation$bin[annotation$bin==2326] <- "2326 g+/l+"
annotation$bin[annotation$bin==1487] <- "1487 +"
annotation$bin[annotation$bin==439] <- "439 -"
annotation$bin[annotation$bin==1851] <- "1851 +"
annotation$bin[annotation$bin==1522] <- "1522 +"
annotation$bin[annotation$bin==2385] <- "2385 g+/l+"
annotation$bin[annotation$bin==2371] <- "2371 -"
annotation$bin[annotation$bin==1374] <- "1374 +"
annotation$bin[annotation$bin==1193] <- "1193 +"
annotation$bin[annotation$bin==606] <- "606 +"
annotation$bin[annotation$bin==2083] <- "2083 +"
annotation$bin[annotation$bin==1337] <- "1337 +"
annotation$bin[annotation$bin==898] <- "898 -"
annotation$bin[annotation$bin==1074] <- "1074 +"
annotation$bin[annotation$bin==2374] <- "2374 g-/l-"
annotation$chr <- as.character((annotation$chr))
annotation$cna <- as.character((annotation$cna))
chromoMap(list(hg19), list(annotation),
          #cna
          data_based_color_map = T, 
          data_type = "categorical",
          data_colors = list(c("#CC0033","#330099","#009999")),
          #labels
          #labels = T, label_font = 0, label_angle = -70,
          #features
          canvas_width = 600, canvas_height = 500,
          chr_color = c("#CCCCCC"),
          #chr_width = 13, chr_length = 6, ch_gap = -6, 
          chr_width = 12, chr_length = 5.5, ch_gap = 12, 
          y_chr_scale = 12, top_margin = 0,
          legend = T, lg_x = 50, lg_y = 150)


legend <- as_ggplot(cowplot::get_legend(ggplot(data.frame(color=c("both","gain","loss"),x=c(1,2,3), y=c(1,2,3))) +
                                                geom_bar(aes(x=x, y=y, fill=factor(color)), stat="identity") + 
                                                scale_fill_manual(values = c("#009999","#CC0033","#330099")) +
                                                labs(colour = NULL, fill = NULL) +
                                                guides(fill = guide_legend(direction = 'horizontal', label.position = "left", label.hjust = 1),
                                                theme(plot.margin = unit(c(t=-100,r=0,b=-100,l=0), "cm"),
                                                      legend.position = "bottom",legend.direction="horizontal",
                                                      legend.title = element_text(size=24, colour='black',face='bold'),
                                                      legend.margin = margin(grid::unit(c(t=-100,r=0,b=-100,l=0),"cm")),
                                                      legend.text = element_text(size=24, colour='black'),
                                                      legend.key.height = grid::unit(0.8,"cm"),
                                                      legend.key.width = grid::unit(1.4,"cm"))  ))) #or 0.8 for square


annotation <- hg19predictors[,c(1:4,6)]
annotation$chr <- as.character((annotation$chr))
annotation$gate <- "0-30%"
annotation$gate[which(annotation$pcSubclonal>0.30)] <- "30-50%"
annotation$gate[which(annotation$pcSubclonal>0.5)] <- "50-70%"
annotation$gate[which(annotation$pcSubclonal>0.7)] <- "70-90%"
annotation$gate[which(annotation$pcSubclonal>0.9)] <- "90-100%"
#annotation$gate <- as.factor(as.character(annotation$gate))
#levels(annotation$gate) <- c("0-30%","90-100%","70-90%","30-50%","50-70%")

annotation <- annotation[-5]

chromoMap(list(hg19), list(annotation),
          #cna
          data_based_color_map = T, 
          data_type = "categorical",
          data_colors = list(c("#00CC66","#6666FF","#660066","#CC66FF","#006600")),
          #labels
          labels = T, label_font = 9, label_angle = -70,
          #features
          canvas_width = 600, canvas_height = 1000,
          chr_color = c("#CCCCCC"),
          chr_width = 12, 
          chr_length = 9, 
          ch_gap = -10, 
          y_chr_scale = -30, 
          top_margin = -20,
          legend = T, 
          lg_x = 220, lg_y = 600
          )

# tile of cna per 27 per patient in training-------------------------------------
x <- car.raw[-c(1:3)]
rownames(x) <- 1:nrow(x)
x <- x[which(rownames(x) %in% hg19predictors$bin), ]
x$bin <- rownames(x)
x$chr <- hg19predictors[match(x$bin, hg19predictors$bin),2]
levels(y$bin) <- c("406","418","439","606","898","1074","1085","1090","1131","1139","1157","1193","1337",
                   "1374","1418","1449","1455","1468","1491","1522","1848","1851","2083","2326","2371","2374","2385" )
y$bin <- as.factor(y$bin)
y$cna <- as.factor(y$cna)

y <- gather(x, key="Sample", value="cna", -bin, -chr)
#y <- y[which(y$cna!=2),]

plot <- ggplot(y, aes(x=bin, y=Sample, fill = as.factor(cna))) +
  geom_tile(color = "black", lwd = 0.1, linetype = 1) +
  scale_fill_manual(values=c("#330099","white","#CC0033")) +
  facet_grid(. ~ chr, scales = "free", space = "free") +
  theme(
    axis.text.x = element_text(size=24, angle=90, vjust = 0.5),
    axis.text.y = element_blank(),
    axis.title = element_text(size=24),strip.text = element_text(size=24),
    panel.spacing = unit(0.1, "lines"),
    legend.position = "none") +
  guides()

jpeg('tempfig.jpeg', width = 1000, height = 1000)
plot
dev.off()

