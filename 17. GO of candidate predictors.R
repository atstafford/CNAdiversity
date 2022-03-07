#NORMAL-------------------------------------------------------------------------
#clusters <- hg19predictors$cluster
#bins <- sig.hclust$bin[which(sig.hclust$cluster %in% clusters)]

# Pull bins that are individually predictive
GO <- rbind(uniReg.out.list[[2]][which(uniReg.out.list[[2]]$sig=='sig'),],
            uniReg.out.list[[3]][which(uniReg.out.list[[3]]$sig=='sig'),])
#GO <- GO[which(GO$bin %in% bins),]

# Gene anno function needs to have cols: bin | chr | start | stop
GO <- merge(GO, car.info$start.stop)
GO <- GO[, c('bin','chr','start','stop','CNA','CNA.freq','coeff','pval')]

#GOlist <- read.delim("~/Documents/CNA/Data/GO_ENTREZ_HGNC_hg38.1.txt")
GOlist <- read.delim("~/Documents/CNA/Data/hg38.1_all_gene_GO_annotations.txt")
GOlist <- GOlist[c(8,11,5,3,4,1)]
GOlist <- GOlist[!duplicated(GOlist),]

gene.anno <- gene_anno(GO, GOlist)
gene.anno$start <- car.info$start.stop[match(gene.anno$bin, car.info$start.stop$bin), 3]
gene.anno$stop <- car.info$start.stop[match(gene.anno$bin, car.info$start.stop$bin), 4]
gene.anno <- gene.anno[c(1,3,7,8,2,6,4,5)]

# Add back in the gain/loss status to the data. 
gene.anno$CNA <- GO[match(gene.anno$bin, GO$bin),5]
gene.anno$CNA.freq <- GO[match(gene.anno$bin, GO$bin),6]

GOgain <- gene.anno$entrezgene[which(gene.anno$CNA=="gain")]
GOgain <- GOgain[!duplicated(GOgain)]
length(GOgain)

GOloss <- gene.anno$entrezgene[which(gene.anno$CNA=="loss")]
GOloss <- GOloss[!duplicated(GOloss)]
length(GOloss)

# KEGG and GO analysis (based on entrezgene)
#kegganno <- kegga(list(Up=GOgain,Down=GOloss))
#topKEGG(kegganno)
GOanno <- goana(list(Up=GOgain,Down=GOloss))
#topgo <- topGO(GOanno, number = 137)

GOanno$p.Up.adj = p.adjust(GOanno$P.Up, method = "fdr")
GOanno$p.Down.adj = p.adjust(GOanno$P.Down, method = "fdr")
GOanno$p.abs <- apply(GOanno[,c(8:9)], 1 , min)
GOanno$adjsig <- 'insig'
GOanno$adjsig[which(GOanno$p.abs<=0.05)] <- 'sig'
GOanno <- GOanno[which(GOanno$adjsig=="sig"),]

x <- data.frame(ID = rownames(GOanno[which(GOanno$p.Up.adj <= 0.05),]))
simMatrix_up <- calculateSimMatrix(x$ID, orgdb="org.Hs.eg.db", ont = c("BP","MF","CC"), method="Rel")
reducedTerms_up <- reduceSimMatrix(simMatrix_up, threshold=0.6, orgdb="org.Hs.eg.db")
jpeg('tempfig.jpeg', width = 1000, height = 1000)
treemapPlot(reducedTerms_up)
dev.off()

x <- data.frame(ID = rownames(GOanno[which(GOanno$p.Down.adj <= 0.05),]))
simMatrix_down <- calculateSimMatrix(x$ID, orgdb="org.Hs.eg.db", ont = c("BP","MF","CC"), method="Rel")
reducedTerms_down <- reduceSimMatrix(simMatrix_down, threshold=0.6, orgdb="org.Hs.eg.db")
jpeg('tempfig.jpeg', width = 1000, height = 1000)
treemapPlot(reducedTerms_down)
dev.off()

x <- data.frame(ID = rownames(GOanno[which(GOanno$p.abs <= 0.05),]))
simMatrix_all <- calculateSimMatrix(x$ID, orgdb="org.Hs.eg.db", ont = c("BP","MF","CC"), method="Rel")
scores <- data.frame(ID = rownames(GOanno[which(GOanno$p.abs <= 0.05),]),  scores = GOanno$p.abs[which(GOanno$p.abs <= 0.05)])
scores <- setNames(-log10(scores$scores), scores$ID)
reducedTerms_all <- reduceSimMatrix(simMatrix_all, threshold=0.6, orgdb="org.Hs.eg.db", scores = scores)
jpeg('tempfig.jpeg', width = 1000, height = 1000)
treemapPlot(reducedTerms_all)
dev.off()

# save
saveRDS(gene.anno, "~/Documents/CNA/Data/gene.anno.rds")
saveRDS(GOanno, "~/Documents/CNA/Data/GOanno.rds")
saveRDS(reducedTerms_all, "~/Documents/CNA/Data/reducedTerms_all.rds")

#bootstrapping 90% of 616-------------------------------------------------------
GO <- rbind(uniReg.out.list[[2]][which(uniReg.out.list[[2]]$sig=='sig'),],
            uniReg.out.list[[3]][which(uniReg.out.list[[3]]$sig=='sig'),])
GOanno.l <- list()
redT.l <- list()
redTscore.l <- list()

foreach(i = 1:100) %do% {
  x <- sample_n(GO, 585)
  x <- merge(x, car.info$start.stop)
  x <- x[, c('bin','chr','start','stop','CNA','CNA.freq','coeff','pval')]
  
  gene.anno <- gene_anno(x, GOlist)
  gene.anno$start <- car.info$start.stop[match(gene.anno$bin, car.info$start.stop$bin), 3]
  gene.anno$stop <- car.info$start.stop[match(gene.anno$bin, car.info$start.stop$bin), 4]
  gene.anno <- gene.anno[c(1,3,7,8,2,6,4,5)]
  
  gene.anno$CNA <- x[match(gene.anno$bin, x$bin),5]
  gene.anno$CNA.freq <- x[match(gene.anno$bin, x$bin),6]
  
  GOgain <- gene.anno$entrezgene[which(gene.anno$CNA=="gain")]
  GOgain <- GOgain[!duplicated(GOgain)]
  GOloss <- gene.anno$entrezgene[which(gene.anno$CNA=="loss")]
  GOloss <- GOloss[!duplicated(GOloss)]
  
  GOanno <- goana(list(Up=GOgain,Down=GOloss))
  GOanno$p.Up.adj = p.adjust(GOanno$P.Up, method = "fdr")
  GOanno$p.Down.adj = p.adjust(GOanno$P.Down, method = "fdr")
  GOanno$p.abs <- apply(GOanno[,c(8:9)], 1 , min)
  GOanno$adjsig <- 'insig'
  GOanno$adjsig[which(GOanno$p.abs<=0.05)] <- 'sig'
  GOanno <- GOanno[which(GOanno$adjsig=="sig"),]
  
  if (nrow(GOanno)==0) {
    next
  }
  
  x <- data.frame(ID = rownames(GOanno[which(GOanno$p.abs <= 0.05),]))
  simMatrix_all <- calculateSimMatrix(x$ID, orgdb="org.Hs.eg.db", ont = c("BP","MF","CC"), method="Rel")
  scores <- data.frame(ID = rownames(GOanno[which(GOanno$p.abs <= 0.05),]),  scores = GOanno$p.abs[which(GOanno$p.abs <= 0.05)])
  scores <- setNames(-log10(scores$scores), scores$ID)
  reducedTerms_all <- reduceSimMatrix(simMatrix_all, threshold=0.6, orgdb="org.Hs.eg.db", scores = scores)
  
  GOanno.l[[i]] <- GOanno
  redT.l[[i]] <- reducedTerms_all
  redTscore.l[[i]] <- aggregate(reducedTerms_all$score, by=list(Category=reducedTerms_all$parentTerm), FUN=sum)
}

saveRDS(redTscore.l, "~/Documents/CNA/Data/redTscore.l2.rds")

x <- do.call("rbind", redTscore.l)
unique(x$Category)


terms <- c("DNA packaging", "nucleosome assembly","DNA replication-dependent nucleosome organization","chromatin assembly",
           #"telomere organization","telomere capping","positive regulation of telomere maintenance via telomerase",
           #"regulation of gene expression, epigenetic","negative regulation of gene expression, epigenetic",
           "chromatin organization involved in negative regulation of transcription","chromatin organization involved in regulation of transcription")
df <- data.frame(i=1:100, chromatin=0)
for (i in 1:length(redTscore.l)) {
  wd <- redTscore.l[[i]]
  wd <- slice_max(wd, order_by = wd$x, n=1)
  if (wd$Category[1] %in% terms) {
    df$chromatin[i] <- 1
  }
}

length(df$chromatin[which(df$chromatin==1)])/nrow(df)

#bootstrapping 616 random bins--------------------------------------------------
GO <- rbind(uniReg.out.list[[2]])
GO <- merge(GO, car.info$start.stop)
GO <- GO[, c('bin','chr','start','stop','CNA','CNA.freq','coeff','pval')]
gene.anno <- gene_anno(GO, GOlist)
gene.anno2 <- gene.anno
gene.anno$start <- car.info$start.stop[match(gene.anno$bin, car.info$start.stop$bin), 3]
gene.anno$stop <- car.info$start.stop[match(gene.anno$bin, car.info$start.stop$bin), 4]
gene.anno <- gene.anno[c(1,3,7,8,2,6,4,5)]
gene.anno$CNA <- GO[match(gene.anno$bin, GO$bin),5]
gene.anno$CNA.freq <- GO[match(gene.anno$bin, GO$bin),6]

GOanno.l <- list()
redT.l <- list()
redTscore.l <- list()

foreach(i = 1:1000) %do% {
  print(i)
  sample <- sample(1:2694, 616)
  x <- gene.anno[which(gene.anno$bin %in% sample),]
  
  GOgain <- gene.anno$entrezgene
  GOgain <- GOgain[!duplicated(GOgain)]
  
  GOanno <- goana(list(Up=GOgain,Down=GOloss))
  GOanno$p.Up.adj = p.adjust(GOanno$P.Up, method = "fdr")
  GOanno$p.Down.adj = p.adjust(GOanno$P.Down, method = "fdr")
  GOanno$p.abs <- apply(GOanno[,c(8:9)], 1 , min)
  GOanno$adjsig <- 'insig'
  GOanno$adjsig[which(GOanno$p.abs<=0.05)] <- 'sig'
  GOanno <- GOanno[which(GOanno$adjsig=="sig"),]
  
  if (nrow(GOanno)==0) {
    next
  }
  
  x <- data.frame(ID = rownames(GOanno[which(GOanno$p.abs <= 0.05),]))
  simMatrix_all <- calculateSimMatrix(x$ID, orgdb="org.Hs.eg.db", ont = c("BP","MF","CC"), method="Rel")
  scores <- data.frame(ID = rownames(GOanno[which(GOanno$p.abs <= 0.05),]),  scores = GOanno$p.abs[which(GOanno$p.abs <= 0.05)])
  scores <- setNames(-log10(scores$scores), scores$ID)
  reducedTerms_all <- reduceSimMatrix(simMatrix_all, threshold=0.6, orgdb="org.Hs.eg.db", scores = scores)
  
  GOanno.l[[i]] <- GOanno
  redT.l[[i]] <- reducedTerms_all
  redTscore.l[[i]] <- aggregate(reducedTerms_all$score, by=list(Category=reducedTerms_all$parentTerm), FUN=sum)
}

saveRDS(redTscore.l, "~/Documents/CNA/Data/redTscore.l_non616.rds")

x <- do.call("rbind", redTscore.l)
unique(x$Category)


terms <- c("DNA packaging", "nucleosome assembly","DNA replication-dependent nucleosome organization","chromatin assembly",
           #"telomere organization","telomere capping","positive regulation of telomere maintenance via telomerase",
           #"regulation of gene expression, epigenetic","negative regulation of gene expression, epigenetic",
           "chromatin organization involved in negative regulation of transcription","chromatin organization involved in regulation of transcription")
df <- data.frame(i=1:164, chromatin=0)
for (i in 1:length(redTscore.l)) {
  wd <- redTscore.l[[i]]
  wd <- slice_max(wd, order_by = wd$x, n=5)
  wd <- wd$Category %in% terms
  if (TRUE %in% wd) {
    df$chromatin[i] <- 1
  }
}

length(df$chromatin[which(df$chromatin==1)])/nrow(df)

#other--------------------------------------------------------------------------
pathview(gene.data = GOgain, 
         species = "hsa", 
         pathway.id = "hsa05203",
         kegg.dir = tempdir())


