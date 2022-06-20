#NORMAL-------------------------------------------------------------------------
beta <- beta_52backrem
summary(beta_52backrem, type = "pearson")
x <- rownames(data.frame(summary(beta)$coefficients))[-1]
hg19predictors <- data.frame(bin=as.numeric(str_extract_all(x, "[0-9]+")), cna=str_sub(x,-4,-1))
hg19predictors <- merge(unique(hg19predictors), car.info$start.stop)
hg19predictors$cluster <- candidate.bins[match(hg19predictors$bin, candidate.bins$bin), 11] 
hg19predictors$pcSubclonal <- candidate.bins[match(hg19predictors$bin, candidate.bins$bin), 10] 
both <- hg19predictors$bin[which(duplicated(hg19predictors$bin)==T)]
hg19predictors$cna[which(hg19predictors$bin %in% hg19predictors$bin[which(duplicated(hg19predictors$bin)==T)] )] <- 'both'

hg19predictors <- hg19predictors[c(1,3:7,2)]
clusters <- unique(hg19predictors$cluster)
bins <- sig.hclust$bin[which(sig.hclust$cluster %in% clusters)]
#bins <- hg19predictors$bin[!duplicated(hg19predictors$bin)]

# Pull bins that are individually predictive
GO <- rbind(uniReg.out.list[[2]][which(uniReg.out.list[[2]]$sig=='sig'),],
            uniReg.out.list[[3]][which(uniReg.out.list[[3]]$sig=='sig'),])
GO <- GO[which(GO$bin %in% bins),]

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
gene.anno$direction <- GO[match(gene.anno$bin, GO$bin), 7]
gene.anno <- gene.anno[c(1,3,7,8,2,6,4,5,9)]

# Add back in the gain/loss status to the data. 
gene.anno$CNA <- GO[match(gene.anno$bin, GO$bin),5]
gene.anno$CNA.freq <- GO[match(gene.anno$bin, GO$bin),6]

GOgainUp <- gene.anno[which(gene.anno$CNA=="gain" & gene.anno$direction > 0),]
GOgainUp <- GOgainUp[!duplicated(GOgainUp),]
nrow(GOgainUp)

GOgainDown <- gene.anno[which(gene.anno$CNA=="gain" & gene.anno$direction < 0),]
GOgainDown <- GOgainDown[!duplicated(GOgainDown),]
nrow(GOgainDown)

GOlossUp <- gene.anno[which(gene.anno$CNA=="loss" & gene.anno$direction > 0),]
GOlossUp <- GOlossUp[!duplicated(GOlossUp),]
nrow(GOlossUp)

GOlossDown <- gene.anno[which(gene.anno$CNA=="loss" & gene.anno$direction < 0),]
GOlossDown <- GOlossDown[!duplicated(GOlossDown),]
nrow(GOlossDown)


# SET UP HALLMARK ======
BiocManager::install("msigdbr")
library(msigdbr)
m_df = msigdbr(species = "Homo sapiens", category = "H")
m_df$gs_name <- gsub('HALLMARK_(\\S+)','\\1',m_df$gs_name)
m_df <- m_df[,c(3:4)]
isc <- read.table('~/Documents/CNA/Data/merlossuarez_2011_ISC_signature.txt')[,1] # Include ISC signature
geneinfo <- fread('~/Documents/CNA/Data/complete_gene_info.txt',data.table=F)
isc_entrez <- geneinfo[which(geneinfo$Symbol %in% isc),]
isc_entrez$gs_name <- "ISC"
isc_entrez <- isc_entrez[,c(12,9)]
colnames(isc_entrez) <- colnames(m_df)
wnt_entrez <- data.frame("WNT", read.table('~/Documents/CNA/Data/wnt_signalling_entrez.txt')[,1]) # Include extra WNT_signalling
colnames(wnt_entrez) <- colnames(m_df)

m_df <- rbind(m_df, isc_entrez, wnt_entrez)

# SEARCH ========
hallmarkGenes <- m_df
chromatinGenes <- read.table('~/Documents/CNA/Data/chromatinGenes', header = T, sep = "\t")
DNAmodGenes <- read.table('~/Documents/CNA/Data/DNAmodGenes', header = T, sep = "\t")
histoneGenes <- read.table('~/Documents/CNA/Data/histoneGenes', header = T, sep = "\t")

x <- GOlist$HGNC.symbol[which(GOlist$NCBI.gene..formerly.Entrezgene..ID %in% GOgainUp$entrezgene)]
x[x %in% chromatinGenes$Symbol]
x[x %in% DNAmodGenes$Symbol]
x[x %in% histoneGenes$Symbol]
y <- data.frame(x[x %in% hallmarkGenes$gene_symbol])
y <- y[order(y),]

