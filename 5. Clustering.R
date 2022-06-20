# GENERATE SIG.HCLUST.rds: ASSIGN BIN TO EACH CLUSTER 

# Use the hiclust rds to generate dataframe showing cluster for each bin

# Extract significant clusters into a dataframe
load("~/Documents/CNA/HPC/hiclust_diploid.aneu.rda")
hiclust <- hiclust_diploid.aneu
sig.hclust <- pvpick(hiclust, alpha=0.95, pv="au", max.only = T)

# For each cluster extract the corresponding bins
list <- list()
for ( i in 1:length(sig.hclust$clusters) ) { 
  list[[i]] <- data.frame(bin = unlist(sig.hclust$clusters[[i]]), cluster = i)
}

# Create a dataframe with columns: bin | cluster
sig.hclust <- do.call('rbind', list)
sig.hclust$bin <- sub('.', '', sig.hclust$bin)

# Convert dataframe to numeric
sig.hclust <- data.frame(apply(sig.hclust, 2, function(x) as.numeric(as.character(x))))

# Add chromosome data based on bin number
sig.hclust$chr <- car.info$start.stop[match(sig.hclust$bin, car.info$start.stop$bin), 2]

# save
saveRDS(sig.hclust, "~/Documents/CNA/Data/sig.hclust.rds")



