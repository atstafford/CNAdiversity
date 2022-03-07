# Function to pull basic info from rawdata eg number of samples/bins/patients
# Requires: chr | start | stop | sample1CN | sample2CN
# Requires sample name in colname, in structure: patientID.sampleNumber
PullDataInfo <- function(rawdata) {
  
  # Dataframe identifying start and stop codon, and chromosome for each bin
  start.stop <- rawdata[,c(1:3)]
  start.stop$bin <- 1:nrow(start.stop)
  start.stop <- start.stop[, c(4, 1, 2, 3)]
  
  # Vector holding sample ID
  sampleIDs <- colnames(rawdata)[-c(1:3)] 
  
  # Vector holding patient identifiers, ie the string before the '.' in sample IDs
  patientIDs <- unique(sub("\\..*", "", colnames(rawdata)))[-c(1:3)] 
  
  # Number of samples
  noSamples <- length(sampleIDs)
  
  # Number of patients
  noPatients <- length(patientIDs)
  
  # list of number of samples per patient
  sampPerPatient <- list()
  for ( i in 1:noPatients ) {
    wd <- rawdata[, which(sub("\\..*", "", colnames(rawdata))==patientIDs[i])]
    sampPerPatient[[i]] <- ncol(wd)
  }
  
  # Number of bins
  noBins <- length(start.stop$bin)
  
  # Visualisation may require a vector to identify of chromosome ends and chromosome midpoints 
  chr.ends <- cumsum(table(start.stop$chr))
  list <- list()
  l <- 1
  for ( i in 1:length(chr.ends) ) {  
    if ( i == 1 ) { 
      list[[l]] <- chr.ends[i]/2 
      l <- l+1
    }
    else { 
      list[[l]] <- chr.ends[i-1] + ((chr.ends[i]-chr.ends[i-1])/2)
      l <- l+1
    }
  }
  chr.mid <- unlist(list)
  chr.ends <- data.frame(start=c(0,chr.ends[-22]), end=(chr.ends), col=c('G','W'))
  
  # Return data
  newData <- list('start.stop'=start.stop, 'sampleIDs'=sampleIDs, 'patientIDs'=patientIDs, 'noSamples'=noSamples, 'noPatients'=noPatients,
                  'sampPerPatient'=sampPerPatient, 'noBins'=noBins,'chr.mid'=chr.mid,  'chr.end'=chr.ends)
  return(newData)
}

# Function to create 5 new data structures relating to clonality
# Requires: chr | start | stop | sample1CN | sample2CN
# dataInfo must be previously generated from rawdata using PullDataInfo
PullDataClonality <- function(rawdata, dataInfo) {
  
  # A dataframe with a bin per row and a patient per column, with values indicating clonality. 
  # 0=notCNA, 1=subclonalCNA, 2=clonalCNA
  clonal.data <- rawdata[,-c(1:3)]
  l <- 1
  clonal <- list()
  for ( k in 1:nrow(clonal.data) ) {
    for ( i in 1:length(dataInfo$patientIDs) ) {
      wd <- clonal.data[k, which(sub("\\..*", "", colnames(clonal.data))==dataInfo$patientIDs[i])]
      wd <- wd[, is.na(wd)!=TRUE]
      
      if ( 1 %in% wd | 3 %in% wd ) { #if one of the samples has a mutation, proceed
        if ( length(unique(t(wd)))==1 ) { #if all the same then clonal
          clonal[[l]] <- 2
          l <- l + 1
        }
        else { #different = subclonal
          clonal[[l]] <- 1
          l <- l + 1
        }
      }
      else if ( length(wd)==0 ) { #all MR for that bit are NA
        clonal[[l]] <- NA
        l <- l + 1
      }
      else { #neither sample has a mutation
        clonal[[l]] <- 0 
        l <- l + 1
      }
    }
  }
  
  clonal.data <- data.frame(t(matrix(unlist(clonal), ncol=dataInfo$noBins)))
  colnames(clonal.data) <- dataInfo$patientIDs
  clonal.data[] <- lapply(clonal.data, factor, levels=unique(unlist(clonal.data)))
  
  
  # Dataframe detailing the counts of gains/losses and whether they are subclonal or clonal
  CNA.clo.counts <- data.frame(bin = dataInfo$start.stop$bin, chr1 = NA, chr2 = NA,
                               clonal.aneu = NA, subclonal.aneu = NA, gain = NA, loss = NA,
                               clonal.gain = NA, clonal.loss = NA, clonal.noCNA = NA, 
                               subclonal.gain = NA, subclonal.loss = NA)
  
  CNA.clo.counts$chr1 <- as.numeric(dataInfo$start.stop[match(CNA.clo.counts$bin, dataInfo$start.stop$bin), 2]) 
  CNA.clo.counts$chr2 <- as.numeric(dataInfo$start.stop[match(CNA.clo.counts$bin, dataInfo$start.stop$bin), 2]) 
  CNA.clo.counts$chr2 <- factor(CNA.clo.counts$chr2,levels = rev(seq(1:22))) 
  
  data <- rawdata[,-c(1:3)]
  for ( k in 1:nrow(data) ) { # for a bin
    clonal.all <- clonal.gain <- clonal.loss <- clonal.noCNA <- subclonal.all <- subclonal.gain <- subclonal.loss <- subclonal.noCNA <- 0
    
    for ( i in 1:dataInfo$noPatients ) { # for a patient
      
      # If its NA
      if ( is.na(clonal.data[k,i]) == TRUE ) {
        next
      }
      
      # If its clonal
      if ( clonal.data[k,i] == 2 ) {
        wd <- data[k, which(sub("\\..*", "", colnames(data))==dataInfo$patientIDs[i])]
        wd <- wd[, is.na(wd)!=TRUE]
        
        if ( 3 %in% wd ) { # if both are gains
          clonal.gain <- clonal.gain + 1
        }
        else if ( 1 %in% wd ) { # if both are losses
          clonal.loss <- clonal.loss + 1
        }
      }
      
      # if its subclonal
      else if ( clonal.data[k,i] == 1 ) {
        wd <- data[k, which(sub("\\..*", "", colnames(data))==dataInfo$patientIDs[i])]
        wd <- wd[, is.na(wd)!=TRUE]
        
        if ( 3 %in% wd ) { # if one is a gain
          subclonal.gain <- subclonal.gain + 1
        }
        if ( 1 %in% wd ) { # if one is a loss
          subclonal.loss <- subclonal.loss + 1
        }
      }
      
      # if its no CNA
      else if ( clonal.data[k,i] == 0 ) {
        clonal.noCNA <- clonal.noCNA + 1
      }
    }
    
    CNA.clo.counts$clonal.gain[k] <- clonal.gain
    CNA.clo.counts$clonal.loss[k] <- clonal.loss
    CNA.clo.counts$clonal.noCNA[k] <- clonal.noCNA
    CNA.clo.counts$subclonal.gain[k] <- subclonal.gain
    CNA.clo.counts$subclonal.loss[k] <- subclonal.loss
  }
  
  CNA.clo.counts$clonal.aneu <- CNA.clo.counts$clonal.gain + CNA.clo.counts$clonal.loss
  CNA.clo.counts$subclonal.aneu <- CNA.clo.counts$subclonal.gain + CNA.clo.counts$subclonal.loss
  CNA.clo.counts$gain <- CNA.clo.counts$clonal.gain + CNA.clo.counts$subclonal.gain
  CNA.clo.counts$loss <- CNA.clo.counts$clonal.loss + CNA.clo.counts$subclonal.loss
  
  # dataframe showing: bin | countGain/NoPatient | countLoss/NoPatient | 
  CloFreq <- cbind(bin=CNA.clo.counts$bin, gain=CNA.clo.counts$gain/dataInfo$noPatients, loss=CNA.clo.counts$loss/dataInfo$noPatients)
  
  # dataframe showing what percent of gain and loss are subclonal. On a patient basis: bin | gain | loss
  pcSubclonal <- data.frame(bin=1:dataInfo$noBins, gain=CNA.clo.counts$subclonal.gain / CNA.clo.counts$gain, loss=CNA.clo.counts$subclonal.loss / CNA.clo.counts$loss)
  
  # Count of noCNA, subclonalCNA, and clonalCNA by patient
  patientClo <- as.data.frame(t(sapply(clonal.data, table))) 
  patientClo <- patientClo[, c('0','1','2')]
  colnames(patientClo) <- c('noCNA','subclonal','clonal')
  patientClo$CNA <- patientClo$subclonal + patientClo$clonal
  patientClo$patient <- rownames(patientClo)
  
  # Return data
  newData <- list('clonal.data'=clonal.data, 'CNA.clo.counts'=CNA.clo.counts, 
                  'CloFreq'=CloFreq, 'pcSubclonal'=pcSubclonal, 'patientClo'=patientClo)
  return(newData)
}

# Function to create 4 new data structures relating to diversity
# Requires: chr | start | stop | sample1CN | sample2CN
# There must only be 2 multiregion samples per patient
# dataInfo must be previously generated from rawdata using PullDataInfo
PullDataDiversity <- function(rawdata, dataInfo) {
  
  # Cannot use averaged raw CN across the MR samples in case of variable number of samples per patient.
  # Will therefore can use a PIC.frac for each patient. 
  # This is a continuous, not binary, measure of clonality
  # Define the maximum pic score depending on the number of samples (up to 13)
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
  
  # A dataframe with a bin per row and a patient per column, with values indicating pic score
  # And a dataframe showing pic.frac per patient
  pic <- list()
  pic.frac <- list()
  for ( i in 1:length(dataInfo$patientIDs)) {
    # Set working data as the cols holding the samples
    wd <- rawdata[, which(sub("\\..*", "", colnames(rawdata))==dataInfo$patientIDs[i])]
    
    # Record the number of samples
    upto <- ncol(wd)
    
    # Use PIC function on wd
    pic[[i]] <- PIC(wd, upto, c(1:upto))
    pic[[i]] <- na.omit(pic[[i]])
    
    # Define the max possible diversity given the number of sample, as maxPIC*number of bins
    max.ITH <- max.pics[[upto]] * length(pic[[i]])
    
    # Store as a dataframe
    pic.frac[[i]] <- data.frame(pic.frac = sum(pic[[i]], na.rm = TRUE)/max.ITH)
  }
  
  pic.data <- as.data.frame(do.call('cbind',pic))
  colnames(pic.data) <- dataInfo$patientIDs
  pic.frac <- do.call('rbind',pic.frac)
  
  # Calulate average PIC per bin as a measure of how diverse that bin tends to be across patients
  ave.pic <- dataInfo$start.stop[,c(1:2)]
  ave.pic$avePic <- rowSums(pic.data)/(ncol(pic.data))
  
  # A dataframe of the propotion of genome gained/lost per sample, alongside ith
  pga <- data.frame(t(rawdata[,-c(1:3)]), check.names = FALSE)
  pga <- data.frame(prop.gain=apply(pga,1,function(x) sum(x == 3, na.rm = TRUE)/ncol(pga)),
                    prop.loss=apply(pga,1,function(x) sum(x == 1, na.rm = TRUE)/ncol(pga)))
  pga$prop.aneu <- pga$prop.gain + pga$prop.loss
  pga <- cbind(as.data.frame(lapply(pic.frac, rep, dataInfo$sampPerPatient)),
               pga)
  
  # Return data
  newData <- list('pic.data'=pic.data, 'pic.frac'=pic.frac, 'ave.pic'=ave.pic, 'pga'=pga)
  return(newData)
}

# Function to generate matrices for dip vs aneu, dip vs gain, dip vs loss, gain vs loss, and dip vs gain vs loss
genMatrices <- function(rawdata) {
  
  # Create list to store correlation matrices
  matrices.list <- rep(list(data.frame((rawdata[,-c(1:3)]))),5)
  names(matrices.list) = c('diploid.aneu', 'diploid.gain', 'diploid.loss', 'loss.gain', 'loss.dip.gain')
  
  # Convert to numeric matrices
  matrices.list <- lapply(matrices.list, function(x) {
    y <- data.frame(apply(x, 2, as.numeric))
    rownames(y) <- rownames(x)
    y
  })
  
  # In the diploid/aueploid matrix make 0=diploid, 1=aneuploid.
  matrices.list$diploid.aneu[matrices.list$diploid.aneu == 1 | matrices.list$diploid.aneu == 3] <- 1
  matrices.list$diploid.aneu[matrices.list$diploid.aneu == 2] <- 0
  
  # In the diploid/gain matrix make 0=diploid, 1=gain.
  matrices.list$diploid.gain[matrices.list$diploid.gain == 1] <- NA
  matrices.list$diploid.gain[matrices.list$diploid.gain == 2] <- 0
  matrices.list$diploid.gain[matrices.list$diploid.gain == 3] <- 1
  
  # In the diploid/loss matrix make 0=diploid, 1=loss. 
  matrices.list$diploid.loss[matrices.list$diploid.loss == 1] <- 1
  matrices.list$diploid.loss[matrices.list$diploid.loss == 2] <- 0
  matrices.list$diploid.loss[matrices.list$diploid.loss == 3] <- NA
  
  # In the loss/gain matrix make 0=loss, 1=gain. 
  matrices.list$loss.gain[matrices.list$loss.gain == 1] <- 0
  matrices.list$loss.gain[matrices.list$loss.gain == 2] <- NA
  matrices.list$loss.gain[matrices.list$loss.gain == 3] <- 1
  
  # # In the loss/dip/gain matrix make 1=loss, 2=diploid, 3=gain. Requires no changes
  
  return(matrices.list)
}

# Function to calc pic score per patient across multiple samples
PIC <- function(data, sample.no, index) { 
  # based on table of unique patient nos
  # returns PIC per row (bin), calc across cols as given by index
  # PIC formula: 1- ((CN1/n)^2 + (CN2/n)^2 + (CN3/n)^2), where CN1 is no. of counts of copy number1
  PIC <- 1 - ((rowSums(data[,index]==1, na.rm = TRUE)/sample.no)^2 + 
                (rowSums(data[,index]==2, na.rm = TRUE)/sample.no)^2 + 
                (rowSums(data[,index]==3, na.rm = TRUE)/sample.no)^2)
  return(PIC)
}

# Function to collect pvalue from a regression
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

# Function to assess ploidy and recentre if average CN is 2.8 or more. skipcol may need to be only 3...
ploidyRecentre <- function (cn.list, skipcol, multi_cn) {
  # Requires list with dataframes per patient 
  # skipcol shows how many columns hold patient info and not cn
  for ( i in 1:length(cn.list) ) {
    wd <- cn.list[[i]][,-c(1:skipcol)]
    
    if ( multi_cn == TRUE ) { # if there is multregion sampling
      if ( mean(colMeans(wd, na.rm = TRUE)) >= 2.8 ) {
        cn.list[[i]] <- cbind(cn.list[[i]][,c(1:skipcol)] ,cn.list[[i]][,-c(1:skipcol)] - 1)
      }
      else {
        cn.list[[i]] <- cn.list[[i]]
      }
    }
    
    else if ( multi_cn == FALSE ){ # if only a single sample per patient
      if ( mean(wd, na.rm = TRUE) >= 2.8 ) {
        cn.list[[i]] <- cbind(cn.list[[i]][,c(1:skipcol)] ,cn.list[[i]][,-c(1:skipcol)] - 1)
      }
      else {
        cn.list[[i]] <- cn.list[[i]]
      }
    }
  }
  return(cn.list)
}

# Function to bin CN data to match training dataset (Cross et al)
alignBins <- function(bins, cn.list) {
  # bins needs to be a dataframe holding: bin | chr | start | stop, for the bins to align to
  # cn.list is the output from ploidyRecentre
  
  # Create dataframe for each patient and put into list
  cnBinned.list <- rep(list(bins), length(cn.list))
  
  # Convert to numeric
  cnBinned.list <- lapply(cnBinned.list, function (x) {
    x[] <- apply(x,2,as.numeric)
    x
  })
  
  # Add empty columns to hold output
  for ( i in 1:length(cnBinned.list) ) {
    cnBinned.list[[i]][,c(5:(1+ncol(cn.list[[i]])))] <- NA
    colnames(cnBinned.list[[i]])[-c(1:4)] <- colnames(cn.list[[i]])[-c(1:3)]
  }
  
  # Bin incoming dataset (cn.list)
  for ( p in 1:length(cnBinned.list) ) {
    for ( b in 1:nrow(cnBinned.list[[p]]) ) {
      chr <- cnBinned.list[[p]]$chr[b]
      start <- cnBinned.list[[p]]$start[b] 
      stop <- cnBinned.list[[p]]$stop[b]
      bin <- cnBinned.list[[p]]$bin[b]
      
      wd <- cn.list[[p]][cn.list[[p]]$chr == chr,]
      wd <- wd[order(wd$start),]
      
      for ( r in 1:nrow(wd) ) {
        
        ## if start is before the earliest row that that chromosome
        if ( start < wd$start[r] ) {
          
          # if stop is in row r of wd or else...
          if ( dplyr::between(stop, wd$start[r], wd$stop[r]) | nrow(wd)==r ) { 
            cnBinned.list[[p]][b,c(5:ncol(cnBinned.list[[p]]))] <- wd[r,c(4:ncol(wd))]
            break
          }
          
          # ...or else stop is beyond row r of wd
          else if ( stop > wd$stop[r] ) { 
            fraction <- list()
            f <- 1
            
            fraction[[f]] <- (wd$stop[r]-start)/(stop-start) # what fraction of the bin is the region in
            f <- f + 1
            
            for ( x in 1:(nrow(wd)-r) ) {
              # stop is within x extra wd rows
              if ( dplyr::between(stop, wd$start[r+x], wd$stop[r+x]) ) { 
                fraction[[f]] <- (wd$stop[r+x]-wd$start[r+x])/(stop-start)
                f <- f + 1
                break
              }
              
              # stop is within a later row of wd
              else if ( stop > wd$stop[r+x] ) { 
                fraction[[f]] <- (wd$stop[r+x]-wd$start[r+x])/(stop-start)
                f <- f + 1
              }
            }
            
            # majority of bin sits in row r+f (q) of wd
            q <- r + (which.max(fraction) - 1)
            cnBinned.list[[p]][b,c(5:ncol(cnBinned.list[[p]]))] <- wd[q,c(4:ncol(wd))]
            break
          }
          
          # stop is also before start of first row of wd
          else if ( stop < wd$start[r] ) {
            
            # if there this r is the first in the chromosome
            if ( r == 1 ) {
              cnBinned.list[[p]][b,c(5:ncol(cnBinned.list[[p]]))] <- NA
              break
            }
            
            # if bin lies between r and r-1 take average
            else if ( r != 1 ) {
              cnBinned.list[[p]][b,c(5:ncol(cnBinned.list[[p]]))] <- NA
              break
            }
          }
        }
        
        ## if start is in row r of wd
        else if ( dplyr::between(start, wd$start[r], wd$stop[r]) ) { 
          
          # if stop is in row r of wd or else...
          if ( dplyr::between(stop, wd$start[r], wd$stop[r]) ) { 
            cnBinned.list[[p]][b,c(5:ncol(cnBinned.list[[p]]))] <- wd[r,c(4:ncol(wd))]
            break
          }
          
          # ...or else stop is beyond row r of wd
          else if ( stop > wd$stop[r] ) { 
            
            # if there are no more rows of wd
            if ( r == nrow(wd) ) {
              cnBinned.list[[p]][b,c(5:ncol(cnBinned.list[[p]]))] <- wd[r,c(4:ncol(wd))]
              break
            }
            
            # if there are more rows of wd
            else if ( r < nrow(wd) ) {
              fraction <- list()
              f <- 1
              
              fraction[[f]] <- (wd$stop[r]-start)/(stop-start) # what fraction of the bin is the region in
              f <- f + 1
              
              for ( x in 1:(nrow(wd)-r) ) {
                # stop is within x extra wd rows
                if ( dplyr::between(stop, wd$start[r+x], wd$stop[r+x]) ) { 
                  fraction[[f]] <- (wd$stop[r+x]-wd$start[r+x])/(stop-start)
                  f <- f + 1
                  break
                }
                
                # stop is within a later row of wd
                else if ( stop > wd$stop[r+x] ) { 
                  fraction[[f]] <- (wd$stop[r+x]-wd$start[r+x])/(stop-start)
                  f <- f + 1
                }
              }
              
              # majority of bin sits in row r+f (q) of wd
              q <- r + (which.max(fraction) - 1)
              cnBinned.list[[p]][b,c(5:ncol(cnBinned.list[[p]]))] <- wd[q,c(4:ncol(wd))]
              break
            }
          }
        }
        
        ## if start is in a later row of wd
        else {
          next
        }
        
      }
    }
  }
  return(cnBinned.list)
}

# a new function to align bins
newAlignBins <- function(bins, cn.list) {
  # bins needs to be a dataframe holding: bin | chr | start | stop, for the bins to align to
  # cn.list is the output from ploidyRecentre with skipcol=3
  
  # Create dataframe for each patient and put into list
  cnBinned.list <- rep(list(bins), length(cn.list))
  
  # Convert to numeric
  cnBinned.list <- lapply(cnBinned.list, function (x) {
    x[] <- apply(x,2,as.numeric)
    x
  })
  
  # Add empty columns to hold output
  for ( i in 1:length(cnBinned.list) ) {
    sampleNo <- ncol(cn.list[[i]])-3 #
    cnBinned.list[[i]] <- do.call("cbind", list(cnBinned.list[[i]], rep(list(NA), sampleNo)))
    colnames(cnBinned.list[[i]])[-c(1:4)] <- colnames(cn.list[[i]])[-c(1:3)]
  }
  
  # Bin incoming dataset (cn.list)
  for ( p in 1:length(cnBinned.list) ) {
    # per bin
    for ( b in 1:nrow(cnBinned.list[[p]]) ) {
      chr <- cnBinned.list[[p]]$chr[b]
      start <- cnBinned.list[[p]]$start[b] 
      stop <- cnBinned.list[[p]]$stop[b]
      bin <- cnBinned.list[[p]]$bin[b]
      
      wd <- cn.list[[p]][cn.list[[p]]$chr == chr,]
      wd <- wd[order(wd$start),]
      
      for ( r in 1:nrow(wd) ) {
        
            ## IF START IS BEFORE ROW R 
            if ( start < wd$start[r] ) {
              
#1.3              # stop is also before start of first row of wd
                  if ( stop < wd$start[r] ) {
                    cnBinned.list[[p]][which(cnBinned.list[[p]]$bin==bin),-c(1:4)] <- "too early"
                    break
                  }
#1.1              # if stop is in row r of wd, or predictor doesnt reach next bin, or this is the last row of wd...
                  else if ( dplyr::between(stop, wd$start[r], wd$stop[r]) | stop < wd$start[r+1] | r == nrow(wd) ) { 
                    cnBinned.list[[p]][which(cnBinned.list[[p]]$bin==bin),-c(1:4)] <- wd[r,-c(1:3)]
                    break
                  }
                
#1.2              # ...or else stop is beyond row r of wd
                  else if ( stop > wd$stop[r] ) {
                    
                    fraction <- list()
                    f <- 1 
                    
                    # find number of extra rows required
                    for ( e in 0:(nrow(wd)-r) ) {
                      
                      # if stop is before row r+e of wd, or stop doesnt reach next bin, or this is the last row of wd...
                      if ( stop <= wd$stop[r+e] | stop < wd$start[r+e+1] | r+e == nrow(wd) ) {
                        fraction[[f]] <- ( min(wd$stop[r+e], stop)-wd$start[r+e])/(stop-start) # we know it always covers wd$start
                        break}
                      
                      else if (stop > wd$stop[r+e]) {
                        fraction[[f]] <- ( min(wd$stop[r+e], stop)-wd$start[r+e])/(stop-start) # take fraction and move onto next
                        f <- f + 1
                        next}
                    }
                    
                    # take weighted mean
                    cnBinned.list[[p]][which(cnBinned.list[[p]]$bin==bin),-c(1:4)] <- colSums( (data.frame(wd[c(r:(r+e)),-c(1:3)]) *unlist(fraction)), na.rm=T) / sum(unlist(fraction))
                      break
                    
                    }
    
            }
              
            ## IF START IS IN ROW R
            else if ( dplyr::between(start, wd$start[r], wd$stop[r]) ) { 
              
#2.1              # if stop is in row r of wd, or predictor doesnt reach next bin, or this is the last row of wd...
                  if ( dplyr::between(stop, wd$start[r], wd$stop[r]) | stop < wd$start[r+1] | r == nrow(wd) ) { 
                    cnBinned.list[[p]][which(cnBinned.list[[p]]$bin==bin),-c(1:4)] <- wd[r,-c(1:3)]
                    break
                  }
                  
#2.2              # ...or else stop is beyond row r of wd
                  else if ( stop > wd$stop[r] ) {
                    
                    fraction <- list()
                    f <- 1 
                    
                    # find number of extra rows required
                    for ( e in 0:(nrow(wd)-r) ) {
                      
                      # if stop is before row r+e of wd, or stop doesnt reach next bin, or this is the last row of wd...
                      if ( stop <= wd$stop[r+e] | stop < wd$start[r+e+1] | r+e == nrow(wd) ) {
                        fraction[[f]] <- ( min(wd$stop[r+e], stop) - max(wd$start[r+e], start))/(stop-start) # we dont know it always covers wd$start
                        break}
                      
                      else if (stop > wd$stop[r+e]) {
                        fraction[[f]] <- ( min(wd$stop[r+e], stop) - max(wd$start[r+e], start))/(stop-start) # take fraction and move onto next
                        f <- f + 1
                        next}
                    }
                    
                    # take weighted mean
                    cnBinned.list[[p]][which(cnBinned.list[[p]]$bin==bin),-c(1:4)] <- colSums( (data.frame(wd[c(r:(r+e)),-c(1:3)]) *unlist(fraction)), na.rm=T) / sum(unlist(fraction))
                    break
                    
                  }
            }
            
            ## IF START IS IN A LATER ROW
            else if ( start > wd$stop[r] ){
              next
            }
        
      }
    }
  }
  return(cnBinned.list)
}

          
          
          
'%!in%' <- function(x,y)!('%in%'(x,y))

# get chromosome lengths and centromere positions
getChrLength <- function(genome = "BSgenome.Hsapiens.UCSC.hg19"){
  g <- getBSgenome(genome, masked=FALSE)
  data.frame(chrom=1:24, length=seqlengths(g)[1:24])
}
.chrAsNum <- function(tbl){
  tbl$chrom <- gsub("chr", "", tbl$chrom)
  tbl$chrom[tbl$chrom=="X"] <- 23
  tbl$chrom[tbl$chrom=="Y"] <- 24
  tbl$chrom <- as.numeric(tbl$chrom)
  tbl[order(tbl$chrom),]
}
getCentromeres <- function( genome="hg19" ){
  mySession <- try(browserSession("UCSC"), silent=TRUE)
  # In case of failure, try another mirror
  if(inherits(mySession, "try-error"))
    mySession <- browserSession("UCSC",
                                url="http://genome-euro.ucsc.edu/cgi-bin/")
  genome(mySession) <- genome
  obj <- ucscTableQuery(mySession, table="gap")
  tbl <- getTable(obj)
  tbl <- tbl[tbl$type=="centromere", c("chrom", "chromStart", "chromEnd")]
  colnames(tbl)[2:3] <- c("centromerStart", "centromerEnd")
  .chrAsNum(tbl)
}
makeHg19 <- function(){
  tbl <- merge(getChrLength(), getCentromeres(), by="chrom")
  cumlen <- c(0, cumsum(as.numeric(tbl$length))[-nrow(tbl)])
  cbind.data.frame(tbl, cumlen=cumlen)
}

gene_anno <- function (top.bins, results) {
  #top.bins data must have bin=col1, chr=col2, start=col3, stop=col4
  #results (gene annotation) must have hgnc_symbol=col1, entrez=col2, chr=col3, start=col4, stop=col5, geneID=col6
  k <- 1
  i <- 1
  l <- 1
  g <- c <- sta <- sto <- id <- b <- NULL
  g.l <- c.l <- sta.l <- sto.l <- id.l <- b.l <- list() #to match genes to bins
  
  for ( k in 1:nrow(top.bins) ) {
    bin <- top.bins[k,1]
    chr <- top.bins[k,2]
    start <- top.bins[k,3]
    stop <- top.bins[k,4]
    
    wd <- results[which(results[,3] == chr),]
    
    for ( i in 1:nrow(wd) ) {
      if ( (chr == wd[i,3]) && ( between(wd[i,4], start, stop) || between(wd[i,5], start, stop)) ) {
        g <- append(g,wd[i,2])
        c <- append(c,wd[i,3])
        sta <- append(sta,wd[i,4])
        sto <- append(sto,wd[i,5])
        id <- append(id,wd[i,6])
        b <- append(b,bin)
      }
    }
    g.l[[l]] <- g
    c.l[[l]] <- c
    sta.l[[l]] <- sta
    sto.l[[l]] <- sto
    id.l[[l]] <- id
    b.l[[l]] <- b
    l <- l+1
    g <- c <- sta <- sto <- id <- b <- NULL
  }
  top.genes <- data.frame(bin=unlist(b.l), entrezgene=unlist(g.l), chr=unlist(c.l), gene.start=unlist(sta.l), gene.stop=unlist(sto.l), geneID=unlist(id.l))
  #top.genes[top.genes==""] <- NA
  top.genes <- na.omit(top.genes) #remove rows with blanks
  #colnames(top.genes) <- c("bin",colnames(wd)[2])
  #top.genes <- merge(top.genes, results, by = colnames(wd)[2]) #add in data from gene annotation data
  #colnames(top.genes)[4] <- "chr"
  return(top.genes)
}

