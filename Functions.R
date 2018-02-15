##########################
#Functions used in analyses reported in "O-GlcNAc transferase recognizes 
#protein substrates using an asparagine ladder in the TPR superhelix"

#######################
###needed libraries
###if any are missing, install them
#######################

library(limma)
library(dplyr)

#######################
###Functions to read in raw data
#######################
ReadInArrays <- function(targets){
  #reads in the arrays listed in targets by their batch, corrects their backgrounds, and outputs
  #a list by batch of array data in maimage format
  Raw.Data <- list()
  spottypes <- readSpotTypes()
  for(i in unique(targets$Batch)){
    Raw.Data[[i]] <- read.maimages(targets$FileName[targets$Batch==i],
                                   source='genepix.median', 
                                   annotation = c('Block', 'Row', 'Column', 'ID', 'Name',
                                                  'Description'),
                                   wt.fun = wtflags(weight=0,cutoff=-50))
    Raw.Data[[i]]$genes$Status <- controlStatus(spottypes, Raw.Data[[i]])
  }
  Data <- list()
  for(i in unique(targets$Batch)){
    Data[[i]] <- backgroundCorrect(Raw.Data[[i]], method = 'normexp')
  }
  Data
}
modifyAccession <- function(genesFrame){
  Names <- sub('^Hs~Ref:([[:alnum:]]+)-[1-2]~CAT.*$','\\1',genesFrame$Name)
  Names <- sub('^Hs~IVGN:(PM_[0-9]{4})~Ext:([[:alnum:]_\\)\\([:blank:]-]+)~.*$','Invitrogen: \\1, \\2', Names)
  Names <- sub('^Hs~IVGN:([[:alnum:]_-]+)~(CAT_[[:alnum:]-]+)~.*$','Invitrogen: \\1, \\2', Names)
  Names <- sub('^Hs~Abcam:([a-zA-Z-]+)~(CAT_[[:alnum:]]+)~.*$','Abcam: \\1, \\2', Names)
  Names <- sub('^Hs~[[:alpha:]/]+:([[:alnum:]_\\.]*)~[[:alpha:]]+.*$','\\1', Names)
  Names <- sub('^([[:alnum:]-]+)~N/A$','\\1', Names)
  Names <- sub('^([[:alnum:]]+)~RFU:[0-9\\.]+$','\\1', Names)
  Names <- sub('^Invitrogen: (NM_[0-9]+), CAT_.*$', '\\1', Names)
  Names
}
#Changes Controls to contain info on which control in the description line
modifyDescription <- function(genesFrame){
  #generates a modified Description list with 'Control' Swapped for 
  #'Control: (Name Of Control)'
  for(i in which(grepl('Control',genesFrame$Description))){
    genesFrame$Description[i] <- paste('Control: ',
                                       sub('^(.*)~.*$','\\1',
                                           genesFrame$Name[i]),
                                       sep = ''
    )
    
  }
  genesFrame$Description
}
addV5Conc <- function(geneList, ContFile){
  #uses the appropriate ContFile to add the V5-biotin control spot concentrations to 
  #the geneList. All Concs are in nM
  V5prots <- grep('V5', ContFile$Control.Proteins, value=T)
  geneList[,'V5conc'] <- 0
  for(i in V5prots){
    geneList[grepl(i,geneList$Name),'V5conc'] <-
      ContFile[ContFile$Control.Proteins == i, 'Concentration..nM.']
  }
  geneList
}

#######################
###Functions to adjust data before analysis (non-normalization)
#######################
ProteinOverlap <- function(uArraySets){
  #takes in a list of RG objects, and references/modifies the Name column of the gene lists
  #to determine the overlapping genes based upon accession. The list is returned with just 
  #the overlapping genes, sorted by Accession number.
  outputArray <- list()
  accession <- list()
  
  for(i in names(uArraySets)){
    currAccession <- modifyAccession(uArraySets[[i]]$genes)
    uArraySets[[i]]$genes$Name <- currAccession
    for(j in names(accession)){#nothing happens in first round
      accession[[j]] <- accession[[j]][accession[[j]] %in% currAccession]
      currAccession <- currAccession[currAccession %in% accession[[j]]]
    }
    accession[[i]] <- currAccession
    
  }
  for(i in names(uArraySets)){
    outputArray[[i]] <- uArraySets[[i]]
    whichProt <- outputArray[[i]]$genes$Name %in% accession[[i]]
    outputArray[[i]]$genes <- outputArray[[i]]$genes[whichProt,]
    outputArray[[i]]$R <- outputArray[[i]]$R[whichProt,]
    outputArray[[i]]$G <- outputArray[[i]]$G[whichProt,]
  }
  outputArray
}
Collapse2Means <- function(uArray, Symbol = FALSE, bySymbol = FALSE, otherCols = NULL){
  #Takes in uArray (with elements $R, $genes$Description) and 
  #collapses to Description (or gene Symbol if the bySymbol argument is TRUE),
  #with first modifying the names to have the controls be sensible.
  #symbol argument is only to be used if HGNC symbols have been added prior to this to the gene
  #data frame of the microarray object
  uArray$Description <- modifyDescription(uArray$genes)
  uArray
  if(Symbol){
    uArray$Symbol <- uArray$genes$Symbol
    if(bySymbol){
      colnamesgenes <- c('Symbol',colnames(uArray$R))
      genesR <- cbind(data.frame(uArray$Symbol),data.frame(uArray$R),data.frame(uArray$genes[,otherCols]))
      colnames(genesR) <- c(colnamesgenes,otherCols)
      byGene <- group_by(genesR, Symbol) %>% summarise_all(funs(mean))
    }else{
      colnamesgenes <- c('Description','Accession','Symbol',colnames(uArray$R))
      genesR <- cbind(data.frame(uArray$Description),data.frame(uArray$genes$Name),data.frame(uArray$Symbol),data.frame(uArray$R),
                      data.frame(uArray$genes[,otherCols]))
      colnames(genesR) <- c(colnamesgenes,otherCols)
      byGene <- group_by(genesR, Accession, Symbol) %>% summarise_all(funs(mean))
      x <- group_by(genesR, Accession, Symbol) %>% summarise('Description' = first(Description))
      byGene$Description <- x$Description
    }
  }else{
    genesR <- cbind(data.frame(uArray$Description),data.frame(uArray$genes$Name),
                    data.frame(uArray$R),data.frame(uArray$genes[,otherCols]))
    colnamesgenes <- c('Description', 'Accession', colnames(uArray$R))
    colnames(genesR) <- c(colnamesgenes,otherCols)
    byGene <- group_by(genesR, Accession) %>% summarise_all(funs(mean)) #this is problematic and should be
    #fixed, but is fine for now...
    #note: gives warnings because you can't take mean of a character vector,
    #but still easy implementation
    x <- group_by(genesR, Accession) %>% summarise('Description' = first(Description))
    byGene$Description <- x$Description
  }
  data.frame(byGene)
}

addHGNC <- function(Data, masterList, masterListDescCol = "Protein.Header",
                    descriptionCol = "Description",
                    knownOGlcNAc = TRUE, indivDB = FALSE){
  #reads through Data "Description" and uses the dataset listed in masterList to add a column of 
  #HGNC symbols. If "knownOGlcNAc" is TRUE, a second column listing if this is in one of the mined databases
  #(in TRUE/FALSE format) is added. if "indivDB" is TRUE, then 2 columns will be added breaking out by 
  #database from with data is obtained
  Data$HGNC <- ''
  if(knownOGlcNAc){
    Data$Known.GlcNAc <- FALSE
  }
  if(indivDB){
    Data$PhosphoSitePlus <- FALSE
    Data$Wang <- FALSE
  }
  Descs <- unique(Data[,descriptionCol])
  for(i in Descs){
    currData <- masterList[(masterList[,masterListDescCol] == i) &
                             !is.na(masterList[,masterListDescCol]),]
    
    Data$HGNC[
      Data[,descriptionCol] == i] <- currData$HGNC[1]
    if(knownOGlcNAc){
      Data$Known.GlcNAc[
        Data[,descriptionCol] == i] <- currData$KnownGlyco[1]
    }
    if(indivDB){
      Data$PhosphoSitePlus[
        Data[,descriptionCol] == i] <- currData$PhosphoSitePlus[1]
      Data$Wang[
        Data[,descriptionCol] == i] <- currData$Wang[1]
    }
  }
  Data
}
medianSignal <- function(arrayDF, listColCategory, nonDataCols = c('Accession','Description','HGNC','Known.GlcNAc',
                                                                   'PhosphoSitePlus','Wang')){
  #gives the by-category log2 median (converted back to non-log2) as columns. Uses names from list
  log2arrayDF <- arrayDF
  log2arrayDF[,!(colnames(log2arrayDF) %in% nonDataCols)] <- log2(log2arrayDF[,!(colnames(log2arrayDF) %in% nonDataCols)])
  CatMedians <- list()
  for(i in names(listColCategory)){
    CatMedians[[paste(i,'Median')]] <- 2^apply(log2arrayDF[,listColCategory[[i]]], 1,median, na.rm = T)
  }
  do.call(cbind.data.frame,CatMedians)
}


#######################
###Normalization Functions for Biotin staining and ERalpha activity normalization
#######################

#normalizes for staining using the biotinylated V5 controls
normStain <- function(dataSet){
  #normalizes data based upon robust linear fit of 
  fitsum <- list()
  for(i in unique(dataSet$genes$Block)){
    fitsum[[i]] <- list()
    for(j in 1:dim(dataSet$R)[2]){
      #now have isolated single block of single array
      x <- dataSet$genes$V5conc[(dataSet$genes$V5conc>0) & (dataSet$genes$Block == i)]
      y <- dataSet$R[(dataSet$genes$V5conc>0) & (dataSet$genes$Block == i),j]
      fit <- rlm(y ~ x, maxit = 200)
      fitsum[[i]][[j]] <- summary(fit)
      #this will allow extraction of slope
    }
  }
  for(i in unique(dataSet$genes$Block)){
    for(j in 1:dim(dataSet$R)[2]){
      tempdata <- dataSet$R[dataSet$genes$Block==i,j]
      if(!is.na(fitsum[[i]][[j]]$coef[1,2]) & (fitsum[[i]][[j]]$coef[2]>0)){
        #If the data is good for a given block, use that fit
        blockintercept <- fitsum[[i]][[j]]$coef[1]
        blockslope <- fitsum[[i]][[j]]$coef[2]
      }else{#if it isn't, use an average of the fits on the array
        blockintercepts <- NULL
        blockslopes <- NULL
        for(k in 1:i){
          blockintercepts <- c(blockintercepts, fitsum[[k]][[j]]$coef[1])
          blockslopes <- c(blockslopes, fitsum[[k]][[j]]$coef[2])
        }
        blockintercept <- mean(blockintercepts)
        blockslope <- mean(blockslopes)
        if(blockslope < 0){
          stop(paste('array: ', as.character(colnames(dataSet$R)[j]), ' block: ', as.character(i), '\n'))
        }
      }
      #ditto slopes
      tempdata <-  tempdata / blockslope
      dataSet$R[dataSet$genes$Block==i,j] <- tempdata
      if(sum(is.na(tempdata))>0){print(c(i,j))}
    }
  }
  dataSet
}
#normalizes for OGT activity based upon relative activity on estrogen receptor alpha
NormERalpha <- function(DataSet,colFactors,relativeActivity,flagInactive = FALSE){
  #set up data as desired into data frame ordered by tidy data principles
  DataSet.df <- data.frame(DataSet$R)
  DataSet.df$Block <- DataSet$genes$Block
  DataSet.df$Description <- DataSet$genes$Description
  DataSet.df$Accession <- DataSet$genes$Name
  DS.tidy <- melt(DataSet.df, id.vars = c('Description','Accession','Block'))
  colnames(DS.tidy)[colnames(DS.tidy)=='variable'] <- 'Array'
  #figure out which arrays are controls, etc.
  contArrays <- make.names(colnames(DataSet$R)[colFactors==0])
  #list OGT-treated arrays by type of array
  OGTArrays <- list()
  for(i in unique(colFactors[colFactors>0])){
    OGTArrays[[i]] <- make.names(colnames(DataSet$R)[colFactors==i])
  }
  #find control median of ERalpha
  contERmedian <- median(DS.tidy$value[(DS.tidy$Accession == 'ERa') &
                                         (DS.tidy$Array %in% contArrays)], na.rm = TRUE)
  #find block ERalpha medians
  blockERmedians <- group_by(DS.tidy[(DS.tidy$Accession == 'ERa') &
                                       !(DS.tidy$Array %in% contArrays),], Array, Block) %>% 
    summarise('ERaMedian' = median(value))
  #find block Activities against ERalpha
  blockERmedians$ERActivity <- blockERmedians$ERaMedian - contERmedian
  #find gene control medians
  geneContMedians <- group_by(DS.tidy[(DS.tidy$Array %in% contArrays),], Accession) %>% 
    summarise('ControlMedian' = median(value))
  #and wt median
  wtERmedian <- median(DS.tidy$value[(DS.tidy$Accession == 'ERa') &
                                       (DS.tidy$Array %in% OGTArrays[[1]])], na.rm = TRUE)
  #now "average" activity for wild-type
  avgActivity <- wtERmedian - contERmedian
  #now go OGT type by OGT type
  for(i in 1:length(OGTArrays)){
    #find by-arrayxBlock accession numbers where there is a sign of activity
    #first join in the control medians, then get gene activities
    geneActivities <- group_by(DS.tidy[(DS.tidy$Array %in% OGTArrays[[i]]),], Array, Block, Accession) %>% 
      summarise('GeneMedian' = median(value)) %>%
      left_join(geneContMedians, by = 'Accession') %>% mutate('Activity' = GeneMedian-ControlMedian) %>%
      left_join(blockERmedians, by = c('Array','Block')) %>%
      mutate('CorrectedActivity' = Activity*avgActivity/ERActivity*relativeActivity[i])
    geneActivities$Activity[geneActivities$Activity<0] <- 0
    geneActivities$CorrectedActivity[geneActivities$CorrectedActivity<0] <- 0
    #now cycle through all the arrays and adjust the single spots, simply subtracting by the Activity
    #and adding the Corrected Activity
    for(j in OGTArrays[[i]]){
      jj <- quo(j)
      toCorrect <- filter(geneActivities, Array == j) %>%
        right_join(DataSet.df[,c('Accession','Block',j)], by = c('Accession','Block')) %>%
        mutate(ValueCorrected = .data[[!!jj]] - .data$Activity + .data$CorrectedActivity)
      DataSet.df[,colnames(DataSet.df)==j] <- toCorrect$ValueCorrected
    }
  }
  DataSet$R <- DataSet.df[,make.names(colnames(DataSet$R))]
  return(DataSet)
}

#######################
###Permutation Test Functions
#######################
eBayesArrays <- function(arrayDF, factors, correctMedian = FALSE, nonDataCols = c('Accession','Description','HGNC','Known.GlcNAc',
                                                                                  'PhosphoSitePlus','Wang')){
  if(correctMedian){
    ExpMatrix <- log2(arrayDF[,!(colnames(arrayDF) %in% nonDataCols)])
    for(f in unique(factors[factors!='control'])){
      ExpMatrix[,factors == f] <- ExpMatrix[,factors == f] - 
        median(unlist(ExpMatrix[,factors == f]),na.rm = T) + 
        median(unlist(ExpMatrix[,factors == 'control']), na.rm = T)
    }
  }else{
    ExpMatrix <- log2(arrayDF[,!(colnames(arrayDF) %in% nonDataCols)])
  }
  factors <- factor(factors)
  designArray <- model.matrix(~ 0 + factors)
  colnames(designArray) <- sub('factors(.*$)', '\\1', colnames(designArray))
  fit <- lmFit(ExpMatrix, design = designArray, method = 'robust', maxit = 400)
  contrast.matrix <- makeContrasts(wt-control, nc5N5A-control, wt-nc5N5A, levels=designArray)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  eBayesOut <- eBayes(fit2, robust = TRUE)
  eBayesOut
}
permuteBackgroundwtcont <- function(p.value,tval,inputData,arrayOrder,toCompare, alpha = 0.1, oneway = TRUE,
                                    N = 100){
  #generates a random permutation distribution of eBayes-derived t tests. Only does two way
  #comparison at a time. Need to give single-column of p value and tval. oneway assumes only
  #looking for greater if used. CANNOT be used with array weights.
  if(oneway){
    workingData <- inputData[tval < quantile(tval, probs = 1-alpha),arrayOrder %in% toCompare]
  }else{
    workingData <- inputData[(tval > quantile(tval, probs = alpha/2))&
                               (tval < quantile(tval, probs = 1-alpha/2)),
                             arrayOrder %in% toCompare]
  }
  newOrder <- arrayOrder[arrayOrder %in% toCompare]
  back.tvals <- list()
  for(i in 1:N){
    workingOrder <- sample(newOrder)
    workingOrder <- factor(workingOrder)
    work.design <- model.matrix(~ 0 + workingOrder )
    colnames(work.design) <- sub('workingOrder(.*$)', '\\1',
                                 colnames(work.design))
    currFit <- lmFit(workingData, design = work.design)
    contrast.matrix <-makeContrasts(wt - control, levels=work.design)
    
    currFit2 <- contrasts.fit(currFit, contrast.matrix)
    currFit2 <- eBayes(currFit2, robust = TRUE)
    back.tvals[[i]] <- currFit2$t
  }
  back <- NULL
  back <- unlist(back.tvals)
  pvalout <- NULL
  if(oneway){
    pvalout <- sapply(tval, function(i,back){sum(back>i)/length(back)},back=back)
  }else{
    pvalout <- sapply(tval, function(i,back){
      if(i>=median(back)){
        2*sum(back>i)/length(back)
      }else{
        2*sum(back<i)/length(back)
      }
    }, back=back)
  }  
  pvalout <- data.frame('p.value' = pvalout)
  rownames(pvalout) <- rownames(inputData)
  output <- list('BackgroundDist' = back, 'pvals' = pvalout)
  output
}
permuteBackground5N5Acont <- function(p.value,tval,inputData,arrayOrder,toCompare, alpha = 0.1, oneway = TRUE,
                                      N = 100){
  #generates a random permutation distribution of eBayes-derived t tests. Only does two way
  #comparison at a time. Need to give single-column of p value and tval. oneway assumes only
  #looking for greater if used. CANNOT be used with array weights.
  if(oneway){
    workingData <- inputData[tval < quantile(tval, probs = 0.9),arrayOrder %in% toCompare]
  }else{
    workingData <- inputData[(tval > quantile(tval, probs = 0.05))&
                               (tval < quantile(tval, probs = 0.95)),
                             arrayOrder %in% toCompare]
  }
  newOrder <- arrayOrder[arrayOrder %in% toCompare]
  back.tvals <- list()
  for(i in 1:N){
    workingOrder <- sample(newOrder)
    workingOrder <- factor(workingOrder)
    work.design <- model.matrix(~ 0 + workingOrder )
    colnames(work.design) <- sub('workingOrder(.*$)', '\\1',
                                 colnames(work.design))
    currFit <- lmFit(workingData, design = work.design)
    contrast.matrix <-makeContrasts(nc5N5A - control, levels=work.design)
    
    currFit2 <- contrasts.fit(currFit, contrast.matrix)
    currFit2 <- eBayes(currFit2, robust = TRUE)
    back.tvals[[i]] <- currFit2$t
  }
  back <- NULL
  back <- unlist(back.tvals)
  pvalout <- NULL
  if(oneway){
    pvalout <- sapply(tval, function(i,back){sum(back>i)/length(back)},back=back)
  }else{
    pvalout <- sapply(tval, function(i,back){
      if(i>=median(back)){
        2*sum(back>i)/length(back)
      }else{
        2*sum(back<i)/length(back)
      }
    }, back=back)
  }  
  pvalout <- data.frame('p.value' = pvalout)
  rownames(pvalout) <- rownames(inputData)
  output <- list('BackgroundDist' = back, 'pvals' = pvalout)
  output
}
permuteBackgroundwtv5N5A <- function(p.value,tval,inputData,arrayOrder,toCompare, alpha = 0.1, oneway = TRUE,
                                     N = 100){
  #generates a random permutation distribution of eBayes-derived t tests. Only does two way
  #comparison at a time. Need to give single-column of p value and tval. oneway assumes only
  #looking for greater if used. CANNOT be used with array weights.
  if(oneway){
    workingData <- inputData[tval < quantile(tval, probs = 0.9),arrayOrder %in% toCompare]
  }else{
    workingData <- inputData[(tval > quantile(tval, probs = 0.05))&
                               (tval < quantile(tval, probs = 0.95)),
                             arrayOrder %in% toCompare]
  }
  newOrder <- arrayOrder[arrayOrder %in% toCompare]
  back.tvals <- list()
  for(i in 1:N){
    workingOrder <- sample(newOrder)
    workingOrder <- factor(workingOrder)
    work.design <- model.matrix(~ 0 + workingOrder )
    colnames(work.design) <- sub('workingOrder(.*$)', '\\1',
                                 colnames(work.design))
    currFit <- lmFit(workingData, design = work.design)
    contrast.matrix <-makeContrasts(wt - nc5N5A, levels=work.design)
    
    currFit2 <- contrasts.fit(currFit, contrast.matrix)
    currFit2 <- eBayes(currFit2, robust = TRUE)
    back.tvals[[i]] <- currFit2$t
  }
  back <- NULL
  back <- unlist(back.tvals)
  pvalout <- NULL
  if(oneway){
    pvalout <- sapply(tval, function(i,back){sum(back>i)/length(back)},back=back)
  }else{
    pvalout <- sapply(tval, function(i,back){
      if(i>=median(back)){
        2*sum(back>i)/length(back)
      }else{
        2*sum(back<i)/length(back)
      }
    }, back=back)
  }  
  pvalout <- data.frame('p.value' = pvalout)
  rownames(pvalout) <- rownames(inputData)
  output <- list('BackgroundDist' = back, 'pvals' = pvalout)
  output
}
medianCorrect <- function(arrayDF, controlCol, listOGTCol, nonDataCols = c('Accession','Description','HGNC','Known.GlcNAc',
                                                                           'PhosphoSitePlus','Wang')){
  #works on data columns as listed and adjusts medians in log2, then converts back.
  log2arrayDF <- arrayDF
  log2arrayDF[,!(colnames(log2arrayDF) %in% nonDataCols)] <- log2(log2arrayDF[,!(colnames(log2arrayDF) %in% nonDataCols)])
  medianControl <- median(unlist(log2arrayDF[,controlCol]), na.rm = TRUE)
  OGTMedians <- list()
  for(i in names(listOGTCol)){
    OGTMedians[[i]] <- median(unlist(log2arrayDF[,listOGTCol[[i]]]), na.rm = TRUE)
    for(j in listOGTCol[[i]]){
      log2arrayDF[,j] <- log2arrayDF[,j] - OGTMedians[[i]] + medianControl
    }
    arrayDF[,listOGTCol[[i]]] <- 2^log2arrayDF[,listOGTCol[[i]]]
  }
  arrayDF
}
