runGenePeakcorr <- function(ATACdf, # SummarizedExperiment object of scATAC data
                            RNAdf, # Paired normalized scRNA-seq data, with gene names as rownames
                            genome, # Must be one of "hg19", "mm10", or "hg38"
                            geneList=NULL, # 2 or more valid gene symbols (if only running on subset of genes)
                            windowPadSize=50000, # base pairs padded on either side of gene TSS
                            normalizeATACmat=TRUE, # Whether or not to normalize scATAC counts (default is yes, assumes raw counts)
                            nCores=4, # Number of cores if parallelization support
                            keepPosCorOnly=TRUE,
                            keepMultiMappingPeaks=FALSE,
                            n_bg=100, # Number of background peaks to use
                            p.cut=NULL # Optional, if specified, will only return sig peak-gene hits
) {
  
  # ATACdf <- ATAC.se
  # RNAdf <- RNAmat@assays$RNA@data
  # geneList=NULL # 2 or more valid gene symbols (if only running on subset of genes)
  # windowPadSize=50000 # base pairs padded on either side of gene TSS
  # normalizeATACmat=TRUE # Whether or not to normalize scATAC counts (default is yes, assumes raw counts)
  # keepPosCorOnly=TRUE
  # keepMultiMappingPeaks=FALSE
  # genome = "hg38"# One of hg19, mm10 or hg38 
  # nCores = 8
  # p.cut = NULL # Set this to NULL and we can filter later
  # n_bg = 100
  
  
  #######################3
  
  stopifnot(inherits(ATACdf,"RangedSummarizedExperiment"))
  stopifnot(inherits(RNAdf,c("Matrix","matrix")))
  
  if(!all.equal(ncol(ATACdf),ncol(RNAdf)))
    stop("Input ATAC and RNA objects must have same number of cells")
  
  message("Assuming paired scATAC/scRNA-seq data ..")
  
  peakRanges.OG <- granges(ATACdf) # Peak ranges in reference input SE (pre-filtering)
  
  # Function needs rownames for both matrices or gives error
  rownames(ATACdf) <- paste0("Peak",1:nrow(ATACdf))
  ATACmat <- assay(ATACdf) # Rownames preserved
  
  # Normalize peak counts
  if(normalizeATACmat)
    ATACmat <- centerCounts(ATACmat) # Rownames preserved
  
  if(is.null(rownames(RNAdf)))
    stop("RNA matrix must have gene names as rownames")
  
  # Check for peaks/genes with 0 accessibility/expression
  
  if(any(Matrix::rowSums(assay(ATACdf))==0)){
    message("Peaks with 0 accessibility across cells exist ..")
    message("Removing these peaks prior to running correlations ..")
    peaksToKeep <- Matrix::rowSums(assay(ATACdf))!=0
    ATACdf <- ATACdf[peaksToKeep,] # Subset ranges
    ATACmat <- ATACmat[peaksToKeep,]
    message("Important: peak indices in returned gene-peak maps are relative to original input SE")
  }
  
  
  peakRanges <- granges(ATACdf) # Peak ranges
  
  if(any(Matrix::rowSums(RNAdf)==0)){
    message("Genes with 0 expression across cells exist ..")
    message("Removing these genes prior to running correlations ..")
    genesToKeep <- Matrix::rowSums(RNAdf)!=0
    RNAdf <- RNAdf[genesToKeep,]
  }
  
  cat("Number of peaks in ATAC data:",nrow(ATACmat),"\n")
  cat("Number of genes in RNA data:",nrow(RNAdf),"\n")
  
  
  if (!genome %in% c("hg19", "hg38", "mm10"))
    stop("You must specify one of hg19, hg38 or mm10 as a genome build for currently supported TSS annotations..\n")
  switch(genome, hg19 = {
    TSSg <- hg19TSSRanges
  }, hg38 = {
    TSSg <- hg38TSSRanges
  }, mm10 = {
    TSSg <- mm10TSSRanges
  })
  
  # Keep genes that have annotation and are in RNA matrix
  names(TSSg) <- as.character(TSSg$gene_name)
  
  if(!is.null(geneList)){
    if(length(geneList)==1)
      stop("Please specify more than 1 valid gene symbol")
    
    if(any(!geneList %in% names(TSSg))){
      cat("One or more of the gene names supplied is not present in the TSS annotation specified: \n")
      cat(geneList[!geneList %in% names(TSSg)], sep = ", ")
      cat("\n")
      stop()
    }
    
    TSSg <- TSSg[geneList]
  }
  
  # Checking in case some genes in RNA don't overlap our TSS annotations
  genesToKeep <- intersect(names(TSSg),rownames(RNAdf))
  cat("\nNum genes overlapping TSS annotation and RNA matrix being considered: ",length(genesToKeep),"\n")
  
  # Match gene order in RNA matrix and TSS ranges
  RNAdf <- RNAdf[genesToKeep,]
  TSSg <- TSSg[genesToKeep]
  
  # Pad TSS by this much *either side*
  TSSflank <- GenomicRanges::flank(TSSg,
                                   width = windowPadSize,
                                   both = TRUE)
  
  # Get peak summit
  cat("\nTaking peak summits from peak windows ..\n")
  peakSummits <- resize(peakRanges,width = 1,fix = "center")
  
  # Find overlap of all peaks to all genes given window
  # Subject is Peaks, query is Gene
  cat("Finding overlapping peak-gene pairs ..\n")
  genePeakOv <- findOverlaps(query = TSSflank,subject = peakSummits)
  numPairs <- length(genePeakOv)
  
  cat("Found ",numPairs,"total gene-peak pairs for given TSS window ..\n")
  
  cat("Number of peak summits that overlap any gene TSS window: ",length(unique(subjectHits(genePeakOv))),"\n")
  cat("Number of gene TSS windows that overlap any peak summit: ",length(unique(queryHits(genePeakOv))),"\n\n")
  
  # For each gene, determine observed correlation of each overlapping peak to its associated gene (gene expression)
  
  # For each of those genes, also determine correlation based on background peaks (run in parallel) and save
  # Do this in parallel, and get data frame of gene-peak-pearson values
  # Fetch background peaks for each peak tested (i.e. that has overlap in window with gene)
  set.seed(123)
  cat("Determining background peaks ..\n")
  
  if(is.null(rowData(ATACdf)$bias)){
    if(genome %in% "hg19")
      myGenome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
    if(genome %in% "mm10")
      myGenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
    if(genome %in% "hg38")
      myGenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    
    ATACdf <- chromVAR::addGCBias(ATACdf,genome=myGenome) }
  
  cat("Using ",n_bg," iterations ..\n\n")
  
  set.seed(123)
  bg <- chromVAR::getBackgroundPeaks(ATACdf,niterations=n_bg)
  
  cat("Computing gene-peak correlations ..\n")
  
  pairsPerChunk <- 500
  
  # This defines the outer (larger chunks)
  largeChunkSize <- 5000
  
  startingPoint <- 1 # If for any reason the workers fail, resume from where it failed by specifying the starting point here
  chunkStarts <- seq(startingPoint, numPairs, largeChunkSize)
  chunkEnds <- chunkStarts + largeChunkSize -1
  chunkEnds[length(chunkEnds)] <- numPairs
  
  library(doParallel)
  
  dorcList <- list()
  for(i in 1:length(chunkStarts)){
    cat("Running pairs: ",chunkStarts[i], "to",chunkEnds[i],"\n")
    # This fill further chunk this up and run in parallel, saving the merged output ObsCor
    ObsCor <- PeakGeneCor(ATAC = ATACmat,
                          RNA = RNAdf,
                          OV = genePeakOv[chunkStarts[i]:chunkEnds[i]],
                          chunkSize = pairsPerChunk,
                          ncores = nCores,
                          bg = bg)
    gc()
    
    dorcList[[i]] <- ObsCor
  }
  
  cat("\nMerging results ..\n")
  dorcTab <- bind_rows(dorcList)
  
  cat("Performing Z-test for correlation significance ..\n")
  permCols <- 4:(ncol(bg)+3)
  
  
  if (keepPosCorOnly){
    # Filter to positive correlations
    cat("Only considering positive correlations ..\n")
    dorcTab <- dorcTab %>% filter(rObs > 0)
  }
  
  if(!keepMultiMappingPeaks){
    # Remove multi-mapping peaks (force 1-1 mapping)
    cat("Keeping max correlation for multi-mapping peaks ..\n")
    dorcTab <- dorcTab %>% group_by(Peak) %>% filter(rObs==max(rObs))
  }
  
  # Swap gene number for gene symbol from TSS annotation lookup
  dorcTab$Gene <- as.character(TSSg$gene_name)[dorcTab$Gene]
  
  # Swap peak numbers to match reference input peak numbers
  # This only changes if some peaks had zero accessibility and were filtered out internally
  # Use rownames from reference matching
  dorcTab$Peak <- as.numeric(splitAndFetch(rownames(ATACmat)[dorcTab$Peak],"Peak",2))
  
  # # Z test pval
  dorcTab$rBgSD <- matrixStats::rowSds(as.matrix(dorcTab[,permCols]))
  dorcTab$rBgMean <- rowMeans(dorcTab[,permCols])
  dorcTab$pvalZ <- 1-stats::pnorm(q = dorcTab$rObs,mean = dorcTab$rBgMean,sd = dorcTab$rBgSD)
  
  
  cat("\nFinished!\n")
  
  if(!is.null(p.cut)){
    cat("Using significance cut-off of ",p.cut," to subset to resulting associations\n")
    dorcTab <- dorcTab[dorcTab$pvalZ <= p.cut,] # Subset to significant correlations only
  }
  
  # Add peak ranges to final data frame output
  dorcTab$PeakRanges <- paste(as.character(seqnames(peakRanges.OG[dorcTab$Peak])),paste(start(peakRanges.OG[dorcTab$Peak]),end(peakRanges.OG[dorcTab$Peak]),sep="-"),sep=":")
  
  return(as.data.frame(dorcTab[,c("Peak","PeakRanges","Gene","rObs","pvalZ")],stringsAsFactors=FALSE))
}
