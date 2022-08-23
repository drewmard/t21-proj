getDORCScores <- function(ATAC.se,
                          dorcTab,
                          normalizeATACmat=TRUE,
                          geneList=NULL,
                          nCores=4){
  
  
  if(!all(c("Peak","Gene") %in% colnames(dorcTab)))
    stop("The provided gene-peak table must have columns named Peak and Gene ..")
  
  if(any(dorcTab$Peak > nrow(ATAC.se)))
    stop("One or more peak indices in the gene-peak table are larger than the total number of peaks in the provided ATAC SE object ..\n Make sure the exact same SummarizedExperiment object is provided here as was for running the runGenePeakcorr function ..\n")
  
  if(!is.null(geneList)){
    if(!(all(geneList %in% as.character(dorcTab$Gene))))
      stop("One or more of the gene names supplied is not present in the gene-peak table provided..\n")
    
    if(length(geneList) > 50){
      message("Running DORC scoring for ",length(geneList)," genes: ",paste(geneList[1:20],collapse=", "),", ... , ... , ... (truncated display)")
    } else {
      message("Running DORC scoring for ",length(geneList)," genes: ",paste(geneList,collapse = "\n"))
    }
    
    cat("........\n")
    
    dorcTab <- dorcTab[dorcTab$Gene %in% geneList,] # Filter only to these genes
    dorcGenes <- sort(as.character(unique(dorcTab$Gene)))
  } else {
    dorcGenes <- sort(as.character(unique(dorcTab$Gene)))
    cat("Running DORC scoring for all genes in annotation! (n = ",length(dorcGenes),")\n",sep="")
  }
  
  
  if(normalizeATACmat){
    # Normalize
    cat("Normalizing scATAC counts ..\n")
    ATAC.mat <- centerCounts(assay(ATAC.se,chunkSize = 5000))
    gc()
  } else {
    cat("Assuming provided scATAC counts are normalized ..\n")
    ATAC.mat <- assay(ATAC.se)
  }
  
  time_elapsed <- Sys.time()
  
  if(Sys.info()['sysname'] %in% "Windows"){
    message("Windows OS detected .. Cannot support parallilzation using mclapply for mc.cores > 1")
    message("Using 1 core instead ..\n")
    nCores <- 1
  }
  
  cat("Computing DORC scores ..\n")
  cat("Running in parallel using ", nCores, "cores ..\n")
  
  dorcMatL <- pbmcapply::pbmclapply(X=dorcGenes,
                                    FUN=function(x) {
                                      
                                      dorcPeaks <- unique(dorcTab$Peak[dorcTab$Gene %in% x])
                                      
                                      if(length(dorcPeaks) > 1) {
                                        dorcCounts <- Matrix::colSums(ATAC.mat[dorcPeaks,])
                                      } else if(length(dorcPeaks==1)) {
                                        dorcCounts <- ATAC.mat[dorcPeaks,]
                                      }
                                    },mc.cores = nCores)
  
  dorcMat <- Matrix::Matrix(do.call('rbind',dorcMatL),sparse=TRUE)
  
  rownames(dorcMat) <- dorcGenes
  
  time_elapsed <- Sys.time() - time_elapsed
  cat(paste("\nTime Elapsed: ",time_elapsed, units(time_elapsed)),"\n\n")
  
  return(dorcMat)
  
}
