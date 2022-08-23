PeakGeneCor <- function(ATAC, # Normalized reads in peaks counts (rownames should  be "Peak1","Peak2" etc.)
                        RNA, # Normalized gene expression counts
                        OV, # Gene TSS - Peak overlap pairs object (Genes: query, Peaks: subject)
                        ncores=4,
                        chunkSize=200,
                        metric="spearman",
                        bg=NULL){
  
  stopifnot(ncol(ATAC)==ncol(RNA))
  
  if(chunkSize > 1000)
    stop("Do not specify very large chunk sizes. Please use chunkSize < 1000")
  
  
  # Number of total gene-peak pairs to chunk up for parallelization
  n <- length(OV)
  starts <- seq(1, n, chunkSize)
  ends <- starts + chunkSize -1
  ends[length(ends)] <- n
  
  OVd <- OV %>% as.data.frame() %>% dplyr::rename("Gene"="queryHits","Peak"="subjectHits")
  
  chunkList <- mapply(c, starts, ends, SIMPLIFY = FALSE)
  
  time_elapsed <- Sys.time()
  
  cat("Running in parallel using ", ncores, "cores ..\n")
  
  cat("Computing observed correlations ..\n")
  
  corList <- pbmcapply::pbmclapply(X=chunkList,
                                   FUN=function(x) {chunkCore(chunk=x,A=ATAC,R=RNA,O=OVd,met=metric)},mc.cores = ncores)
  
  
  if(any(unlist(sapply(corList,is.null)))){
    message("One or more of the chunk processes failed unexpectedly (returned NULL) ..")
    message("Please check to see you have enough cores/memory allocated")
    message("Also make sure you have filtered down to non-zero peaks/genes")
  }
  
  OVd$rObs <- unlist(corList)
  
  
  cat("Finished!\n")
  
  time_elapsed <- Sys.time() - time_elapsed
  cat(paste("\nTime Elapsed: ",time_elapsed, units(time_elapsed)),"\n\n")
  
  if(!is.null(bg)){
    n_iter <- ncol(bg)
    cat("Computing background correlations ..\n")
    
    time_elapsed <- Sys.time()
    
    bgCor <- foreach(i=1:n_iter,.combine = 'cbind') %do% {
                       OVdBg <- OVd[,1:2] # Initialize gene-peak pairing to observed
                       OVdBg$Peak <- bg[OVdBg$Peak,i] # Swap actual peaks with bg peaks for given iteration in pairing
                       bgCorList <- pbmcapply::pbmclapply(X=chunkList,
                                                          FUN=function(x) {chunkCore(chunk=x,A=ATAC,R=RNA,O=OVdBg,met=metric)},mc.cores = ncores)
                       unlist(bgCorList) # Vector of background permuted correlation values for that set of background peaks
                     }
    
    
    if(sum(is.null(bgCor))!=0 | sum(is.na(bgCor))!=0)
      stop("One or more of the chunk processes failed unexpectedly (returned NULL) .. Please check to see you have enough cores/m
           emory allocated")
    
    time_elapsed <- Sys.time() - time_elapsed
    cat(paste("\nTime Elapsed: ",time_elapsed, units(time_elapsed)),"\n\n")
    
    colnames(bgCor) <- paste0("rBg",1:ncol(bgCor))
    OVd <- cbind(OVd,bgCor)
  }
  
  return(OVd)
}
