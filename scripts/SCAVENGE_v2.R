# library(SCAVENGE)
# library(chromVAR)
# library(gchromVAR)
# library(BuenColors)
# library(SummarizedExperiment)
# library(data.table)
# library(BiocParallel)
# # library(BSgenome.Hsapiens.UCSC.hg19)
# library(BSgenome.Hsapiens.UCSC.hg38)
# library(dplyr)
# library(igraph)
# library(viridis)

library(Seurat)
library(Signac)
dfseurat = readRDS("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h_v2/Multiome.RNA_ATAC.h.ChromVAR.rds")

# subset to run things faster
Ncell = 5000
dfseurat = dfseurat[,sample(1:ncol(dfseurat),Ncell,replace = F)]
dfseurat$disease = "H"

dfseurat2 = readRDS("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_ds_v2/Multiome.RNA_ATAC.ds.ChromVAR.rds")

# subset to run things faster
Ncell = 5000
dfseurat2 = dfseurat2[,sample(1:ncol(dfseurat2),Ncell,replace = F)]
dfseurat2$disease = "T21"

dfseurat = merge(dfseurat,dfseurat2)

# numFrag = apply(dfseurat@assays$ATAC@counts,1,sum)

rowcount <- rowSums(dfseurat@assays$ATAC@counts > 0)
rowKeep = rowcount > 0
dfseurat@assays$ATAC@counts = (dfseurat@assays$ATAC@counts[rowKeep, ])
dfseurat@assays$ATAC@ranges = (dfseurat@assays$ATAC@ranges[rowKeep,])
library(SummarizedExperiment)
SE_data <- SummarizedExperiment(assays = list(counts = dfseurat@assays$ATAC@counts),
                             rowData = dfseurat@assays$ATAC@ranges, 
                             colData = DataFrame(names = colnames(dfseurat@assays$ATAC),dataset=dfseurat@meta.data$dataset))
library(BSgenome.Hsapiens.UCSC.hg38)
library(chromVAR)
library(gchromVAR)
SE_data <- addGCBias(SE_data, genome = BSgenome.Hsapiens.UCSC.hg38)
SE_data_bg <- getBackgroundPeaks(SE_data, niterations=200)

# add trait
# use monocytes from example file:
library(SCAVENGE)
# trait_file <- paste0(system.file('extdata', package='SCAVENGE'), "/mono.PP001.bed") # this is hg19
# trait_file = "/home/amarder/R/x86_64-pc-linux-gnu-library/4.1/SCAVENGE/extdata/mono.PP001.hg38.bed"
trait_file = "/oak/stanford/groups/smontgom/amarder/t21-proj/scripts/ALL_finemap.hg38.bed"

trait_import <- importBedScore(rowRanges(SE_data), trait_file, colidx=5)
SE_data_DEV <- computeWeightedDeviations(SE_data, trait_import, background_peaks =
                                             SE_data_bg)
library(dplyr)
z_score_mat <- data.frame(colData(SE_data), z_score=t(assays(SE_data_DEV)[["z"]]) %>% c)

seed_idx <- seedindex(z_score_mat$z_score, 0.05)
scale_factor <- cal_scalefactor(z_score=z_score_mat$z_score, 0.01)

peak_by_cell_mat <- assay(SE_data)
tfidf_mat <- tfidf(bmat=peak_by_cell_mat, mat_binary=TRUE, TF=TRUE, log_TF=TRUE)
lsi_mat <- do_lsi(tfidf_mat, dims=30)
library(harmony)
harmony_mat <- HarmonyMatrix(lsi_mat, SE_data@colData, "dataset", do_pca=FALSE)
# mutualknn30 <- getmutualknn(lsi_mat, 30)
mutualknn30 <- getmutualknn(harmony_mat, 30)

# Network propagation
np_score <- randomWalk_sparse(intM=mutualknn30, rownames(mutualknn30)[seed_idx], gamma=0.05)

omit_idx <- np_score==0
sum(omit_idx)

mutualknn30 <- mutualknn30[!omit_idx, !omit_idx]
np_score <- np_score[!omit_idx]
TRS <- np_score %>% capOutlierQuantile(., 0.95) %>% max_min_scale
TRS <- TRS * scale_factor
trait_mat <- data.frame(z_score_mat[!omit_idx, ], seed_idx[!omit_idx], np_score, TRS)
trait_mat = merge(data.frame(names=rownames(dfseurat@meta.data),clust=dfseurat@meta.data$subclust_v6),trait_mat,by='names')
aggregate(trait_mat$TRS,by=list(trait_mat$clust),median)

trait_permu <- get_sigcell_simple(knn_sparse_mat=mutualknn30, seed_idx=trait_mat$seed_idx,
                                 topseed_npscore=trait_mat$np_score, permutation_times=1000,
                                 true_cell_significance=0.05, rda_output=F, mycores=8, rw_gamma=0.05)
trait_mat2 <- data.frame(trait_mat, trait_permu)

trait_mat2 = merge(trait_mat2,data.frame(names=rownames(dfseurat@meta.data),disease=dfseurat@meta.data$disease),by='names')

trait_mat2 %>%
  group_by(clust) %>%
  summarise(enriched_cell=mean(true_cell_top_idx))

trait_mat2$rev_true_cell_top_idx <- !trait_mat2$true_cell_top_idx
trait_mat2 %>%  
  group_by(clust) %>%
  summarise(depleted_cell=mean(rev_true_cell_top_idx))

trait_mat2 %>%
  group_by(clust,disease) %>%
  summarise(enriched_cell=mean(true_cell_top_idx)) %>% print(n=Inf)



