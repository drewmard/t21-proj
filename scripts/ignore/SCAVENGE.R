library(SCAVENGE)
library(chromVAR)
library(gchromVAR)
library(BuenColors)
library(SummarizedExperiment)
library(data.table)
library(BiocParallel)
library(BSgenome.Hsapiens.UCSC.hg19)
# library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(igraph)
library(viridis)

set.seed(9527)

trait_file <- paste0(system.file('extdata', package='SCAVENGE'), "/mono.PP001.bed")
pbmc5krda <- paste0(system.file('rda', package='SCAVENGE'), "/pbmc5k_SE.rda")
load(pbmc5krda)

SE_pbmc5k <- addGCBias(SE_pbmc5k, genome = BSgenome.Hsapiens.UCSC.hg19)
SE_pbmc5k_bg <- getBackgroundPeaks(SE_pbmc5k, niterations=200)
trait_import <- importBedScore(rowRanges(SE_pbmc5k), trait_file, colidx=5)
SE_pbmc5k_DEV <- computeWeightedDeviations(SE_pbmc5k, trait_import, background_peaks =
                                             SE_pbmc5k_bg)
z_score_mat <- data.frame(colData(SE_pbmc5k), z_score=t(assays(SE_pbmc5k_DEV)[["z"]]) %>% c)
head(z_score_mat)

seed_idx <- seedindex(z_score_mat$z_score, 0.05)
scale_factor <- cal_scalefactor(z_score=z_score_mat$z_score, 0.01)

##

peak_by_cell_mat <- assay(SE_pbmc5k)
tfidf_mat <- tfidf(bmat=peak_by_cell_mat, mat_binary=TRUE, TF=TRUE, log_TF=TRUE)
lsi_mat <- do_lsi(tfidf_mat, dims=30)
mutualknn30 <- getmutualknn(lsi_mat, 30)

library(uwot)
umap_mat = uwot::umap(lsi_mat[,-1],
                min_dist=0.3,
                n_neighbors = 30,
                nn_method = "annoy",
                metric="cosine"
                ) # remove first dim in scATAC data

##
# Network propagation

np_score <- randomWalk_sparse(intM=mutualknn30, rownames(mutualknn30)[seed_idx], gamma=0.05)

omit_idx <- np_score==0
sum(omit_idx)

mutualknn30 <- mutualknn30[!omit_idx, !omit_idx]
np_score <- np_score[!omit_idx]
TRS <- np_score %>% capOutlierQuantile(., 0.95) %>% max_min_scale
TRS <- TRS * scale_factor
mono_mat <- data.frame(z_score_mat[!omit_idx, ], seed_idx[!omit_idx], np_score, TRS)
head(mono_mat)

# UMAP plots of cell type annotation and cell-to-cell graph
p <- ggplot(data=mono_mat, aes(x, y, color=color)) + geom_point(size=1, na.rm = TRUE) +
  pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
p

# Visualize cell-to-cell graph if you have low-dimensional coordinates such as UMAP1 and UMAP2
mutualknn30_graph <- graph_from_adjacency_matrix(mutualknn30, mode = "undirected", diag = F)
plot.igraph(mutualknn30_graph, vertex.size=0.8, vertex.label=NA,
            vertex.color=adjustcolor("#c7ce3d", alpha.f = 1), vertex.frame.color=NA,
            edge.color=adjustcolor("#443dce", alpha.f = 1), edge.width=0.3, edge.curved=.5,
            layout=as.matrix(data.frame(mono_mat$x, mono_mat$y)))

# SCAVENGE TRS
viridis = c("#440154FF", "#472D7BFF", "#3B528BFF", "#2C728EFF", "#21908CFF", "#27AD81FF",
            "#5DC863FF", "#AADC32FF", "#FDE725FF")
p2 <- ggplot(data=mono_mat, aes(x, y, color=TRS)) + 
  geom_point(size=1, na.rm = TRUE, alpha = 0.6) +
  scale_color_gradientn(colors = viridis) + 
  scale_alpha()+
  pretty_plot() + theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
p2

pp2 <- ggplot(data=mono_mat,  aes(x=color, y=TRS))  +
  geom_boxplot(aes(fill=color, color=color), outlier.shape=NA) +
  guides(fill=FALSE) + pretty_plot(fontsize = 10) +
  stat_summary(geom = "crossbar", width=0.65, fatten=0, color="black", 
               fun.data = function(x) { return(c(y=median(x), ymin=median(x), ymax=median(x))) }) + theme(legend.position ="none")
pp2


mono_permu <- get_sigcell_simple(knn_sparse_mat=mutualknn30, seed_idx=mono_mat$seed_idx,
                                 topseed_npscore=mono_mat$np_score, permutation_times=1000,
                                 true_cell_significance=0.05, rda_output=F, mycores=8, rw_gamma=0.05)
mono_mat2 <- data.frame(mono_mat, mono_permu)


# enriched cells
mono_mat2 %>%
  group_by(color) %>%
  summarise(enriched_cell=mean(true_cell_top_idx)) %>%
  ggplot(aes(x=color, y=enriched_cell, fill=color)) + geom_bar(stat="identity") +
  theme_classic()

# depleted cells
mono_mat2$rev_true_cell_top_idx <- !mono_mat2$true_cell_top_idx
mono_mat2 %>%  group_by(color) %>%
  summarise(depleted_cell=mean(rev_true_cell_top_idx)) %>% 
  ggplot(aes(x=color, y=depleted_cell, fill=color)) + 
  geom_bar(stat="identity") +
  theme_classic()




