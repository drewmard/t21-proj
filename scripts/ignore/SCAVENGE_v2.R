# srun --account=default --partition=interactive --time=24:00:00 --mem=128G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash
# module load R/4.1.2

library(Seurat)
library(Signac)
library(SummarizedExperiment)
library(BSgenome.Hsapiens.UCSC.hg38)
library(chromVAR)
library(gchromVAR)
library(SCAVENGE)
library(harmony)
library(dplyr)
library(uwot)
library(parallel)
library(ggplot2)

numThreads = detectCores()/2
viridis = c("#440154FF", "#472D7BFF", "#3B528BFF", "#2C728EFF", "#21908CFF", "#27AD81FF",
                      "#5DC863FF", "#AADC32FF", "#FDE725FF")

########################################3
# Read in data: 

# Healthy data
dfseurat = readRDS("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h_v2/Multiome.RNA_ATAC.h.ChromVAR.rds")
# subset to run things faster
# Ncell = 5000
# dfseurat = dfseurat[,sample(1:ncol(dfseurat),Ncell,replace = F)]
dfseurat$disease = "H"

dfseurat2 = readRDS("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_ds_v2/Multiome.RNA_ATAC.ds.ChromVAR.rds")
# subset to run things faster
# Ncell = 5000
# dfseurat2 = dfseurat2[,sample(1:ncol(dfseurat2),Ncell,replace = F)]
dfseurat2$disease = "T21"

########################################3

# combine H and T21 data:
dfseurat = merge(dfseurat,dfseurat2)
dfseurat = dfseurat[,!(dfseurat@meta.data$subclust_v6 %in% c("Unknown","No markers"))]

# remove empty peaks if necessary
rowcount <- rowSums(dfseurat@assays$ATAC@counts > 0)
rowKeep = rowcount > 0
dfseurat@assays$ATAC@counts = (dfseurat@assays$ATAC@counts[rowKeep, ])
dfseurat@assays$ATAC@ranges = (dfseurat@assays$ATAC@ranges[rowKeep,])

########################

saveRDS(ATAC_df,"/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_ds_v2/Multiome.RNA_ATAC.all.scavenge_input.rds")


########################################3

# transform seurat to summarizedexperiment
SE_data <- SummarizedExperiment(assays = list(counts = dfseurat@assays$ATAC@counts),
                             rowData = dfseurat@assays$ATAC@ranges, 
                             colData = DataFrame(names = colnames(dfseurat@assays$ATAC),dataset=dfseurat@meta.data$dataset))

SE_data <- addGCBias(SE_data, genome = BSgenome.Hsapiens.UCSC.hg38)

# preprocess gchromvar
SE_data_bg <- getBackgroundPeaks(SE_data, niterations=200)

# preprocess for dimensional reduction, visualization, kNN, and network propagation
peak_by_cell_mat <- assay(SE_data)
peak_by_cell_mat[peak_by_cell_mat > 1] = 1
num_cell_detected = rowSums(peak_by_cell_mat)
sum(num_cell_detected > 200)
num_cell_detected[num_cell_detected > 200]

print("Binarizing matrix...")
peak_by_cell_mat@x[peak_by_cell_mat@x >= 1] <- 1

print("Identifying number of cells in each peak...")
featurecounts <- rowSums(x = peak_by_cell_mat)
min.cutoff=200
peaks_to_keep = featurecounts >= min.cutoff
print(paste0("Keeping ",sum(peaks_to_keep)," peaks and removing ",sum(!peaks_to_keep)," peaks for LSI, UMAP, and KNN analyses..."))
peak_by_cell_mat = peak_by_cell_mat[peaks_to_keep,]
tfidf_mat <- tfidf(bmat=peak_by_cell_mat, mat_binary=FALSE, TF=TRUE, log_TF=TRUE) # mat_binary does: bmat@x[bmat@x >= 1] <- 1 to binarize the matrix
lsi_mat <- do_lsi(tfidf_mat, dims=30)
harmony_mat <- HarmonyMatrix(lsi_mat, SE_data@colData, "dataset", do_pca=FALSE)
set.seed(03191995)
umap_mat = uwot::umap(harmony_mat[,-1],
                min_dist=0.3,
                metric="cosine",
                n_neighbors = 30,
                nn_method = "annoy",
                n_threads = numThreads,
                n_trees = 50
                ) # remove first dim in scATAC data
tmp = data.frame(names=rownames(dfseurat@meta.data),
                 cell=dfseurat@meta.data$cell,
                 disease=dfseurat@meta.data$disease,
                 subclust_v6=dfseurat@meta.data$subclust_v6,
                 UMAP_1=umap_mat[,1],
                 UMAP_2=umap_mat[,2])
fwrite(tmp,"/oak/stanford/groups/smontgom/amarder/tmp/umap.txt",
       row.names = F,col.names = T,na = "NA",sep = "\t",quote = F)
# mutualknn30 <- getmutualknn(lsi_mat, 30)
mutualknn30 <- getmutualknn(harmony_mat[,-1], 30)
saveRDS(ATAC_df,"/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_ds_v2/Multiome.RNA_ATAC.all.scavenge_input.rds")
saveRDS(mutualknn30,file="/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_ds_v2/Multiome.RNA_ATAC.all.scavenge_input.knn.rds")

ATAC_df = readRDS("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_ds_v2/Multiome.RNA_ATAC.all.scavenge_input.rds")
mutualknn30 = readRDS("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_ds_v2/Multiome.RNA_ATAC.all.scavenge_input.knn.rds")
ATAC_df@meta.data$disease2 = as.numeric(ATAC_df@meta.data$disease=="T21")
ATAC_df@meta.data$knn_disease_score = unlist(mclapply(1:ncol(mutualknn30),function(i) mean(ATAC_df@meta.data$disease2[mutualknn30[,i]==1]),mc.cores = numThreads))

aggregate(ATAC_df@meta.data$knn_disease_score,by=list(ATAC_df@meta.data$disease),mean)
saveRDS(ATAC_df,"/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_ds_v2/Multiome.RNA_ATAC.all.scavenge_input.add_knn.rds")

ATAC_df = readRDS("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_ds_v2/Multiome.RNA_ATAC.all.scavenge_input.add_knn.rds")

umap_mat = ATAC_df@reductions$umap@cell.embeddings
tmp = data.frame(names=rownames(ATAC_df@meta.data),
                 cell=ATAC_df@meta.data$cell,
                 dataset=ATAC_df@meta.data$dataset,
                 disease=ATAC_df@meta.data$disease,
                 # subclust_v6=dfseurat@meta.data$subclust_v6,
                 UMAP_1=umap_mat[,1],
                 UMAP_2=umap_mat[,2],
                 knn=ATAC_df@meta.data$knn_disease_score)
fwrite(tmp,"/oak/stanford/groups/smontgom/amarder/tmp/knn.txt",
       row.names = F,col.names = T,na = "NA",sep = "\t",quote = F)

library(viridis)
p2 <- ggplot(data=tmp[sample(1:nrow(tmp),nrow(tmp),replace = F),], aes(UMAP_1, UMAP_2)) + 
  geom_point(aes(color=knn),size=rel(0.5), na.rm = TRUE, alpha = 0.6) +
  scale_alpha()+
  theme_bw() + 
  theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2") + 
  scale_color_viridis(begin=1,end=0) +
  geom_richtext(fill = NA, label.color = NA, # remove background and outline
                    label.padding = grid::unit(rep(0, 4), "pt")) # remove padding

f.plot = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scavenge_plots/umap.knn.pdf")
pdf(f.plot,width=9,height=7)
print(p2)
dev.off()

unlist(lapply(1:ncol(mutualknn30),function(i) mean(ATAC_df@meta.data$disease2[mutualknn30[,i]])))

# 
p2 <- ggplot(data=tmp, aes(UMAP_1, UMAP_2)) + 
  geom_point(aes(color=disease),size=1, na.rm = TRUE, alpha = 0.6) +
  scale_alpha()+
  theme_bw() + 
  theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")

f.plot = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scavenge_plots/umap.disease.2.pdf")
pdf(f.plot,width=9,height=7)
print(p2)
dev.off()


# add trait
# use monocytes from example file:
# trait_file = "/home/amarder/R/x86_64-pc-linux-gnu-library/4.1/SCAVENGE/extdata/mono.PP001.hg38.bed"
# trait_file = "/oak/stanford/groups/smontgom/amarder/t21-proj/scripts/ALL_finemap.hg38.bed"

fileDir="/oak/stanford/groups/smontgom/amarder/t21_multiome/output/finemap/finemappedtraits_hg19/hg38"
flst = list.files(fileDir)
flst = list.files(fileDir,pattern = "hg38",include.dirs = FALSE)

trait_mat.save.all = list()
for (k in 1:length(flst)) {
  
  f = flst[k]
  trait_file = paste0(fileDir,"/",f)
  traitName = gsub("\\..*","",f)
  
  print(k)
  print(traitName)
  
  trait_import <- importBedScore(rowRanges(SE_data), trait_file, colidx=5)
  SE_data_DEV <- computeWeightedDeviations(SE_data, trait_import, background_peaks =
                                             SE_data_bg)
  z_score_mat <- data.frame(colData(SE_data), z_score=t(assays(SE_data_DEV)[["z"]]) %>% c)
  
  seed_idx <- seedindex(z_score_mat$z_score, 0.05)
  scale_factor <- cal_scalefactor(z_score=z_score_mat$z_score, 0.01)
  
  # Network propagation
  np_score <- randomWalk_sparse(intM=mutualknn30, rownames(mutualknn30)[seed_idx], gamma=0.05)
  
  # remove cells not connected to other cells
  omit_idx <- np_score==0
  sum(omit_idx)
  
  mutualknn30.2 <- mutualknn30[!omit_idx, !omit_idx]
  np_score <- np_score[!omit_idx]
  TRS <- np_score %>% capOutlierQuantile(., 0.95) %>% max_min_scale
  TRS <- TRS * scale_factor
  trait_mat <- data.frame(z_score_mat[!omit_idx, ], seed_idx[!omit_idx], np_score, TRS)
  trait_mat = merge(data.frame(names=rownames(dfseurat@meta.data),
                               cell=dfseurat@meta.data$cell,
                               disease=dfseurat@meta.data$disease,
                               subclust_v6=dfseurat@meta.data$subclust_v6,
                               UMAP_1=umap_mat[,1],
                               UMAP_2=umap_mat[,2]),
                    trait_mat,
                    by='names')
  
  # traitName = "mono"
  f.out = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scavenge/",traitName,".txt")
  trait_mat.save = trait_mat[,c("names","cell","dataset","disease","subclust_v6","UMAP_1","UMAP_2","TRS")]
  library(data.table)
  fwrite(trait_mat.save,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
  
  p2 <- ggplot(data=trait_mat, aes(UMAP_1, UMAP_2, color=TRS)) + 
    geom_point(size=1, na.rm = TRUE, alpha = 0.6) +
    scale_color_gradientn(colors = viridis) + 
    scale_alpha()+
    theme_bw() + 
    theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
  f.plot = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scavenge_plots/",traitName,".pdf")
  pdf(f.plot,width=7,height=7)
  print(p2)
  dev.off()
  
  # colnames(trait_mat.save)[colnames(trait_mat.save)=="TRS"] <- traitName
  # 
  # trait_mat.save.all[[k]] = trait_mat.save

}



##########################

aggregate(trait_mat$TRS,by=list(trait_mat$disease,trait_mat$subclust_v6),median)
wilcox.test(trait_mat$TRS[trait_mat$subclust_v6=="Granulocyte progenitors" & trait_mat$disease=="H"],
            trait_mat$TRS[trait_mat$subclust_v6=="Granulocyte progenitors" & trait_mat$disease=="T21"])
wilcox.test(trait_mat$TRS[trait_mat$disease=="H"],
            trait_mat$TRS[trait_mat$disease=="T21"])
aggregate(trait_mat$TRS,by=list(trait_mat$disease),mean)

wilcox.test(trait_mat$TRS[trait_mat$subclust_v6=="Neutrophils" & trait_mat$disease=="H"],
            trait_mat$TRS[trait_mat$subclust_v6=="Neutrophils" & trait_mat$disease=="T21"])

length(trait_mat$TRS[trait_mat$clust=="Granulocyte progenitors" & trait_mat$disease=="H"])
length(trait_mat$TRS[trait_mat$clust=="Granulocyte progenitors" & trait_mat$disease=="T21"])

trait_permu <- get_sigcell_simple(knn_sparse_mat=mutualknn30.2, seed_idx=trait_mat$seed_idx,
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

trait_mat2 %>%
  group_by(clust,disease) %>%
  summarise(enriched_cell=mean(true_cell_top_idx),TRS=median(TRS)) %>% print(n=Inf)

traitName = "mono"
f.out = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scavenge/",traitName,".txt")
trait_mat.save = trait_mat[,c("names","cell","dataset","disease","subclust_v6","TRS")]
library(data.table)
fwrite(trait_mat,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

