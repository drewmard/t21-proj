library(Seurat)
library(Signac)
library(data.table)
library(ggplot2)
library(dplyr)
library(FNN)
library(chromVAR)
library(doParallel)
library(BuenColors)
library(Matrix)
library(pbmcapply)
library(SummarizedExperiment)

NumBGPerm1 = 100
NumBGPerm2 = 50
dir="/oak/stanford/groups/smontgom/amarder/neuro-variants"
DATASET="DS_Multiome_h"

dir.create(paste0(dir,"/output/data/",DATASET,"/FigR/"))

print("Reading RNA...")

dfrna <- readRDS(paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET,"/RNA_FindClusters.rds"))
meta2 = fread(paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET,"/RNA_meta_v3.txt"),data.table = F,stringsAsFactors = F)
rownames(meta2) <- meta2[,1]; meta2 <- meta2[,-1]
dfrna@meta.data <- meta2; #rm(meta); rm(meta2)
# dfatac <- readRDS(paste0(dir,"/output/data/",DATASET,"/round2_FindClusters.rds"))
rnacluster <- dfrna$Final
Idents(dfrna) <- "Final"

print("Reading ATAC...")

# f.sce="/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h/round2_FindClusters.sce.rds"
f.sce="/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h/RNA_cluster.ATAC.sce.rds"
dfatac <- readRDS(f.sce)

print("ATAC granges...")

y = strsplit(rownames(dfatac),"-")
chr=unlist(lapply(y,function(x) x[[1]]))
start=unlist(lapply(y,function(x) x[[2]]))
end=unlist(lapply(y,function(x) x[[3]]))
grangedf <- data.frame(chr,start,end)
grangedf <- makeGRangesFromDataFrame(grangedf)
rowRanges(dfatac) <- grangedf

print("data subset based on QC and HSCs...")

dfrna <- dfrna[,QC_cellsFailATAC]

# set.seed(123)
# cellsToKeep <- sample(colnames(dfatac),size = 3000,replace = FALSE)
cellsToKeep <- colnames(dfrna)[dfrna$Final=="HSC/MPPs"]

# which(colnames(dfatac)!=colnames(dfrna))
ATAC.se <- dfatac[,cellsToKeep]
RNAmat <- dfrna[,cellsToKeep]

# Remove genes with zero expression across all cells
RNAmat <- RNAmat[Matrix::rowSums(RNAmat)!=0,]

# Derive cell kNN using this
print("derive cell kNN...")
set.seed(123)
cellkNN <- get.knn(reducedDim(ATAC.se,"HARMONY"),k = 30)$nn.index
rownames(cellkNN) <- cellsToKeep


# source("/oak/stanford/groups/smontgom/amarder/bin/FigR/R/")
source("/oak/stanford/groups/smontgom/amarder/bin/FigR/R/utils.R")
source("/oak/stanford/groups/smontgom/amarder/bin/FigR/R/FigR.R")
source("/oak/stanford/groups/smontgom/amarder/bin/FigR/R/DORCs.R")
source("/oak/stanford/groups/smontgom/amarder/bin/FigR/R/cellPairing.R")

source("/oak/stanford/groups/smontgom/amarder/t21-proj/scripts/FigR/runGenePeakcorr.R")
source("/oak/stanford/groups/smontgom/amarder/t21-proj/scripts/FigR/PeakGeneCor.R")
source("/oak/stanford/groups/smontgom/amarder/t21-proj/scripts/FigR/getDORCScores.R")
source("/oak/stanford/groups/smontgom/amarder/t21-proj/scripts/FigR/runFigRGRN.R")

print("hg38TSSRanges...")

hg38TSSRanges = readRDS("/oak/stanford/groups/smontgom/amarder/bin/FigR/data/hg38TSSRanges.RDS")

# ATAC.se <- ATAC.se[1:10000,]
print("runGenePeakcorr ...")
cisCorr <- runGenePeakcorr(ATACdf = ATAC.se,
                                 RNAdf = RNAmat@assays$RNA@data,
                                 genome = "hg38", # One of hg19, mm10 or hg38 
                                 nCores = 8,
                                 p.cut = NULL, # Set this to NULL and we can filter later
                                 n_bg = NumBGPerm)

f.out=paste0(dir,"/output/data/",DATASET,"/FigR/HSC_MPPs.runGenePeakcorr.rds")
saveRDS(cisCorr,file=f.out)
head(cisCorr)

cisCorr.filt <- cisCorr %>% filter(pvalZ <= 0.05)

f.out=paste0(dir,"/output/data/",DATASET,"/FigR/HSC_MPPs.dorcJPlot.pdf")
pdf(f.out,width = 6,height=6)
print(dorcGenes <- dorcJPlot(dorcTab = cisCorr.filt,
                       cutoff = 10, # No. sig peaks needed to be called a DORC
                       labelTop = 20,
                       returnGeneList = TRUE, # Set this to FALSE for just the plot
                       force=2))
dev.off()
print(f.out)

# Unfiltered
numDorcs <- cisCorr.filt %>% group_by(Gene) %>% tally() %>% arrange(desc(n))
numDorcs

print("getDORCScores ...")
dorcMat <- getDORCScores(ATAC.se = ATAC.se, # Has to be same SE as used in previous step
                         dorcTab = cisCorr.filt,
                         geneList = dorcGenes,
                         nCores = 4)

f.out=paste0(dir,"/output/data/",DATASET,"/FigR/HSC_MPPs.getDORCScores.rds")
saveRDS(dorcMat,file=f.out)

# Smooth dorc scores using cell KNNs (k=30)
print("dorcMat.s ...")
dorcMat.s <- smoothScoresNN(NNmat = cellkNN[,1:30],mat = dorcMat,nCores = 4)
f.out=paste0(dir,"/output/data/",DATASET,"/FigR/HSC_MPPs.dorcMat.s.rds")
saveRDS(dorcMat.s,file=f.out)

# Smooth RNA using cell KNNs
# This takes longer since it's all genes
print("RNAmat.s ...")
RNAmat.s <- smoothScoresNN(NNmat = cellkNN[,1:30],mat = RNAmat@assays$RNA@data,nCores = 16)
f.out=paste0(dir,"/output/data/",DATASET,"/FigR/HSC_MPPs.RNAmat.s.rds")
saveRDS(RNAmat.s,file=f.out)

# # Visualize on pre-computed UMAP
# umap.d <- as.data.frame(colData(ATAC.se)[,c("UMAP1","UMAP2")])
# 
# # DORC score for Dlx3
# dorcg <- plotMarker2D(umap.d,dorcMat.s,markers = c("Dlx3"),maxCutoff = "q0.99",colorPalette = "brewer_heat") + ggtitle("Dlx3 DORC")
# 
# # RNA for Dlx3
# rnag <- plotMarker2D(umap.d,RNAmat.s,markers = c("Dlx3"),maxCutoff = "q0.99",colorPalette = "brewer_purple") + ggtitle("Dlx3 RNA")
# 
# library(patchwork)
# dorcg + rnag

print("figR.d ...")
numCores_to_use <- parallel::detectCores()
print("Using # cores:")
print(numCores_to_use)
figR.d <- runFigRGRN(ATAC.se = ATAC.se, # Must be the same input as used in runGenePeakcorr()
                     dorcTab = cisCorr.filt, # Filtered peak-gene associations
                     dorcK=min(30,nrow(dorcMat.s)-1),
                     genome = "hg38",
                     dorcMat = dorcMat.s,
                     rnaMat = RNAmat.s, 
                     nCores = numCores_to_use,
                     n_bg = NumBGPerm)
print("Saving figR.d ...")

f.out=paste0(dir,"/output/data/",DATASET,"/FigR/HSC_MPPs.figR.rds")
saveRDS(figR.d,file=f.out)

print("Visualizations! ...")
print("#1 ...")
f.out=paste0(dir,"/output/data/",DATASET,"/FigR/HSC_MPPs.TF_DORC_scatter.pdf")
pdf(f.out,width = 6,height=6)
figR.d %>% 
  ggplot(aes(Corr.log10P,Enrichment.log10P,color=Score)) + 
  ggrastr::geom_point_rast(size=0.01,shape=16) + 
  theme_classic() + 
  scale_color_gradientn(colours = jdb_palette("solar_extra"),limits=c(-3,3),oob = scales::squish,breaks=scales::breaks_pretty(n=3))
dev.off()

print("#2 ...")
f.out=paste0(dir,"/output/data/",DATASET,"/FigR/HSC_MPPs.TF_DORC_rank.pdf")
pdf(f.out,width = 6,height=6)
rankDrivers(figR.d,rankBy = "meanScore",interactive = FALSE)
dev.off()




