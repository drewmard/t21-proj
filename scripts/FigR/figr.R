library(Seurat)
library(Signac)
library(data.table)
library(ggplot2)
library(dplyr)
library(FNN)
library(chromVAR)
library(doParallel)
library(BuenColors)


dir="/oak/stanford/groups/smontgom/amarder/neuro-variants"
DATASET="DS_Multiome_h"

dfrna <- readRDS(paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET,"/RNA_FindClusters.rds"))
meta2 = fread(paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET,"/RNA_meta_v3.txt"),data.table = F,stringsAsFactors = F)
rownames(meta2) <- meta2[,1]; meta2 <- meta2[,-1]
dfrna@meta.data <- meta2; #rm(meta); rm(meta2)
# dfatac <- readRDS(paste0(dir,"/output/data/",DATASET,"/round2_FindClusters.rds"))
rnacluster <- dfrna$Final
Idents(dfrna) <- "Final"

f.sce="/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h/round2_FindClusters.sce.rds"
dfatac <- readRDS(f.sce)

y = strsplit(rownames(dfatac),"-")
chr=unlist(lapply(y,function(x) x[[1]]))
start=unlist(lapply(y,function(x) x[[2]]))
end=unlist(lapply(y,function(x) x[[3]]))
grangedf <- data.frame(chr,start,end)
grangedf <- makeGRangesFromDataFrame(grangedf)
rowRanges(dfatac) <- grangedf

dfrna <- dfrna[,QC_cellsFailATAC]

# set.seed(123)
# cellsToKeep <- sample(colnames(dfatac),size = 3000,replace = FALSE)
cellsToKeep <- colnames(dfrna)[dfrna$Final=="HSC/MPPs"]

# which(colnames(dfatac)!=colnames(dfrna))
ATAC.se <- dfatac[,cellsToKeep]
RNAmat <- dfrna[,cellsToKeep]

# Remove genes with zero expression across all cells
RNAmat <- RNAmat[Matrix::rowSums(RNAmat)!=0,]

# source("/oak/stanford/groups/smontgom/amarder/bin/FigR/R/")
source("/oak/stanford/groups/smontgom/amarder/bin/FigR/R/utils.R")
source("/oak/stanford/groups/smontgom/amarder/bin/FigR/R/FigR.R")
source("/oak/stanford/groups/smontgom/amarder/bin/FigR/R/DORCs.R")
source("/oak/stanford/groups/smontgom/amarder/bin/FigR/R/cellPairing.R")

source("/oak/stanford/groups/smontgom/amarder/t21-proj/scripts/FigR/runGenePeakcorr.R")
source("/oak/stanford/groups/smontgom/amarder/t21-proj/scripts/FigR/PeakGeneCor.R")

hg38TSSRanges = readRDS("/oak/stanford/groups/smontgom/amarder/bin/FigR/data/hg38TSSRanges.RDS")

cisCorr <- runGenePeakcorr(ATACdf = ATAC.se,
                                 RNAdf = RNAmat@assays$RNA@data,
                                 genome = "hg38", # One of hg19, mm10 or hg38 
                                 nCores = 8,
                                 p.cut = NULL, # Set this to NULL and we can filter later
                                 n_bg = 10)
saveRDS(cisCorr,file="/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h/HSC_MPPs.runGenePeakcorr.rds")
# head(cisCorr)
# 
# cisCorr.filt <- cisCorr %>% filter(pvalZ <= 0.05)
# 
# dorcGenes <- dorcJPlot(dorcTab = cisCorr.filt,
#                        cutoff = 10, # No. sig peaks needed to be called a DORC
#                        labelTop = 20,
#                        returnGeneList = TRUE, # Set this to FALSE for just the plot
#                        force=2)
# 
# # Unfiltered
# numDorcs <- cisCorr.filt %>% group_by(Gene) %>% tally() %>% arrange(desc(n))
# numDorcs
# 
# 
# dorcMat <- getDORCScores(ATAC.se = ATAC.se, # Has to be same SE as used in previous step
#                          dorcTab = cisCorr.filt,
#                          geneList = dorcGenes,
#                          nCores = 4)
# 
# 
# 
# 
# 
# 
# 
