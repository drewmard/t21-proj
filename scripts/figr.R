library(Seurat)
library(Signac)

f.sce="/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h/round2_FindClusters.sce.rds"

library(Seurat)
# library(Signac)
library(data.table)
library(ggplot2)

dir="/oak/stanford/groups/smontgom/amarder/neuro-variants"
DATASET="DS_Multiome_h"

dfrna <- readRDS("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h/RNA_FindClusters.rds")

f <- "/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h/RNA_meta_v1.txt"
meta <- read.table(f,stringsAsFactors = F,header=T,sep='\t',row.names=1)
f=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/scripts/subcluster/rna_h_subclust_v1_mapping_v3.csv")
mapping <- fread(f,stringsAsFactors = F,data.table = F)
mapping$ClustName <- paste0(mapping$Number,":",mapping$Name,",",mapping$Mclust)
meta2 <- meta
meta2$id  <- 1:nrow(meta2)
meta2 <- merge(meta2,mapping[,c("Number","ClustName","Final")],by.x="subclust_v1",by.y="Number")
meta2 <- meta2[order(meta2$id), ]
meta2 <- meta2[,colnames(meta2)!="id"]
rownames(meta2) <- rownames(meta)
# meta2 <- meta2[,-1]
meta2$subclust_v1_name <- meta2$ClustName
meta2 <- meta2[,colnames(meta2)!="ClustName"]
meta2$Final[is.na(meta2$Final)] <- "No markers"

s="HSC/MPPs, Progenitors, MEMPs, Early erythroid, Late erythroid, Mast cells, Megakaryocytes, Granulocyte progenitors, Neutrophils, Monocyte progenitors, Inflammatory macrophages, Kupffer cells, NK cells, Lymphoid progenitors, Pre pro B cells, Pro B cells, B cells, pDCs, cDCs, Hepatocytes, LSECs, No markers"
cellOrder=unlist(strsplit(s,", "))
cellOrder=rev(cellOrder)

y=as.character(unique(meta2$Final))
cellLevels=c(cellOrder[cellOrder %in% y],y[!(y %in% cellOrder)])
cellLevels <- rev(cellLevels)
meta2$Final <- factor(as.character(meta2$Final),levels=cellLevels)

dfrna@meta.data <- meta2

f.sce="/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h/round2_FindClusters.sce.rds"
dfatac <- readRDS(f.sce)

dfrna <- dfrna[,colnames(dfrna) %in% colnames(dfatac)]

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
hg38TSSRanges = readRDS("/oak/stanford/groups/smontgom/amarder/bin/FigR/data/hg38TSSRanges.RDS")

cisCorr <- runGenePeakcorr(ATAC.se = ATAC.se,
                                 RNAmat = RNAmat@assays$RNA@data,
                                 genome = "hg38", # One of hg19, mm10 or hg38 
                                 nCores = 8,
                                 p.cut = NULL, # Set this to NULL and we can filter later
                                 n_bg = 100)
saveRDS(cisCorr,file="/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h/HSC_MPPs.runGenePeakcorr.rds")
head(cisCorr)

cisCorr.filt <- cisCorr %>% filter(pvalZ <= 0.05)

dorcGenes <- dorcJPlot(dorcTab = cisCorr.filt,
                       cutoff = 10, # No. sig peaks needed to be called a DORC
                       labelTop = 20,
                       returnGeneList = TRUE, # Set this to FALSE for just the plot
                       force=2)

# Unfiltered
numDorcs <- cisCorr.filt %>% group_by(Gene) %>% tally() %>% arrange(desc(n))
numDorcs


dorcMat <- getDORCScores(ATAC.se = ATAC.se, # Has to be same SE as used in previous step
                         dorcTab = cisCorr.filt,
                         geneList = dorcGenes,
                         nCores = 4)







