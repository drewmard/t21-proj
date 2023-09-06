# srun --account=default --partition=interactive --time=24:00:00 --mem=128G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash

library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(Seurat)
library(Signac)
library(GenomicRanges)
library(future)
library(data.table)
library(harmony)
library("JASPAR2020")
library("TFBSTools")
library(parallel)

DATASET="DS_Multiome_h"
dir="/oak/stanford/groups/smontgom/amarder/neuro-variants"

dfrna <- readRDS(paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET,"/RNA_FindClusters.rds"))
dfatac <- readRDS(paste0(dir,"/output/data/",DATASET,"/round2_FindClusters.GeneActivity.rds"))

if (DATASET=="DS_Multiome_ds") {
  mapping = fread("/oak/stanford/groups/smontgom/amarder/t21-proj/scripts/subcluster/rna_t21_subclust_v1.csv",data.table = F,stringsAsFactors = F)
  
  cellNames = rownames(dfrna@meta.data)
  dfrna@meta.data$i = 1:nrow(dfrna@meta.data)
  meta <- dfrna@meta.data
  meta = merge(meta,mapping,by.x="seurat_clusters",by.y="ClustNum")
  meta <- meta[order(meta$i),]
  rownames(meta) <- cellNames
  meta[is.na(meta$Name),"Name"] <- "No markers"
  y=as.character(unique(meta$Name))
  meta$Name <- factor(meta$Name,levels=c(cellOrder[cellOrder %in% y],y[!(y %in% cellOrder)]))
  dfrna@meta.data <- meta
  
} else if (DATASET=="DS_Multiome_h") {
  f <- "/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h/RNA_meta_v1.txt"
  meta <- read.table(f,stringsAsFactors = F,header=T,sep='\t',row.names=1)
  f=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/scripts/subcluster/rna_subclust_v1_mapping_v2.csv")
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
  meta2$Name <- as.character(meta2$Final)
  dfrna@meta.data <- meta2
  fwrite(meta2,"/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h/RNA_meta_v2.txt",quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
}

rnacluster <- dfrna$Name
dfatac <- AddMetaData(dfatac,metadata = rnacluster,col.name="rnacluster")
dfatac

dataf = dfatac@meta.data[,c("seurat_clusters","rnacluster")]
tab = aggregate(data.frame(value=dataf$rnacluster),dataf,length)

# Stacked + percent
f.out <- paste0(dir,"/output/data/",DATASET,"/stacked_barplot.x_RNA.fill_ATAC.pdf")
pdf(f.out,width=14,height=9)
g <- ggplot(tab, aes(fill=seurat_clusters, y=value, x=rnacluster)) + 
  geom_bar(position="fill", stat="identity") +
  labs(y = "Percentage (%)",x="RNA Label",fill="ATAC Cluster") +
  theme_bw() +
  theme(panel.grid=element_blank(),axis.text.x = element_text(angle=30,hjust=1))
print(g)
dev.off()

f.out <- paste0(dir,"/output/data/",DATASET,"/stacked_barplot.x_ATAC.fill_RNA.pdf")
pdf(f.out,width=14,height=9)
g <- ggplot(tab, aes(x=seurat_clusters, y=value, fill=rnacluster)) + 
  geom_bar(position="fill", stat="identity") +
  labs(y = "Percentage (%)",fill="RNA Label",x="ATAC Cluster") +
  theme_bw() +
  theme(panel.grid=element_blank(),axis.text.x = element_text(angle=30,hjust=1,vjust=0.5))
print(g)
dev.off()


