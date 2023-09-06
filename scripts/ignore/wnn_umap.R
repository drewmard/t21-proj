library(ggplot2)
library(Seurat)
library(Signac)
library(data.table)
library(harmony)
library("JASPAR2020")
library("TFBSTools")
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)

####################

# args = commandArgs(trailingOnly=TRUE)
# dir=args[1]
# # start=args[2]
# # end=args[3]
# DATASET=args[4]

dir="/oak/stanford/groups/smontgom/amarder/neuro-variants"
# DATASET="DS_Multiome_ds"
DATASET="DS_Multiome_h"

####################################################################

print("Reading ATAC data...")
# f <- paste0(dir,"/output/data/",DATASET,"/round2_FindClusters.rds")
f <- paste0(dir,"/output/data/","DS_Multiome_all","/round2_FindClusters.rds")
ATAC_df <- readRDS(file = f)
y = strsplit(rownames(ATAC_df@meta.data),"_")
ATAC_df@meta.data$cell <- unlist(lapply(y,function(x) x[1]))
ATAC_df@meta.data$dataNum <- unlist(lapply(y,function(x) x[2]))
ATAC_df@meta.data$cell_dataset <- paste(ATAC_df@meta.data$cell,ATAC_df@meta.data$dataset)

f <- paste0(dir,"/output/data/",DATASET,"_v2","/RNA_FindClusters.rds")
print(paste0("Loading RNA data: ",f,"..."))
RNA_df <- readRDS(file = f)
smallRNA_meta.path=paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET,"_v2","/meta.091522_subclust.txt")
meta <- read.table(smallRNA_meta.path,stringsAsFactors = F,header=T,sep='\t',row.names=1)
RNA_df@meta.data <- meta

RNA_df@meta.data$cell_dataset <- paste(RNA_df@meta.data$cell,RNA_df@meta.data$dataset)

ATAC_df.meta.sub = subset(ATAC_df@meta.data,cell_dataset %in% RNA_df$cell_dataset)
RNA_df.meta.sub = subset(RNA_df@meta.data,cell_dataset %in% ATAC_df$cell_dataset)
which(ATAC_df.meta.sub$cell_dataset != RNA_df.meta.sub$cell_dataset)

ATAC_df.sub = subset(ATAC_df,cell_dataset %in% RNA_df$cell_dataset)
RNA_df.sub = subset(RNA_df,cell_dataset %in% ATAC_df$cell_dataset)
which(ATAC_df.sub$cell_dataset != RNA_df.sub$cell_dataset)
mean(ATAC_df.sub$cell_dataset == RNA_df.sub$cell_dataset)

# colnames(ATAC_df@meta.data)[colnames(ATAC_df@meta.data)=="seurat_clusters"] <- "ATAC_clusters"
# colnames(RNA_df@meta.data)[colnames(RNA_df@meta.data)=="seurat_clusters"] <- "RNA_clusters"

if (DATASET=="DS_Multiome_ds") {
  # colnames(ATAC_df.sub) <- paste(ATAC_df.sub@meta.data$cell,as.numeric(ATAC_df.sub@meta.data$dataNum) - 6,sep="_")
  # rownames(ATAC_df.sub@meta.data) <- paste(ATAC_df.sub@meta.data$cell,as.numeric(ATAC_df.sub@meta.data$dataNum) - 6,sep="_")
  ATAC_df.sub = RenameCells(ATAC_df.sub,old.names=Cells(ATAC_df.sub),new.names=paste(ATAC_df.sub@meta.data$cell,as.numeric(ATAC_df.sub@meta.data$dataNum) - 6,sep="_"))
}

# 
dfcombined <- ATAC_df.sub
dfcombined[["RNA"]] <- RNA_df.sub[["RNA"]] # CreateAssayObject(data=x)
dfcombined[["RNA_harmony"]] <- RNA_df.sub[["harmony"]]
# which(rownames(RNA_df@meta.data)!=rownames(ATAC_df@meta.data))
dfcombined@meta.data <- cbind(dfcombined@meta.data,
                              RNA_df.sub@meta.data[,c("nCount_RNA","nFeature_RNA","percent.mt","subclust_v5")]) # "nCount_RNA","nFeature_RNA" already in ATAC_df

print("Creating UMAP with seurat clusters...")
# p1 <- DimPlot(dfcombined, group.by = "seurat_clusters")
# DATASET="DS_Multiome_h"
f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/ATAC_umap.RNA_clusters.png")
png(filename=f.out,width = 5000,height=4000,res=500)
p1 <- DimPlot(dfcombined, group.by = "subclust_v5",reduction="umap") & theme(plot.title = element_text(hjust = 0.5)) 
print(LabelClusters(plot = p1, id = "subclust_v5"))
dev.off()

dfcombined <- FindMultiModalNeighbors(dfcombined, reduction.list = list("harmony", "RNA_harmony"), dims.list = list(2:50, 1:50))
dfcombined <- RunUMAP(dfcombined, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

print("Creating UMAP with seurat clusters...")
# p1 <- DimPlot(dfcombined, group.by = "seurat_clusters")
f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/WNN_umap.RNA_clusters.png")
png(filename=f.out,width = 5000,height=4000,res=500)
p1 <- DimPlot(dfcombined, group.by = "subclust_v5",reduction="wnn.umap") & theme(plot.title = element_text(hjust = 0.5)) 
print(LabelClusters(plot = p1, id = "subclust_v5"))
dev.off()

# Save data:
f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/WNN_ATAC_RNA.rds")
print(paste0("Saving data: ",f.out,"..."))
saveRDS(dfcombined,file = f.out)

# Stacked + percent

dataf = dfcombined@meta.data[,c("seurat_clusters","subclust_v5")]
tab = aggregate(data.frame(value=dataf$subclust_v5),dataf,length)

f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/stacked_barplot.x_RNA.fill_ATAC.pdf")
pdf(f.out,width=14,height=9)
g <- ggplot(tab, aes(fill=seurat_clusters, y=value, x=subclust_v5)) + 
  geom_bar(position="fill", stat="identity",col='black') +
  labs(y = "Percentage (%)",x="RNA Label",fill="ATAC Cluster") +
  theme_bw() +
  theme(panel.grid=element_blank(),axis.text.x = element_text(angle=30,hjust=1))
print(g)
dev.off()

f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/stacked_barplot.x_ATAC.fill_RNA.pdf")
pdf(f.out,width=14,height=9)
g <- ggplot(tab, aes(x=seurat_clusters, y=value, fill=subclust_v5)) + 
  geom_bar(position="fill", stat="identity",col='black') +
  labs(y = "Percentage (%)",fill="RNA Label",x="ATAC Cluster") +
  theme_bw() +
  theme(panel.grid=element_blank(),axis.text.x = element_text(angle=30,hjust=1))
print(g)
dev.off()





