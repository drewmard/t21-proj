library(Seurat)
df = readRDS("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_all/RNA_FindClusters.rds")

# library(data.table)
# f.out = "/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_all/HVG.txt"
# fwrite(data.frame(VariableFeatures(df)),f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
# mean(VariableFeatures(df) %in% VariableFeatures(dfrna2))

# check batch correction:
print("Creating split-UMAP with dataset-based coloring (round 2)...")
p1 <- DimPlot(df, split.by = "dataset",ncol = 3)
f.out <- paste0(dir,"/output/data/",DATASET,"/RNA_umap.dataset_integration.split.png")
png(filename=f.out,width = 7000,height=7000,res=400)
print(p1 & theme(plot.title = element_text(hjust = 0.5)) + NoLegend())
dev.off()

library(data.table)
DATASET="DS_Multiome_h"
meta2 = fread(paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET,"/RNA_meta_v3.txt"),data.table = F,stringsAsFactors = F)
rownames(meta2) <- meta2[,1]
meta2 <- meta2[,-1]

DATASET="DS_Multiome_ds"
meta3 = fread(paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET,"/RNA.meta.clust_v3.txt"),data.table = F,stringsAsFactors = F)
meta3$Name[meta3$Name=="Neutrophils,2"] <- "Granulocyte progenitors"
meta3$Name[meta3$Name=="Neutrophils,1"] <- "Neutrophils"
meta3$Name[meta3$Name=="Pro B cells,1"] <- "pre pro B cells"
meta3$Name[meta3$Name=="Pro B cells,2"] <- "pro B cells"
rownames(meta3) <- meta3[,1]
meta3 <- meta3[,-1]

meta <- df@meta.data
meta$cell <- unlist(lapply(strsplit(rownames(meta),"_"),function(x) x[1]))
meta2$cell <- unlist(lapply(strsplit(rownames(meta2),"_"),function(x) x[1]))
meta3$cell <- unlist(lapply(strsplit(rownames(meta3),"_"),function(x) x[1]))

meta2$Disease <- "H"
meta3$Disease <- "T21"

tmp1 <- meta2[,c("cell","dataset","Disease","Final")]
tmp2 <- meta3[,c("cell","dataset","Disease","Name")]
colnames(tmp2)[4] <- "Final"
meta.split <- as.data.frame(rbind(tmp1,tmp2))

meta$i <- 1:nrow(meta)
meta.mg <- merge(meta,meta.split,by=c("cell","dataset"))
meta.mg <- meta.mg[order(meta.mg$i),]

meta.mg$Final[meta.mg$Final=="Pro B cells"] <- "pro B cells"
meta.mg$Final[meta.mg$Final=="Pre pro B cells"] <- "pre pro B cells"
meta.mg$Final[meta.mg$Final=="Late erythroid cells"] <- "Late erythroid"
meta.mg$Final[meta.mg$Final=="HSC/MPPs"] <- "HSCs"
meta.mg$Final_V2 <- paste0(meta.mg$Disease,":",meta.mg$Final)

rownames(meta.mg) <- rownames(meta)
meta.mg <- meta.mg[,!(colnames(meta.mg) %in% c("i"))]
df@meta.data <- meta.mg

library(ggplot2)
dir="/oak/stanford/groups/smontgom/amarder/neuro-variants"
DATASET="DS_Multiome_all"
print("Creating UMAP with seurat clusters...")
f.out <- paste0(dir,"/output/data/",DATASET,"/RNA_umap.final.png")
png(filename=f.out,width = 5000,height=4000,res=500)
p1 <- DimPlot(df, group.by = "Final_V2") & theme(plot.title = element_text(hjust = 0.5)) 
print(LabelClusters(plot = p1, id = "Final_V2"))
dev.off()

dir="/oak/stanford/groups/smontgom/amarder/neuro-variants"
DATASET="DS_Multiome_all"
print("Creating UMAP with seurat clusters...")
f.out <- paste0(dir,"/output/data/",DATASET,"/RNA_umap.final.png")
png(filename=f.out,width = 5000,height=4000,res=500)
p1 <- DimPlot(df, group.by = "Final_V2") & theme(plot.title = element_text(hjust = 0.5)) 
print(LabelClusters(plot = p1, id = "Final_V2"))
dev.off()

cellType="HSCs"
for (cellType in unique(df$Final)) {
  
  print(cellType)
  cellNames1 <- rownames(df@meta.data[df$Final_V2 == paste0("T21:",cellType),])
  cellNames2 <- rownames(df@meta.data[df$Final_V2 == paste0("H:",cellType),])
  f.out <- paste0(dir,"/output/data/",DATASET,"/RNA_umap.final.",cellType,".png")
  png(filename=f.out,width = 5000,height=4000,res=500)
  print(DimPlot(df, label=T,  cells.highlight= list(cellNames1,cellNames2), cols.highlight = c("blue", "red"), cols= "grey") & theme(plot.title = element_text(hjust = 0.5)))
  # print(LabelClusters(plot = p1, id = "Final_V2"))
  dev.off()
  
  f.out <- paste0(dir,"/output/data/",DATASET,"/RNA_umap.final.",cellType,".split.png")
  png(filename=f.out,width = 8000,height=4000,res=500)
  print(DimPlot(df, label=T,  split.by="Disease",cells.highlight= list(cellNames1,cellNames2), cols.highlight = c("blue", "red"), cols= "grey") & theme(plot.title = element_text(hjust = 0.5))) 
  dev.off()
  
}

geneName="CD34"
f.out <- paste0(dir,"/output/data/",DATASET,"/RNA_umap.final.gene_",geneName,".split.png")
png(filename=f.out,width = 8000,height=4000,res=500)
print(FeaturePlot(df, split.by="Disease",features = geneName) & theme(plot.title = element_text(hjust = 0.5))) 
dev.off()

geneName="percent.mt"
f.out <- paste0(dir,"/output/data/",DATASET,"/RNA_umap.final.",geneName,".split.png")
png(filename=f.out,width = 8000,height=4000,res=500)
print(FeaturePlot(df, split.by="Disease",features = geneName) & theme(plot.title = element_text(hjust = 0.5))) 
dev.off()

##########

tmp <- sapply(unique(meta.mg$Final), function(x) as.integer(x == df$Final))
df@meta.data <- cbind(meta.mg,tmp)
Idents(df) <- "seurat_clusters"
print("Creating dotplot with clusters...")
f.out <- paste0(dir,"/output/data/",DATASET,"/final_vs_seurat.dotplot.pdf")
pdf(f.out,width = 8,height=12)
print(
  DotPlot(
    df,
    features=unique(meta.mg$Final),
    col.min	= 0,
    col.max = 1,
    cluster.idents = FALSE
  ) + 
    theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1))
)
dev.off()

tmp <- sapply(sort(unique(meta.mg$seurat_clusters)), function(x) as.integer(x == df$seurat_clusters))
colnames(tmp) <- paste0("seurat_",sort(unique(meta.mg$seurat_clusters)))
df@meta.data <- cbind(meta.mg,tmp)
Idents(df) <- "Final"
print("Creating dotplot with clusters...")
f.out <- paste0(dir,"/output/data/",DATASET,"/final_vs_seurat.dotplot.v2.pdf")
pdf(f.out,width = 12,height=12)
print(
  DotPlot(
    df,
    features=paste0("seurat_",sort(unique(meta.mg$seurat_clusters))),
    col.min	= 0,
    col.max = 1,
    cluster.idents = FALSE
  ) + 
    theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1))
)
dev.off()

######
# 

clusterNum=13
Idents(df) <- "Final"
print(cellType)
cellNames1 <- rownames(df@meta.data[df$seurat_clusters == clusterNum,])
f.out <- paste0(dir,"/output/data/",DATASET,"/cluster",clusterNum,".dotplot.pdf")
pdf(f.out,width = 12,height=12)
print(
  DotPlot(
    df[,cellNames1],
    features=c("CD34","SPINK2","MLLT3","TESPA1","GATA2","GATA1","HDC","CTSB"),
    col.min	= 0,
    col.max = 1,
    cluster.idents = FALSE
  ) + 
    theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1))
)
dev.off()

tab <- data.frame(table(df[,cellNames1]$Final))
f.out <- paste0(dir,"/output/data/",DATASET,"/cluster",clusterNum,".cell_count.pdf")
pdf(f.out,width = 9,height=5)
print(
  ggplot(tab,aes(x=reorder(Var1,Freq,median),y=Freq)) + 
    geom_bar(stat='identity') + theme_bw() +
    theme(plot.title = element_text(hjust=0.5),axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x="Cell type",y='Cell count',title=paste0('Cluster ',clusterNum))
)
dev.off()

tab <- table(df$Disease,df$Final)
tab
tab[1,] <- tab[1,]/sum(tab[1,])
tab[2,] <- tab[2,]/sum(tab[2,])
tab

idVar="dataset"
Idents(df) <- idVar
f.out <- paste0(dir,"/output/data/",DATASET,"/RNA.QC.Violin.",idVar,".pdf")
pdf(f.out,width=17,height=9)
print(VlnPlot(df, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) + NoLegend())
dev.off()


f.out <- paste0(dir,"/output/data/",DATASET,"/RNA_umap.final.","all_cells",".split.png")
png(filename=f.out,width = 8000,height=4000,res=500)
print(DimPlot(df, label=T,  split.by="Disease", group.by = "Final") & 
        theme(plot.title = element_text(hjust = 0.5)) + NoLegend()) 
dev.off()

