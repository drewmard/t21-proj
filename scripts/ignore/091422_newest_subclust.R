library(Seurat)
library(data.table)
library(ggplot2)
dir="/oak/stanford/groups/smontgom/amarder/neuro-variants"
# system(paste0("rm -r ",dir,"/output/data/",DATASET,"_v2","/091422_clust"))

DATASET="DS_Multiome_ds"
# DATASET="DS_Multiome_h"
df = readRDS(paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET,"_v2","/RNA_FindClusters.rds"))
f = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET,"_v2","/meta.091322_subclust.txt")
meta = fread(f,data.table = F,stringsAsFactors = F)
rownames(meta) <- meta[,1]; meta <- meta[,-1]


if (DATASET=="DS_Multiome_h") { 
  meta$subclust_v4 <- meta$subclust_v3
  meta$subclust_v4[meta$subclust_v3 %in% c("Progenitors","HSCs")] <- "HSCs"
  meta$subclust_v4[meta$subclust_v3 %in% c("pre pro B cells","pro B cells","B cells")] <- "B cells"
  meta$subclust_v4[meta$subclust_v3 %in% c("MEMPs","?")] <- "MEMPs"

} 
if (DATASET=="DS_Multiome_ds") {
  meta_pre_swap = fread("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_ds_v2/meta.090922_subclust.pre_swap.txt",data.table = F,stringsAsFactors = F)
  meta$subclust_v2 <- meta_pre_swap$subclust_v2
  meta$subclust_v2_short <- meta_pre_swap$subclust_v2_short
  i = meta_pre_swap$subclust_v2_short == "Megakaryocytes"
  j = meta_pre_swap$subclust_v2_short == "MEMPs"
  meta$subclust_v4 <- meta$subclust_v3
  meta$subclust_v4[i] <- "MEMPs"
  meta$subclust_v4[j] <- "Megakaryocytes"
  meta$subclust_v4[meta$subclust_v4 %in% c("pre pro B cells","pro B cells","B cells")] <- "B cells"
  
  # # Mast cells.3 – split in 2
  cellName = "Mast cells.3"; kclust = 2
  i = meta$subclust_v3 %in% c(cellName)
  dftmp = df[,i]
  kmeans_results <- kmeans(Embeddings(dftmp,reduction="harmony"),kclust)
  kmeans_results1 <- data.frame(kmeans_results$cluster)
  meta$subclust_v4[i] <- paste0(cellName,".",kmeans_results1[,1])

  #
  system(paste0("mkdir -p ",dir,"/output/data/",DATASET,"_v2","/091422_clust"))
  cellType = "Mast cells.3.1"
  cellNames1 <- rownames(df@meta.data[meta$subclust_v4 == cellType,])
  cellType = gsub("/","_",cellType)
  f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/091422_clust/RNA_umap.",cellType,".png")
  png(filename=f.out,width = 5000,height=4000,res=500)
  print(DimPlot(df, label=T,  cells.highlight= list(cellNames1), cols.highlight = c("red"), cols= "grey") & theme(plot.title = element_text(hjust = 0.5)))
  # print(LabelClusters(plot = p1, id = "Final_V2"))
  dev.off()

  cellType = "Mast cells.3.2"
  cellNames1 <- rownames(df@meta.data[meta$subclust_v4 == cellType,])
  cellType = gsub("/","_",cellType)
  f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/091422_clust/RNA_umap.",cellType,".png")
  png(filename=f.out,width = 5000,height=4000,res=500)
  print(DimPlot(df, label=T,  cells.highlight= list(cellNames1), cols.highlight = c("red"), cols= "grey") & theme(plot.title = element_text(hjust = 0.5)))
  # print(LabelClusters(plot = p1, id = "Final_V2"))
  dev.off()
}
df@meta.data <- meta

f.out = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET,"_v2","/meta.091422_subclust.txt")
fwrite(df@meta.data,
       f.out,
       quote=F,na = "NA",sep = "\t",row.names = T,col.names = T)

###################################################################

system(paste0("mkdir -p ",dir,"/output/data/",DATASET,"_v2","/091422_clust"))

###################################################################

Idents(df) <- "subclust_v4"
geneOrder=read.table(text=
                       "CD34, SPINK2, MLLT3
                     GATA1, KLF1, TESPA1
                     AHSP, ALAS2, HBA1, GYPA
                     GATA2, HDC, CPA3
                     ITGA2B, GP9
                     MPO, AZU1, SPI1, LYZ
                     CD14, CD68, S100A9, MNDA, FCN1, VCAN, IL1B
                     CD163, MS4A7, CTSB, STAB1
                     NKG7, PRF1, GZMA
                     IL7R, DHFR, PAX5, MME, IGLL1, IGHM, CD79A, CD19
                     JCHAIN, IRF8, CLEC4C, ILR3A
                     CD1C, CLEC4A, CLEC10A
                     KDR, CDH5, COL1A1, ALB, ACTA2
                     MKI67"
                     ,sep='\n')
geneOrder1=unlist(strsplit(gsub("\\ ","",geneOrder[,1]),","))
s="CD34, PROM1, MLLT3, TESPA1, GATA2, GATA1, BLVRB, ALAS2, HBA1, HDC, ITGA2B, MPO, AZU1, FCN1, VCAN, CTSB, STAB1, IL2RB, DHFR, PAX5, IGHM, IRF8, CLEC4C, IL3RA, ALB, CDH5, MKI67"
geneOrder2=unlist(strsplit(s,", "))
geneOrder = unique(geneOrder1,geneOrder2)
genes.sub <- unique(geneOrder[geneOrder %in% rownames(df)])

rng = seq(1,length(genes.sub));
for (k in rng) {
  f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/091422_clust/Violin.",geneOrder[k],".pdf")
  print(paste0(k,": ",f.out))
  pdf(f.out,width=17,height=9)
  print(VlnPlot(df,features=genes.sub[k],sort=TRUE) + NoLegend() + theme(axis.text.x=element_text(angle=90)))
  dev.off()
}

print("Creating UMAP with seurat clusters...")
f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/091422_clust/RNA_umap.png")
png(filename=f.out,width = 5000,height=4000,res=500)
p1 <- DimPlot(df, group.by = "subclust_v4") & theme(plot.title = element_text(hjust = 0.5)) 
print(LabelClusters(plot = p1, id = "subclust_v4"))
dev.off()

cellType="HSCs"
for (cellType in unique(df$subclust_v4)) {
  print(cellType)
  cellNames1 <- rownames(df@meta.data[df$subclust_v4 == cellType,])
  cellType = gsub("/","_",cellType)
  f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/091422_clust/RNA_umap.",cellType,".png")
  png(filename=f.out,width = 5000,height=4000,res=500)
  print(DimPlot(df, label=T,  cells.highlight= list(cellNames1), cols.highlight = c("red"), cols= "grey") & theme(plot.title = element_text(hjust = 0.5)))
  # print(LabelClusters(plot = p1, id = "Final_V2"))
  dev.off()
}

s = "Neutrophils, Monocyte progenitors, Inflammatory macrophages, Kupffer cells, NK cells, B cells, pDCs, cDCs, Hepatocytes, LSECs, No markers"
cellOrder = c("HSCs",
      "MEMPs",
      "Early erythroid",
      "Late erythroid",
      paste0("Mast cells.",c(1,2,4)),
      paste0("Mast cells.3.",c(1,2)),
      "Megakaryocytes", "Granulocyte progenitors",
      unlist(strsplit(s,", ")))
# cellOrder=rev(cellOrder)
y=as.character(unique(df@meta.data$subclust_v4))
cellLevels=c(cellOrder[cellOrder %in% y],y[!(y %in% cellOrder)])
cellLevels <- rev(cellLevels)

df@meta.data$subclust_v4b <- factor(as.character(df@meta.data$subclust_v4),levels=cellLevels)
Idents(df) <- "subclust_v4b"

print("Creating dotplot with clusters...")
f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/091422_clust/RNA_dotplot.pdf")
pdf(f.out,width = 20,height=10)
print(
  DotPlot(
    df,
    features=genes.sub,
    col.min	= 0,
    col.max = 1,
    cluster.idents = FALSE
  ) + 
    theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1))
)
dev.off()



# print("FindAllMarkers...")
# res <- FindAllMarkers(df,logfc.threshold=0.25,min.pct = 0.1)
# f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/091422_clust/FindAllMarkers.RNA.txt")
# print(paste0("Saving to file: ",f.out))
# print("...")
# fwrite(res,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)



