library(Seurat)
library(data.table)
library(ggplot2)
dir="/oak/stanford/groups/smontgom/amarder/neuro-variants"
# system(paste0("rm -r ",dir,"/output/data/",DATASET,"_v2","/091522_clust"))

# DATASET="DS_Multiome_ds"
DATASET="DS_Multiome_h"
df = readRDS(paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET,"_v2","/RNA_FindClusters.rds"))
f = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET,"_v2","/meta.091422_subclust.txt")
meta = fread(f,data.table = F,stringsAsFactors = F)
rownames(meta) <- meta[,1]; meta <- meta[,-1]


if (DATASET=="DS_Multiome_h") { 
  meta$subclust_v5 <- meta$subclust_v4
  meta$subclust_v5[meta$subclust_v4 %in% c("Proinflammatory macrophages")] <- "Pro-inflammatory macrophages"
  meta$subclust_v5[meta$subclust_v4 %in% c("Mast cell")] <- "Mast cells"
} 
if (DATASET=="DS_Multiome_ds") {
  meta$subclust_v5 <- meta$subclust_v4
  meta$subclust_v5[meta$subclust_v4 %in% c("HSCs")] <- "HSCs1"
  meta$subclust_v5[meta$subclust_v4 %in% c("?(8).1","?(8).8","?(8).6","?(26).2")] <- "HSCs2"
  meta$subclust_v5[meta$subclust_v4 %in% paste0("Mast cells.",c("1","2","3.1","3.2","4"))] <- "Mast cells"
  meta$subclust_v5[meta$subclust_v4 %in% c("Megakaryocytes","?(26).1")] <- "Megakaryocytes"
  meta$subclust_v5[meta$subclust_v4 %in% c("Neutrophils","?(8).2")] <- "Neutrophils"
  meta$subclust_v5[meta$subclust_v4 %in% c("Proinflammatory macrophages","Proinflammatory macrophages 2")] <- "Pro-inflammatory macrophages"
  meta$subclust_v5[meta$subclust_v4 %in% c("Pre pro B cells","pro B cells","B cells")] <- "B cells"
  meta$subclust_v5[meta$subclust_v4 %in% c("?(8).3","?(8).4","?(8).5","?(8).7")] <- "Unknown"
  meta$subclust_v5[meta$subclust_v4 %in% c("Mix-up/endothelial","Hepatocytes")] <- "Stroma"
  meta$subclust_v5[meta$subclust_v4 %in% c("?(12)")] <- "No marker"
}
df@meta.data <- meta

f.out = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET,"_v2","/meta.091522_subclust.txt")
fwrite(df@meta.data,
       f.out,
       quote=F,na = "NA",sep = "\t",row.names = T,col.names = T)

###################################################################

system(paste0("mkdir -p ",dir,"/output/data/",DATASET,"_v2","/091522_clust"))

###################################################################

Idents(df) <- "subclust_v5"

s="CD34, SPINK2, PROM1, MLLT3, TESPA1, GATA2, GATA1, ALAS2, HBA1, HDC, ITGA2B, MPO, AZU1, SPI1, FCN1,VCAN,CTSB,NKG7, PRF1, IL2RB, IL7R, PAX5, IGHM, IRF8, CLEC4C, CLEC10A, MKI67"
geneOrder=unlist(strsplit(s,", |,"))
genes.sub <- unique(geneOrder[geneOrder %in% rownames(df)])
# geneOrder[!(geneOrder %in% genes.sub)]

# rng = seq(1,length(genes.sub));
# for (k in rng) {
#   f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/091522_clust/Violin.",geneOrder[k],".pdf")
#   print(paste0(k,": ",f.out))
#   pdf(f.out,width=17,height=9)
#   print(VlnPlot(df,features=genes.sub[k],sort=TRUE) + NoLegend() + theme(axis.text.x=element_text(angle=90)))
#   dev.off()
# }

print("Creating UMAP with seurat clusters...")
f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/091522_clust/RNA_umap.png")
png(filename=f.out,width = 5000,height=4000,res=500)
p1 <- DimPlot(df, group.by = "subclust_v5") & theme(plot.title = element_text(hjust = 0.5)) 
print(LabelClusters(plot = p1, id = "subclust_v5"))
dev.off()

# cellType="HSCs"
for (cellType in unique(df$subclust_v5)) {
  print(cellType)
  cellNames1 <- rownames(df@meta.data[df$subclust_v5 == cellType,])
  cellType = gsub("/","_",cellType)
  f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/091522_clust/RNA_umap.",cellType,".png")
  png(filename=f.out,width = 5000,height=4000,res=500)
  print(DimPlot(df, label=T,  cells.highlight= list(cellNames1), cols.highlight = c("red"), cols= "grey") & theme(plot.title = element_text(hjust = 0.5)))
  # print(LabelClusters(plot = p1, id = "Final_V2"))
  dev.off()
}

if (DATASET=="DS_Multiome_ds") {
  cellOrder = c("HSCs1",
                "HSCs2",
                "MEMPs",
                "Early erythroid",
                "Late erythroid",
                "Cycling erythroid",
                "Mast cells",
                "Megakaryocytes",
                "Granulocyte progenitors",
                "Neutrophils",
                "Pro-inflammatory macrophages",
                "Kupffer cells",
                "NK cells",
                "T cells",
                "B cells",
                "pDCs",
                "cDCs",
                "Stroma",
                "Unknown",
                "No marker")
  # cellOrder=rev(cellOrder)
}
if (DATASET=="DS_Multiome_h") {
  cellOrder = c("HSCs",
                "MEMPs",
                "Early erythroid",
                "Late erythroid",
                "Mast cells",
                "Megakaryocytes",
                "Granulocyte progenitors",
                "Neutrophils",
                "Pro-inflammatory macrophages",
                "Kupffer cells",
                "NK cells",
                "B cells",
                "pDCs",
                "Stroma",
                "No marker")
}
y=as.character(unique(df@meta.data$subclust_v5))
cellLevels=c(cellOrder[cellOrder %in% y],y[!(y %in% cellOrder)])
cellLevels <- rev(cellLevels)
df@meta.data$subclust_v5b <- factor(as.character(df@meta.data$subclust_v5),levels=cellLevels)
Idents(df) <- "subclust_v5b"

print("Creating dotplot with clusters...")
f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/091522_clust/RNA_dotplot.pdf")
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

