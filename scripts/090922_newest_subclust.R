library(Seurat)
library(data.table)
library(ggplot2)
dir="/oak/stanford/groups/smontgom/amarder/neuro-variants"

DATASET="DS_Multiome_ds"
df = readRDS(paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET,"_v2","/RNA_FindClusters.rds"))
f = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET,"_v2","/meta.090822_subclust.txt")
meta = fread(f,data.table = F,stringsAsFactors = F)
rownames(meta) <- meta[,1]; meta <- meta[,-1]

if (DATASET=="DS_Multiome_h") { 
  meta$subclust_v2 <- meta$subclust
  meta$subclust_v2[meta$subclust %in% c(paste0("Maybe HSCs/Prog:2.",c(1,3)))] <- "Granulocyte progenitors:2.1,3"
  meta$subclust_v2[meta$subclust %in% c(paste0("pro B cells:4.",1:4),"?:18.2")] <- "pro B cells:4.1,2,3,4+18.2"
  meta$subclust_v2[meta$subclust %in% c(paste0("?:18.",c(1,3,4)))] <- "pre pro B cells:18.1,3,4"
  meta$subclust_v2[meta$subclust %in% c(paste0("Maybe HSCs/Prog:2.",c(2,4)))] <- "Progenitors:2.2,4"
  meta$subclust_v2_short <- unlist(lapply(strsplit(meta$subclust_v2,":"),function(x) x[[1]]))
} 
if (DATASET=="DS_Multiome_ds") {
  meta$subclust_v2 <- meta$subclust
  meta$subclust_v2[meta$subclust %in% c(paste0("Megakaryocytes:11,14.",c(1,2,3,4)),paste0("MEMPs:3,5,7.",c(1)))] <- "MEMPs"
  meta$subclust_v2[meta$subclust %in% c(paste0("MEMPs:3,5,7.",c(2,3,4)),paste0("Early erythroid:0,2,10,16,17,19.",c(2,4)))] <- "Megakaryocytes"
  meta$subclust_v2[meta$subclust %in% c(paste0("Early erythroid:0,2,10,16,17,19.",c(1,3)))] <- "Early erythroid"
  meta$subclust_v2[meta$subclust %in% c(paste0("Neutrophils:21.",2))] <- "Neutrophils"
  meta$subclust_v2[meta$subclust %in% c(paste0("Neutrophils:21.",c(1,3,4)))] <- "Granulocyte progenitors"
  meta$subclust_v2[meta$subclust %in% c(paste0("Proinflammatory macrophages:15.",c(3,4)))] <- "Proinflammatory macrophages"
  meta$subclust_v2[meta$subclust %in% c(paste0("Proinflammatory macrophages:15.",c(1)),paste0("Kupffer cells/cDCs:6,29.",c(1,2,3)))] <- "cDCs"
  meta$subclust_v2[meta$subclust %in% c(paste0("Kupffer cells/cDCs:6,29.",c(4)))] <- "Kupffer cells"
  meta$subclust_v2[meta$subclust %in% c(paste0("B lineage:22.",c(3)))] <- "Pre pro B cells"
  meta$subclust_v2[meta$subclust %in% c(paste0("B lineage:22.",c(1,2,4)))] <- "pro B cells"
  meta$subclust_v2[meta$subclust %in% c(paste0("?:28"))] <- "Hepatocytes"
  meta$subclust_v2[meta$subclust %in% c(paste0("?:27"))] <- "Endothelial cells"
  meta$subclust_v2[meta$subclust %in% c(paste0("Proinflammatory macrophages:15.2"))] <- "Proinflammatory macrophages 2"
  meta$subclust_v2[meta$subclust %in% c(paste0("?:12.",c(1,2,3,4)))] <- "?(12)"
  meta$subclust_v2[meta$subclust %in% c(paste0("?:26.",c(1,2,3,4)))] <- "?(26)"
  meta$subclust_v2[meta$subclust %in% c(paste0("?:8.",c(1,2,3,4)))] <- "?(8)"
  meta$subclust_v2_short <- unlist(lapply(strsplit(meta$subclust_v2,":"),function(x) x[[1]]))
  # Mast cells – split in 4
  i = meta$subclust_v2_short %in% c("Mast cells")
  dftmp = df[,i]
  kmeans_results <- kmeans(Embeddings(dftmp,reduction="harmony"),4)
  kmeans_results1 <- data.frame(kmeans_results$cluster)
  meta$subclust_v2_short[i] <- paste0("Mast cells.",kmeans_results1[,1])
}
df@meta.data <- meta

df@meta.data$subclust_v2[df@meta.data$subclust_v2=="MEMPs"] <- "Megakaryocytes"
df@meta.data$subclust_v2[df@meta.data$subclust_v2=="Megakaryocytes"] <- "MEMPs"
df@meta.data$subclust_v2_short[df@meta.data$subclust_v2_short=="MEMPs"] <- "Megakaryocytes"
df@meta.data$subclust_v2_short[df@meta.data$subclust_v2_short=="Megakaryocytes"] <- "MEMPs"

f.out = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET,"_v2","/meta.090922_subclust.txt")
fwrite(df@meta.data,
       f.out,
       quote=F,na = "NA",sep = "\t",row.names = T,col.names = T)

###################################################################

Idents(df) <- "subclust_v2_short"
system(paste0("rm -r ",dir,"/output/data/",DATASET,"_v2","/090922_subclust"))
system(paste0("mkdir -p ",dir,"/output/data/",DATASET,"_v2","/090922_subclust"))
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

system(paste0("rm -r ",dir,"/output/data/",DATASET,"_v2","/090922_clust"))
system(paste0("mkdir -p ",dir,"/output/data/",DATASET,"_v2","/090922_clust"))
rng = seq(1,length(genes.sub));
for (k in rng) {
  f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/090922_clust/Violin.",geneOrder[k],".pdf")
  print(paste0(k,": ",f.out))
  pdf(f.out,width=17,height=9)
  print(VlnPlot(df,features=genes.sub[k],sort=TRUE) + NoLegend() + theme(axis.text.x=element_text(angle=90)))
  dev.off()
}

print("Creating UMAP with seurat clusters...")
f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/090922_clust/RNA_umap.png")
png(filename=f.out,width = 5000,height=4000,res=500)
p1 <- DimPlot(df, group.by = "subclust_v2_short") & theme(plot.title = element_text(hjust = 0.5)) 
print(LabelClusters(plot = p1, id = "subclust_v2_short"))
dev.off()

cellType="HSCs"
for (cellType in unique(df$subclust_v2_short)) {
  
  print(cellType)
  cellNames1 <- rownames(df@meta.data[df$subclust_v2_short == cellType,])
  f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/090922_clust/RNA_umap.",cellType,".png")
  png(filename=f.out,width = 5000,height=4000,res=500)
  print(DimPlot(df, label=T,  cells.highlight= list(cellNames1), cols.highlight = c("red"), cols= "grey") & theme(plot.title = element_text(hjust = 0.5)))
  # print(LabelClusters(plot = p1, id = "Final_V2"))
  dev.off()
  
}

print("FindAllMarkers...")
res <- FindAllMarkers(df,logfc.threshold=0.25,min.pct = 0.1)
f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/090922_clust/FindAllMarkers.RNA.txt")
print(paste0("Saving to file: ",f.out))
print("...")
fwrite(res,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)



