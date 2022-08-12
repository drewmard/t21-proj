# module load R/4.1.2
library(Seurat)
library(Signac)
library(data.table)
library(ggplot2)

dir="/oak/stanford/groups/smontgom/amarder/neuro-variants"
DATASET="DS_Multiome_h"

dfrna <- readRDS("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h/RNA_FindClusters.rds")

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
# f.out <- "/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h/RNA_meta_v1b.txt"
# fwrite(meta2,f.out,quote = F,na = "NA",sep = '\t',row.names = T,col.names = T)
dfrna@meta.data <- meta2; rm(meta); rm(meta2)

Idents(dfrna) <- "Final"

clusterData <- fread("/oak/stanford/groups/smontgom/amarder/t21-proj/scripts/subcluster/rna_subclust_v1.csv",data.table = F,stringsAsFactors = F)
genes <- unlist(lapply(strsplit(unique(clusterData$Genes),","),function(x) gsub(" ","",x)))
genes.sub <- genes[genes %in% rownames(dfrna)]
print("Creating dotplot with clusters...")
# dev.off()
f.out <- paste0(dir,"/output/data/",DATASET,"/RNA_dotplot.final.png")
png(filename=f.out,width = 15000,height=4000,res=500)

f.out <- paste0(dir,"/output/data/",DATASET,"/RNA_dotplot.final.pdf")
pdf(f.out,width = 17,height=5)
print(
  DotPlot(
    dfrna[,!is.na(dfrna$Final)],
    features=genes.sub,
    col.min	= 0,
    col.max = 1,
    cluster.idents = TRUE
  ) + 
    theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1))
)
dev.off()

system(paste0("mkdir -p ",dir,"/output/data/",DATASET,"/subclust_v1_v2"))
rng = seq(1,length(genes.sub));
for (k in rng) {
  # for (k in seq(1,length(genes))) {
  print(k)
  f.out <- paste0(dir,"/output/data/",DATASET,"/subclust_v1_v2/Violin.",genes.sub[k],".pdf")
  # png(filename=f.out,width = 5000,height=4000,res=500)
  pdf(f.out,width=17,height=9)
  # print(RidgePlot(dfrna,features=genes[k]) + NoLegend())
  print(VlnPlot(dfrna,features=genes.sub[k],sort=TRUE) + NoLegend())
  dev.off()
}

