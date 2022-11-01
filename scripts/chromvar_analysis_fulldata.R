# df = subset(dfcombined1,subclust_v6=="HSCs")
print("FindTopFeatures...")
dfcombined1 <- FindTopFeatures(dfcombined1,min.cutoff = 200,assay="ATAC")
# df <- FindTopFeatures(df,min.cutoff = "q75",assay="ATAC")
print("RunTFIDF...")
dfcombined1 <- RunTFIDF(dfcombined1)
print("RunSVD...")
dfcombined1 <- RunSVD(dfcombined1)
print("RunHarmony...")
lambda_val=1; 
dfcombined1 <- RunHarmony(dfcombined1,"dataset",reduction="lsi",assay.use="ATAC",project.dim=FALSE,lambda=lambda_val)
print("RunUMAP...")
dfcombined1 <- RunUMAP(dfcombined1, reduction = "harmony", dims = 2:30,reduction.name = "atac_umap")

for (gene_to_use in c("REL","RELB","NFKB2","GATA1")) { 
  
  motif_to_use = subset(tmp,gene_name==gene_to_use)$motif_name
  print(gene_to_use)
  # print(AverageExpression(df,features = c(motif_to_use),assay="chromvar",group.by="ATAC_snn_res.1"))
  # # motif_to_use = subset(tmp,gene_name=="GATA1")$motif_name
  # 
  # DefaultAssay(df) <- 'chromvar'
  # print("Creating UMAP with seurat clusters (round 2)...")
  # f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/ChromVAR_HSCs_vln.",gene_to_use,".png")
  # png(filename=f.out,width = 5000,height=4000,res=500)
  # p1 = VlnPlot(df, features = c(motif_to_use),assay="chromvar",group.by="ATAC_snn_res.1")
  # print(p1)
  # dev.off()
  print("Creating UMAP with seurat clusters (round 2)...")
  f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/ChromVAR_all_umap.",gene_to_use,".png")
  png(filename=f.out,width = 5000,height=4000,res=500)
  p1 <- FeaturePlot(
    object = dfcombined1, # df,
    features = c(motif_to_use),
    reduction="atac_umap",
    pt.size = 0.1,
    min.cutoff = 0,
    max.cutoff = 'q95'
  ); print(p1)
  dev.off()
}

library(ggplot2)
print("Creating UMAP with seurat clusters (round 2)...")
f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/ATAC_all_umap.seurat_clusters.png")
png(filename=f.out,width = 5000,height=4000,res=500)
p1 <- DimPlot(dfcombined1, group.by = "subclust_v6",reduction="atac_umap") & theme(plot.title = element_text(hjust = 0.5)) 
print(LabelClusters(plot = p1, id = "subclust_v6"))
dev.off()

dfcombined1 <- FindVariableFeatures(dfcombined1, selection.method = "vst", assay = "chromvar")
top10 <- head(VariableFeatures(dfcombined1), 10)
hv = data.frame(i=1:length(VariableFeatures(dfcombined1)),motif_name=VariableFeatures(dfcombined1))
hv = merge(hv,tmp,by="motif_name")
hv = hv[order(hv$i),]
library(data.table)
f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/chromvar.highlyvariable.h.txt")
fwrite(hv,f.out,quote = F,na = "NA",sep = "\t",row.names = F,col.names = T)
# plot variable features with and without labels
f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/chromvar.highly_variable.png")
png(filename=f.out,width = 5000,height=4000,res=500)
plot1 <- VariableFeaturePlot(dfcombined1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
dev.off()



