library(data.table)
library(Seurat)
library(Matrix.utils) # for pseudobulking
library(BiocParallel) # for parallel
library(variancePartition) # for limma voom

cluster_label="combi_annot"
disease_status="Healthy"
sampletype="Femur"
subset_column="sample"
cell_type_filename="Cycling HSC"
dfseurat=list()
i=0; for (sampletype in c("Femur","Liver")) {
  for (disease_status in c("Healthy","DownSyndrome")) {
    i = i + 1; print(i)
    f = paste0("~/Documents/Research/t21-proj/out/full/downsample_induced/input/",disease_status,".",sampletype,".",cell_type_filename,".",cluster_label,".rds")
    dfseurat[[i]] = readRDS(f)
    dfseurat[[i]]@meta.data$dataset = paste0("Integ: ",disease_status," ",sampletype)
  }
}

i = i + 1
disease_status="DownSyndrome"
sampletype="Liver"
cell_type = "Cycling HSCs/MPPs"
cluster_label = "leiden_names"
cell_type_filename = gsub("/","_",cell_type)
f = paste0("~/Documents/Research/t21-proj/out/full/downsample_induced/input/",disease_status,".",sampletype,".",cell_type_filename,".",cluster_label,".rds")
dfseurat[[i]] = readRDS(f)
dfseurat[[i]]@meta.data$dataset = paste0("Separate: ",disease_status," ",sampletype)

i = i + 1
disease_status="Healthy"
sampletype="Liver"
cell_type = "Cycling HSCs/MPPs"
cluster_label = "Predicted"
cell_type_filename = gsub("/","_",cell_type)
f = paste0("~/Documents/Research/t21-proj/out/full/downsample_induced/input/",disease_status,".",sampletype,".",cell_type_filename,".",cluster_label,".rds")
dfseurat[[i]] = readRDS(f)
dfseurat[[i]]@meta.data$dataset = paste0("T21 Transfer: ",disease_status," ",sampletype)

print("Merging...")
dfcombined <- merge(dfseurat[[1]],
                    dfseurat[2:length(dfseurat)],
                    add.cell.ids=NULL)
Idents(dfcombined) <- "dataset"
dfcombined <- NormalizeData(dfcombined,assay="RNA")
dfcombined <- ScaleData(dfcombined,assay="RNA")

# k=1
# genes_lst = c("CD34", "SPINK2", "MLLT3","MPO", "LYZ", "SPI1", "IL7R", "IGLL1", "IGHM", "CD79A")
# # genes_lst = c("IGLL1", "IGHM", "CD79A")
# dir.create("~/Documents/Research/t21-proj/out/full/cycling_plots_v2/")
# for (k in 1:length(genes_lst)) {
#   f.out <- paste0("~/Documents/Research/t21-proj/out/full/cycling_plots_v2/Violin.",genes_lst[k],".pdf")
#   print(paste0(k,": ",f.out))
#   pdf(f.out,width=17,height=7)
#   print(VlnPlot(dfcombined,features=genes_lst[k],sort=FALSE) + NoLegend())
#   dev.off()
# }

k=1
genes_lst = c("CD34", "SPINK2", "MLLT3","MPO", "LYZ", "SPI1", "IL7R", "IGLL1", "IGHM", "CD79A")
dir.create("~/Documents/Research/t21-proj/out/full/cycling_plots_v2/")
# genes_lst = c("IGLL1", "IGHM", "CD79A")
for (k in 1:length(genes_lst)) {
  f.out <- paste0("~/Documents/Research/t21-proj/out/full/cycling_plots_v2/Violin.",genes_lst[k],".pdf")
  print(paste0(k,": ",f.out))
  pdf(f.out,width=17,height=7)
  print(VlnPlot(dfcombined,features=genes_lst[k],sort=FALSE,slot="data") + NoLegend() + geom_boxplot(width=0.1, fill="white"))
  dev.off()
  
  f.out <- paste0("~/Documents/Research/t21-proj/out/full/cycling_plots_v2/Ridge.",genes_lst[k],".ridge.pdf")
  print(paste0(k,": ",f.out))
  pdf(f.out,width=17,height=7)
  print(RidgePlot(dfcombined, features = genes_lst[k], ncol = 1,slot="data"))
  dev.off()
}


# Single cell heatmap of feature expression
f.out <- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/cycling_plots/doheatmap.pdf")
print(paste0(k,": ",f.out))
pdf(f.out,width=17,height=7)
print(DoHeatmap(subset(dfcombined, downsample = 100), features = genes_lst, size = 3,slot = "scale.data"))
# DoHeatmap(subset(dfcombined, downsample = 50), features = genes_lst, size = 3,slot = "data")
dev.off()

f.out <- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/cycling_plots/dotplot.pdf")
print(paste0(k,": ",f.out))
pdf(f.out,width=17,height=7)
print(DotPlot(dfcombined, features = genes_lst) + RotatedAxis())
# DoHeatmap(subset(dfcombined, downsample = 50), features = genes_lst, size = 3,slot = "data")
dev.off()



