library(data.table)
library(Seurat)
library(Matrix.utils) # for pseudobulking
library(BiocParallel) # for parallel
library(variancePartition) # for limma voom

##################

integrated_labels = fread("/oak/stanford/groups/smontgom/amarder/t21-proj/out/a3_tab.obs.withannot.csv",data.table = F,stringsAsFactors = F)
integrated_labels$cellname=unlist(lapply(strsplit(integrated_labels$barcodekey,"-"),function(x) paste(x[1:2],collapse = "-")))
integrated_labels$sample=unlist(lapply(strsplit(as.character(integrated_labels$patient_sample)," "),function(x) paste(x[2],collapse = "-")))
integrated_labels = integrated_labels[!duplicated(integrated_labels[,c("cellname","sample","combi_annot")]),] # why??

cluster_label="combi_annot"
disease_status="Healthy"
sampletype="Liver"
subset_column="sample"
f<- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/10X_",disease_status,'_',sampletype,'.umap2d.cells_removed.rds')
dfseurat=readRDS(f)
dfseurat@assays$RNA@key = "rna_"
y=grep("leiden_v",colnames(dfseurat@meta.data))
y=colnames(dfseurat@meta.data)[y[length(y)]]
dfseurat@meta.data$leiden_names = dfseurat@meta.data[,y]

# if (cluster_label!="leiden_names") { 
seurat_meta = dfseurat@meta.data 
rows_to_save = rownames(seurat_meta)
seurat_meta$i = 1:nrow(seurat_meta)
seurat_meta$cellname=unlist(lapply(strsplit(rownames(seurat_meta),"-"),function(x) paste(x[1:2],collapse = "-")))
df.mg = merge(seurat_meta,integrated_labels[,c("cellname","sample","combi_annot")],by=c("cellname","sample"),all.x=TRUE)
df.mg = df.mg[order(df.mg$i),]
print(nrow(seurat_meta))
print(nrow(df.mg))
rownames(df.mg)=rows_to_save
dfseurat@meta.data=df.mg
rm(df.mg)
rm(seurat_meta)

cell_type="Cycling HSC"
cell_type_filename = gsub("/","_",cell_type)
f.out = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/downsample_induced/input/",disease_status,".",sampletype,".",cell_type_filename,".",cluster_label,".rds")
dfsub = dfseurat[,dfseurat@meta.data[,cluster_label]==cell_type]
saveRDS(dfsub,f.out)


##############

sampletype="Femur"
for (disease_status in c("Healthy","DownSyndrome")) { 
  cluster_label="combi_annot"
  subset_column="sample"
  f<- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/10X_",disease_status,'_',sampletype,'.umap2d.cells_removed.rds')
  dfseurat=readRDS(f)
  dfseurat@assays$RNA@key = "rna_"
  y=grep("leiden_v",colnames(dfseurat@meta.data))
  y=colnames(dfseurat@meta.data)[y[length(y)]]
  dfseurat@meta.data$leiden_names = dfseurat@meta.data[,y]
  
  # if (cluster_label!="leiden_names") { 
  seurat_meta = dfseurat@meta.data 
  rows_to_save = rownames(seurat_meta)
  seurat_meta$i = 1:nrow(seurat_meta)
  seurat_meta$cellname=unlist(lapply(strsplit(rownames(seurat_meta),"-"),function(x) paste(x[1:2],collapse = "-")))
  df.mg = merge(seurat_meta,integrated_labels[,c("cellname","sample","combi_annot")],by=c("cellname","sample"),all.x=TRUE)
  df.mg = df.mg[order(df.mg$i),]
  print(nrow(seurat_meta))
  print(nrow(df.mg))
  rownames(df.mg)=rows_to_save
  dfseurat@meta.data=df.mg
  rm(df.mg)
  rm(seurat_meta)
  
  cell_type="Cycling HSC"
  cell_type_filename = gsub("/","_",cell_type)
  f.out = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/downsample_induced/input/",disease_status,".",sampletype,".",cell_type_filename,".",cluster_label,".rds")
  dfsub = dfseurat[,dfseurat@meta.data[,cluster_label]==cell_type]
  saveRDS(dfsub,f.out)
}
