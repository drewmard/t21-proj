library(data.table)
library(Seurat)
library(Matrix.utils) # for pseudobulking

##################

integrated_labels = fread("/oak/stanford/groups/smontgom/amarder/t21-proj/out/Ts21_reference_scArches_Healthy_Liver_metadata.csv",data.table = F,stringsAsFactors = F)
integrated_labels$cellname=unlist(lapply(strsplit(integrated_labels$barcode,"-"),function(x) paste(x[1:2],collapse = "-")))
integrated_labels = integrated_labels[!duplicated(integrated_labels[,c("cellname","sample","Predicted")]),] # why??

cluster_label="Predicted"
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
df.mg = merge(seurat_meta,integrated_labels[,c("cellname","sample","Predicted")],by=c("cellname","sample"),all.x=TRUE)
df.mg = df.mg[order(df.mg$i),]
print(nrow(seurat_meta))
print(nrow(df.mg))
rownames(df.mg)=rows_to_save
dfseurat@meta.data=df.mg
rm(df.mg)
rm(seurat_meta)


# Subset matrix down to cell type of interest
cluster_label="Predicted"
cell_type="HSCs/MPPs"
cell_type_filename = gsub("/","_",cell_type)
f.out = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/downsample_induced/input/",disease_status,".",sampletype,".",cell_type_filename,".",cluster_label,".rds")
dfsub = dfseurat[,dfseurat@meta.data[,cluster_label]==cell_type]
saveRDS(dfsub,f.out)
cell_type="Cycling HSCs/MPPs"
cell_type_filename = gsub("/","_",cell_type)
f.out = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/downsample_induced/input/",disease_status,".",sampletype,".",cell_type_filename,".",cluster_label,".rds")
dfsub = dfseurat[,dfseurat@meta.data[,cluster_label]==cell_type]
saveRDS(dfsub,f.out)

