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
disease_status="DownSyndrome"
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
# }

# Subset matrix down to cell type of interest
cluster_label="combi_annot"
cell_type="HSCs"
cell_type_filename = gsub("/","_",cell_type)
f.out = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/downsample_induced/input/",disease_status,".",sampletype,".",cell_type_filename,".",cluster_label,".rds")
dfsub = dfseurat[,dfseurat@meta.data[,cluster_label]==cell_type]
saveRDS(dfsub,f.out)
cell_type="Cycling HSC"
cell_type_filename = gsub("/","_",cell_type)
f.out = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/downsample_induced/input/",disease_status,".",sampletype,".",cell_type_filename,".",cluster_label,".rds")
dfsub = dfseurat[,dfseurat@meta.data[,cluster_label]==cell_type]
saveRDS(dfsub,f.out)

cluster_label="leiden_names"
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

################################






disease_status="Healthy"
sampletype="Liver"
subset_column="sample"
print(paste0("Reading data: ",disease_status," ",subset_column,"..."))
# cell_type = "HSCs/MPPs"
cell_type = "Cycling HSC"
cell_type_filename = gsub("/","_",cell_type)
print("Reading metadata2...")
f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/cellComp/10X_",disease_status,"_",sampletype,".cellComp.csv")
meta2.full<-fread(f,data.table = F,stringsAsFactors = F)
colnames(meta2.full)[6] <- "leiden_names"
meta2.full$cellname = unlist(lapply(strsplit(meta2.full$V1,"-"),function(x) paste(x[1:2],collapse = "-")))
integrated_labels = fread("/oak/stanford/groups/smontgom/amarder/t21-proj/out/a3_tab.obs.withannot.csv",data.table = F,stringsAsFactors = F)
integrated_labels$cellname=unlist(lapply(strsplit(integrated_labels$barcodekey,"-"),function(x) paste(x[1:2],collapse = "-")))
integrated_labels$sample=unlist(lapply(strsplit(as.character(integrated_labels$patient_sample)," "),function(x) paste(x[2],collapse = "-")))
integrated_labels = integrated_labels[!duplicated(integrated_labels[,c("cellname","sample","combi_annot")]),] # why??
nrow(meta2.full)
meta2.full = merge(meta2.full,integrated_labels,by=c("sample","cellname"),all.x=TRUE)
nrow(meta2.full)

celltype1="Cycling HSC"
meta1 = subset(meta2.full,combi_annot==celltype1)[,c("patient","sample","sorting","leiden_latest","combi_annot")]
tab = aggregate(meta1$sample,by=meta1[,c("sample","patient","sorting")],length)
colnames(tab)[ncol(tab)] = "num_cells"
tab = tab[,-1] # remove sample
tab = subset(tab,num_cells >= 10)
dir.create("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/downsample_induced")
dir.create("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/downsample_induced/input")
f.out = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/downsample_induced/input/meta1.txt")
fwrite(tab,f.out,quote = F,na = "NA",sep = "\t",row.names = F,col.names = T)

celltype1="HSCs"
meta1 = subset(meta2.full,combi_annot==celltype1)[,c("patient","sample","sorting","leiden_latest","combi_annot")]
tab = aggregate(meta1$sample,by=meta1[,c("sample","patient","sorting")],length)
colnames(tab)[ncol(tab)] = "num_cells"
tab = tab[,-1] # remove sample
tab = subset(tab,num_cells >= 10)
dir.create("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/downsample_induced")
dir.create("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/downsample_induced/input")
f.out = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/downsample_induced/input/meta2.txt")
fwrite(tab,f.out,quote = F,na = "NA",sep = "\t",row.names = F,col.names = T)



library(data.table)
# library(Seurat)
library(Matrix.utils) # for pseudobulking
library(BiocParallel) # for parallel
library(variancePartition) # for limma voom

disease_status="DownSyndrome"
sampletype="Liver"
subset_column="sample"
print(paste0("Reading data: ",disease_status," ",subset_column,"..."))
# cell_type = "HSCs/MPPs"
cell_type = "Cycling HSC"
cell_type_filename = gsub("/","_",cell_type)
print("Reading metadata2...")
f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/cellComp/10X_",disease_status,"_",sampletype,".cellComp.csv")
meta2.full<-fread(f,data.table = F,stringsAsFactors = F)
colnames(meta2.full)[6] <- "leiden_names"
meta2.full$cellname = unlist(lapply(strsplit(meta2.full$V1,"-"),function(x) paste(x[1:2],collapse = "-")))
integrated_labels = fread("/oak/stanford/groups/smontgom/amarder/t21-proj/out/a3_tab.obs.withannot.csv",data.table = F,stringsAsFactors = F)
integrated_labels$cellname=unlist(lapply(strsplit(integrated_labels$barcodekey,"-"),function(x) paste(x[1:2],collapse = "-")))
integrated_labels$sample=unlist(lapply(strsplit(as.character(integrated_labels$patient_sample)," "),function(x) paste(x[2],collapse = "-")))
integrated_labels = integrated_labels[!duplicated(integrated_labels[,c("cellname","sample","combi_annot")]),] # why??
nrow(meta2.full)
meta2.full = merge(meta2.full,integrated_labels,by=c("sample","cellname"),all.x=TRUE)
nrow(meta2.full)

celltype1="Cycling HSC"
meta1.input = subset(meta2.full,combi_annot==celltype1)

celltype1="HSCs"
meta2.input = subset(meta2.full,combi_annot==celltype1)

[,c("patient","sample","sorting","leiden_latest","combi_annot")]


