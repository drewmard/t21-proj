library(data.table)
# library(Seurat)
library(Matrix.utils) # for pseudobulking
library(BiocParallel) # for parallel
library(variancePartition) # for limma voom

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


