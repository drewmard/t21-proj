# srun --account=smontgom --partition=batch --time=24:00:00 --mem=256G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash 

# conda activate /oak/stanford/groups/smontgom/amarder/micromamba/envs/seurat5v2

library(data.table)

dir="/oak/stanford/groups/smontgom/amarder/Trisomy"
DATASET="Trisomy18_Trisomy21_Ctrl"
SUFFIX=".PRE_CB"

####################################################################################

genes_to_use = c("CD34","SPINK2","PROM1","MLLT3",
                 "GATA1","KLF1","TESPA1","AHSP",
                 "ALAS2","HBA1","GYPA","GATA2","HDC",
                 "CPA3","ITGA2B","GP9",
                 "MPO","AZU1","SPI1","LYZ","CD14","CD68",
                 "S100A9","MNDA","FCN1","CD163","MS4A7","VCAN","VCAM1","MARCO",
                 "C1QA","C1QB","C1QC",
                 "CTSB","NKG7","PRF1",
                 "GZMA",
                 "TIGIT","TRAC","CD3D","RORC", # T CELLS
                 "IL2RB",
                 "IL7R","DHFR","PAX5","MME","IGLL1","IGHM",
                 "CD79A","CD19","JCHAIN","IRF8","CLEC4C","IL3RA",
                 "CD1C","CLEC4A",
                 "CLEC10A", #cdc1
                 "CLEC9A" , "THBD","XCR1", "BATF3", 
                 "MKI67")

####################################################################################

meta=fread("/oak/stanford/groups/smontgom/amarder/Trisomy/output/data/Trisomy18_Trisomy21_Ctrl/annotations.meta_data.v3.txt",data.table = F,stringsAsFactors = F)
meta = meta[,-1]
rownames(meta) = meta[,1]
meta = meta[,-1]

meta$WNN_Ts18_Labels[meta$WNN_Ts18_Labels%in%c("Granulocyte progenitors","Granulocytes")]="Granulocytic lineage"
meta$WNN_Ts18_Labels[meta$WNN_Ts18_Labels%in%c("Granulocyte/monocyte progenitors","Monocytes")]="Monocytic lineage"

meta$WNN_Ts21_Labels[meta$WNN_Ts21_Labels%in%c("Granulocytes")]="Granulocytic lineage"
meta$WNN_Ts21_Labels[meta$WNN_Ts21_Labels%in%c("Monocytes")]="Monocytic lineage"

meta$WNN_Ctrl_Labels[meta$WNN_Ctrl_Labels%in%c("Granulocyte progenitors","Granulocytes")]="Granulocytic lineage"
meta$WNN_Ctrl_Labels[meta$WNN_Ctrl_Labels%in%c("Granulocyte/monocyte progenitors","Monocytes")]="Monocytic lineage"

meta$RNA_Ts18_Labels[meta$RNA_Ts18_Labels=="31"] = "unknown"

meta$final_integrated_labels[meta$final_integrated_labels%in%c("Granulocyte progenitors","Granulocytes")] = "Granulocytic lineage"
meta$final_integrated_labels[meta$final_integrated_labels%in%c("Granulocyte/monocyte progenitors","Monocytes")] = "Monocytic lineage"

####################################################################################

# Read Ts18, Ts21, Ctrl data into R, generate UMAP
for (dataset in c("Ts18","Ts21","Ctrl")) {
  
  if (dataset=="Ts18") {
    f="/oak/stanford/groups/smontgom/amarder/Trisomy/output/data/Trisomy18/ATAC_RNA_combined.wnn_cluster_v2.rds"
  } else if (dataset=="Ts21") {
    f="/oak/stanford/groups/smontgom/amarder/Trisomy/output/data/Trisomy21/ATAC_RNA_combined.GeneActivity.DimRed.rds"
  } else if (dataset=="Ctrl") {
    f="/oak/stanford/groups/smontgom/amarder/Trisomy/output/data/Ctrl/ATAC_RNA_combined.GeneActivity.DimRed.rds"
  }
  
  df = readRDS(f)
  
  tmp = df@meta.data
  y = rownames(tmp)
  tmp$i = 1:nrow(tmp)
  tmp = merge(tmp[,c("dataset","cell","i")],meta,by=c("dataset","cell"))
  tmp = tmp[order(tmp$i),]
  rownames(tmp) = y
  
  df@meta.data = tmp
  
  colName = paste0("WNN_",dataset,"_Labels")
  
  DefaultAssay(df) = "RNA"; Idents(df) = colName
  f.out <- paste0(dir,"/output/data/Trisomy18_Trisomy21_Ctrl/RNA_dotplot.",dataset,SUFFIX,".png")
  png(filename=f.out,width = 9500,height=5000,res=500)
  print(DotPlot(object = df, features = genes_to_use) + RotatedAxis())
  dev.off()
  
}

f="/oak/stanford/groups/smontgom/amarder/Trisomy/output/data/Trisomy18_Trisomy21_Ctrl/WNN.GeneActivity.clusters.rds"
# add in meta data
rownames(meta) = meta[,1]
meta = meta[,-1]
dfcombined@meta.data = meta


df = readRDS(f)

# Read combined data into R, generate UMAP


# Pull out MEMPs, mast cells, erythroid - see if I can disentangle which boundary got screwed up...