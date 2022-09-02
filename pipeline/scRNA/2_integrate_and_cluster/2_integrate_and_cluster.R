##########

library(Seurat)
library(harmony)
library(ggplot2)
library(data.table)

###########

# PARAM:

NUM_HVG=10000

#############

# dir=args[1]
# start=args[2]
# end=args[3]
# DATASET=args[4]

dir="/oak/stanford/groups/smontgom/amarder/t21-proj"
start=1
end=99
DATASET="DS_Multiome_combined"

###############

dir.out <- paste0(dir,"/out/combined/2_integrate_and_cluster/",DATASET)
dir.create(dir.out)

f.samp_lst=paste0(dir,"/pipeline/scRNA/2_integrate_and_cluster/input/",DATASET,"/samp_lst")
f.samp_ids=paste0(dir,"/pipeline/scRNA/2_integrate_and_cluster/input/",DATASET,"/samp_ids")
samp_lst <- fread(f.samp_lst,data.table = F,stringsAsFactors = F,header=F,sep='\t')[,1]
samp_ids <- fread(f.samp_ids,data.table = F,stringsAsFactors = F,header=F,sep='\t')[,1]
f.metadata <- paste0(dir,"/pipeline/scRNA/2_integrate_and_cluster/input/",DATASET,"/scRNA.metadata.txt")
metadata <- fread(f.metadata,data.table = F,stringsAsFactors = F,header=T,sep='\t')

#################

# Subset:
ind <- which(metadata$ID1 %in% c(15633,15593))
samp_lst <- samp_lst[ind]
samp_ids <- samp_ids[ind]
metadata <- metadata[ind,]

#################

for (i in 1:length(samp_lst)) {
  samp <- samp_lst[[i]]
  # print(paste0("Reading data from ",i,"/",length(samp_lst),": ",samp," ..."))
  f=paste0(samp,"/seurat_obj.rds")
  if (!dir.exists(samp)) {print(paste0(i,": ",samp))}
}

if (start==1) {
  
  print("read_10X_and_saveRNAcounts...")
  dfseurat <- list()
  i=1
  for (i in 1:length(samp_lst)) {
    samp <- samp_lst[[i]]
    print(paste0("Reading data from ",i,"/",length(samp_lst),": ",samp," ..."))
    f=paste0(samp,"/seurat_obj.rds")
    dfseurat[[i]] <- readRDS(file = f)
    
    dfseurat[[i]] <- dfseurat[[i]][,1:min(ncol(dfseurat[[i]]),2000)]

    dfseurat[[i]]$i <- i
    dfseurat[[i]]$cell <- rownames(dfseurat[[i]]@meta.data)
    dfseurat[[i]]@meta.data <- cbind(dfseurat[[i]]@meta.data,subset(metadata,ID2==samp_ids[i]))
    dfseurat[[i]]$dataset <- dfseurat[[i]]$ID2
    
  }
  
  print("Merging...")
  for (i in 1:length(dfseurat)) { dfseurat[[i]]$dataset <- samp_ids[i]}
  dfcombined <- merge(dfseurat[[1]],
                      dfseurat[2:length(dfseurat)])#, # cell_names <- unlist(lapply(1:length(dfseurat),function(i) paste0("data",i,'_',rownames(dfseurat[[i]]@meta.data))))
  
  print(paste0(ncol(dfcombined)," cells total..."))
  
  print("NormalizeData + FindVariableFeatures + ScaleData + RunPCA...")
  dfcombined <- NormalizeData(dfcombined,normalization.method = "LogNormalize", scale.factor = 10000) %>% 
    FindVariableFeatures(selection.method = "vst", nfeatures = NUM_HVG) %>% 
    ScaleData() %>% # change to all genes using , features=rownames(pbmc) if using heatmap
    RunPCA()
  
  print("RunHarmony...")
  lambda_val=1; dfcombined <- RunHarmony(dfcombined,"dataset",lambda=lambda_val)
  
  print("RunUMAP...")
  dfcombined <- RunUMAP(dfcombined, reduction = "harmony", dims = 1:30)
  
  print("FindNeighbors...")
  dfcombined <- FindNeighbors(object = dfcombined, reduction = 'harmony', dims = 1:30)
  
  f.out <- paste0(dir,"/out/combined/2_integrate_and_cluster/",DATASET,"/RNA_FindNeighbors.rds")
  print(paste0("Saving data: ",f.out,"..."))
  saveRDS(dfcombined,file = f.out)
  
  print("Creating UMAP with dataset-based coloring (round 2)...")
  dir.create(paste0(dir,"/out/combined/2_integrate_and_cluster/",DATASET))
  f.out <- paste0(dir,"/out/combined/2_integrate_and_cluster/",DATASET,"/RNA_umap.dataset_integration.png")
  png(filename=f.out,width = 5000,height=2600,res=500)
  p1 <- DimPlot(dfcombined, group.by = "dataset")
  print(p1 & theme(plot.title = element_text(hjust = 0.5)))
  dev.off()
  
  # check batch correction:
  print("Creating split-UMAP with dataset-based coloring (round 2)...")
  p1 <- DimPlot(dfcombined, split.by = "dataset",ncol = 2)
  f.out <- paste0(dir,"/out/combined/2_integrate_and_cluster/",DATASET,"/RNA_umap.dataset_integration.split.png")
  png(filename=f.out,width = 5000,height=5000,res=400)
  print(p1 & theme(plot.title = element_text(hjust = 0.5)) + NoLegend())
  dev.off()
  
}

if (start == 2) {
  f <- paste0(dir,"/out/combined/2_integrate_and_cluster/",DATASET,"/RNA_FindNeighbors.rds")
  print(paste0("Loading data: ",f,"..."))
  dfcombined <- readRDS(file = f)
}

if (end > 1 & start <= 2) {
  
  print("FindClusters...")
  # resolution_val=0.1; algnum=4; dfcombined <- FindClusters(object = dfcombined, verbose = FALSE, algorithm = algnum,resolution=resolution_val)
  resolution_val=1.5; algnum=1; dfcombined <- FindClusters(object = dfcombined, verbose = FALSE, algorithm = algnum,resolution=resolution_val)
  
  print("Creating UMAP with seurat clusters...")
  # p1 <- DimPlot(dfcombined, group.by = "seurat_clusters")
  f.out <- paste0(dir,"/out/combined/2_integrate_and_cluster/",DATASET,"/RNA_umap.seurat_clusters.png")
  png(filename=f.out,width = 5000,height=4000,res=500)
  p1 <- DimPlot(dfcombined, group.by = "seurat_clusters") & theme(plot.title = element_text(hjust = 0.5)) 
  print(LabelClusters(plot = p1, id = "seurat_clusters"))
  dev.off()
  
  # Save data:
  f.out <- paste0(dir,"/out/combined/2_integrate_and_cluster/",DATASET,"/RNA_FindClusters.rds")
  print(paste0("Saving data: ",f.out,"..."))
  saveRDS(dfcombined,file = f.out)
  
}


