# module load R/4.1.2
library(Seurat)
library(Signac)
library(data.table)
library(ggplot2)
library(GenomicRanges)
library(future)
library(parallel)
library("JASPAR2020")
library("TFBSTools")
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)
library(dplyr)

dir="/oak/stanford/groups/smontgom/amarder/neuro-variants"
DATASET="DS_Multiome_h"
args = commandArgs(trailingOnly=TRUE)
dir=args[1]
start=args[2]
end=args[3]
DATASET=args[4]

########################################################################

source(paste0(dir,"/scripts/single_cell/helper_func/","integration_function_sheet.R"))
f.samp_lst=paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/scripts/single_cell/input/",DATASET,"/samp_lst")
f.samp_ids=paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/scripts/single_cell/input/",DATASET,"/samp_ids")
samp_lst <- fread(f.samp_lst,data.table = F,stringsAsFactors = F,header=F,sep='\t')[,1]
samp_ids <- fread(f.samp_ids,data.table = F,stringsAsFactors = F,header=F,sep='\t')[,1]

########################################################################

system(paste0("mkdir -p ",dir,"/output/data/",DATASET))

if (start <= 1) {
  print("read_10X_and_acquire_cell_names...")
  read_10X_and_acquire_cell_names(samp_lst)
}

if (start <= 1 & end > 0) {
  macs2_path <- '/home/amarder/anaconda3/envs/colocalization/bin/macs2'
  
  # 1: Call peaks within each sample, and combine into a single peak set
  print("identify_cellnames...")
  cellnames <- identify_cellnames(samp_lst)
  print("create_fragment_object...")
  frags <- create_fragment_object(samp_lst)
}

########################################################################

dfrna <- readRDS(paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET,"/RNA_FindClusters.rds"))
meta2 = fread(paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET,"/RNA_meta_v3.txt"),data.table = F,stringsAsFactors = F)
dfrna@meta.data <- meta2; #rm(meta); rm(meta2)
dfatac <- readRDS(paste0(dir,"/output/data/",DATASET,"/round2_FindClusters.rds"))

rnacluster <- dfrna$Final
Idents(dfrna) <- "Final"

if (start <= 0 & end > 0) {
  
  dfatac <- AddMetaData(dfatac,metadata = rnacluster,col.name="rnacluster")
  
  print("cluster_peaks: CallPeaks... (cluster-specific)")
  macs2_path <- '/home/amarder/anaconda3/envs/colocalization/bin/macs2'
  dfatac <- dfatac[,dfatac$rnacluster!="No markers"]
  # table(dfatac$rnacluster)
  cluster_peaks <- CallPeaks(
    object = dfatac,
    macs2.path = macs2_path,
    group.by = "rna"
  )
  
  grange.use <- seqnames(cluster_peaks) %in% standardChromosomes(cluster_peaks)
  cluster_peaks <- cluster_peaks[as.vector(grange.use), ]
  
  saveRDS(cluster_peaks,paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET,"/RNA_cluster.peaks.rds"))
  
}

if (start <= 1) {
  
  cluster_peaks <- readRDS(paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/",DATASET,"/RNA_cluster.peaks.rds"))

}

if (start<=1 & end >= 1) {
  print("create_feature_matrix...")
  dfcts <- create_feature_matrix(frags,cellnames,cluster_peaks)
  print("create_chrom_assay...")
  dfassay <- create_chrom_assay(dfcts,frags)
  print("create_seurat_object...")
  dfseurat <- create_seurat_object(dfassay)
  print("set_dataset...")
  dfseurat <- set_dataset(dfseurat,samp_ids)
  
  print("Saving seurat...")
  save_seurat_object(samp_lst,dfseurat,file_name="RNA_cluster.peaks.seurat.rds")
}

if (start==2 & end >= 2) {
  print("Reading data at start point 3...")
  dfseurat <- read_seurat_object(samp_lst,file_name="RNA_cluster.peaks.seurat.rds")
}

if (start<=2 & end >= 2) {
  # 3: combine into single dataset and batch correct. save peak calls within each cluster
  print("Merging...")
  dfcombined <- merge(dfseurat[[1]],
                      dfseurat[2:length(dfseurat)],
                      add.cell.ids=NULL)
  
  # print("Performing QC...")
  # dfcombined <- TSSEnrichment(dfcombined)
  # dfcombined <- subset(x=dfcombined,subset = TSS.enrichment > 3 & nCount_ATAC >= (1e3))
  
  print("FindTopFeatures...")
  dfcombined <- FindTopFeatures(dfcombined,min.cutoff = 10,assay="ATAC")
  print("RunTFIDF...")
  dfcombined <- RunTFIDF(dfcombined)
  print("RunSVD...")
  dfcombined <- RunSVD(dfcombined)
  print("RunHarmony...")
  lambda_val=1; dfcombined <- RunHarmony(dfcombined,"dataset",reduction="lsi",assay.use="ATAC",project.dim=FALSE,lambda=lambda_val)
  
  # 4: cluster
  print("RunUMAP...")
  dfcombined <- RunUMAP(dfcombined, reduction = "harmony", dims = 2:30)
  
  print("Creating UMAP with dataset-based coloring...")
  f.out <- paste0(dir,"/output/data/",DATASET,"/RNA_cluster.umap.dataset_integration.png")
  png(filename=f.out,width = 5000,height=2600,res=500)
  p1 <- DimPlot(dfcombined, group.by = "dataset")
  print(p1 & theme(plot.title = element_text(hjust = 0.5)))
  dev.off()
  
  # check batch correction:
  print("Creating split-UMAP with dataset-based coloring...")
  p1 <- DimPlot(dfcombined, split.by = "dataset",ncol = 2)
  f.out <- paste0(dir,"/output/data/",DATASET,"/RNA_cluster.umap.dataset_integration.split.png")
  png(filename=f.out,width = 5000,height=5000,res=400)
  print(p1 & theme(plot.title = element_text(hjust = 0.5)) + NoLegend())
  dev.off()
  
  dfcombined <- AddMetaData(dfcombined,metadata = rnacluster,col.name="rnacluster")
  
  print("Creating UMAP with RNA clusters...")
  f.out <- paste0(dir,"/output/data/",DATASET,"/RNA_cluster.umap.RNA_clusters.png")
  png(filename=f.out,width = 5000,height=4000,res=500)
  p1 <- DimPlot(dfcombined, group.by = "rnacluster") & theme(plot.title = element_text(hjust = 0.5)) 
  print(LabelClusters(plot = p1, id = "rnacluster"))
  dev.off()
  
  print("Saving data before start point 4...")
  f.out <- paste0(dir,"/output/data/",DATASET,"/RNA_cluster.ATAC.rds")
  saveRDS(dfcombined,file = f.out)
  
  
}

if (start==3) {
  print("Reading data at start point 4...")
  f <- paste0(dir,"/output/data/",DATASET,"/RNA_cluster.ATAC.rds")
  dfcombined <- readRDS(file = f)
}

# if (start<=4 & end>= 4) {
#   print('Calculating GeneActivity...')
#   dfcombined[['GeneActivity']] <- CreateAssayObject(counts=GeneActivity(dfcombined))
#   print("Normalizing data...")
#   dfcombined <- NormalizeData(
#     object = dfcombined,
#     assay = 'GeneActivity',
#     normalization.method = 'LogNormalize',
#     scale.factor = median(dfcombined$nCount_GeneActivity)
#   )
#   DefaultAssay(dfcombined) <- 'GeneActivity'
#   print("FindAllMarkers (for GeneActivity)...")
#   res.all <- FindAllMarkers(dfcombined)
#   f.out <- paste0(dir,"/output/data/",DATASET,"/FindAllMarkers.GeneActivity.txt")
#   print(paste0("Saving to file: ",f.out))
#   print("...")
#   fwrite(res.all,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
#   
#   print("Saving round 2 results + GeneActivity...")
#   f.out <- paste0(dir,"/output/data/",DATASET,"/round2_FindClusters.GeneActivity.rds")
#   saveRDS(dfcombined,file = f.out)
#   
# }
# 
# if (start==9) {
#   print("Loading round 2 results + GeneActivity...")
#   f <- paste0(dir,"/output/data/",DATASET,"/round2_FindClusters.GeneActivity.rds")
#   dfcombined <- readRDS(file = f)
# }
# 
# if (start<=9 & end>= 9) {
#   
#   DefaultAssay(dfcombined) <- 'ATAC'
#   
#   print("Beginning chromvar analysis...")
#   print("Get a list of motif position frequency matrices from the JASPAR database...")
#   pfm <- getMatrixSet(
#     x = JASPAR2020,
#     opts = list(species = 9606, all_versions = FALSE)
#   )
#   print("add motif information...")
#   dfcombined <- AddMotifs(
#     object = dfcombined,
#     genome = BSgenome.Hsapiens.UCSC.hg38,
#     pfm = pfm
#   )
#   print("RunChromVAR...")
#   DefaultAssay(dfcombined) <- 'ATAC'
#   dfcombined <- RunChromVAR(
#     object = dfcombined,
#     genome = BSgenome.Hsapiens.UCSC.hg38
#   )
#   
#   print("Saving round 2 results + GeneActivity + ChromVAR...")
#   f.out <- paste0(dir,"/output/data/",DATASET,"/round2_FindClusters.GeneActivity_ChromVAR.rds")
#   saveRDS(dfcombined,file = f.out)
#   
#   print("Normalizing data...")
#   dfcombined <- NormalizeData(
#     object = dfcombined,
#     assay = 'chromvar',
#     normalization.method = 'LogNormalize',
#     scale.factor = median(dfcombined$nCount_chromvar)
#   )
#   
#   print("FindAllMarkers (for ChromVAR)...")
#   DefaultAssay(dfcombined) <- 'chromvar'
#   res.all <- FindAllMarkers(dfcombined)
#   f.out <- paste0(dir,"/output/data/",DATASET,"/FindAllMarkers.ChromVAR.v1.txt")
#   print(paste0("Saving to file: ",f.out))
#   print("...")
#   fwrite(res.all,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
#   
#   print("getMatrixByID...")
#   genelst <- getMatrixByID(JASPAR2020, ID = unique(res.all$gene))
#   
#   print("motif_name...")
#   motif_name <- as.character(
#     unlist(
#       mclapply(
#         genelst,
#         function(x) ID(x),
#         mc.cores = 4)
#     )
#   )
#   
#   print("gene_name...")
#   gene_name <- as.character(
#     unlist(
#       mclapply(
#         genelst,
#         function(x) name(x),
#         mc.cores = 4)
#     )
#   )
#   
#   print("motif-to-gene linking...")
#   res.all$i <- 1:nrow(res.all)
#   res.all.2 <- merge(res.all,
#                      data.frame(motif_name,gene_name),
#                      by.x='gene',
#                      by.y='motif_name')
#   
#   print("reordering...")
#   res.all.2 <- res.all.2[order(res.all.2$i),]
#   
#   print(paste0("Saving to file: ",f.out))
#   f.out <- paste0(dir,"/output/data/",DATASET,"/FindAllMarkers.ChromVAR.txt")
#   fwrite(res.all.2,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
#   
#   print("Saving round 2 results + GeneActivity + ChromVAR...")
#   f.out <- paste0(dir,"/output/data/",DATASET,"/round2_FindClusters.GeneActivity_ChromVAR.rds")
#   saveRDS(dfcombined,file = f.out)
#   
# }


