print("Reading ATAC data...")
# f <- paste0(dir,"/output/data/",DATASET,"/round2_FindClusters.rds")
dir="/oak/stanford/groups/smontgom/amarder/neuro-variants"
f <- paste0(dir,"/output/data/","DS_Multiome_all","/round2_FindClusters.rds")
ATAC_df <- readRDS(file = f)

unique(dfseurat@meta.data[,c("dataset","disease")])

disease_label <- function(dataset) {
  dplyr::case_when(
    dataset == "15836 A" ~ "H",
    dataset == "15836 B" ~ "H",
    dataset == "16171" ~ "H",
    dataset == "16216 A" ~ "H",
    dataset == "16216 B" ~ "H",
    dataset == "16216 C" ~ "H",
    dataset == "15828 C" ~ "T21",
    dataset == "15828 D" ~ "T21",
    dataset == "15582 nuclei A" ~ "T21",
    dataset == "15582 nuclei B" ~ "T21",
    dataset == "15669_A" ~ "T21",
    dataset == "15669_B" ~ "T21",
    TRUE ~ NA_character_
  )
}

ATAC_df@meta.data$disease <- disease_label(ATAC_df@meta.data$dataset)

print("Creating UMAP with dataset-based coloring (round 1)...")
p1 <- DimPlot(ATAC_df, group.by = "disease")
f.plot = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scavenge_plots/umap.disease.3.pdf")
pdf(f.plot,width=9,height=7)
print(p1)
dev.off()

ATAC_df@meta.data$cell = unlist(lapply(strsplit(rownames(ATAC_df@meta.data),"_"),function(x) x[[1]]))
ATAC_df@meta.data$cell_dataset = paste(ATAC_df@meta.data$cell,ATAC_df@meta.data$dataset)
ind = which(ATAC_df@meta.data$cell_dataset %in% dfseurat@meta.data$cell_dataset)
ATAC_df = ATAC_df[,ind]

print("Creating UMAP with dataset-based coloring (round 1)...")
p1 <- DimPlot(ATAC_df, group.by = "disease")
f.plot = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scavenge_plots/umap.disease.4.pdf")
pdf(f.plot,width=9,height=7)
print(p1)
dev.off()

###########3

print("RunTFIDF...")
ATAC_df <- RunTFIDF(ATAC_df)
print("FindTopFeatures...")
ATAC_df <- FindTopFeatures(ATAC_df,min.cutoff = 'q0',assay="ATAC")
print("RunSVD...")
ATAC_df <- RunSVD(ATAC_df)
print("RunHarmony...")
lambda_val=1; ATAC_df <- RunHarmony(ATAC_df,"dataset",reduction="lsi",assay.use="ATAC",project.dim=FALSE,lambda=lambda_val,dims.use=2:30)

# 4: cluster
print("RunUMAP...")
ATAC_df <- RunUMAP(ATAC_df, reduction = "harmony",dims=1:50)

##########

print("Creating UMAP with dataset-based coloring (round 1)...")
p1 <- DimPlot(ATAC_df, group.by = "disease")
f.plot = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scavenge_plots/umap.disease.5.pdf")
pdf(f.plot,width=9,height=7)
print(p1)
dev.off()

umap_mat = ATAC_df@reductions$umap@cell.embeddings
tmp = data.frame(names=rownames(dfseurat@meta.data),
                 cell=dfseurat@meta.data$cell,
                 disease=dfseurat@meta.data$disease,
                 subclust_v6=dfseurat@meta.data$subclust_v6,
                 UMAP_1=umap_mat[,1],
                 UMAP_2=umap_mat[,2])
library(data.table)
fwrite(tmp,"/oak/stanford/groups/smontgom/amarder/tmp/umap.txt",
       row.names = F,col.names = T,na = "NA",sep = "\t",quote = F)

ATAC_df@reductions$harmony@
# mutualknn30 <- getmutualknn(lsi_mat, 30)
mutualknn30 <- getmutualknn(ATAC_df@reductions$harmony@cell.embeddings, 30)
mutualknn30.v2 <- getmutualknn(ATAC_df@reductions$harmony@cell.embeddings, 30)

# transform seurat to summarizedexperiment
SE_data <- SummarizedExperiment(assays = list(counts = ATAC_df@assays$ATAC@counts),
                                rowData = ATAC_df@assays$ATAC@ranges, 
                                colData = DataFrame(names = colnames(dfseurat@assays$ATAC),dataset=dfseurat@meta.data$dataset))

SE_data <- addGCBias(SE_data, genome = BSgenome.Hsapiens.UCSC.hg38)

# preprocess gchromvar
SE_data_bg <- getBackgroundPeaks(SE_data, niterations=200)




