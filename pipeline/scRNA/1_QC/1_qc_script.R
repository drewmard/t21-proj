# Script based on Mike Nelson's 1_mergingAndQualityControl.py script.
# Exception: doublets were identified and removed prior to merging

RETICULATE_PYTHON.path="/home/amarder/anaconda3/envs/minimal_env/bin/python"
RETICULATE_PYTHON.path="/scg/apps/software/r/4.1.2/bin/python"
Sys.setenv(RETICULATE_PYTHON=RETICULATE_PYTHON.path)
library(reticulate)
# use_condaenv(condaenv="minimal_env", required=TRUE)
scrub <<- reticulate::import(module = "scrublet")
library(Seurat)
library(data.table)

seqDataType="scRNA-seq"
print("read_10X_and_saveRNAcounts...")
# dflst <- list()

fpath <- "/oak/stanford/groups/smontgom/amarder/data/t21/Cellranger"
flist <- list.files(fpath)
flist=read.table(text="cellranger310_count_H2NVKDSX2_and_HFVYKDSX2_L15593_E_GRCh38-3_0_0
cellranger310_count_H2NVKDSX2_and_HFVYKDSX2_L15593_F_GRCh38-3_0_0
cellranger310_count_H2NVKDSX2_and_HFVYKDSX2_L15633_A_GRCh38-3_0_0
cellranger310_count_H2NVKDSX2_and_HFVYKDSX2_L15636_B_GRCh38-3_0_0
cellranger310_count_H2NVKDSX2_and_HFVYKDSX2_L15636_C_GRCh38-3_0_0
cellranger310_count_H2NVKDSX2_and_HFVYKDSX2_L15636_D_GRCh38-3_0_0
cellranger310_count_H2NVKDSX2_and_HFVYKDSX2_T21_15619_A_GRCh38-3_0_0
cellranger310_count_H2NVKDSX2_and_HFVYKDSX2_T21_15619_B_GRCh38-3_0_0
cellranger310_count_H2NVKDSX2_and_HFVYKDSX2_T21_15619_C_GRCh38-3_0_0
cellranger310_count_H2NVKDSX2_and_HFVYKDSX2_T21_15619_E_GRCh38-3_0_0
cellranger310_count_H3FGTDSX2_and_HFYGWDSX2_15646A_GRCh38-3_0_0
cellranger310_count_H3FGTDSX2_and_HFYGWDSX2_15646B_GRCh38-3_0_0
cellranger310_count_H3FGTDSX2_and_HFYGWDSX2_15656D_GRCh38-3_0_0
cellranger310_count_H3FGTDSX2_and_HFYGWDSX2_15667E_GRCh38-3_0_0
cellranger310_count_H3FGTDSX2_and_HFYGWDSX2_15667F_GRCh38-3_0_0
cellranger310_count_H3FGTDSX2_and_HFYGWDSX2_15669G_GRCh38-3_0_0
cellranger310_count_H3FGTDSX2_and_HFYGWDSX2_15669H_GRCh38-3_0_0
cellranger310_count_H3FGTDSX2_and_HFYGWDSX2_15712_A_GRCh38-3_0_0
cellranger310_count_H3FGTDSX2_and_HFYGWDSX2_15712_B_GRCh38-3_0_0
cellranger310_count_H3FGTDSX2_and_HFYGWDSX2_F15633_C_GRCh38-3_0_0
cellranger310_count_H3FGTDSX2_and_HFYGWDSX2_F15669F_GRCh38-3_0_0
cellranger310_count_H5KT2DSX2_and_HJJ3FDRXY_F15657G_GRCh38-3_0_0
cellranger310_count_HCCVYDSX2_and_HJJ3FDRXY_F15781B_GRCh38-3_0_0
cellranger310_count_HFVYKDSX2_and_37323_T21_15582_E_GRCh38-3_0_0
cellranger310_count_HFVYKDSX2_and_37323_T21_15582_F_GRCh38-3_0_0
")[,1]

samp_lst <- paste0(fpath,"/",flist)

if (seqDataType=="scRNA-seq") {
  topDir="outputEmptyDrops"
} else if (seqDataType=="10X Multiome") {
  topDir="filtered_feature_bc_matrix"
}

# need to do 40, 45, 52, 53

for (i in c(1:length(samp_lst))) {
  samp <- samp_lst[[i]]
  f.out=paste0(samp,"/",topDir,"/barcodes.tsv.gz")
  if (!file.exists(f.out)) {print(i)}
}

for (i in c(1:length(samp_lst))) {
# for (i in c(85,86,88)) {
  tryCatch({
  samp <- samp_lst[[i]]
  print(paste0("Reading data from ",i,"/",length(samp_lst),": ",samp," ..."))
  f=paste0(samp,"/",topDir)
  dflst <- Read10X(f)
  if (seqDataType=="scRNA-seq") {
    dfsample <- CreateSeuratObject(counts = dflst)
  } else if (seqDataType=="10X Multiome") {
    dfsample <- CreateSeuratObject(counts = dflst[["Gene Expression"]])
  }
  
  print("Estimating mito %...")
  dfsample[["percent.mt"]] <- PercentageFeatureSet(dfsample, pattern = "^MT-")
  
  print("Performing QC...")
  nFeature_RNA.limits=c(250,8500)
  nCount_RNA.limits=c(750,110000)
  percent.mt.limits=c(-Inf,20)
  QC_data = data.frame(
    nFeature_RNA.lower = sum(dfsample$nFeature_RNA < nFeature_RNA.limits[1]),
    nFeature_RNA.upper = sum(dfsample$nFeature_RNA > nFeature_RNA.limits[2]),
    nCount_RNA.lower = sum(dfsample$nCount_RNA < nCount_RNA.limits[1]),
    nCount_RNA.upper = sum(dfsample$nCount_RNA > nCount_RNA.limits[2]),
    percent.mt.lower = sum(dfsample$percent.mt < percent.mt.limits[1]),
    percent.mt.upper = sum(dfsample$percent.mt > percent.mt.limits[2])
  )
  QC_data2 = QC_data / nrow(dfsample)
  f.out=paste0(samp,"/QC_sum.txt")
  fwrite(QC_data,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
  f.out=paste0(samp,"/QC_prop.txt")
  fwrite(QC_data2,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
  
  dfsample <- subset(dfsample, subset = (nFeature_RNA > nFeature_RNA.limits[1] & nFeature_RNA < nFeature_RNA.limits[2]) &
                       (nCount_RNA > nCount_RNA.limits[1] & nCount_RNA < nCount_RNA.limits[2]) & 
                       percent.mt > percent.mt.limits[1] & percent.mt < percent.mt.limits[2])
  
  
  ############
  
  print("Doublets...")
  expected_doublet_rate=0.06
  scr <- eval(rlang::expr(scrub$Scrublet(
    counts_matrix = reticulate::r_to_py(Seurat::GetAssayData(object = dfsample, slot = "counts"))$T$tocsc(),
    expected_doublet_rate = expected_doublet_rate)))
  
  doublet_results <- 
    scr$scrub_doublets(
      min_counts=2,
      min_cells=3,
      min_gene_variability_pctl=85,
      n_prin_comps=as.integer(30)
    )
  
  doublet_score <- doublet_results[[1]]
  names(doublet_score) <- colnames(dfsample)
  doublet_prediction <- doublet_results[[2]]
  names(doublet_prediction) <- colnames(dfsample)
  
  dfsample = dfsample[,names(doublet_prediction)[!(doublet_prediction)]]
  
  print("Saving seurat object...")
  f.out=paste0(samp,"/seurat_obj.rds")
  saveRDS(dfsample,file = f.out)
  
  },error=function(e) {print(paste0("FAILURE: ",i))})
}
