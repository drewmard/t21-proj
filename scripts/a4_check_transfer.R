rm(list = ls())

# setwd("~/Projects/DownSyndrome/downsyndrome/MultiomeAnalysis")

cmd.args <- commandArgs(trailingOnly = TRUE)
# bigRNA.path <- cmd.args[1]
bigRNA.path="/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/10X_Healthy_Liver.umap2d.cells_removed.rds"

# bigRNA.path <- "../../../DownSyndrome/Multiome_results/test-data/10X_Healthy_Liver_all.rds"
# transfered.smallRNA.path <- cmd.args[2]
multiomePath="/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h/"
output.prefix <- paste0(multiomePath,"multiome_transfer")
transfered.smallRNA.path <- paste0(output.prefix, ".RNA.rds")

# transfered.smallRNA.path <- "../../../../Google Drive/Shared drives/Cvejic-group/Data/Andrew/multiome_transfer.RNA.rds"
# fig.path <- cmd.args[3]
fig.path <- "/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h/bridging_result_analysis_test.pdf"
# fig.path <- "../../Multiome_results/single_scale/bridging_result_analysis_test.pdf"

library(Seurat)
library(magrittr)
library(stringr)

# big RNA-seq data with raw counts
bigRNA.res <- readRDS(bigRNA.path) 
bigRNA.res@assays$RNA@key <- "rna_"
bigRNA.res <- bigRNA.res %>%
  NormalizeData(normalization.method = "LogNormalize", 
                scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

# queried small RNA-seq data with normalised and scalled counts
transfered.res <- readRDS(transfered.smallRNA.path)

# retrieve original counts and re-do normalisation and HVG searching
merged.res <- merge(x = bigRNA.res, y = transfered.res,
                    merge.data = FALSE, 
                    add.cell.ids = c("big", "transfered.small"))
merged.res@meta.data$protocol <- str_extract(rownames(merged.res@meta.data), 
                                             pattern = "big|transfered.small")

merged.res$celltype <- NA
merged.res$celltype[merged.res$protocol == "big"] <- merged.res$leiden_v7[merged.res$protocol == "big"]
merged.res$celltype[merged.res$protocol == "small"] <- merged.res$predicted.id[merged.res$protocol == "transferred.small"]

res.list <- SplitObject(merged.res, split.by = "protocol")
res.list <- lapply(X = res.list, FUN = function(x) {
    x <- NormalizeData(x, normalization.method = "LogNormalize", 
                       scale.factor = 10000)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# integrate two datasets
anchors <- FindIntegrationAnchors(object.list = res.list, dims = 1:30)
integrated.res <- IntegrateData(anchorset = anchors, dims = 1:30)

# process integrated counts and plot UMAPs
DefaultAssay(integrated.res) <- "integrated"
integrated.res <- ScaleData(integrated.res) %>%
  RunPCA() %>%
  RunUMAP(dims = 1:30)

integrated.res$celltype <- NA
integrated.res$celltype[integrated.res$protocol == "big"] <- integrated.res$leiden_v7[integrated.res$protocol == "big"]
integrated.res$celltype[integrated.res$protocol == "transfered.small"] <- integrated.res$predicted.id[integrated.res$protocol == "transfered.small"]

pdf(file = fig.path, height = 9, width = 22)
plot1 <- DimPlot(integrated.res, reduction = "umap", group.by = "protocol")
plot2 <- DimPlot(integrated.res, reduction = "umap", group.by = "celltype")
plot1 + plot2
dev.off()
