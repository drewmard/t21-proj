# Script based on Mike Nelson's script.

# Specify a CRAN mirror
cranMirror <- "https://cran.ma.imperial.ac.uk/"
options(warn=-1)

library(knitr)

opts_chunk$set(error=FALSE, message=FALSE, warning=FALSE)
opts_chunk$set(dpi=300, dev="png", dev.args=list(pointsize=15))
options(bitmapType="cairo")

library(DropletUtils)
library(data.table)

# Local and global paths, where the local input is passed as an argument 
globalPath = "/oak/stanford/groups/smontgom/amarder/data/t21/Cellranger/"
pltPath <- "/oak/stanford/groups/smontgom/amarder/t21-proj/out/combined/0_emptyDrop/"
dir.create(file.path(pltPath)) 

# flist <- list.files(globalPath)
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

args <- commandArgs(TRUE)
i=(args[1])
i=as.numeric(i)
print("running...")
sampleName = flist[i]
print(i)
print(sampleName)
# for (sampleName in flist[1]) { 
samplePath <- paste(sampleName, "/raw_feature_bc_matrix/", sep="")
fname <- paste(globalPath, samplePath, sep="")

# Read from the 10X counts matrix
print("Read from the 10X counts matrix")
sce   <- read10xCounts(fname, col.names=TRUE)
set.seed(100)

# Cell barcode ranking
print("Cell barcode ranking")
bcrank = NULL
# Alternative ranking if a lower bound is specified
if (grepl("LowerBound", sampleName, fixed=TRUE)) {
  print("Analysing lower bound inputs")
  bcrank = barcodeRanks(counts(sce), lower=500)
} else {
  bcrank = barcodeRanks(counts(sce))
}

# Running EmptyDrops
print("Running EmptyDrops")
e.out  <- emptyDrops(counts(sce))
# Keep cell rows of the matix with low false discovery rate
is.cell = (e.out$FDR <= 0.001)

print("Plotting...")
# Only showing unique points for plotting speed.
uniq <- duplicated(bcrank$rank)
plt1 <- paste(pltPath, "Rank_vs_UMI_", sampleName, ".png", sep="")
png(file=plt1)
par(mar=c(5,4,2,1), bty="n")
plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy", xlab="Rank", ylab="Total UMI count", cex=0.5, cex.lab=1.2)
abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2)
abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2)
legend("left", legend=c("Inflection", "Knee"), bty="n", col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)
dev.off()

# Plotting after cell preselection
plt2 <- paste(pltPath, "UMI_vs_LogProb_", sampleName, ".png", sep="")
png(file=plt2)
par(mar=c(5,4,1,1), mfrow=c(1,2), bty="n")
plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"), xlab="Total UMI count", ylab="-Log Probability", cex=0.5)
abline(v = metadata(bcrank)$inflection, col="darkgreen")
abline(v = metadata(bcrank)$knee, col="dodgerblue")
legend("bottomright", legend=c("Inflection", "Knee"), bty="n", col=c("darkgreen", "dodgerblue"), lty=1, cex=1.2)
dev.off()

print("Subset...")
w2kp = which(is.cell & e.out$Total >= metadata(bcrank)$inflection)
sce = sce[,w2kp]

subdirPath <- paste(sampleName, "/outputEmptyDrops", sep="")
fname2 <- paste(globalPath, subdirPath, sep="")

if (dir.exists(fname2)) {
  print("Folder /outputEmptyDrops already exist!")
} else {
  print("Creating folder /outputEmptyDrops...")
  dir.create(file.path(fname2)) 
}

print("Saving...")
write10xCounts(fname2,
               sce@assays@data$counts,
               barcodes    = colData(sce)$Barcode,
               gene.id     = rowData(sce)$ID,
               gene.symbol = rowData(sce)$Symbol,
               version     = "3",
               overwrite   = TRUE)

# }
