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

flist <- list.files(globalPath)
args <- commandArgs(TRUE)
sampleName = flist[args[1]]
# for (sampleName in flist[1]) { 
samplePath <- paste(sampleName, "/raw_feature_bc_matrix/", sep="")
fname <- paste(globalPath, samplePath, sep="")

# Read from the 10X counts matrix
sce   <- read10xCounts(fname, col.names=TRUE)
set.seed(100)

# Cell barcode ranking
bcrank = NULL
# Alternative ranking if a lower bound is specified
if (grepl("LowerBound", sampleName, fixed=TRUE)) {
  print("Analysing lower bound inputs")
  bcrank = barcodeRanks(counts(sce), lower=500)
} else {
  bcrank = barcodeRanks(counts(sce))
}

# Running EmptyDrops
e.out  <- emptyDrops(counts(sce))
# Keep cell rows of the matix with low false discovery rate
is.cell = (e.out$FDR <= 0.001)

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

write10xCounts(fname2,
               sce@assays@data$counts,
               barcodes    = colData(sce)$Barcode,
               gene.id     = rowData(sce)$ID,
               gene.symbol = rowData(sce)$Symbol,
               version     = "3",
               overwrite   = TRUE)

# }
