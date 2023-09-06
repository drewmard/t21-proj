library(Seurat)
library(harmony)
library(ggplot2)
library(data.table)

###########

# PARAM:

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

i=0
dflst = list()
for (disease_status in c("Healthy","DownSyndrome")) {
  for (organ in c("Femur","Liver")) {
    i = i +1
    # f = (paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/meta.10X_",disease_status,"_",organ,".umap2d.cells_removed.txt"))
    # real file:
    f = (paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/meta.10X_",disease_status,"_",organ,".umap2d.cells_removed.v2.txt"))
    dflst[[i]] = fread(f,data.table = F,stringsAsFactors = F)
    colName1 <- paste0("leiden_v",max(as.numeric(substring(colnames(dflst[[i]])[grep("leiden_v",colnames(dflst[[i]]))],nchar("leiden_v")+1)),na.rm=T))
    dflst[[i]][,"leiden_names"] <- dflst[[i]][,colName1]
    dflst[[i]] <- dflst[[i]][,c("V1","sample","cell_type_groups","leiden_names")]
  }
}

meta = do.call(rbind,dflst)
meta$cell = unlist(lapply(strsplit(meta$V1,"-"),function(x) paste0(x[1],"-",x[2])))

sum(duplicated(paste(meta$cell,meta$dataset)))

meta$DUP = duplicated(paste(meta$cell,meta$dataset))
table(meta$sample,meta$DUP)

meta.orig <- dfcombined@meta.data
meta.orig$i <- 1:nrow(meta.orig)

meta.mg <- merge(meta.orig,meta,by.x=c("cell","dataset"),by.y=c("cell","sample"))

which(duplicated(paste(meta.mg$cell,meta.mg$dataset)))[1:3]

subset(meta.mg,cell=="AAACCCAAGCGTGTCC-1" & dataset=="F15593P")
subset(meta,cell=="AAACCCAAGCGTGTCC-1" & sample=="F15593P")
subset(meta.orig,cell=="AAACCCAAGCGTGTCC-1" & dataset=="F15593P")

dflst[[i]]$cell = unlist(lapply(strsplit(dflst[[i]]$V1,"-"),function(x) paste0(x[1],"-",x[2])))
subset(dflst[[i]],cell=="AAACCCAAGCGTGTCC-1" & sample=="F15593P")

dflst[[i]]

meta.mg[order(meta.m)]

f <- paste0(dir,"/out/combined/2_integrate_and_cluster/",DATASET,"/RNA_FindClusters.rds")
print(paste0("Reading data: ",f,"..."))
dfcombined = readRDS(file = f)

