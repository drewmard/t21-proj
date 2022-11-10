library(data.table)
dir="/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/cpdb_Results"
# disease_status="Healthy"
disease_status="DownSyndrome"
envir="Liver"
print(paste0(dir,"/",disease_status,"_",envir))
dfmu = fread(paste0(dir,"/",disease_status,"_",envir,"/means.txt"),data.table = F,stringsAsFactors = F)
dfp = fread(paste0(dir,"/",disease_status,"_",envir,"/pvalues.txt"),data.table = F,stringsAsFactors = F)

# Get the number of cell types
# N_cellTypes = sqrt(ncol(dfmu)-11) # may be wrong

keep <- !(grepl("HLA",dfmu$interacting_pair) | grepl("complex",dfmu$interacting_pair) | dfmu$gene_a=="" | dfmu$gene_b=="" | apply(dfmu[,c("receptor_a","receptor_b")],1,sum) != 1)
# removing interacting_pairs like: 12oxoLeukotrieneB4_byPTGR1_LTB4R
# sometimes both genes not annotated as a receptor: SFRP5 or WNT6
dfmu <- dfmu[keep,]

# keep <- !(grepl("HLA",dfp$interacting_pair) | grepl("complex",dfp$interacting_pair)) # same as dfmu
dfp <- dfp[keep,]

rng <- 12:ncol(dfp)

dfmu.melt <- reshape2::melt(dfmu[,c(5,6,8,9,rng)],id.vars=c("gene_a","gene_b","receptor_a","receptor_b"))
x = strsplit(as.character(dfmu.melt$variable),"\\|")
dfmu.melt$cell1 = as.character(unlist(lapply(x,function(y){y[1]})))
dfmu.melt$cell2 = as.character(unlist(lapply(x,function(y){y[2]})))
dfmu.melt = subset(dfmu.melt,!(cell1==cell2))
dfmu.melt$Ligand = ifelse(!dfmu.melt$receptor_b,dfmu.melt$gene_b,dfmu.melt$gene_a)
dfmu.melt$Receptor = ifelse(dfmu.melt$receptor_b,dfmu.melt$gene_b,dfmu.melt$gene_a)
dfmu.melt$Sender = ifelse(!dfmu.melt$receptor_b,dfmu.melt$cell2,dfmu.melt$cell1)
dfmu.melt$Receiver = ifelse(dfmu.melt$receptor_b,dfmu.melt$cell2,dfmu.melt$cell1)
dfmu.melt = dfmu.melt[,c("Ligand","Receptor","Sender","Receiver","value")]
colnames(dfmu.melt)[5] <- "MeanExpr"

# # check to make sure it worked:
# subset(dfmu.melt,Ligand %in% c("CCL4","SLC7A1") & Receptor %in% c("CCL4","SLC7A1") &
#          Sender%in%c("Activated stellate cells","B cells") & Receiver%in%c("Activated stellate cells","B cells") )

# subset(dfp,gene_a=="CCL3L1" & gene_b=="CCR1")
dfp.melt <- reshape2::melt(dfp[,c(5,6,8,9,rng)],id.vars=c("gene_a","gene_b","receptor_a","receptor_b"))
# subset(dfp.melt,gene_a=="CCL3L1" & gene_b=="CCR1" & variable=="Inflammatory macrophages|Kupffer cells")
# subset(dfp.melt,gene_a=="CCR1" & gene_b=="CCL3L1" & variable=="Inflammatory macrophages|Kupffer cells")
x = strsplit(as.character(dfp.melt$variable),"\\|")
dfp.melt$cell1 = as.character(unlist(lapply(x,function(y){y[1]})))
dfp.melt$cell2 = as.character(unlist(lapply(x,function(y){y[2]})))
dfp.melt = subset(dfp.melt,!(cell1==cell2))
dfp.melt$Ligand = ifelse(!dfp.melt$receptor_b,dfp.melt$gene_b,dfp.melt$gene_a)
dfp.melt$Receptor = ifelse(dfp.melt$receptor_b,dfp.melt$gene_b,dfp.melt$gene_a)
dfp.melt$Sender = ifelse(!dfp.melt$receptor_b,dfp.melt$cell2,dfp.melt$cell1)
dfp.melt$Receiver = ifelse(dfp.melt$receptor_b,dfp.melt$cell2,dfp.melt$cell1)
dfp.melt = dfp.melt[,c("Ligand","Receptor","Sender","Receiver","value")]
colnames(dfp.melt)[5] <- "P"

# subset(dfp.melt,Ligand=="CCL3L1" & Receptor=="CCR1" & Sender=="Inflammatory macrophages" & Receiver=="Kupffer cells")

df.mg <- merge(dfmu.melt,dfp.melt,by=c("Ligand","Receptor","Sender","Receiver"))
df.mg <- df.mg[order(df.mg$P,-1*df.mg$MeanExpr),]
df.mg <- subset(df.mg,P < 0.05)

# there are duplicated interactions due to curation from multiple annotation strategies. 
# for example, "annotation_strategy" is "curated" and "I2D,InnateDB-All"
# it got duplicated into 4 in the final df.mg output because i needed to merge a mean expr and p value table which each had 2 rows, and merges on all possibilities so 2x2 = 4 rows
df.mg = df.mg[!duplicated(df.mg),] # 

f.out = paste0(dir,"/",disease_status,"_",envir,"/",disease_status,"_",envir,".processed_output.txt")
fwrite(df.mg,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)


