sampletype="Femur"
cell_type="Late erythroid cells"
subset_column="sample"

print(sampletype)

# 
iter = 0
df.aggre.stats <- list()
for (sampletype in c("Liver","Femur")) {
  
  print("Reading metadata1...")
  disease_status="Healthy"
  f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/cellComp/10X_",disease_status,"_",sampletype,".cellComp.csv")
  meta1.full<-fread(f,data.table = F,stringsAsFactors = F)
  cells1 <- unique(meta1.full[,6])
  meta1 <- unique(meta1.full[,c("patient","sample","sorting")])
  meta1$environment <- disease_status
  
  #
  print("Reading metadata2...")
  disease_status="DownSyndrome"
  f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/cellComp/10X_",disease_status,"_",sampletype,".cellComp.csv")
  meta2.full<-fread(f,data.table = F,stringsAsFactors = F)
  cells2 <- unique(meta2.full[,6])
  meta2 <- unique(meta2.full[,c("patient","sample","sorting")])
  meta2$environment <- disease_status
  
  colnames(meta1.full)[6] <- "leiden_names"
  colnames(meta2.full)[6] <- "leiden_names"
  meta.all.full <- rbind(meta1.full,meta2.full)
  
  #
  print("Merge metadata...")
  x <- rbind(meta1,meta2)
  rownames(x) <- x[,subset_column]
  x$sorting[!(x$sorting %in% c("CD235a-","CD45+"))] <- "Other"
  
  
  print("cell types of interest...")
  clusters_for_DE <- cells1[cells1 %in% cells2]
  P <- length(clusters_for_DE)
  
  cell_type="HSCs/MPPs"
  
  for (cell_type in clusters_for_DE) { 
    #
    iter = iter + 1
    
    cell_type_filename = gsub("/","_",cell_type)
    
    f=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data_pb_leiden/",sampletype,".pb.",cell_type_filename,".",subset_column,".txt")
    df.aggre <- fread(f,data.table = F,stringsAsFactors = F)
    rownames(df.aggre) <- df.aggre[,1]; df.aggre <- df.aggre[,-1]
    
    metadata_to_use <- x[rownames(x) %in% colnames(df.aggre),]
    to_merge <- names(table(metadata_to_use$sorting)[table(metadata_to_use$sorting) < 5])
    df.aggre <- as.matrix(df.aggre[,match(rownames(metadata_to_use),colnames(df.aggre))])
    
    tab <- table(meta.all.full[meta.all.full[,"leiden_names"]==cell_type,'sample'])
    samples_to_keep <- names(tab)[tab >= 10]
    samples_to_keep <- samples_to_keep[samples_to_keep %in% rownames(metadata_to_use)]
    df.aggre <- df.aggre[,samples_to_keep]
    metadata_to_use <- metadata_to_use[samples_to_keep,]
    
    df.aggre.stats[[iter]] <- data.frame(sampletype,cell_type,
                                         N_sample=ncol(df.aggre),N_cell=sum(tab[tab>=10]),avg_reads=mean(apply(df.aggre,2,sum)))
    
    
  }
  
}

df.aggre.stats.all <- as.data.frame(do.call(rbind,df.aggre.stats))

sampletype <- 'Liver'
res.df.all.p <- readRDS(paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".p.rds"))
res.df.all.lfc <- readRDS(paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".lfc.rds"))
cell_type_groups <- colnames(res.df.all.p[[1]])[(which(colnames(res.df.all.p[[1]])=="names")+1):(which(colnames(res.df.all.p[[1]])=="chromosome_name")-1)]

for (cell_type in cell_type_groups) {
  keep_genes <- res.df.all.p[[2]][!is.na(res.df.all.p[[2]][,cell_type]),'names']
  res.df.all.p[[1]][!(res.df.all.p[[1]][,'names']%in%keep_genes),cell_type] <- NA
  res.df.all.lfc[[1]][!(res.df.all.lfc[[1]][,'names']%in%keep_genes),cell_type] <- NA
}

tab <- lapply(res.df.all.p,function(x) {aggregate(x[,cell_type_groups]<0.1,by=list(x$chr21),mean,na.rm=T)})
tab <- lapply(tab,function(x) {rownames(x) <- x[,1]; x <- x[,-1]; x <- as.data.frame(t(x))})
tab <- do.call(cbind,tab)
tab <- tab[,c(3,4)]
tab$sampletype <- sampletype
tab$cell_type <- rownames(tab)
rownames(tab) <- NULL
tab1 <- tab

sampletype <- 'Femur'
res.df.all.p <- readRDS(paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".p.rds"))
res.df.all.lfc <- readRDS(paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".lfc.rds"))
cell_type_groups <- colnames(res.df.all.p[[1]])[(which(colnames(res.df.all.p[[1]])=="names")+1):(which(colnames(res.df.all.p[[1]])=="chromosome_name")-1)]

for (cell_type in cell_type_groups) {
  keep_genes <- res.df.all.p[[2]][!is.na(res.df.all.p[[2]][,cell_type]),'names']
  res.df.all.p[[1]][!(res.df.all.p[[1]][,'names']%in%keep_genes),cell_type] <- NA
  res.df.all.lfc[[1]][!(res.df.all.lfc[[1]][,'names']%in%keep_genes),cell_type] <- NA
}

tab <- lapply(res.df.all.p,function(x) {aggregate(x[,cell_type_groups]<0.1,by=list(x$chr21),mean,na.rm=T)})
tab <- lapply(tab,function(x) {rownames(x) <- x[,1]; x <- x[,-1]; x <- as.data.frame(t(x))})
tab <- do.call(cbind,tab)
tab <- tab[,c(3,4)]
tab$sampletype <- sampletype
tab$cell_type <- rownames(tab)
rownames(tab) <- NULL
tab2 <- tab

df.aggre.stats.all$cell_type <- gsub("/","_",df.aggre.stats.all$cell_type)

df.aggre.stats.all.mg <- merge(df.aggre.stats.all,rbind(tab1,tab2),by=c("sampletype","cell_type"))
subset(df.aggre.stats.all.mg,cell_type=="HSCs_MPPs")
subset(df.aggre.stats.all.mg,cell_type=="Late erythroid cells")

cor(df.aggre.stats.all.mg[,-c(1:2)])

df.aggre.stats.all.mg.sub <- subset(df.aggre.stats.all.mg,sampletype=="Femur")
cor(df.aggre.stats.all.mg.sub[,-c(1:2)])
cor(df.aggre.stats.all.mg.sub[,"Chr 21"],df.aggre.stats.all.mg.sub[,"avg_reads"]*df.aggre.stats.all.mg.sub[,"N_sample"])

df.aggre.stats.all.mg.sub[,c("Chr 21","N_sample","cell_type")]
summary(lm(`Chr 21`~N_sample + sampletype,df.aggre.stats.all.mg))

# cor.test(df.aggre.stats.all.mg.sub[,'N_sample'],df.aggre.stats.all.mg.sub[,'Not Chr 21'])
cor.test(df.aggre.stats.all.mg.sub[,'N_sample'],df.aggre.stats.all.mg.sub[,'Chr 21'])

cor.test(df.aggre.stats.all.mg[,'N_sample'],df.aggre.stats.all.mg[,'Chr 21'])

df.aggre


