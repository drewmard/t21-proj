sampletype="Liver"
res.df.all.p <- readRDS(paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".p.rds"))
res.df.all.lfc <- readRDS(paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".lfc.rds"))
cell_type_groups <- colnames(res.df.all.p[[1]])[(which(colnames(res.df.all.p[[1]])=="names")+1):(which(colnames(res.df.all.p[[1]])=="chromosome_name")-1)]

for (cell_type in cell_type_groups) {
  keep_genes <- res.df.all.p[[2]][!is.na(res.df.all.p[[2]][,cell_type]),'names']
  res.df.all.p[[1]][!(res.df.all.p[[1]][,'names']%in%keep_genes),cell_type] <- NA
  res.df.all.lfc[[1]][!(res.df.all.lfc[[1]][,'names']%in%keep_genes),cell_type] <- NA
}

res.df.all.lfc.liver <- res.df.all.lfc[[2]]
res.df.all.p.liver <- res.df.all.p[[2]]
cell_type_groups.liver <- cell_type_groups

sampletype="Femur"
res.df.all.p <- readRDS(paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".p.rds"))
res.df.all.lfc <- readRDS(paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".lfc.rds"))
cell_type_groups <- colnames(res.df.all.p[[1]])[(which(colnames(res.df.all.p[[1]])=="names")+1):(which(colnames(res.df.all.p[[1]])=="chromosome_name")-1)]

for (cell_type in cell_type_groups) {
  keep_genes <- res.df.all.p[[2]][!is.na(res.df.all.p[[2]][,cell_type]),'names']
  res.df.all.p[[1]][!(res.df.all.p[[1]][,'names']%in%keep_genes),cell_type] <- NA
  res.df.all.lfc[[1]][!(res.df.all.lfc[[1]][,'names']%in%keep_genes),cell_type] <- NA
}

res.df.all.lfc.femur <- res.df.all.lfc[[2]]
res.df.all.p.femur <- res.df.all.p[[2]]
cell_type_groups.femur <- cell_type_groups

cell_type_groups <- cell_type_groups.femur[cell_type_groups.femur %in% cell_type_groups.liver]

library(data.table)
fdir="/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_int/"
clusters_for_DE <- c("cDC2","Early erythroid cells"   ,
           "HSCs/MPPs","Inflammatory macrophages",
           "Late erythroid cells","Mast cells"    ,          
           "Megakaryocytes","MEMPs"     ,              
           "NK cells","pDCs"        ,            
           "Pre pro B cells","Pro B cells")

cell_type = clusters_for_DE[11]
for (cell_type in clusters_for_DE) { 
  cell_type_filename = gsub("/","_",cell_type)
  
  # cell_type="Early erythroid cells"
  # cell_type="Pro B cells"
  res.df.mg <- fread(paste0(fdir,cell_type_filename,'.txt'),data.table = F,stringsAsFactors = F)
  print(cell_type)
  df1 <- res.df.all.lfc.liver[,c("names",cell_type_filename)]; df1$P.liver <- res.df.all.p.liver[,c(cell_type_filename)]
  df2 <- res.df.all.lfc.femur[,c("names","chromosome_name",cell_type_filename)]; df2$P.femur <- res.df.all.p.femur[,c(cell_type_filename)]
  colnames(df1)[2] <- "Liver"
  colnames(df2)[3] <- "Femur"
  
  df.mg <- merge(df1,df2,by='names')
  
  ###########
  
  # res.df.mg <- df
  # adj.P.Val.x
  res.df.mg2 <- merge(df.mg,res.df.mg,by='names')
  
  ##########
  
  fdr_thres = 0.5
  # tmp <- subset(res.df.mg2,(!is.na(Liver) | !is.na(Femur)))
  tmp <- subset(res.df.mg2,(!is.na(Liver) & !is.na(Femur)))
  tmp <- subset(tmp,adj.P.Val.x < 0.05)
  tmp <- tmp[order(tmp$adj.P.Val.x),]
  tmp$P.liver[is.na(tmp$Liver)] <- 1
  tmp$P.femur[is.na(tmp$Femur)] <- 1
  tmp$Liver[is.na(tmp$Liver)] <- 0
  tmp$Femur[is.na(tmp$Femur)] <- 0
  tmp$int_type <- 'unknown'
  tmp$int_type[sign(tmp$Liver) != sign(tmp$Femur) & tmp$P.liver < fdr_thres & tmp$P.femur < fdr_thres] <- 'switch'
  tmp$int_type[tmp$int_type!='switch' & sign(tmp$logFC.x) == sign(tmp$logFC.y) & tmp$adj.P.Val.y < fdr_thres] <- 'amplified'
  tmp$int_type[tmp$int_type!='switch' & sign(tmp$logFC.x) != sign(tmp$logFC.y) & tmp$adj.P.Val.y < fdr_thres] <- 'reduced'
  tmp$int_type[tmp$int_type!='switch' & sign(tmp$logFC.x) == 1 & tmp$adj.P.Val.y >= fdr_thres] <- 'amplified'
  tmp$int_type[tmp$int_type!='switch' & sign(tmp$logFC.x) == -1 & tmp$adj.P.Val.y >= fdr_thres] <- 'reduced'
  
  print(cell_type)
  print(table(tmp$int_type))
  print(tmp[1,])
}
# # tmp$int_type[tmp$int_type!='switch' & abs(tmp$Liver) > abs(tmp$Femur)] <- 'amplified'
# tmp$int_type[tmp$int_type!='switch' & sign(tmp$logFC.x) == sign(tmp$logFC.y) & abs(tmp$Liver) > abs(tmp$Femur)] <- 'amplified'
# # tmp$int_type[tmp$int_type!='switch' & abs(tmp$Liver) < abs(tmp$Femur)] <- 'reduced'
# tmp$int_type[tmp$int_type!='switch' & abs(tmp$Liver) < abs(tmp$Femur) & sign(tmp$logFC.x) != sign(tmp$logFC.y)] <- 'reduced'
# tmp$int_type[tmp$int_type!='switch' & abs(tmp$Liver) < abs(tmp$Femur) & sign(tmp$logFC.x) != sign(tmp$logFC.y)] <- 'reduced'
# tmp$int_type[tmp$int_type=='unknown' & sign(tmp$logFC.x) != sign(tmp$logFC.y)] <- 'reduced'
# tmp$int_type[tmp$int_type=='unknown' & sign(tmp$logFC.x) == sign(tmp$logFC.y)] <- 'amplified'

table(tmp$int_type)

score <- sum(tmp$int_type=="amplified")/(sum(tmp$int_type=="amplified") + sum(tmp$int_type=="reduced"))
sample(1:nrow(tmp),size = )



#############


res.df.mg2.sub <- subset(res.df.mg2,sign(Femur)==sign(Liver))
sum(res.df.mg2.sub$P.Value.x < 0.05)
sum(res.df.mg2.sub$adj.P.Val.x < 0.05)

tmp <- subset(res.df.mg2.sub,adj.P.Val.x < 0.05)
tmp <- subset(res.df.mg2.sub,P.Value.x < 0.05)
mean(sign(tmp$logFC.x) == sign(tmp$Liver))

tmp <- subset(res.df.mg2,adj.P.Val.x < 0.05)
tmp$main_effect <- ifelse(!is.na(tmp$Liver),sign(tmp$Liver),sign(tmp$Femur))
table(sign(tmp$logFC.x) == tmp$main_effect)
ma

# tmp <- subset(res.df.mg2,adj.P.Val.x < 0.05)
fdr_thres = 0.25
tmp <- subset(res.df.mg2,(!is.na(Liver) | !is.na(Femur)))
tmp <- subset(tmp,adj.P.Val.x < 0.05)
tmp$P.liver[is.na(tmp$Liver)] <- 1
tmp$P.femur[is.na(tmp$Femur)] <- 1
tmp$Liver[is.na(tmp$Liver)] <- 0
tmp$Femur[is.na(tmp$Femur)] <- 0
tmp$int_type <- 'unknown'
tmp$int_type[sign(tmp$Liver) != sign(tmp$Femur) & tmp$P.liver < fdr_thres & tmp$P.femur < fdr_thres] <- 'switch'
tmp$int_type[tmp$int_type!='switch' & abs(tmp$Liver) > abs(tmp$Femur)] <- 'amplified'
tmp$int_type[tmp$int_type!='switch' & abs(tmp$Liver) < abs(tmp$Femur)] <- 'reduced'
table(tmp$int_type)

sum(abs(tmp$Liver) > sum(tmp$Femur))


  sign(tmp$logFC.x) == sign(tmp$Liver)






