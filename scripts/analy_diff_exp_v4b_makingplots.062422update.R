library(data.table)
library(ggplot2)
library(cowplot)
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

dir <- '/oak/stanford/groups/smontgom/amarder'
# dir <- '~/Documents/Research'

sampletype <- "Liver"
sampletype <- "Femur"
datadir="DE_pb_leiden_names";

#########33

for (sampletype in c("Liver","Femur")) { 
  res.df.all <- list()
  
  print(iter)
  
  dir <- '/oak/stanford/groups/smontgom/amarder'
  pathtodir <- paste0(dir,"/t21-proj/out/full/",datadir)
  flist <- list.files(pathtodir)
  flist <- flist[grep(sampletype,flist)]
  cell_type_groups <- flist[!grepl("res",flist)]
  cell_type_groups <- substring(cell_type_groups,nchar(sampletype)+2,nchar(cell_type_groups)-4)
  cell_type_groups <- cell_type_groups[!grepl("sample",cell_type_groups)]
  
  # for (i in 1:21) {
  for (i in 1:length(cell_type_groups)) {
    print(i)
    cell_type <- cell_type_groups[i]
    f<- paste0(dir,"/t21-proj/out/full/",datadir,"/",sampletype,".",cell_type,".sample.txt")
    df <- fread(f,data.table = F,stringsAsFactors = F)
    tmp.lfc <- df[,c("names","logFC")]
    tmp.p <- df[,c("names","P.Value")]
    tmp.fdr <- df[,c("names","adj.P.Val")]
    colnames(tmp.lfc)[2] <- cell_type
    colnames(tmp.p)[2] <- cell_type
    colnames(tmp.fdr)[2] <- cell_type
    if (i==1) {
      res.lfc <- tmp.lfc
      res.p <- tmp.p
      res.fdr <- tmp.fdr
    } else {
      res.lfc <- merge(res.lfc,tmp.lfc,by='names',all = TRUE)
      res.p <- merge(res.p,tmp.p,by='names',all = TRUE)
      res.fdr <- merge(res.fdr,tmp.fdr,by='names',all = TRUE)
    }
  }
  annot <- getBM(attributes = c('hgnc_symbol', 'chromosome_name',
                                'start_position', 'end_position'),
                 filters = 'hgnc_symbol', 
                 values = res.lfc$names, 
                 mart = ensembl)
  fileOut <- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/data_small/gene_annot.txt")
  geneAnnot <- fread(fileOut,data.table = F,stringsAsFactors = F)
  geneAnnot <- geneAnnot[!duplicated(geneAnnot$hgnc_symbol),]
  df1.lfc <- merge(res.lfc,geneAnnot[,c("hgnc_symbol","chromosome_name","start_position")],by.x="names",by.y="hgnc_symbol")
  df1.p <- merge(res.p,geneAnnot[,c("hgnc_symbol","chromosome_name","start_position")],by.x="names",by.y="hgnc_symbol")
  df1.fdr <- merge(res.fdr,geneAnnot[,c("hgnc_symbol","chromosome_name","start_position")],by.x="names",by.y="hgnc_symbol")
  f.out <- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".FE.lfc.txt")
  fwrite(df1.lfc,f.out,quote = F,na = "NA",sep = "\t",row.names = F,col.names = T)
  f.out <- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".FE.p.txt")
  fwrite(df1.p,f.out,quote = F,na = "NA",sep = "\t",row.names = F,col.names = T)
  f.out <- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".FE.fdr.txt")
  fwrite(df1.fdr,f.out,quote = F,na = "NA",sep = "\t",row.names = F,col.names = T)
  
}

#############

# saveRDS(res.df.all.p,file=f.out)



##############################

median_lfc_across_cell_types <- function(sampletype) {
  res.df.all.p <- fread(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".FE.fdr.txt"),data.table = F,stringsAsFactors=F)
  res.df.all.lfc <- fread(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".FE.lfc.txt"),data.table = F,stringsAsFactors = F)
  cell_type_groups <- colnames(res.df.all.p)[(which(colnames(res.df.all.p)=="names")+1):(which(colnames(res.df.all.p)=="chromosome_name")-1)]
  
  x <- aggregate(res.df.all.lfc[,cell_type_groups],by=list(chr21=res.df.all.lfc$chromosome_name==21),median,na.rm=T)
  rownames(x) <- x[,1]; x <- x[,-1]; x <- as.data.frame(t(x))
  colnames(x) <- c("Not Chr 21","Chr 21")
  x$cell <- rownames(x)
  tab <- reshape2::melt(x,id.vars="cell")
  
  g <- ggplot(tab,aes(x=cell,fill=variable,y=value)) + 
    geom_bar(stat = 'identity',col='black',position = position_dodge2(preserve = 'single')) + 
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle=55,hjust = 1)
    ) +
    labs(x="Cell type",y="Median LFC (T21 vs Healthy)",fill='Chromosome') 
  print(g)
}

median_lfc_across_cell_types("Liver")
median_lfc_across_cell_types("Femur")

####################

# lfc_for_gene_across_cell_types <- function(gene_to_use,sampletype) { 
#   
#   res.df.all.p <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".p.rds"))
#   res.df.all.lfc <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".lfc.rds"))
#   cell_type_groups <- colnames(res.df.all.p[[1]])[(which(colnames(res.df.all.p[[1]])=="names")+1):(which(colnames(res.df.all.p[[1]])=="chromosome_name")-1)]
#   
#   # look at only same genes for both methods
#   for (cell_type in cell_type_groups) {
#     keep_genes <- res.df.all.p[[2]][!is.na(res.df.all.p[[2]][,cell_type]),'names']
#     res.df.all.p[[1]][!(res.df.all.p[[1]][,'names']%in%keep_genes),cell_type] <- NA
#     res.df.all.lfc[[1]][!(res.df.all.lfc[[1]][,'names']%in%keep_genes),cell_type] <- NA
#   }
#   
#   iter=1
#   tab <- subset(res.df.all.lfc[[iter]],names==gene_to_use)
#   tab <- reshape2::melt(tab[,cell_type_groups])
#   iter=2
#   tab2 <- subset(res.df.all.lfc[[iter]],names==gene_to_use)
#   tab2 <- reshape2::melt(tab2[,cell_type_groups])
#   tab$sample <- "Fetus"
#   tab2$sample <- "Sample"
#   tab <- rbind(tab,tab2)
#   
#   ggplot(tab,aes(x=variable,fill=sample,y=value)) + 
#     geom_bar(stat = 'identity',col='black',position = position_dodge2(preserve = 'single')) + 
#     theme_bw() +
#     theme(
#       panel.grid = element_blank(),
#       axis.text.x = element_text(angle=55,hjust = 1),
#       plot.title = element_text(hjust=0.5)
#     ) +
#     labs(x="Cell type",y="LFC (T21 vs Healthy)",fill='Pseudobulk',title="SOD1")
# }
# lfc_for_gene_across_cell_types("SOD1","Femur")
# lfc_for_gene_across_cell_types("SOD1","Liver")
# 
# #####################
# 
# library(cowplot)
# library(ggplot2)
# DEG_comparison <- function(sampletype) {
#   
#   res.df.all.p <- fread(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".FE.fdr.txt"),data.table = F,stringsAsFactors=F)
#   res.df.all.lfc <- fread(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".FE.lfc.txt"),data.table = F,stringsAsFactors = F)
#   cell_type_groups <- colnames(res.df.all.p)[(which(colnames(res.df.all.p)=="names")+1):(which(colnames(res.df.all.p)=="chromosome_name")-1)]
#   
#   x <- res.df.all.p
#   x <- aggregate(res.df.all.p[,cell_type_groups]<0.1,by=list(chr21=res.df.all.lfc$chromosome_name==21),mean,na.rm=T)
#   rownames(x) <- x[,1]; x <- x[,-1]; x <- as.data.frame(t(x))
# 
#   g1 <- ggplot(tab,aes(x=100*`Chr 21 1`,y=100*`Chr 21 2`)) + 
#     geom_point() + 
#     theme_bw() + 
#     theme(panel.grid = element_blank(),
#           plot.title=element_text(hjust=0.5)) +
#     geom_smooth(method='loess',span=3,se=F,col='purple') + 
#     geom_abline(slope=1,intercept = 0,col='red',lty='dashed') +
#     labs(x="% DEG per cell type (fetus)",y="% DEG per cell type (sample)",title="Chr 21")
#   g2 <- ggplot(tab,aes(x=100*`Not Chr 21 1`,y=100*`Not Chr 21 2`)) + 
#     geom_point() + 
#     theme_bw() + 
#     theme(panel.grid = element_blank(),
#           plot.title=element_text(hjust=0.5)) +
#     geom_smooth(method='loess',span=3,se=F,col='purple') + 
#     geom_abline(slope=1,intercept = 0,col='red',lty='dashed') +
#     labs(x="% DEG per cell type (fetus)",y="% DEG per cell type (sample)",title="Not Chr 21")
#   
#   plot_grid(g1,g2,ncol=2)
#   
# }
# 
# DEG_comparison("Liver")
# DEG_comparison("Femur")
# sampletype="Femur"


DEG_comparison_chr21_v_notchr21 <- function(sampletype) { 
  res.df.all.p <- fread(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".FE.fdr.txt"),data.table = F,stringsAsFactors=F)
  res.df.all.lfc <- fread(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".FE.lfc.txt"),data.table = F,stringsAsFactors = F)
  cell_type_groups <- colnames(res.df.all.p)[(which(colnames(res.df.all.p)=="names")+1):(which(colnames(res.df.all.p)=="chromosome_name")-1)]
  
  x <- res.df.all.p
  # x <- aggregate(res.df.all.p[,cell_type_groups]<0.1,by=list(chr21=res.df.all.lfc$chromosome_name==21),mean,na.rm=T)
  x <- aggregate(res.df.all.p[,cell_type_groups]<0.05,by=list(chr21=res.df.all.lfc$chromosome_name==21),mean,na.rm=T)
  rownames(x) <- x[,1]; x <- x[,-1]; x <- as.data.frame(t(x))
  colnames(x) <- c("Not Chr 21","Chr 21")
  
  ggplot(x,aes(x=100*`Chr 21`,y=100*`Not Chr 21`)) + 
    geom_point() + 
    theme_bw() + 
    theme(panel.grid = element_blank(),
          plot.title=element_text(hjust=0.5)) +
    geom_smooth(method='loess',span=3,se=F,col='purple') + 
    labs(x="% DEG per cell type (chr 21)",y="% DEG per cell type (not chr 21)",title=paste0(sampletype))
}

DEG_comparison_chr21_v_notchr21("Liver")
DEG_comparison_chr21_v_notchr21("Femur")

DEG_comparison_chrX_v_notchrX <- function(sampletype,chrNum) { 
  res.df.all.p <- fread(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".FE.fdr.txt"),data.table = F,stringsAsFactors=F)
  # res.df.all.lfc <- fread(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".FE.lfc.txt"),data.table = F,stringsAsFactors = F)
  cell_type_groups <- colnames(res.df.all.p)[(which(colnames(res.df.all.p)=="names")+1):(which(colnames(res.df.all.p)=="chromosome_name")-1)]
  
  # if (chrNum!=21) {res.df.all.p <- subset(res.df.all.p,chromosome_name != 21)}
  # x <- aggregate(res.df.all.p[,cell_type_groups]<0.1,by=list(chr21=res.df.all.lfc$chromosome_name==21),mean,na.rm=T)
  x <- aggregate(res.df.all.p[,cell_type_groups]<0.05,by=list(chr21=res.df.all.p$chromosome_name==chrNum),mean,na.rm=T)
  rownames(x) <- x[,1]; x <- x[,-1]; x <- as.data.frame(t(x))
  colnames(x) <- c("Not Chr 21","Chr 21")
  
  return(ggplot(x,aes(x=100*`Chr 21`,y=100*`Not Chr 21`)) + 
           geom_point() + 
           theme_bw() + 
           theme(panel.grid = element_blank(),
                 plot.title=element_text(hjust=0.5)) +
           geom_smooth(method='loess',span=3,se=F,col='purple') + 
           # labs(x=paste0("% DEG per cell type (chr ",chrNum,")"),y=paste0("% DEG per cell type (not chr ",chrNum,")"),title=paste0(sampletype," chr",chrNum)))
           labs(x=paste0("chr",chrNum),y=paste0("non-chr",chrNum,""),title=paste0(sampletype," chr",chrNum)))
}

glst <- list()
for (i in 1:22) {
  glst[[i]] <- DEG_comparison_chrX_v_notchrX("Liver",i)
}
plot_grid(plotlist=glst,ncol=4)
# DEG_comparison_chr1_v_notchr1("Femur",1)


##########################

draw_logfc_plots <- function(sampletype,iter) { 
  
  res.df.all.p <- fread(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".FE.fdr.txt"),data.table = F,stringsAsFactors=F)
  res.df.all.lfc <- fread(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".FE.lfc.txt"),data.table = F,stringsAsFactors = F)
  cell_type_groups <- colnames(res.df.all.p)[(which(colnames(res.df.all.p)=="names")+1):(which(colnames(res.df.all.p)=="chromosome_name")-1)]

  df1.lfc.perc <- res.df.all.lfc
  df1.lfc.perc[,cell_type_groups] <- apply(apply(df1.lfc.perc[,cell_type_groups],2,rank,na.last="keep"),2,function(x) x/sum(!is.na(x)))
  df1.lfc.perc$chr21 <- ifelse(df1.lfc.perc$chromosome_name==21,"Chr 21","Not Chr 21")
  # df1.lfc.perc[,cell_type_groups] <- lapply(res.df.all.lfc,function(x) {(apply(apply(x[,cell_type_groups],2,rank,na.last="keep"),2,function(x) x/sum(!is.na(x))))})[[iter]]
  
  # tab <- lapply(res.df.all.lfc,function(x) {aggregate(apply(apply(x[,cell_type_groups],2,rank),2,function(x) x/sum(!is.na(x))),by=list(x$chr21),median,na.rm=T)})
  # tab <- lapply(tab,function(x) {rownames(x) <- x[,1]; x <- x[,-1]; x <- as.data.frame(t(x))})
  # tab <- do.call(cbind,tab)
  
  
  # df1.lfc.perc[,cell_type_groups] <- apply(df1.lfc.perc[,cell_type_groups],2,rank)/nrow(df1.lfc.perc)
  df1.lfc.perc.melt <- reshape2::melt(df1.lfc.perc[,c(cell_type_groups,"chr21")],id.vars=c("chr21"))
  
  plot_order <- list()
  plot_order[["Femur"]] <- c("HSCs_MPPs","MEMPs","Granulocyte progenitors",
                             "Early erythroid cells","Late erythroid cells",
                             "Mast cells",
                             "Megakaryocytes",
                             "Neutrophils","Cycling neutrophils","Inflammatory macrophages","Tolerogenic macrophages","Collagen+ myeloid cells","Osteoclasts","pDCs","cDC2",
                             "NK cells",
                             "Pre pro B cells","Pro B cells","B cells",
                             "CAR cells","Fibroblasts","LEPR+ CAR cells","Osteoprogenitors","Pericytes","Schwann cells","Sinusoidal endothelial cells","Transitioning endothelial cells","Vascular endothelial cells")
  plot_order[["Liver"]] <- c("HSCs_MPPs","MEMPs","Granulocyte progenitors",
                             "Early erythroid cells","Late erythroid cells",
                             "Mast cells",
                             "Megakaryocytes",
                             "Monocyte progenitors","Inflammatory macrophages","Kupffer cells","pDCs","cDC2",
                             "NK cells","NK progenitors",
                             "Pre pro B cells","Pro B cells",
                             "Activated stellate cells","Hepatic stellate cells","Hepatocytes","LSECs","Vascular endothelial cells")
  
  df1.lfc.perc.melt$variable <- factor(df1.lfc.perc.melt$variable,plot_order[[sampletype]])
  
  ggplot(df1.lfc.perc.melt,aes(x=variable,fill=chr21,y=value)) + 
    geom_boxplot() + 
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle=55,hjust = 1)
    ) +
    labs(x="Cell type",y="LFC (Percentile) (T21 vs Healthy)",fill='') 
}

draw_logfc_plots("Liver",2)
draw_logfc_plots("Femur",2)

sampletype="Liver"
res.df.all.lfc <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".lfc.rds"))
res.df.all.p <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".p.rds"))
cell_type_groups <- colnames(res.df.all.lfc[[1]])[(which(colnames(res.df.all.lfc[[1]])=="names")+1):(which(colnames(res.df.all.lfc[[1]])=="chromosome_name")-1)]
res.df.p <- res.df.all.p[[2]]
res.df.p.med <- data.frame(names=res.df.p$names,chromosome_name=res.df.p$chromosome_name,med_padj=apply(res.df.p[,cell_type_groups],1,median,na.rm=T),m=apply(res.df.p[,cell_type_groups],1,function(x) {sum(!is.na(x))}))
res.df.p.med <- subset(res.df.p.med,m>=10)
aggregate(res.df.p.med$med_padj,by=list(res.df.p.med$chromosome_name==21),median)
tmp <- res.df.p.med[order(res.df.p.med$med_padj,decreasing = F),]
gene_lst <- subset(tmp,med_padj < 0.1 & chromosome_name!=21)$names
# gene_lst <- subset(tmp,med_padj < 0.1)$names # & chromosome_name==21)$names
library("enrichR")
dbs <- c("GO_Molecular_Function_2018", "GO_Biological_Process_2018","ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")
enriched <- enrichr(gene_lst, dbs)
y <- enriched[[dbs[1]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); head(subset(y,geneCt>3),3)
y <- enriched[[dbs[2]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); head(subset(y,geneCt>3 & Adjusted.P.value<0.1),25)$Term
y <- enriched[[dbs[3]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); head(subset(y,geneCt>3 & Adjusted.P.value<0.1),20)$Term

enriched[[dbs[2]]][1:10,]
enriched[[dbs[3]]][1:3,]

subset(res.df.all.p[[2]],names=="MAX")
subset(res.df.all.p[[2]],names=="TAF7")
subset(res.df.all.p[[2]],names=="MYC")

enriched <- enrichr(gene_lst, dbs)


subset(tmp,chromosome_name!=21)[1:5,]
subset(res.df.all.p[[2]],names=="IGHM")
subset(res.df.all.lfc[[2]],names=="IGHM")

subset(res.df.lfc,names=="IL13")

aggregate(res.df.lfc$`Late erythroid cells`,by=list(res.df.lfc$chr21),median,na.rm=T)
aggregate(rank(res.df.lfc$`Late erythroid cells`,na.last = "keep"),by=list(res.df.lfc$chr21),median,na.rm=T)


###############

draw_heatmap <- function(sampletype,iter) { 
  res.df.all.p <- fread(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".FE.fdr.txt"),data.table = F,stringsAsFactors=F)
  res.df.all.lfc <- fread(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".FE.lfc.txt"),data.table = F,stringsAsFactors = F)
  cell_type_groups <- colnames(res.df.all.p)[(which(colnames(res.df.all.p)=="names")+1):(which(colnames(res.df.all.p)=="chromosome_name")-1)]
  
  df1.p <- res.df.all.p
  df1.lfc <- res.df.all.lfc
  df1.p <- subset(df1.p,chromosome_name %in% seq(1,22))
  
  aggregate(df1.p$HSCs_MPPs,by=list(df1.p$chromosome_name==21),function(x) {mean(x < 0.05,na.rm=T)})
  
  # df1.p.melt <- reshape2::melt(aggregate(df1.p[,cell_type_groups]<0.1,by=list(df1.p$chromosome_name),mean,na.rm=T),id.vars="Group.1")
  df1.p.melt <- reshape2::melt(aggregate(df1.p[,cell_type_groups]<0.05,by=list(df1.p$chromosome_name),mean,na.rm=T),id.vars="Group.1")
  # df1.p.melt <- reshape2::melt(aggregate(df1.p[,cell_type_groups]<0.1 & df1.lfc[,cell_type_groups] > 0,by=list(df1.p$chromosome_name),mean,na.rm=T),id.vars="Group.1")
  
  plot_order <- list()
  plot_order[["Femur"]] <- c("HSCs_MPPs","MEMPs","Granulocyte progenitors",
                             "Early erythroid cells","Late erythroid cells",
                             "Mast cells",
                             "Megakaryocytes",
                             "Neutrophils","Cycling neutrophils","Inflammatory macrophages","Tolerogenic macrophages","Collagen+ myeloid cells","Osteoclasts","pDCs","cDC2",
                             "NK cells",
                             "Pre pro B cells","Pro B cells","B cells",
                             "CAR cells","Fibroblasts","LEPR+ CAR cells","Osteoprogenitors","Pericytes","Schwann cells","Sinusoidal endothelial cells","Transitioning endothelial cells","Vascular endothelial cells")
  plot_order[["Liver"]] <- c("HSCs_MPPs","MEMPs","Granulocyte progenitors",
                             "Early erythroid cells","Late erythroid cells",
                             "Mast cells",
                             "Megakaryocytes",
                             "Monocyte progenitors","Inflammatory macrophages","Kupffer cells","pDCs","cDC2",
                             "NK cells","NK progenitors",
                             "Pre pro B cells","Pro B cells",
                             "Activated stellate cells","Hepatic stellate cells","Hepatocytes","LSECs","Vascular endothelial cells")
  
  df1.p.melt$variable <- factor(df1.p.melt$variable,plot_order[[sampletype]])
  df1.p.melt$Group.1 <- factor(df1.p.melt$Group.1,levels = seq(1,22))
  g <- ggplot(df1.p.melt,aes(x=Group.1,y=variable,fill=100*value)) + 
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
    theme_minimal()+
    theme(
      plot.title=element_text(hjust=0.5)) + 
    # guides(fill=F) +
    coord_fixed() +
    labs(title='% DEG',x='Chromosome',y='Cell type',fill="% sig")
  print(g)
}

# draw_heatmap("Femur",1)
draw_heatmap("Femur",2)
draw_heatmap("Liver",1)
# draw_heatmap("Liver",2)

sampletype="Femur"
res.df.all.p <- fread(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".FE.fdr.txt"),data.table = F,stringsAsFactors=F)
df1.p <- subset(res.df.all.p,chromosome_name %in% seq(1,22))
aggregate(df1.p$HSCs_MPPs,by=list(df1.p$chromosome_name==21),function(x) {mean(x < 0.05,na.rm=T)})
sampletype="Liver"
res.df.all.p <- fread(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".FE.fdr.txt"),data.table = F,stringsAsFactors=F)
df1.p <- subset(res.df.all.p,chromosome_name %in% seq(1,22))
aggregate(df1.p$HSCs_MPPs,by=list(df1.p$chromosome_name==21),function(x) {mean(x < 0.05,na.rm=T)})

###################3

sampletype="Femur"; iter = 2
df1.p <- fread(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".FE.fdr.txt"),data.table = F,stringsAsFactors=F)
df1.p.melt <- reshape2::melt(aggregate(df1.p[,c("Pre pro B cells","Pro B cells","B cells")]<0.05,by=list(df1.p$chromosome_name==21),mean,na.rm=T),id.vars="Group.1")
ggplot(df1.p.melt,aes(x=variable,y=100*value,fill=Group.1)) + 
  geom_bar(stat='identity',position = 'dodge',col='black') +
  scale_fill_brewer(palette = 'Set2',labels=c("Not Chr 21","Chr 21")) + 
  theme_bw() + 
  theme(panel.grid = element_blank()) +
  labs(y="% significant DEG",x="Cell type",fill='') 


#

sampletype="Liver"
rng <- c(1e-1,1e-2,1e-3,1e-4,1e-5)
dftmp = list()
iter=0
for (sampletype in c("Liver","Femur")) {
  res.df.all.p <- fread(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".FE.p.txt"),data.table = F,stringsAsFactors=F)
  for (i in 1:length(rng)) {
    iter = iter + 1
    thres = rng[i]
    dftmp[[iter]] <- data.frame(sampletype,thres,sig=mean(res.df.all.p$HSCs_MPPs < thres,na.rm=T))
  }
}
dftmp <- do.call(rbind,dftmp)
ggplot(dftmp,aes(x=-log10(thres),y=sig,col=sampletype)) + geom_line() + geom_point() + theme_bw() + theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5)) +
  labs(x="-log10 P-value threshold",y="Proportion of genes with nominal P-value below threshold",col="Environment",title="T21 vs Healthy") + scale_color_brewer(palette = "Set1")

sampletype="Liver"
rng <- c(0.25,0.2,0.15,0.1,0.05,0.01)
dftmp = list()
iter=0
for (sampletype in c("Liver","Femur")) {
  res.df.all.p <- fread(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".FE.fdr.txt"),data.table = F,stringsAsFactors=F)
  for (i in 1:length(rng)) {
    iter = iter + 1
    thres = rng[i]
    dftmp[[iter]] <- data.frame(sampletype,thres,sig=mean(res.df.all.p$HSCs_MPPs < thres,na.rm=T))
  }
}
dftmp <- do.call(rbind,dftmp)
ggplot(dftmp,aes(x=-log10(thres),y=sig,col=sampletype)) + geom_line() + geom_point() + theme_bw() + theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5)) +
  labs(x="-log10 FDR threshold",y="Proportion of genes with FDR below threshold",col="Environment",title="T21 vs Healthy") + scale_color_brewer(palette = "Set1")

