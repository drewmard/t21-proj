median_lfc_across_cell_types <- function(sampletype) {
  res.df.all.p <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".p.rds"))
  res.df.all.lfc <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".lfc.rds"))
  cell_type_groups <- colnames(res.df.all.p[[1]])[(which(colnames(res.df.all.p[[1]])=="names")+1):(which(colnames(res.df.all.p[[1]])=="chromosome_name")-1)]
  
  # look at only same genes for both methods
  for (cell_type in cell_type_groups) {
    keep_genes <- res.df.all.p[[2]][!is.na(res.df.all.p[[2]][,cell_type]),'names']
    res.df.all.p[[1]][!(res.df.all.p[[1]][,'names']%in%keep_genes),cell_type] <- NA
    res.df.all.lfc[[1]][!(res.df.all.lfc[[1]][,'names']%in%keep_genes),cell_type] <- NA
  }
  
  tab <- lapply(res.df.all.lfc,function(x) {aggregate(x[,cell_type_groups],by=list(x$chr21),median,na.rm=T)})
  tab <- lapply(tab,function(x) {rownames(x) <- x[,1]; x <- x[,-1]; x <- as.data.frame(t(x))})
  tab <- do.call(cbind,tab)
  tab <- tab[,c(2,4)]
  colnames(tab) <- c("Fetus","Sample")
  tab$cell <- rownames(tab)
  tab <- reshape2::melt(tab,id.vars="cell")
  
  g <- ggplot(tab,aes(x=cell,fill=variable,y=value)) + 
    geom_bar(stat = 'identity',col='black',position = position_dodge2(preserve = 'single')) + 
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle=55,hjust = 1)
    ) +
    labs(x="Cell type",y="Median LFC (T21 vs Healthy)",fill='Pseudobulk') 
  print(g)
}

median_lfc_across_cell_types("Liver")
median_lfc_across_cell_types("Femur")

####################

lfc_for_gene_across_cell_types <- function(gene_to_use,sampletype) { 
  
  res.df.all.p <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".p.rds"))
  res.df.all.lfc <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".lfc.rds"))
  cell_type_groups <- colnames(res.df.all.p[[1]])[(which(colnames(res.df.all.p[[1]])=="names")+1):(which(colnames(res.df.all.p[[1]])=="chromosome_name")-1)]
  
  # look at only same genes for both methods
  for (cell_type in cell_type_groups) {
    keep_genes <- res.df.all.p[[2]][!is.na(res.df.all.p[[2]][,cell_type]),'names']
    res.df.all.p[[1]][!(res.df.all.p[[1]][,'names']%in%keep_genes),cell_type] <- NA
    res.df.all.lfc[[1]][!(res.df.all.lfc[[1]][,'names']%in%keep_genes),cell_type] <- NA
  }
  
  iter=1
  tab <- subset(res.df.all.lfc[[iter]],names==gene_to_use)
  tab <- reshape2::melt(tab[,cell_type_groups])
  iter=2
  tab2 <- subset(res.df.all.lfc[[iter]],names==gene_to_use)
  tab2 <- reshape2::melt(tab2[,cell_type_groups])
  tab$sample <- "Fetus"
  tab2$sample <- "Sample"
  tab <- rbind(tab,tab2)
  
  ggplot(tab,aes(x=variable,fill=sample,y=value)) + 
    geom_bar(stat = 'identity',col='black',position = position_dodge2(preserve = 'single')) + 
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle=55,hjust = 1),
      plot.title = element_text(hjust=0.5)
    ) +
    labs(x="Cell type",y="LFC (T21 vs Healthy)",fill='Pseudobulk',title="SOD1")
}
lfc_for_gene_across_cell_types("SOD1","Femur")
lfc_for_gene_across_cell_types("SOD1","Liver")

#####################

library(cowplot)
library(ggplot2)
DEG_comparison <- function(sampletype) {
  
  res.df.all.p <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".p.rds"))
  res.df.all.lfc <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".lfc.rds"))
  cell_type_groups <- colnames(res.df.all.p[[1]])[(which(colnames(res.df.all.p[[1]])=="names")+1):(which(colnames(res.df.all.p[[1]])=="chromosome_name")-1)]
  
  for (cell_type in cell_type_groups) {
    keep_genes <- res.df.all.p[[2]][!is.na(res.df.all.p[[2]][,cell_type]),'names']
    res.df.all.p[[1]][!(res.df.all.p[[1]][,'names']%in%keep_genes),cell_type] <- NA
    res.df.all.lfc[[1]][!(res.df.all.lfc[[1]][,'names']%in%keep_genes),cell_type] <- NA
  }
  
  tab <- lapply(res.df.all.p,function(x) {aggregate(x[,cell_type_groups]<0.1,by=list(x$chr21),mean,na.rm=T)})
  tab <- lapply(tab,function(x) {rownames(x) <- x[,1]; x <- x[,-1]; x <- as.data.frame(t(x))})
  tab <- do.call(cbind,tab)
  
  colnames(tab)[1:2] <- paste0(colnames(tab)[1:2]," 1")
  colnames(tab)[3:4] <- paste0(colnames(tab)[3:4]," 2")
  colnames(tab)[5:6] <- paste0(colnames(tab)[5:6]," 3")
  
  g1 <- ggplot(tab,aes(x=100*`Chr 21 1`,y=100*`Chr 21 2`)) + 
    geom_point() + 
    theme_bw() + 
    theme(panel.grid = element_blank(),
          plot.title=element_text(hjust=0.5)) +
    geom_smooth(method='loess',span=3,se=F,col='purple') + 
    geom_abline(slope=1,intercept = 0,col='red',lty='dashed') +
    labs(x="% DEG per cell type (fetus)",y="% DEG per cell type (sample)",title="Chr 21")
  g2 <- ggplot(tab,aes(x=100*`Not Chr 21 1`,y=100*`Not Chr 21 2`)) + 
    geom_point() + 
    theme_bw() + 
    theme(panel.grid = element_blank(),
          plot.title=element_text(hjust=0.5)) +
    geom_smooth(method='loess',span=3,se=F,col='purple') + 
    geom_abline(slope=1,intercept = 0,col='red',lty='dashed') +
    labs(x="% DEG per cell type (fetus)",y="% DEG per cell type (sample)",title="Not Chr 21")
  
  plot_grid(g1,g2,ncol=2)
  
}

DEG_comparison("Liver")
DEG_comparison("Femur")
sampletype="Femur"


DEG_comparison_chr21_v_notchr21 <- function(sampletype) { 
  res.df.all.p <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".p.rds"))
  res.df.all.lfc <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".lfc.rds"))
  cell_type_groups <- colnames(res.df.all.p[[1]])[(which(colnames(res.df.all.p[[1]])=="names")+1):(which(colnames(res.df.all.p[[1]])=="chromosome_name")-1)]
  
  for (cell_type in cell_type_groups) {
    keep_genes <- res.df.all.p[[2]][!is.na(res.df.all.p[[2]][,cell_type]),'names']
    res.df.all.p[[1]][!(res.df.all.p[[1]][,'names']%in%keep_genes),cell_type] <- NA
    res.df.all.lfc[[1]][!(res.df.all.lfc[[1]][,'names']%in%keep_genes),cell_type] <- NA
  }
  
  tab <- lapply(res.df.all.p,function(x) {aggregate(x[,cell_type_groups]<0.1,by=list(x$chr21),mean,na.rm=T)})
  tab <- lapply(tab,function(x) {rownames(x) <- x[,1]; x <- x[,-1]; x <- as.data.frame(t(x))})
  tab <- do.call(cbind,tab)
  
  colnames(tab)[1:2] <- paste0(colnames(tab)[1:2]," 1")
  colnames(tab)[3:4] <- paste0(colnames(tab)[3:4]," 2")
  colnames(tab)[5:6] <- paste0(colnames(tab)[5:6]," 3")
  
  ggplot(tab,aes(x=100*`Chr 21 2`,y=100*`Not Chr 21 2`)) + 
    geom_point() + 
    theme_bw() + 
    theme(panel.grid = element_blank(),
          plot.title=element_text(hjust=0.5)) +
    geom_smooth(method='loess',span=3,se=F,col='purple') + 
    labs(x="% DEG per cell type (chr 21)",y="% DEG per cell type (not chr 21)",title=paste0(sampletype))
}

DEG_comparison_chr21_v_notchr21("Liver")
DEG_comparison_chr21_v_notchr21("Femur")

##########################

draw_logfc_plots <- function(sampletype,iter) { 
  
  res.df.all.p <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".p.rds"))
  res.df.all.lfc <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".lfc.rds"))
  cell_type_groups <- colnames(res.df.all.p[[1]])[(which(colnames(res.df.all.p[[1]])=="names")+1):(which(colnames(res.df.all.p[[1]])=="chromosome_name")-1)]
  
  for (cell_type in cell_type_groups) {
    keep_genes <- res.df.all.p[[2]][!is.na(res.df.all.p[[2]][,cell_type]),'names']
    res.df.all.p[[1]][!(res.df.all.p[[1]][,'names']%in%keep_genes),cell_type] <- NA
    res.df.all.lfc[[1]][!(res.df.all.lfc[[1]][,'names']%in%keep_genes),cell_type] <- NA
  }
  
  iter=2
  df1.lfc.perc <- res.df.all.lfc[[iter]]
  df1.lfc.perc[,cell_type_groups] <- lapply(res.df.all.lfc,function(x) {(apply(apply(x[,cell_type_groups],2,rank,na.last="keep"),2,function(x) x/sum(!is.na(x))))})[[iter]]
  
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
  res.df.all.p <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".p.rds"))
  res.df.all.lfc <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".lfc.rds"))
  cell_type_groups <- colnames(res.df.all.p[[1]])[(which(colnames(res.df.all.p[[1]])=="names")+1):(which(colnames(res.df.all.p[[1]])=="chromosome_name")-1)]
  
  for (cell_type in cell_type_groups) {
    keep_genes <- res.df.all.p[[2]][!is.na(res.df.all.p[[2]][,cell_type]),'names']
    res.df.all.p[[1]][!(res.df.all.p[[1]][,'names']%in%keep_genes),cell_type] <- NA
    res.df.all.lfc[[1]][!(res.df.all.lfc[[1]][,'names']%in%keep_genes),cell_type] <- NA
  }
  
  
  df1.p <- res.df.all.p[[iter]]
  df1.lfc <- res.df.all.lfc[[iter]]
  
  df1.p.melt <- reshape2::melt(aggregate(df1.p[,cell_type_groups]<0.1,by=list(df1.p$chromosome_name),mean,na.rm=T),id.vars="Group.1")
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

draw_heatmap("Femur",1)
draw_heatmap("Femur",2)
draw_heatmap("Liver",1)
draw_heatmap("Liver",2)

###################3

sampletype="Femur"; iter = 2
res.df.all.p <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".p.rds"))
res.df.all.lfc <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".lfc.rds"))
cell_type_groups <- colnames(res.df.all.p[[1]])[(which(colnames(res.df.all.p[[1]])=="names")+1):(which(colnames(res.df.all.p[[1]])=="chromosome_name")-1)]

for (cell_type in cell_type_groups) {
  keep_genes <- res.df.all.p[[2]][!is.na(res.df.all.p[[2]][,cell_type]),'names']
  res.df.all.p[[1]][!(res.df.all.p[[1]][,'names']%in%keep_genes),cell_type] <- NA
  res.df.all.lfc[[1]][!(res.df.all.lfc[[1]][,'names']%in%keep_genes),cell_type] <- NA
}


df1.p <- res.df.all.p[[iter]]
df1.lfc <- res.df.all.lfc[[iter]]

df1.p.melt <- reshape2::melt(aggregate(df1.p[,c("Pre pro B cells","Pro B cells","B cells")]<0.1,by=list(df1.p$chromosome_name==21),mean,na.rm=T),id.vars="Group.1")
ggplot(df1.p.melt,aes(x=variable,y=100*value,fill=Group.1)) + 
  geom_bar(stat='identity',position = 'dodge',col='black') +
  scale_fill_brewer(palette = 'Set2',labels=c("Not Chr 21","Chr 21")) + 
  theme_bw() + 
  theme(panel.grid = element_blank()) +
  labs(y="% significant DEG",x="Cell type",fill='') #+
  scale_fill_discrete(labels=c("Not Chr 21","Chr 21"))

