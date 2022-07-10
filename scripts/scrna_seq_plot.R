library(data.table)
library(ggplot2)
library(cowplot)
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
      axis.text.x = element_text(angle=60,hjust = 1,color='black'),
      panel.border = element_rect(colour = "black", fill=NA, size=1.5)
    ) +
    labs(x="Cell type",y="Median LFC (T21 vs Healthy)",fill='Chromosome') +
    scale_fill_manual(values=c("#E41A1C","#377EB8"))
  print(g)
  return(g)
}

g1=median_lfc_across_cell_types("Liver")
g2=median_lfc_across_cell_types("Femur")
system("mkdir -p ~/Documents/Research/t21-proj/out/full/figures")
x=1.3; pdf("~/Documents/Research/t21-proj/out/full/figures/median_lfc.pdf",width = 11*x,height=3*x)
plot_grid(g1,g2,nrow=1)
dev.off()

#####

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
                 plot.title=element_text(hjust=0.5),
                 axis.text=element_text(color='black'),
                 panel.border = element_rect(colour = "black", fill=NA, size=1.5)) +
           geom_smooth(method='loess',span=3,se=F,col='purple') + 
           labs(x=paste0("% DEG per cell type (chr ",chrNum,")"),y=paste0("% DEG per cell type (not chr ",chrNum,")"),title=paste0(sampletype," (chr",chrNum,")")))
           # labs(x=paste0("chr",chrNum),y=paste0("non-chr",chrNum,""),title=paste0(sampletype," chr",chrNum)))
}

glst <- list()
for (i in 1:22) {
  glst[[i]] <- DEG_comparison_chrX_v_notchrX("Liver",i)
}
plot_grid(plotlist=glst,ncol=4)

x=1.3; pdf("~/Documents/Research/t21-proj/out/full/figures/deg_chrX.liver.pdf",width = 6*x,height=3*x)
plot_grid(plotlist=glst[c(21,1)],nrow=1)
dev.off()

glst <- list()
for (i in 1:22) {
  glst[[i]] <- DEG_comparison_chrX_v_notchrX("Femur",i)
}
x=1.3; pdf("~/Documents/Research/t21-proj/out/full/figures/deg_chrX.femur.pdf",width = 6*x,height=3*x)
plot_grid(plotlist=glst[c(21,1)],nrow=1)
dev.off()

#########

df <- fread("~/Documents/Research/t21-proj/out/full/data/CellCycle.stats.txt",data.table = F,stringsAsFactors = F)
for (env in c("Femur","Liver")) {
  df.sub <- subset(df,sampletype==env)
  df.sub$Cycling <- 1 - df.sub$G1
  cell_types_to_use = subset(df.sub,disease_status=="Healthy")$cell[subset(df.sub,disease_status=="Healthy")$cell %in% subset(df.sub,disease_status=="DownSyndrome")$cell]
  # df.sub <- subset(df.sub,
  #                  cell %in% 
  #                    cell_types_to_use)
  df.sub1 <- subset(df.sub[,c("cell","Cycling")],df.sub$disease_status=="Healthy"); colnames(df.sub1)[2] <- "Healthy"
  df.sub2 <- subset(df.sub[,c("cell","Cycling")],df.sub$disease_status=="DownSyndrome"); colnames(df.sub2)[2] <- "DownSyndrome"
  df.mg <- merge(df.sub1,df.sub2,by='cell')
  print(wilcox.test(df.mg$DownSyndrome,df.mg$Healthy,paired = T))
  df.mg$cellLabel <- NA
  df.mg$cellLabel[df.mg$cell=="Pre pro B cells"] <- "Pre pro B cells"
  g=ggplot(df.mg,aes(x=100*Healthy,y=100*DownSyndrome)) + 
    geom_point(col='red3') +
    theme_bw() + 
    theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5),
          axis.text=element_text(color='black'),
          panel.border = element_rect(colour = "black", fill=NA, size=1.5)) +
    # geom_smooth(method='lm',se=F,col='red3') +
    geom_abline(slope=1,intercept=0,lty='dashed') +
    labs(x="% Cycling (per cell type) - Healthy",y="% Cycling (per cell type) - T21",title=env) 
  g1 = g +
    geom_text(aes(label=cell),vjust=2)
  g2 = g #+ geom_text(aes(label=cellLabel),vjust=2)
  pdf(paste0("~/Documents/Research/t21-proj/out/full/figures/",env,"_cycling.pdf"),width=8,height=4)
  print(plot_grid(g1,g2,nrow=1))
  dev.off()
}

############

library(Seurat)
library(stringr)
source("~/Documents/Research/Useful_scripts/rntransform.R")

cell_type_filename="HSCs_MPPs"
f.out=paste0("~/Documents/Research/t21-proj/out/full/data/cell_subset/",cell_type_filename,".all_sig_genes.rds")
res <- fread("~/Documents/Research/t21-proj/out/full/DEG_list/liver_v_femur.txt",data.table = F,stringsAsFactors = F)
dfcombined=readRDS(f.out)
upreg_t21_induced <- apply(apply(dfcombined[["RNA"]][subset(res,class=="t21-induced" & logFC.t21 > 0)$names,],1,scale),1,mean)
env <- factor(dfcombined$environment,levels = c("Healthy","Down Syndrome"))
glst = list()
df.mg = list()
env <- factor(dfcombined$environment,levels = c("Healthy","Down Syndrome"))
# for (class_to_use in c("environment-independent","liver-induced")) {
for (class_to_use in c("environment-driven","t21-induced")) {
  
  input=subset(res,class==class_to_use & logFC.t21 > 0)$names
  
  y <- apply(apply(dfcombined[["RNA"]][input,],1,scale),1,mean)
  y <- rntransform(y)
  summary(lm(y~env*organ,dftmp <- data.frame(y=y,env=env,organ=dfcombined$organ)))
  
  df.mg[[class_to_use]] = merge(aggregate(y,by=list(dfcombined$environment,dfcombined$organ),function(x) sd(x)/sqrt(length(input))),
                                aggregate(y,by=list(dfcombined$environment,dfcombined$organ),mean),
                                by=c("Group.1","Group.2"))
  
  glst[[class_to_use]] <- ggplot(df.mg[[class_to_use]],aes(x=as.numeric(as.factor(Group.2)),
                                                           y=x.y,
                                                           col=Group.1,
                                                           ymin=x.y-2*x.x,
                                                           ymax=x.y+2*x.x)) +
    geom_line() + 
    geom_point(size=rel(2)) +
    geom_errorbar(aes(ymin=x.y-2*x.x,ymax=x.y+2*x.x),width=.2) + #,position=position_dodge(0.05)) +
    theme_bw() + 
    theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5),
          axis.text=element_text(color='black'),
          panel.border = element_rect(colour = "black", fill=NA, size=1.5)) +
    labs(x="Environment",y="Average normalized expression (+/- 2 SE)",col="Disease status",title=str_to_title(class_to_use)) +
    scale_x_continuous(labels=c("Femur","Liver"),breaks=c(1,2)) +
    scale_colour_manual(values=c("#E41A1C","#377EB8"))
}

pdf(paste0("~/Documents/Research/t21-proj/out/full/figures/t21_induced_gxe.pdf"),width=10,height=4)
plot_grid(plotlist = glst)
dev.off()


#
library(data.table)
library("enrichR")
dbs <- c("GO_Molecular_Function_2021", "GO_Biological_Process_2021","ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")

df1.full.mg = fread("~/Documents/Research/t21-proj/out/full/DEG_list/liver_v_femur.txt",data.table = F,stringsAsFactors = F)

gene_lst <- subset(df1.full.mg,class=="t21-induced" & logFC.t21 > 0)$names
enriched <- enrichr(gene_lst, dbs)
y <- enriched[[dbs[1]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); y = (subset(y,geneCt>=3 & Adjusted.P.value<0.1)); y$Set <- dbs[1]; ykeep <- y
y <- enriched[[dbs[2]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); ytmp = (subset(y,geneCt>=3 & Adjusted.P.value<0.1)); ytmp$Set <- dbs[2]; ykeep<-rbind(ykeep,ytmp)
y <- enriched[[dbs[3]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); ytmp = (subset(y,geneCt>=3 & Adjusted.P.value<0.1)); ytmp$Set <- dbs[3]; ykeep<-rbind(ykeep,ytmp)

g <- ggplot(ykeep,aes(x=Set,y=-log10(Adjusted.P.value),col=Set)) +
  geom_jitter(width=0.1) + 
  theme_bw()  +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 60, hjust=1,color='black'), #,colour = axis_colors),
        axis.text=element_text(color='black'),
        panel.border = element_rect(colour = "black", fill=NA, size=1.5),
        legend.position = 'none'
  ) +
  labs(x='Set',y= expression(-log[10](italic(FDR)))) + scale_color_brewer(palette="Set2") +
  scale_x_discrete(labels=c("ENCODE + ChEA TFs","GO BP","GO MF"))
# coord_flip()

pdf(paste0("~/Documents/Research/t21-proj/out/full/figures/t21_induced_gxe.pdf"),width=5,height=5)
print(g)
dev.off()

fwrite(ykeep,"~/Documents/Research/t21-proj/out/full/figures/t21_induced_gxe.csv",row.names = F)

