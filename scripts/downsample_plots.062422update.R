library(data.table)
library(ggplot2)
library(cowplot)
library(reshape2)
df1.full.p <- fread("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names_downsamp/Liver.HSCs_MPPs.sample.permALL.p.txt",data.table = F,stringsAsFactors = F)
df1.full.fdr <- fread("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names_downsamp/Liver.HSCs_MPPs.sample.permALL.fdr.txt",data.table = F,stringsAsFactors = F)
df1.full.lfc <- fread("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names_downsamp/Liver.HSCs_MPPs.sample.permALL.lfc.txt",data.table = F,stringsAsFactors = F)

# aggregate(df1.full.fdr[,c(8:12)],by=list(chr21=df1.full.fdr$chromosome_name==21),function(x) {mean(x<0.05)})

# M=100
# df1.full.fdr$medFDR = apply(df1.full.fdr[,paste0("sim",1:M)],1,median)
# aggregate(df1.full.fdr$medFDR,by=list(chr21=df1.full.fdr$chromosome_name==21),function(x) {mean(x<0.05)})
# df1.full.fdr[order(df1.full.fdr$medFDR)[1],]

M=100
df1.full.p$medP <- apply(df1.full.p[,paste0("sim",1:M)],1,median)
M=75
df1.full.p$medP.75 <- apply(df1.full.p[,paste0("sim",1:M)],1,median)
M=25
df1.full.p$medP.25 <- apply(df1.full.p[,paste0("sim",1:M)],1,median)
M=10
df1.full.p$medP.10 <- apply(df1.full.p[,paste0("sim",1:M)],1,median)
cor(df1.full.p$medP.75,df1.full.p$medP)
g1=ggplot(df1.full.p,aes(x=sim1,y=medP)) + 
  geom_point() + geom_abline(slope=1,intercept=0,col='red',lty='dashed') + 
  theme_bw() + theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5)) +
  labs(x="Median P-value after 1 downsampling",y="Median P-value after 100 downsamplings")
g2=ggplot(df1.full.p,aes(x=medP.10,y=medP))  + 
  geom_point() + geom_abline(slope=1,intercept=0,col='red',lty='dashed') + 
  theme_bw() + theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5)) +
  labs(x="Median P-value after 10 downsamplings",y="Median P-value after 100 downsamplings")
g2b=ggplot(df1.full.p,aes(x=medP.25,y=medP))  + 
  geom_point() + geom_abline(slope=1,intercept=0,col='red',lty='dashed') + 
  theme_bw() + theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5)) +
  labs(x="Median P-value after 25 downsamplings",y="Median P-value after 100 downsamplings")
g3=ggplot(df1.full.p,aes(x=medP.75,y=medP))   + 
  geom_point() + geom_abline(slope=1,intercept=0,col='red',lty='dashed') + 
  theme_bw() + theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5)) +
  labs(x="Median P-value after 75 downsamplings",y="Median P-value after 100 downsamplings")
plot_grid(g1,g2,g2b,g3,nrow=1)

df1.full.p$adj.P.Val.liv <- df1.full.fdr$adj.P.Val.liv
df1.full.p$adj.P.Val.fem <- df1.full.fdr$adj.P.Val.fem
df1.full.p$class <- "none"
df1.full.p$class[df1.full.p$adj.P.Val.liv < 0.05 | df1.full.p$adj.P.Val.fem < 0.05] <- "unknown"
df1.full.p$class[df1.full.p$adj.P.Val.liv < 0.05 & df1.full.p$P.Value.fem < 0.05] <- "environment-independent"
df1.full.p$class[df1.full.p$adj.P.Val.fem < 0.05 & df1.full.p$P.Value.liv < 0.05] <- "environment-independent"
df1.full.p$class[df1.full.p$adj.P.Val.liv < 0.05 & df1.full.p$medP < 0.05 & df1.full.p$P.Value.fem >= 0.05] <- "liver-induced"
df1.full.p$class[df1.full.p$adj.P.Val.fem < 0.05 & df1.full.p$P.Value.liv >= 0.05] <- "femur-induced"
table(df1.full.p$class)
tab <- table(subset(df1.full.p,class != "none")$class)
pct <- round(tab/sum(tab)*100,1)
pie(tab,labels=paste0(names(tab)," (",unname(tab),', ',unname(pct),"%)"),radius = 1,col = c("purple", "violetred1", "green3","cornsilk"))

library(enrichR)
# df1.full.mg <- merge(df1.full.p,df1.full.lfc[,c("names","logFC.t21")],by='names')
df1.full.mg <- merge(df1.full.p,df1.full.lfc[,c("names","logFC.liv","logFC.fem")],by='names')
dbs <- c("GO_Molecular_Function_2021", "GO_Biological_Process_2021","ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")
# gene_lst <- subset(df1.full.p,class=="environment-driven")$names# 
# gene_lst <- subset(df1.full.mg,class=="liver-induced" & logFC.liv > 0)$names
gene_lst <- subset(df1.full.mg,class=="liver-induced" & logFC.liv < 0)$names
# gene_lst <- subset(df1.full.mg,class=="unknown" & logFC.liv < 0)$names

# gene_lst <- subset(df1.full.p,class=="t21-induced")$names
# gene_lst <- subset(df1.full.p,class=="t21-induced" & chromosome_name != 21)$names
# subset(df1.full.p,class=="t21-reverted")$names
enriched <- enrichr(gene_lst, dbs)
y <- enriched[[dbs[1]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); (subset(y,geneCt>=3 & Adjusted.P.value<0.1))$Term
y <- enriched[[dbs[2]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); (subset(y,geneCt>=3 & Adjusted.P.value<0.1))$Term
y <- enriched[[dbs[3]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); (subset(y,geneCt>=3 & Adjusted.P.value<0.1))$Term

subset(df1.full.mg,names=="HLA-DMA")
subset(df1.full.mg,names=="HLA-DRA")
subset(df1.full.mg,names=="HLA-DRB1")
subset(df1.full.mg,names=="CIITA")
# subset(df1.full.lfc,names=="HLA-DRB1")
subset(df1.full.lfc,names=="CIITA")


gene_names="GATA1"
gene_names="APOC1"
gene_names=c("GATA1","APOC1","MYL4","DANT2","CERS5")
glist <- list(); dflst <- list()
for (j in 1:length(gene_names)) {
  gene_of_interest=gene_names[j]
  # k=which(df1.full.p$names=="GATA1")
  # k=which(df1.full.p$names=="APOC1")
  k=which(df1.full.p$names==gene_of_interest)
  rng=seq(0,100,by=2)[-1]
  medP = rep(NA,length(rng))
  for (i in 1:length(rng)) {
    M=rng[i]
    medP[i] = apply(df1.full.p[k,paste0("sim",1:M)],1,median)
    # medP[i] = apply(df1.full.p[k,paste0("sim",1:M)],1,mean)
  }
  dflst[[j]] <- data.frame(rng,medP)
  glist[[j]] <- ggplot(dflst[[j]],aes(x=rng,y=medP)) + 
    geom_line() + geom_point() + 
    theme_bw() + 
    theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5)) +
    labs(x="Downsampling iteration",y='Median P-value',title=bquote(~italic(.(gene_of_interest)))) +
    ylim(0,1)
}

plot_grid(glist[[1]],glist[[2]],glist[[3]],glist[[4]],glist[[5]],nrow = 1)

ggplot(df1.full.p,aes(x=-log10(P.Value.liv),y=-log10(medP))) + 
  geom_point(col='red3') + 
  geom_smooth(method = 'lm',se=F,col='black') + 
  theme_bw() +
  labs(x="-log10 P (liver)",y="-log10 P (median downsamp)")
cor(-log10(df1.full.p$P.Value.liv),-log10(df1.full.p$medP))
cor(-log10(df1.full.p$sim1),-log10(df1.full.p$sim2))
cor(-log10(df1.full.p$sim1),-log10(df1.full.p$sim3))
cor(-log10(df1.full.p$sim1),-log10(df1.full.p$sim34))

cor(df1.full.p$P.Value.liv,df1.full.p$medP)

fwrite(df1.full.mg[order(df1.full.mg$P.Value.liv),c("names","chromosome_name","logFC.fem","logFC.liv","P.Value.fem","adj.P.Val.fem","P.Value.liv","adj.P.Val.liv","medP","class")],"~/Downloads/t21_v_healthy.txt",quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

subset(df1.full.mg,)


