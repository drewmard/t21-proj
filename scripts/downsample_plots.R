library(data.table)
library(ggplot2)
library(cowplot)
library(reshape2)
df1.full.p <- fread("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names_downsamp/DownSyndrome.HSCs_MPPs.sample.permALL.p.txt",data.table = F,stringsAsFactors = F)
df1.full.fdr <- fread("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names_downsamp/DownSyndrome.HSCs_MPPs.sample.permALL.fdr.txt",data.table = F,stringsAsFactors = F)
df1.full.p$sim1[df1.full.p$names=="KLF1"] <- 0.04323491
df1.full.lfc <- fread("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names_downsamp/DownSyndrome.HSCs_MPPs.sample.permALL.lfc.txt",data.table = F,stringsAsFactors = F)


# df1.full.p[1500,]

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

df1.full.p$adj.P.Val.t21 <- df1.full.fdr$adj.P.Val.t21
ggplot(df1.full.p,aes(x=-log10(P.Value.t21),y=-log10(medP))) + 
  geom_point(col='red3') + 
  geom_smooth(method = 'lm',se=F,col='black') + 
  theme_bw() +
  labs(x="logFC (t21)",y="logFC (median downsamp)")
cor(df1.full.p$P.Value.t21,df1.full.p$medP)
cor(df1.full.lfc$logFC.t21,df1.full.lfc$sim1)
cor(df1.full.p$sim1,df1.full.p$sim2)
cor(df1.full.p$sim1,df1.full.p$sim3)
cor(df1.full.p$sim2,df1.full.p$sim3)
cor(df1.full.lfc$sim1,df1.full.lfc$sim2)
cor(df1.full.lfc$sim1,df1.full.lfc$sim3)
cor(df1.full.lfc$sim2,df1.full.lfc$sim3)

cor(df1.full.p$P.Value.h,df1.full.p$medP)
cor(df1.full.p$P.Value.h,df1.full.p$P.Value.t21)

ggplot(df1.full.lfc,aes(x=(logFC.t21),y=(sim1))) + 
  geom_point(col='red3') + 
  geom_smooth(method = 'lm',se=F,col='black') + 
  theme_bw() +
  labs(x="logFC (t21)",y="logFC (median downsamp)")

M=100
df1.full.p$medP <- apply(df1.full.p[,paste0("sim",1:M)],1,median)
df1.full.p$adj.P.Val.t21 <- df1.full.fdr$adj.P.Val.t21
df1.full.p$adj.P.Val.h <- df1.full.fdr$adj.P.Val.h
df1.full.p$class <- "none"
df1.full.p$class[df1.full.p$adj.P.Val.t21 < 0.05 | df1.full.p$adj.P.Val.h < 0.05] <- "unknown"
df1.full.p$class[df1.full.p$adj.P.Val.t21 < 0.05 & df1.full.p$P.Value.h < 0.05] <- "environment-driven"
df1.full.p$class[df1.full.p$adj.P.Val.h < 0.05 & df1.full.p$P.Value.t21 < 0.05] <- "environment-driven"
df1.full.p$class[df1.full.p$adj.P.Val.t21 < 0.05 & df1.full.p$medP < 0.05 & df1.full.p$P.Value.h >= 0.05] <- "t21-induced"
df1.full.p$class[df1.full.p$adj.P.Val.h < 0.05 & df1.full.p$P.Value.t21 >= 0.05] <- "t21-reverted"
table(df1.full.p$class)
table(subset(df1.full.p,adj.P.Val.t21<0.05)$class)
table(subset(df1.full.p,adj.P.Val.h<0.05)$class)
tab <- table(subset(df1.full.p,class != "none")$class)
pct <- round(tab/sum(tab)*100,1)
pie(tab,labels=paste0(names(tab)," (",unname(tab),', ',unname(pct),"%)"),radius = 1,col = c("purple", "violetred1", "green3","cornsilk"))

table(subset(df1.full.p,chromosome_name==21)$class)
subset(df1.full.p,chromosome_name==21 & class=="t21-induced")$names

library(enrichR)
library("enrichR")

# df1.full.mg <- merge(df1.full.p,df1.full.lfc[,c("names","logFC.t21")],by='names')
df1.full.mg <- merge(df1.full.p,df1.full.lfc[,c("names","logFC.t21","logFC.h")],by='names')
dbs <- c("GO_Molecular_Function_2021", "GO_Biological_Process_2021","ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")
# gene_lst <- subset(df1.full.p,class=="environment-driven")$names# 
gene_lst <- subset(df1.full.mg,class=="t21-induced" & logFC.t21 > 0)$names
# gene_lst <- subset(df1.full.mg,class=="t21-induced" & logFC.t21 < 0)$names
# gene_lst <- subset(df1.full.p,class=="t21-induced")$names
# gene_lst <- subset(df1.full.p,class=="t21-induced" & chromosome_name != 21)$names
# subset(df1.full.p,class=="t21-reverted")$names
enriched <- enrichr(gene_lst, dbs)
y <- enriched[[dbs[1]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); (subset(y,geneCt>3 & Adjusted.P.value<0.1))$Term
y <- enriched[[dbs[2]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); (subset(y,geneCt>3 & Adjusted.P.value<0.1))$Term
y <- enriched[[dbs[3]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); (subset(y,geneCt>3 & Adjusted.P.value<0.1))$Term

subset(y,Term=="myeloid cell differentiation (GO:0030099)")

subset(df1.full.p,names=="BTG3")
subset(df1.full.p,names=="NFE2L2")

subset(df1.full.lfc,names=="BTG3")

df1.full.p[order(df1.full.p$P.Value.t21)[1:5],c("names","P.Value.h","P.Value.t21","adj.P.Val.t21","medP","class")]
create_permutation_panel_plotII <- function(gene_of_interest,M=100) {
  tmp3 <- subset(df1.full.p,names==gene_of_interest)
  tmp4 <- reshape2::melt(tmp3[,c('names',paste0('sim',1:M))],id.vars="names")
  mean(tmp4$value < tmp3$P.Value.h)
  g2 <- ggplot(tmp4,aes(x=value)) + geom_histogram(fill='bisque2',col='black',bins=20) + 
    theme_bw() + theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5)) +
    geom_vline(xintercept=tmp3$P.Value.h,col='red') + #xlim(min(tmp4$value)*0.9,1.1*tmp3$P.Value.h) +
    geom_vline(xintercept=tmp3$medP,col='blue') +
    geom_vline(xintercept=0.05,col='black',lty='dashed') +
    labs(x='P-value',title=bquote(~italic(.(gene_of_interest))))
  return(print(g2))
}
g1 <- create_permutation_panel_plotII("MYL4")
g2 <- create_permutation_panel_plotII("GATA1")
g3 <- create_permutation_panel_plotII("APOC1")
g4 <- create_permutation_panel_plotII("AZU1")
g5 <- create_permutation_panel_plotII("DLGAP5")

plot_grid(g1,g2,g3,g4,g5,nrow=1)
create_permutation_panel_plotII("CASP1")
create_permutation_panel_plotII("DLGAP5")
create_permutation_panel_plotII("APOC1")

fwrite(df1.full.mg[order(df1.full.mg$P.Value.t21),c("names","chromosome_name","logFC.h","logFC.t21","P.Value.h","adj.P.Val.h","P.Value.t21","adj.P.Val.t21","medP","class")],"~/Downloads/liver_v_femur.txt",quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
# fwrite(df1.full.p[order(df1.full.p$P.Value.t21),c("names","chromosome_name","P.Value.h","adj.P.Val.h","P.Value.t21","adj.P.Val.t21","medP","class")],"~/Downloads/liver_v_femur.txt",quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)


df1.full.p[df1.full.p$P.Value.t21 >= 0.05 & df1.full.p$medP < 0.05,][1,]


