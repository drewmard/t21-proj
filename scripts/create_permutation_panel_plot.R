library(data.table)
library(ggplot2)
library(cowplot)
library(reshape2)
df1.full.lfc <- readRDS("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.Femur.lfc.rds")
subset(df1.full.lfc[[2]],names=="GATA1")
subset(df1.full.lfc[[2]],names=="GATA1")

df1.full.lfc <- fread("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names_downsamp/DownSyndrome.HSCs_MPPs.sample.permALL.lfc.txt",data.table = F,stringsAsFactors = F)
df1.full.p <- fread("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names_downsamp/DownSyndrome.HSCs_MPPs.sample.permALL.p.txt",data.table = F,stringsAsFactors = F)

gene_of_interest="GATA1"

create_permutation_panel_plot <- function(gene_of_interest) {
  tmp <- subset(df1.full.lfc,names==gene_of_interest)
  tmp2 <- reshape2::melt(tmp[,c('names',paste0('sim',1:50))],id.vars="names")
  g1 <- ggplot(tmp2,aes(x=value)) + geom_histogram(fill='orange',col='black',bins=20) + theme_bw() +
    theme_bw() + theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5)) +
    geom_vline(xintercept=tmp$logFC.h,col='red') + #xlim(tmp$logFC.h*0.9,max(tmp2$value)*1.1)  +
    geom_vline(xintercept=tmp$logFC.t21,col='blue') +
    labs(x='Log-fold change',title=bquote(~italic(.(gene_of_interest))))
  tmp3 <- subset(df1.full.p,names==gene_of_interest)
  tmp4 <- reshape2::melt(tmp3[,c('names',paste0('sim',1:50))],id.vars="names")
  mean(tmp4$value < tmp3$P.Value.h)
  g2 <- ggplot(tmp4,aes(x=value)) + geom_histogram(fill='bisque2',col='black',bins=20) + 
    theme_bw() + theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5)) +
    geom_vline(xintercept=tmp3$P.Value.h,col='red') + #xlim(min(tmp4$value)*0.9,1.1*tmp3$P.Value.h) +
    labs(x='P-value',title=bquote(~italic(.(gene_of_interest))))
  print(plot_grid(g1,g2,ncol = 2))
}

tmp3 <- subset(df1.full.p,names==gene_of_interest)
p.val <- rep(NA,nrow(df1.full.p))
for (i in 1:nrow(df1.full.p)) { 
  p.val[i] <- mean(df1.full.p$P.Value.h[i] < df1.full.p[i,paste0("sim",1:50)])
}
sum(df1.full.p$P.Value.h < 0.05 & df1.full.p$P.Value.t21 < 0.05)
ggplot(data.frame(p.val),aes(p.val)) + geom_histogram(col='black',fill='orange',bins=40) + 
  theme_bw() + labs(x="% of t21 downsamples where P-value is less than Healthy P-value")
ggplot(data.frame(p.val)[df1.full.p$P.Value.t21 < 0.05,,drop=FALSE],aes(p.val)) + geom_histogram(col='black',fill='orange',bins=40) + 
  theme_bw() + labs(x="% of t21 downsamples where P-value is less than Healthy P-value")

df1.full.p[df1.full.p$P.Value.h < 0.05 & df1.full.p$P.Value.t21 < 0.05,][1,]

gene_of_interest="GATA1"; create_permutation_panel_plot(gene_of_interest); subset(df1.full.lfc,names==gene_of_interest)
gene_of_interest="APOC1"; create_permutation_panel_plot(gene_of_interest); subset(df1.full.lfc,names==gene_of_interest)
gene_of_interest="MYL4"; create_permutation_panel_plot(gene_of_interest); subset(df1.full.lfc,names==gene_of_interest)

df1.full.fdr <- fread("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names_downsamp/DownSyndrome.HSCs_MPPs.sample.permALL.fdr.txt",data.table = F,stringsAsFactors = F)
mean(df1.full.fdr$adj.P.Val.h < 0.1)
mean(df1.full.fdr$adj.P.Val.t21 < 0.1)
tmp <- aggregate(df1.full.fdr[,-c(1:2)],by=list(df1.full.fdr$chromosome_name==21),function(x) {mean(x<0.1)})
# tmp <- aggregate(df1.full.p[,-c(1:2)],by=list(df1.full.fdr$chromosome_name==21),function(x) {mean(x<0.1)})
tmp2 <- reshape2::melt(tmp[,c('Group.1',paste0('sim',1:100))],id.vars="Group.1")
tmp3 <- subset(tmp2,Group.1)
g1 <- ggplot(tmp3,aes(x=value)) + geom_histogram(fill='orange',col='black',bins=30) + theme_bw() +
  theme_bw() + theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5)) +
  geom_vline(xintercept=subset(tmp,Group.1)$adj.P.Val.h,col='red') + #xlim(tmp$logFC.h*0.9,max(tmp2$value)*1.1)  +
  labs(x='% DEG',title="Chr21")
g1
tmp4 <- subset(tmp2,!Group.1)
g2 <- ggplot(tmp4,aes(x=value)) + geom_histogram(fill='orange',col='black',bins=30) + theme_bw() +
  theme_bw() + theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5)) +
  geom_vline(xintercept=subset(tmp,!Group.1)$adj.P.Val.h,col='red') + #xlim(tmp$logFC.h*0.9,max(tmp2$value)*1.1)  +
  labs(x='% DEG',title="Non-Chr21")
g2
plot_grid(g1,g2)
mean(tmp3$value > subset(tmp,Group.1)$adj.P.Val.h)
mean(tmp4$value > subset(tmp,!Group.1)$adj.P.Val.h)
gene_of_interest="GATA1"; subset(df1.full.p,names==gene_of_interest)

tmp <- aggregate(df1.full.p[,-c(1:2)],by=list(df1.full.fdr$chromosome_name==21),function(x) {mean(x<0.05)})
tmp2 <- reshape2::melt(tmp[,c('Group.1',paste0('sim',1:100))],id.vars="Group.1")
tmp3 <- subset(tmp2,Group.1)
g1 <- ggplot(tmp3,aes(x=value)) + geom_histogram(fill='orange',col='black',bins=30) + theme_bw() +
  theme_bw() + theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5)) +
  geom_vline(xintercept=subset(tmp,Group.1)$P.Value.h,col='red') + #xlim(tmp$logFC.h*0.9,max(tmp2$value)*1.1)  +
  labs(x='% genes with nominal P-value < 0.05',title="Chr21")
g1
tmp4 <- subset(tmp2,!Group.1)
g2 <- ggplot(tmp4,aes(x=value)) + geom_histogram(fill='orange',col='black',bins=30) + theme_bw() +
  theme_bw() + theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5)) +
  geom_vline(xintercept=subset(tmp,!Group.1)$P.Value.h,col='red') + #xlim(tmp$logFC.h*0.9,max(tmp2$value)*1.1)  +
  labs(x='% genes with nominal P-value < 0.05',title="Non-Chr21")
g2
plot_grid(g1,g2)

mean(tmp3$value > subset(tmp,Group.1)$P.Value.h)
mean(tmp4$value > subset(tmp,!Group.1)$P.Value.h)

cor(df1.full.lfc$logFC.h,df1.full.lfc$logFC.t21)
cor(df1.full.lfc$logFC.t21,df1.full.lfc$sim1)
cor(df1.full.lfc$logFC.t21,df1.full.lfc$sim2)
cor(df1.full.lfc$logFC.t21,df1.full.lfc$sim3)
cor(df1.full.lfc$logFC.h,df1.full.lfc$sim1)
cor(df1.full.lfc$logFC.h,df1.full.lfc$sim2)
cor(df1.full.lfc$logFC.h,df1.full.lfc$sim3)

ggplot(df1.full.lfc,aes(x=logFC.h,y=logFC.t21)) + 
  geom_point(col='red3') + 
  geom_smooth(method = 'lm',se=F,col='black') + 
  theme_bw() +
  labs(x="logFC (Healthy)",y="logFC (t21)")
ggplot(df1.full.lfc,aes(x=logFC.t21,y=sim7)) +
# ggplot(df1.full.lfc,aes(x=logFC.h,y=sim1)) + 
  geom_point(col='red3') + 
  geom_smooth(method = 'lm',se=F,col='black') + 
  geom_abline(slope = 1,intercept = 0,col='steelblue',lty='dashed') +
  theme_bw() +
  labs(x="logFC (t21)",y="logFC (downsamp)")
ggplot(df1.full.lfc,aes(x=sim7,y=logFC.t21)) +
  # ggplot(df1.full.lfc,aes(x=logFC.h,y=sim1)) + 
  geom_point(col='red3') + 
  geom_smooth(method = 'lm',se=F,col='black') + 
  geom_abline(slope = 1,intercept = 0,col='steelblue',lty='dashed') +
  theme_bw() +
  labs(x="logFC (t21)",y="logFC (downsamp)")

plot(df1.full.lfc$logFC.h,df1.full.lfc$logFC.t21)
plot(df1.full.lfc$logFC.t21,df1.full.lfc$sim1)
plot(df1.full.lfc$logFC.t21,df1.full.lfc$sim2)
plot(df1.full.lfc$logFC.t21,df1.full.lfc$sim3)
plot(df1.full.lfc$logFC.h,df1.full.lfc$sim1)
plot(df1.full.lfc$logFC.h,df1.full.lfc$sim2)
plot(df1.full.lfc$logFC.h,df1.full.lfc$sim3)

