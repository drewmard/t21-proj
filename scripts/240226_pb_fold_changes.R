library(data.table)
library(ggplot2)
res.save=list(); i =0; sampletype="DownSyndrome"
# res.save=list(); i =0; for (sampletype in c("DownSyndrome","Healthy")) { 
i = i+1
res.mg = fread(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/pseudobulks/DE.",sampletype,".txt"),data.table = F,stringsAsFactors = F)
# res.mg = fread("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/pseudobulks/DE.DownSyndrome.txt",data.table = F,stringsAsFactors = F)
# ggplot(res.mg,aes(x=logFC.sample.FALSE,y=logFC.patient.FALSE)) + geom_point() + geom_abline(slope=1,intercept = 0,lty='dashed',col='red') + ggpubr::theme_pubr() + labs(x="LFC (sample-level pseudobulk)",y="LFC (fetus-level pseudobulk)")
# ggplot(res.mg,aes(x=logFC.sample.TRUE,y=logFC.patient.TRUE)) + geom_point() + geom_abline(slope=1,intercept = 0,lty='dashed',col='red')

colnames(res.mg)[2:5] = c("Fetus (+ age)","Fetus","Sample (+ age)","Sample")
cor.mat = cor(res.mg[,-1],use='na.or.complete')
cor.melt = melt(cor.mat)

g <- ggplot(cor.melt,aes(x=Var1,y=Var2,fill=value)) +
  geom_tile(col='black') +
  theme_bw()  +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 60, hjust=1),
        panel.border = element_blank(),
        legend.position = 'none',
        plot.title = element_text(hjust=0.5)
  ) +
  scale_fill_gradient2(low = "blue", mid='white',high = "red",name='Power') +
  geom_text(aes(label=round(value,3))) + 
  labs(x='Differential expression analysis',y='Differential expression analysis',title=sampletype); g

pdf(paste0("~/Documents/Research/t21-proj/out/full/pseudobulks/lfc_correlation.",sampletype,".pdf"),width=6,height=6)
print(g)
dev.off()

genes = fread("~/Documents/Research/data/grch38/genes.extended.txt",data.table = F,stringsAsFactors = F)
res.mg2 = merge(res.mg,genes[,c("Gene name","Chromosome/scaffold name")],by.x="names",by.y = "Gene name")
res.mg2 = subset(res.mg2,`Chromosome/scaffold name` != 21)[,colnames(res.mg2)!="Chromosome/scaffold name"]
cor.mat = cor(res.mg2[,-1],use='na.or.complete')
cor.melt = melt(cor.mat)

g <- ggplot(cor.melt,aes(x=Var1,y=Var2,fill=value)) +
  geom_tile(col='black') +
  theme_bw()  +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 60, hjust=1),
        panel.border = element_blank(),
        legend.position = 'none',
        plot.title = element_text(hjust=0.5)
  ) +
  scale_fill_gradient2(low = "blue", mid='white',high = "red",name='Power') +
  geom_text(aes(label=round(value,3))) + 
  labs(x='Differential expression analysis',y='Differential expression analysis',title=sampletype); g
g <- g + scale_x_discrete(limits = rev(levels(cor.melt$Var1)))  # Flip the x-axis and place tick labels at the bottom
pdf(paste0("~/Documents/Research/t21-proj/out/full/pseudobulks/lfc_correlation.",sampletype,".no_chr21.pdf"),width=6,height=6)
print(g)
dev.off()

de_res = fread("~/Documents/Research/t21-proj/out/full/pseudobulks/sample/DE/DownSyndrome_HSCs_MPPs.de.txt",data.table = F,stringsAsFactors = F)
res.mg3 = subset(res.mg2,names %in% subset(de_res,P.Value < 0.05)$names)
cor.mat = cor(res.mg3[,-1],use='na.or.complete')
cor.melt = melt(cor.mat)
g <- ggplot(cor.melt,aes(x=Var1,y=Var2,fill=value)) +
  geom_tile(col='black') +
  theme_bw()  +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 60, hjust=1),
        panel.border = element_blank(),
        legend.position = 'none',
        plot.title = element_text(hjust=0.5)
  ) +
  scale_fill_gradient2(low = "blue", mid='white',high = "red",name='Power') +
  geom_text(aes(label=round(value,3))) + 
  labs(x='Differential expression analysis',y='Differential expression analysis',title=sampletype); g
g <- g + scale_x_discrete(limits = rev(levels(cor.melt$Var1)))  # Flip the x-axis and place tick labels at the bottom
pdf(paste0("~/Documents/Research/t21-proj/out/full/pseudobulks/lfc_correlation.",sampletype,".no_chr21.pdf"),width=6,height=6)
print(g)
dev.off()

res.mg3 = subset(res.mg2,names %in% subset(de_res,adj.P.Val < 0.05)$names)
cor.mat = cor(res.mg3[,-1],use='na.or.complete')
cor.melt = melt(cor.mat)
g <- ggplot(cor.melt,aes(x=Var1,y=Var2,fill=value)) +
  geom_tile(col='black') +
  theme_bw()  +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 60, hjust=1),
        panel.border = element_blank(),
        legend.position = 'none',
        plot.title = element_text(hjust=0.5)
  ) +
  scale_fill_gradient2(low = "blue", mid='white',high = "red",name='Power') +
  geom_text(aes(label=round(value,3))) + 
  labs(x='Differential expression analysis',y='Differential expression analysis',title=sampletype); g

res.mg3 = subset(res.mg2,names %in% head(de_res,200)$names)
cor.mat = cor(res.mg3[,-1],use='na.or.complete')
cor.melt = melt(cor.mat)
g <- ggplot(cor.melt,aes(x=Var1,y=Var2,fill=value)) +
  geom_tile(col='black') +
  theme_bw()  +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 60, hjust=1),
        panel.border = element_blank(),
        legend.position = 'none',
        plot.title = element_text(hjust=0.5)
  ) +
  scale_fill_gradient2(low = "blue", mid='white',high = "red",name='Power') +
  geom_text(aes(label=round(value,3))) + 
  labs(x='Differential expression analysis',y='Differential expression analysis',title=sampletype); g
plot(res.mg3[,2:5])

cb = fread("~/Downloads/cellbender/Down Syndrome_HSCs_MPPs.de.txt",data.table = F,stringsAsFactors = F)
cb = cb[,c("names","logFC")]
colnames(cb) = c("names","CellBender")
df.mg = merge(res.mg2,cb,by="names")
# df = fread("~/Downloads/DownSyndrome_HSC_PeakGeneSets/liver_v_femur.txt",data.table = F,stringsAsFactors = F)
# df.mg = merge(df.mg,df[,c("names","logFC.t21")],by="names")
cor(df.mg[,-1])
res.mg3 = subset(df.mg,names %in% subset(de_res,P.Value < 0.05)$names)
plot(res.mg3$`Fetus (+ age)`,res.mg3$CellBender)
cor(res.mg3[,-1])
res.mg3 = subset(df.mg,names %in% subset(de_res,adj.P.Val < 0.05)$names)
cor(res.mg3[,-1])
res.mg3 = subset(df.mg,names %in% head(de_res,200)$names)
cor(res.mg3[,-1])







res.save = merge(res.save[[1]],res.save[[2]],by='names')
cor.mat=cor(res.save[,-1],use='na.or.complete')
cor.melt = melt(cor.mat)
g <- ggplot(cor.melt,aes(x=Var1,y=Var2,fill=scale(value))) +
  geom_tile(col='black') +
  theme_bw()  +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 60, hjust=1),
        panel.border = element_blank(),
        legend.position = 'none',
        plot.title = element_text(hjust=0.5)
  ) +
  scale_fill_gradient2(low = "orange", mid='red',high = "red",name='Power') +
  geom_text(aes(label=round(value,3))) + 
  labs(x='Differential expression analysis',y='Differential expression analysis',"Liver vs Femur"); g
g <- g + scale_x_discrete(limits = rev(levels(cor.melt$Var1)))  # Flip the x-axis and place tick labels at the bottom
pdf(paste0("~/Documents/Research/t21-proj/out/full/pseudobulks/lfc_correlation.","All",".pdf"),width=6,height=6)
print(g)
dev.off()

pdf("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/pseudobulks/t21.pairs.pdf",width=7,height = 7)
print(pairs(res.mg3[,-1]))
dev.off()

