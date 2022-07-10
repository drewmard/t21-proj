library(data.table)
library(ggplot2)
library(cowplot)
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


df.melt <- reshape2::melt(df.sub[,c("G1","cell","disease_status")],id.vars=c("disease_status","cell"))
reshape2::dcast(G1~cell+disease_status,df.sub[,c("G1","cell","disease_status")])
rownames(df.melt) <- df.melt$cell
