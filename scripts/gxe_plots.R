library(Seurat)
library(stringr)
library(data.table)
library(ggplot2)
library(cowplot)
source("~/Documents/Research/Useful_scripts/rntransform.R")

cell_type_filename="HSCs_MPPs"
f.out=paste0("~/Documents/Research/t21-proj/out/full/data/cell_subset/",cell_type_filename,".all_sig_genes.rds")
res <- fread("~/Documents/Research/t21-proj/out/full/DEG_list/liver_v_femur.txt",data.table = F,stringsAsFactors = F)
dfcombined=readRDS(f.out)

# res <- fread("~/Documents/Research/t21-proj/out/full/DEG_list/t21_v_healthy.txt",data.table = F,stringsAsFactors = F)
# f.out=paste0("~/Documents/Research/t21-proj/out/full/data/cell_subset/",cell_type_filename,".all_sig_genes.t21.rds")
# dfcombined=readRDS(f.out)

upreg_t21_induced <- apply(apply(dfcombined[["RNA"]][subset(res,class=="t21-induced" & logFC.t21 > 0)$names,],1,scale),1,mean)
aggregate(upreg_t21_induced,by=list(dfcombined$environment,dfcombined$organ),mean)
aggregate(upreg_t21_induced,by=list(dfcombined$environment,dfcombined$organ),median)

env <- factor(dfcombined$environment,levels = c("Healthy","Down Syndrome"))
summary(lm(upreg_t21_induced~env*dfcombined$organ))


# y <- apply(apply(dfcombined[["RNA"]][subset(res,class=="t21-induced" & logFC.t21 > 0)$names,],1,rntransform),1,mean)
# input=subset(res,class=="t21-induced" & logFC.t21 > 0)$names

glst = list()
df.mg = list()
env <- factor(dfcombined$environment,levels = c("Healthy","Down Syndrome"))
# for (class_to_use in c("environment-independent","liver-induced")) {
for (class_to_use in c("environment-driven","t21-induced")) {
  
  # input=subset(res,class==class_to_use & logFC.t21 < 0)$names
  # input=subset(res,class==class_to_use & logFC.liv > 0)$names
  input=subset(res,class==class_to_use & logFC.t21 > 0)$names
  
  
  # input=subset(res,class=="environment-driven" & logFC.t21 > 0)$names
  # input=subset(res,class==class_to_use & logFC.t21 < 0)$names
  # input=subset(res,class==class_to_use & logFC.t21 < 0)$names
  y <- apply(apply(dfcombined[["RNA"]][input,],1,scale),1,mean)
  # y <- apply(apply(dfcombined[["RNA"]]["GATA1",],1,scale),1,mean)
  
  # y <- apply(apply(dfcombined[["RNA"]][subset(res,class=="environment-driven" & logFC.t21 < 0)$names,],1,scale),1,mean)
  # aggregate(y,by=list(dfcombined$environment,dfcombined$organ),mean)
  # aggregate(y,by=list(dfcombined$environment,dfcombined$organ),median)
  y <- rntransform(y)
  summary(lm(y~env*organ,dftmp <- data.frame(y=y,env=env,organ=dfcombined$organ)))
  
  # ggplot(dftmp,aes(x=env,y=rntransform(y),col=organ)) + geom_boxplot()
  # 
  df.mg[[class_to_use]] = merge(aggregate(y,by=list(dfcombined$environment,dfcombined$organ),function(x) sd(x)/sqrt(length(input))),
                aggregate(y,by=list(dfcombined$environment,dfcombined$organ),mean),
                by=c("Group.1","Group.2"))
  
  glst[[class_to_use]] <- ggplot(df.mg[[class_to_use]],aes(x=as.numeric(as.factor(Group.2)),
                   y=x.y,
                   col=Group.1,
                   ymin=x.y-2*x.x,
                   ymax=x.y+2*x.x)) +
    geom_line() + 
    # geom_pointrange(position=position_dodge(0.05)) +
    # geom_linerange(position=position_dodge(0.05)) +
    # geom_pointrange(stat="summary",fun.ymin=df.mg[[class_to_use]]$x.y-2*df.mg[[class_to_use]]$x.x,fun.ymax=df.mg[[class_to_use]]$x.y+2*df.mg[[class_to_use]]$x.x,fun.y=median) + 
    geom_point(size=rel(2)) +
    geom_errorbar(aes(ymin=x.y-2*x.x,ymax=x.y+2*x.x),width=.2) + #,position=position_dodge(0.05)) +
    theme_bw() + theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5)) +
      labs(x="Environment",y="Average normalized expression (+/- 2 SE)",col="Disease status",title=str_to_title(class_to_use)) +
    scale_x_continuous(labels=c("Femur","Liver"),breaks=c(1,2)) +
    scale_color_brewer(palette = "Set1")

  # glst[[class_to_use]] <- ggplot(dftmp,aes(x=organ,y=rntransform(y),col=env)) + geom_boxplot() +
  #   theme_bw() + theme(panel.grid = element_blank()) +
  #   labs(x="Environment",y="Normalized expression",col="Disease status",title=str_to_title(class_to_use))
  
}
plot_grid(plotlist = glst)

# ggplot(aggregate(y,by=list(dfcombined$environment,dfcombined$organ),mean),
#        aes(x=as.numeric(as.factor(Group.2)),y=x,col=Group.1)) + geom_line()

