# ST
library(data.table)
df1=fread("~/Documents/Research/t21-proj/out/full/cpdb_Results/DownSyndrome_Liver/DownSyndrome_Liver.processed_output.txt",data.table = F,stringsAsFactors = F)
df2=fread("~/Documents/Research/t21-proj/out/full/cpdb_Results/Healthy_Liver/Healthy_Liver.processed_output.txt",data.table = F,stringsAsFactors = F)
df1 <- df1[df1$Receiver!="X",] ; df2 <- df2[df2$Receiver!="X",] 
df1.sub = subset(df1,(Sender=="HSCs/MPPs" | Receiver=="HSCs/MPPs") & P < 0.001)
df2.sub = subset(df2,(Sender=="HSCs/MPPs" | Receiver=="HSCs/MPPs") & P < 0.001)
df1.sub$Data <- "T21"
df2.sub$Data <- "H"
df.sub <- rbind(df1.sub,df2.sub)
df.sub$Other <- ifelse(df.sub$Sender=="HSCs/MPPs",df.sub$Receiver,df.sub$Sender)
df.sub$All = paste(df.sub$Ligand,df.sub$Receptor,df.sub$Sender,df.sub$Receiver)
df.sub <- subset(df.sub,!(All %in% df.sub$All[duplicated(df.sub$All)]))

DEgenes = fread("~/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/Liver.HSCs_MPPs.sample.permALL.p.txt",data.table = F,stringsAsFactors = F)
DEgenes.sub = DEgenes$names[DEgenes$P.Value.liv<0.05]
DE = lapply(1:nrow(df.sub),function(i) {ifelse(df.sub$Sender[i]=="HSCs/MPPs",df.sub$Ligand[i] %in% DEgenes.sub,df.sub$Receptor[i] %in% DEgenes.sub)})
df.sub$DE = unlist(DE)

tab <- aggregate(df.sub$Sender,by=list(Cell=df.sub$Other,Data=df.sub$Data),length)
df.mg1=merge(
  subset(tab,Data=="H")[,-2],
  subset(tab,Data=="T21")[,-2],
  by="Cell",
  all=T
)
tab <- aggregate(df.sub$Sender,by=list(Cell=df.sub$Other,Data=df.sub$Data,DE=df.sub$DE),length)
tab <- subset(tab,DE)
df.mg2=merge(
  subset(tab,Data=="H" & DE)[,-c(2:3)],
  subset(tab,Data=="T21" & DE)[,-c(2:3)],
  by="Cell",
  all=T
)
df.mg=merge(df.mg1,df.mg2,by="Cell")
subset(df.sub,DE)
df.sub$Pair=paste0(df.sub$Ligand,"_",df.sub$Receptor)

df.sub2 = subset(df.sub,DE & Sender=="HSCs/MPPs" & Data=="H")
library(ggplot2)
ggplot(df.sub2,aes(x=Receiver,
                  y=Pair,
                  fill=MeanExpr)) +
  geom_tile(col='black') +
  theme_bw()  +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 60, hjust=1), #,colour = axis_colors),
        panel.border = element_blank(),
        legend.position = 'none'
  ) +
  scale_fill_gradient2(low = "blue", mid='white',high = "red",name='Power') +
  labs(x='Receiver',y='L-R Pair',title="Healthy")

df.sub2 = subset(df.sub,DE & Sender=="HSCs/MPPs" & Data=="T21")
library(ggplot2)
ggplot(df.sub2,aes(x=Receiver,
                   y=Pair,
                   fill=MeanExpr)) +
  geom_tile(col='black') +
  theme_bw()  +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 60, hjust=1), #,colour = axis_colors),
        panel.border = element_blank(),
        legend.position = 'none'
  ) +
  scale_fill_gradient2(low = "blue", mid='white',high = "red",name='Power') +
  labs(x='Receiver',y='L-R Pair',title="T21")

df.sub2 = subset(df.sub,DE & Receiver=="HSCs/MPPs" & Data=="H")
library(ggplot2)
ggplot(df.sub2,aes(x=Sender,
                   y=Pair,
                   fill=MeanExpr)) +
  geom_tile(col='black') +
  theme_bw()  +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 60, hjust=1), #,colour = axis_colors),
        panel.border = element_blank(),
        legend.position = 'none'
  ) +
  scale_fill_gradient2(low = "blue", mid='white',high = "red",name='Power') +
  labs(x='Sender',y='L-R Pair',title="Healthy")

df.sub2 = subset(df.sub,DE & Receiver=="HSCs/MPPs" & Data=="T21")
library(ggplot2)
ggplot(df.sub2,aes(x=Sender,
                   y=Pair,
                   fill=MeanExpr)) +
  geom_tile(col='black') +
  theme_bw()  +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 60, hjust=1), #,colour = axis_colors),
        panel.border = element_blank(),
        legend.position = 'none',
        
  ) +
  scale_fill_gradient2(low = "blue", mid='white',high = "red",name='Power') +
  labs(x='Sender',y='L-R Pair',title="T21") 


ggplot(aes(x=Receiver,y=Pair,fill=P)) 


merge(
  ,
  merge(
    subset(tab,Data=="H" & DE)[,-2],
    subset(tab,Data=="T21" & DE)[,-2],
    by="Cell",
    all=T
  ),
  by="Cell"
)

subset(df.sub,Other=="Early erythroid cells")


unique(c(df1.sub$Sender,df1.sub$Receiver)) 
unique(c(df1$Sender,df1$Receiver))
unique(c(df2.sub$Sender,df2.sub$Receiver))
unique(c(df2$Sender,df2$Receiver))


