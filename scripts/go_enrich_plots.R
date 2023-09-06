library(data.table)
library("enrichR")
dbs <- c("GO_Molecular_Function_2021", "GO_Biological_Process_2021","ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")

df1.full.mg = fread("~/Documents/Research/t21-proj/out/full/DEG_list/liver_v_femur.txt",data.table = F,stringsAsFactors = F)

gene_lst <- subset(df1.full.mg,class=="t21-induced" & logFC.t21 > 0)$names
enriched <- enrichr(gene_lst, dbs)
y <- enriched[[dbs[1]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); y = (subset(y,geneCt>=3 & Adjusted.P.value<0.1)); y$Set <- dbs[1]; ykeep <- y
y <- enriched[[dbs[2]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); ytmp = (subset(y,geneCt>=3 & Adjusted.P.value<0.1)); ytmp$Set <- dbs[2]; ykeep<-rbind(ykeep,ytmp)
y <- enriched[[dbs[3]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); ytmp = (subset(y,geneCt>=3 & Adjusted.P.value<0.1)); ytmp$Set <- dbs[3]; ykeep<-rbind(ykeep,ytmp)
ykeep[,1]
# ggplot(ykeep,aes(x=Set,
#                  y=reorder(Term,as.numeric(as.factor(Set))),
ggplot(ykeep2,aes(x=Set,
                 y=Term,
                 fill=-log10(Adjusted.P.value))) +
  geom_tile(col='black') +
  theme_bw()  +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 60, hjust=1), #,colour = axis_colors),
        panel.border = element_blank(),
        legend.position = 'none'
  ) +
  scale_fill_gradient2(low = "blue", mid='white',high = "red",name='Power') +
  labs(x='Set',y='Ter')

require(dplyr)
require(forcats)
ykeep2 = ykeep
ykeep2 = ykeep2 %>% 
  mutate(ordering = -as.numeric(as.factor(Set)) + Adjusted.P.value,
         Term = fct_reorder(Term, ordering, .desc = T))
  
ggplot(ykeep2,aes(fill=Set,
                 x=Term,
                 # x=reorder(Term,as.numeric(as.factor(Set))),
                 y=-log10(Adjusted.P.value))) +
  geom_bar(stat='identity') +
  theme_bw()  +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 60, hjust=1), #,colour = axis_colors),
        panel.border = element_blank()
  ) +
  labs(x='Term',y='-log10 FDR')  + scale_fill_brewer(palette="Set2")


ggplot(ykeep,aes(x=Set,y=-log10(Adjusted.P.value),col=Set)) +
  geom_jitter(width=0.1) + 
  theme_bw()  +
  theme(panel.grid = element_blank(),
        # axis.text.x = element_text(angle = 60, hjust=1,color='black'), #,colour = axis_colors),
        axis.text=element_text(color='black'),
        panel.border = element_rect(colour = "black", fill=NA, size=1.5),
        legend.position = 'none'
  ) +
  labs(x='Set',y='-log10 FDR') + scale_color_brewer(palette="Set2") +
  coord_flip()

ykeep$Set
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

print(g)




