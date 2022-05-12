library(data.table)
library(ggplot2)

# df <- fread("~/Documents/Research/t21-proj/out/full/DE_cell_type_groups/HSC_Progenitors.txt",data.table = F,stringsAsFactors = F)
# df <- fread("~/Documents/Research/t21-proj/out/full/DE_cell_type_groups/Erythroid.txt",data.table = F,stringsAsFactors = F)
# df <- fread("~/Documents/Research/t21-proj/out/full/DE_leiden_names/Early erythroid cells.txt",data.table = F,stringsAsFactors = F)
df <- fread("~/Documents/Research/t21-proj/out/full/DE_cell_type_groups/Liver.Erythroid.txt",data.table = F,stringsAsFactors = F)
df_pb <- fread("~/Documents/Research/t21-proj/out/full/DE_pb_cell_type_groups/Liver.Erythroid.txt",data.table = F,stringsAsFactors = F)
df.mg <- merge(df,df_pb,by.x="names",by.y="gene")
cor(df.mg$logfoldchanges,df.mg$logFC)
cor(-log10(df.mg$pvals+1e-216),-log10(df.mg$PValue+1e-216),use="complete.obs")

gene_info <- fread("~/Documents/Research/data/grch38/protein_coding_genes.txt",data.table = F,stringsAsFactors = F)

df <- merge(df.mg,gene_info[,c("Chromosome/scaffold name","Gene name")],by.x="names",by.y="Gene name")

df <- df[df$`Chromosome/scaffold name` %in% seq(1,22),]
df$`Chromosome/scaffold name` <- factor(df$`Chromosome/scaffold name`,levels=seq(1,22))
# ggplot(df,aes(x=`Chromosome/scaffold name`,y=logfoldchanges)) + geom_boxplot()

# aggregate(df$logfoldchanges,by=list(df$`Chromosome/scaffold name`),median)

# ggplot(df,aes(x=`Chromosome/scaffold name`==21,y=logfoldchanges)) + geom_boxplot()

source("~/Documents/Research/Useful_scripts/rntransform.R")
aggregate(df$logfoldchanges,by=list(df$`Chromosome/scaffold name`==21),median)

df$chr21 <- factor(ifelse(df$`Chromosome/scaffold name`==21,'Chr 21','Not Chr 21'),levels=c('Not Chr 21','Chr 21'))
ggplot(df,aes(x=chr21,y=rank(logfoldchanges)/nrow(df),fill=chr21)) + 
  geom_boxplot() + 
  theme_bw() +
  labs(x="Chromosome",y='Log Fold Change Percentile (T21 vs Healthy)') +
  scale_fill_brewer(palette = "Set2") + guides(fill=F)
ggplot(df,aes(x=chr21,y=rank(logFC)/nrow(df),fill=chr21)) + 
  geom_boxplot() + 
  theme_bw() +
  labs(x="Chromosome",y='Log Fold Change Percentile (T21 vs Healthy)') +
  scale_fill_brewer(palette = "Set2") + guides(fill=F)

ggplot(df,aes(x=`Chromosome/scaffold name`,y=logfoldchanges)) + geom_boxplot()
ggplot(df,aes(x=`Chromosome/scaffold name`,y=logFC)) + geom_boxplot()

aggregate(df$pvals_adj<0.01,by=list(df$`Chromosome/scaffold name`),mean)



