library(data.table)
library(ggplot2)

# df <- fread("~/Documents/Research/t21-proj/out/full/DE_int/MEMPs.gene_CEP55.txt",data.table = F,stringsAsFactors = F)
# df <- fread("~/Documents/Research/t21-proj/out/full/DE_int/Late erythroid cells.gene_SRSF8.txt",data.table = F,stringsAsFactors = F)
df <- fread("~/Documents/Research/t21-proj/out/full/DE_int/Pre pro B cells.gene_BIRC5.txt",data.table = F,stringsAsFactors = F)
# df <- fread("~/Documents/Research/t21-proj/out/full/DE_int/Pre pro B cells.gene_CCS.txt",data.table = F,stringsAsFactors = F)

df$environment <- factor(df$environment,levels=c("Healthy","DownSyndrome"))
# ggplot(df,aes(x=sampletype,fill=environment,y=geneExpr)) + geom_boxplot() + theme_bw()
ggplot(df,aes(x=environment,y=geneExpr,fill=environment)) + geom_boxplot() + theme_bw() + facet_wrap(~sampletype) + scale_fill_brewer(palette = "Set3")

 # table(df$environment,df$sampletype)

