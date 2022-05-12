library(data.table)
library(biomaRt)

cell_type_lst <- c("B cells","Erythroid","HSC_Progenitors","Mast cells","Megakaryocytes","Myeloid","NK_T cells")
iter=0; for (cell_type in cell_type_lst) {
  iter=iter+1
  f=paste0("~/Documents/Research/t21-proj/out/full/DE_cell_type_groups/Liver.",cell_type,".txt")
  df <- fread(f,data.table = F,stringsAsFactors = F)
  df$logfoldchanges <- rank(df$logfoldchanges)/nrow(df)
  x <- df[,c('names','logfoldchanges')]
  colnames(x)[2] <- cell_type
  if (iter==1) {
    df.save <- x
  } else {
    df.save <- merge(df.save,x,by='names')
  }
}

cor(df.save[,-1])

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
annot <- getBM(attributes = c('hgnc_symbol', 'chromosome_name',
                              'start_position', 'end_position'),
               filters = 'hgnc_symbol', 
               values = df.save$names, 
               mart = ensembl)
df1 <- merge(df.save,annot,by.x="names",by.y="hgnc_symbol")
df1 <- df1[df1$chromosome_name %in% seq(1,22),]
df1$chromosome_name <- factor(df1$chromosome_name,levels=seq(1,22))
df1$chr21 <- factor(ifelse(df1$chromosome_name==21,'Chr 21','Not Chr 21'),levels=c('Not Chr 21','Chr 21'))

df1.melt <- melt(df1[,c(cell_type_lst,"chr21")],id.vars="chr21")
ggplot(df1.melt,aes(x=variable,y=value,fill=chr21)) + geom_boxplot() + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle=60,hjust=1)) +
  labs(x="Chromosome",y='Log Fold Change Percentile (T21 vs Healthy)') +
  scale_fill_brewer(palette = "Set2") #+ guides(fill=F)



