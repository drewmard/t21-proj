library(data.table)
# df1=fread("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/pseudobulks/sample/DE/DownSyndrome_cyc_vs_hsc.de.txt",data.table = F,stringsAsFactors = F)
# df2=fread("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/pseudobulks/sample/DE/DownSyndrome_cyc_vs_hsc.de.integ1.txt",data.table = F,stringsAsFactors = F)
# df3=fread("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/pseudobulks/sample/DE/Healthy_cyc_vs_hsc.de.integ1.txt",data.table = F,stringsAsFactors = F)
df1=fread("~/Documents/Research/t21-proj/out/full/pseudobulks/sample/DE/DownSyndrome_cyc_vs_hsc.de.txt",data.table = F,stringsAsFactors = F)
df2=fread("~/Documents/Research/t21-proj/out/full/pseudobulks/sample/DE/DownSyndrome_cyc_vs_hsc.de.integ1.txt",data.table = F,stringsAsFactors = F)
df3=fread("~/Documents/Research/t21-proj/out/full/pseudobulks/sample/DE/Healthy_cyc_vs_hsc.de.integ1.txt",data.table = F,stringsAsFactors = F)

df.mg = merge(merge(df1,df2,by="names"),df3,by='names')
cor(df.mg$logFC.x,df.mg$logFC.y)
cor(df.mg$logFC.x,df.mg$logFC)
cor(df.mg$logFC.y,df.mg$logFC)

library(data.table)
library("enrichR")
dbs <- c("GO_Molecular_Function_2021", "GO_Biological_Process_2021","ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")

prefix="original_labels"
gene_lst <- subset(df1,adj.P.Val < 0.05 & logFC > 0)$names
enriched <- enrichr(gene_lst, dbs)
y <- enriched[[dbs[1]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); y = (subset(y,geneCt>=3 & Adjusted.P.value<0.1)); y$Set <- dbs[1]; 
f.out = paste0("~/Documents/Research/t21-proj/out/full/pseudobulks/sample/DE/gsea_cyc_v_hsc/",prefix,".upreg.",dbs[1],".txt")
fwrite(y,f.out,quote = F,na = NA,sep = "\t",row.names = F,col.names = T)
y <- enriched[[dbs[2]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); y = (subset(y,geneCt>=3 & Adjusted.P.value<0.1)); y$Set <- dbs[2]; 
f.out = paste0("~/Documents/Research/t21-proj/out/full/pseudobulks/sample/DE/gsea_cyc_v_hsc/",prefix,".upreg.",dbs[2],".txt")
fwrite(y,f.out,quote = F,na = NA,sep = "\t",row.names = F,col.names = T)
y <- enriched[[dbs[3]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); y = (subset(y,geneCt>=3 & Adjusted.P.value<0.1)); y$Set <- dbs[3]; 
f.out = paste0("~/Documents/Research/t21-proj/out/full/pseudobulks/sample/DE/gsea_cyc_v_hsc/",prefix,".upreg.",dbs[3],".txt")
fwrite(y,f.out,quote = F,na = NA,sep = "\t",row.names = F,col.names = T)
gene_lst <- subset(df1,adj.P.Val < 0.05 & logFC < 0)$names
enriched <- enrichr(gene_lst, dbs)
y <- enriched[[dbs[1]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); y = (subset(y,geneCt>=3 & Adjusted.P.value<0.1)); y$Set <- dbs[1]; 
f.out = paste0("~/Documents/Research/t21-proj/out/full/pseudobulks/sample/DE/gsea_cyc_v_hsc/",prefix,".downreg.",dbs[1],".txt")
fwrite(y,f.out,quote = F,na = NA,sep = "\t",row.names = F,col.names = T)
y <- enriched[[dbs[2]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); y = (subset(y,geneCt>=3 & Adjusted.P.value<0.1)); y$Set <- dbs[2]; 
f.out = paste0("~/Documents/Research/t21-proj/out/full/pseudobulks/sample/DE/gsea_cyc_v_hsc/",prefix,".downreg.",dbs[2],".txt")
fwrite(y,f.out,quote = F,na = NA,sep = "\t",row.names = F,col.names = T)
y <- enriched[[dbs[3]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); y = (subset(y,geneCt>=3 & Adjusted.P.value<0.1)); y$Set <- dbs[3]; 
f.out = paste0("~/Documents/Research/t21-proj/out/full/pseudobulks/sample/DE/gsea_cyc_v_hsc/",prefix,".downreg.",dbs[3],".txt")
fwrite(y,f.out,quote = F,na = NA,sep = "\t",row.names = F,col.names = T)


prefix="integrated_t21"
gene_lst <- subset(df2,adj.P.Val < 0.05 & logFC > 0)$names
enriched <- enrichr(gene_lst, dbs)
y <- enriched[[dbs[1]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); y = (subset(y,geneCt>=3 & Adjusted.P.value<0.1)); y$Set <- dbs[1]; 
f.out = paste0("~/Documents/Research/t21-proj/out/full/pseudobulks/sample/DE/gsea_cyc_v_hsc/",prefix,".upreg.",dbs[1],".txt")
fwrite(y,f.out,quote = F,na = NA,sep = "\t",row.names = F,col.names = T)
y <- enriched[[dbs[2]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); y = (subset(y,geneCt>=3 & Adjusted.P.value<0.1)); y$Set <- dbs[2]; 
f.out = paste0("~/Documents/Research/t21-proj/out/full/pseudobulks/sample/DE/gsea_cyc_v_hsc/",prefix,".upreg.",dbs[2],".txt")
fwrite(y,f.out,quote = F,na = NA,sep = "\t",row.names = F,col.names = T)
y <- enriched[[dbs[3]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); y = (subset(y,geneCt>=3 & Adjusted.P.value<0.1)); y$Set <- dbs[3]; 
f.out = paste0("~/Documents/Research/t21-proj/out/full/pseudobulks/sample/DE/gsea_cyc_v_hsc/",prefix,".upreg.",dbs[3],".txt")
fwrite(y,f.out,quote = F,na = NA,sep = "\t",row.names = F,col.names = T)
gene_lst <- subset(df2,adj.P.Val < 0.05 & logFC < 0)$names
enriched <- enrichr(gene_lst, dbs)
y <- enriched[[dbs[1]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); y = (subset(y,geneCt>=3 & Adjusted.P.value<0.1)); y$Set <- dbs[1]; 
f.out = paste0("~/Documents/Research/t21-proj/out/full/pseudobulks/sample/DE/gsea_cyc_v_hsc/",prefix,".downreg.",dbs[1],".txt")
fwrite(y,f.out,quote = F,na = NA,sep = "\t",row.names = F,col.names = T)
y <- enriched[[dbs[2]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); y = (subset(y,geneCt>=3 & Adjusted.P.value<0.1)); y$Set <- dbs[2]; 
f.out = paste0("~/Documents/Research/t21-proj/out/full/pseudobulks/sample/DE/gsea_cyc_v_hsc/",prefix,".downreg.",dbs[2],".txt")
fwrite(y,f.out,quote = F,na = NA,sep = "\t",row.names = F,col.names = T)
y <- enriched[[dbs[3]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); y = (subset(y,geneCt>=3 & Adjusted.P.value<0.1)); y$Set <- dbs[3]; 
f.out = paste0("~/Documents/Research/t21-proj/out/full/pseudobulks/sample/DE/gsea_cyc_v_hsc/",prefix,".downreg.",dbs[3],".txt")
fwrite(y,f.out,quote = F,na = NA,sep = "\t",row.names = F,col.names = T)

prefix="integrated_healthy"
gene_lst <- subset(df3,adj.P.Val < 0.05 & logFC > 0)$names
enriched <- enrichr(gene_lst, dbs)
y <- enriched[[dbs[1]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); y = (subset(y,geneCt>=3 & Adjusted.P.value<0.1)); y$Set <- dbs[1]; 
f.out = paste0("~/Documents/Research/t21-proj/out/full/pseudobulks/sample/DE/gsea_cyc_v_hsc/",prefix,".upreg.",dbs[1],".txt")
fwrite(y,f.out,quote = F,na = NA,sep = "\t",row.names = F,col.names = T)
y <- enriched[[dbs[2]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); y = (subset(y,geneCt>=3 & Adjusted.P.value<0.1)); y$Set <- dbs[2]; 
f.out = paste0("~/Documents/Research/t21-proj/out/full/pseudobulks/sample/DE/gsea_cyc_v_hsc/",prefix,".upreg.",dbs[2],".txt")
fwrite(y,f.out,quote = F,na = NA,sep = "\t",row.names = F,col.names = T)
y <- enriched[[dbs[3]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); y = (subset(y,geneCt>=3 & Adjusted.P.value<0.1)); y$Set <- dbs[3]; 
f.out = paste0("~/Documents/Research/t21-proj/out/full/pseudobulks/sample/DE/gsea_cyc_v_hsc/",prefix,".upreg.",dbs[3],".txt")
fwrite(y,f.out,quote = F,na = NA,sep = "\t",row.names = F,col.names = T)
gene_lst <- subset(df3,adj.P.Val < 0.05 & logFC < 0)$names
enriched <- enrichr(gene_lst, dbs)
y <- enriched[[dbs[1]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); y = (subset(y,geneCt>=3 & Adjusted.P.value<0.1)); y$Set <- dbs[1]; 
f.out = paste0("~/Documents/Research/t21-proj/out/full/pseudobulks/sample/DE/gsea_cyc_v_hsc/",prefix,".downreg.",dbs[1],".txt")
fwrite(y,f.out,quote = F,na = NA,sep = "\t",row.names = F,col.names = T)
y <- enriched[[dbs[2]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); y = (subset(y,geneCt>=3 & Adjusted.P.value<0.1)); y$Set <- dbs[2]; 
f.out = paste0("~/Documents/Research/t21-proj/out/full/pseudobulks/sample/DE/gsea_cyc_v_hsc/",prefix,".downreg.",dbs[2],".txt")
fwrite(y,f.out,quote = F,na = NA,sep = "\t",row.names = F,col.names = T)
y <- enriched[[dbs[3]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); y = (subset(y,geneCt>=3 & Adjusted.P.value<0.1)); y$Set <- dbs[3]; 
f.out = paste0("~/Documents/Research/t21-proj/out/full/pseudobulks/sample/DE/gsea_cyc_v_hsc/",prefix,".downreg.",dbs[3],".txt")
fwrite(y,f.out,quote = F,na = NA,sep = "\t",row.names = F,col.names = T)


ykeep <- y
y <- enriched[[dbs[2]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); ytmp = (subset(y,geneCt>=3 & Adjusted.P.value<0.1)); ytmp$Set <- dbs[2]; ykeep<-rbind(ykeep,ytmp)
y <- enriched[[dbs[3]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); ytmp = (subset(y,geneCt>=3 & Adjusted.P.value<0.1)); ytmp$Set <- dbs[3]; ykeep<-rbind(ykeep,ytmp)


# ggplot(ykeep,aes(x=Set,
#                   y=Term,
#                   fill=-log10(Adjusted.P.value))) +
#   geom_tile(col='black') +
#   theme_bw()  +
#   theme(panel.grid = element_blank(),
#         axis.text.x = element_text(angle = 60, hjust=1), #,colour = axis_colors),
#         panel.border = element_blank(),
#         legend.position = 'none'
#   ) +
#   scale_fill_gradient2(low = "blue", mid='white',high = "red",name='Power') +
#   labs(x='Set',y='Ter')


# ABHD16B for 15582E cycling has 0 counts in new, but 136 counts in old
# df.aggre["ABHD16B","15582E"]
# 
# f.out=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data_pb_leiden/Liver.pb.",cell_type_filename,".txt")
# tmp = fread(f.out,data.table = F,stringsAsFactors = F)
# subset(tmp,tmp[,1]=="ABHD16B")[,"15582E.cyc"]
# 
# df.aggre1["ABHD16B","15582E.cyc"]
# subset(df.aggre1,df.aggre1[,1]=="ABHD16B")[,"15582E.cyc"]
# df.aggre["ABHD16B","15582E"]


# it turns out, the error was in merging two cells of different types which required a new command in my script
# which had a bug!! and threw off the gene expression counts for > 90% of genes...
