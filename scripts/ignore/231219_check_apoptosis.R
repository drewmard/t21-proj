library(data.table)
f = "~/Downloads/DownSyndrome_HSC_PeakGeneSets/Liver.HSCs_MPPs.sample.txt"
# f = "~/Downloads/DownSyndrome_HSC_PeakGeneSets/liver_v_femur.txt"
df = fread(f,data.table = F,stringsAsFactors = F)

gene_lst=c("BAX","BAK","BID","CASP3",'CASP8','CASP9','CYCS',"APAF1") # pro apoptosis
subset(df,names %in% gene_lst)

gene_lst=c("BCL2","BCL2L1","MCL1","XIAP",'BIRC5','FLIP',"TNFRSF1A") # anti apoptosis
subset(df,names %in% gene_lst)

gene_set = fread("~/Downloads/Apoptosis_from_KEGG_2021_Human.gmt",data.table = F,stringsAsFactors = F)
gene_set = unname(unlist(gene_set[,3:ncol(gene_set)]))
gene_set

df$apoptosis = df$names %in% gene_set
fisher.test(df$adj.P.Val < 0.05,df$apoptosis)
fisher.test(df$logFC < 0 & df$adj.P.Val < 0.1,df$apoptosis)

subset(df,df$adj.P.Val < 0.1 & df$logFC > 0 & df$apoptosis)
subset(df,df$adj.P.Val < 0.1 & df$logFC < 0 & df$apoptosis)


fisher.test(df$class == "t21-induced",df$apoptosis)
fisher.test(df$class %in% c("t21-induced","unknown","environment-driven","t21-reverted") & df$logFC.t21 > 0,df$apoptosis)
fisher.test(df$class %in% c("t21-induced","unknown","environment-driven") & df$logFC.t21 > 0,df$apoptosis)
fisher.test(df$class %in% c("t21-induced","unknown") & df$logFC.t21 > 0,df$apoptosis)
fisher.test(df$class %in% c("t21-induced") & df$logFC.t21 > 0,df$apoptosis)

subset(df,df$class %in% c("t21-induced") & df$logFC.t21 > 0 & df$apoptosis)

fisher.test(df$logFC < 0 & df$adj.P.Val < 0.1,df$apoptosis)



