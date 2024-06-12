library(data.table)
df1 = fread("~/Downloads/DownSyndrome_HSC_PeakGeneSets/Liver.HSCs_MPPs.sample.txt",data.table = F,stringsAsFactors = F)
df2 = fread("~/Downloads/Liver_HSCs_MPPs.de.txt",data.table = F,stringsAsFactors = F)
df.mg = merge(df1,df2,by='names',all.x = TRUE)
df.mg = df.mg[order(df.mg$P.Value.x),]
df.mg = df.mg[,c("names","logFC.x","P.Value.x","adj.P.Val.x","logFC.y","P.Value.y","adj.P.Val.y")]
colnames(df.mg) = c("names","logFC.pre_cellbender","P.Value.pre_cellbender","adj.P.Val.pre_cellbender","logFC.post_cellbender","P.Value.post_cellbender","adj.P.Val.post_cellbender")
fwrite(df.mg,"~/Downloads/Fetal_Liver_HSC.Ts21_vs_Disomy.cellbender.txt",quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)