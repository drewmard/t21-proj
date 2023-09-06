library(data.table)
# f="~/Documents/Research/t21-proj/out/full/data/meta.10X_Healthy_Liver.umap2d.cells_removed.txt"
# sep=fread(f,data.table = F,stringsAsFactors = F)
# sep$cell
# sep[,c()]
f="~/Downloads/a3_tab.obs.withannot.csv"
integrated_labels=fread(f,data.table = F,stringsAsFactors = F)
integrated_labels$cellname=unlist(lapply(strsplit(integrated_labels$barcodekey,"-"),function(x) paste(x[1:2],collapse = "-")))
integrated_labels$sample=unlist(lapply(strsplit(as.character(integrated_labels$patient_sample)," "),function(x) paste(x[2],collapse = "-")))
integrated_labels = integrated_labels[!duplicated(integrated_labels[,c("cellname","sample")]),] # why??
df.mg = subset(integrated_labels,environment=="Healthy_Liver")[,c("cellname","sample","leiden_latest","combi_annot")]
colnames(df.mg)[3:4] = c("Separate","Integrated")
f="~/Downloads/scArches_Healthy_Liver_metadata.csv"
popescu=fread(f,data.table = F,stringsAsFactors = F)
popescu$cellname=unlist(lapply(strsplit(popescu$barcode,"-"),function(x) paste(x[1:2],collapse = "-")))
popescu = popescu[!duplicated(popescu[,c("cellname","sample","Predicted")]),] # why??
popescu = popescu[,c("cellname","sample","Predicted")]
colnames(popescu)[3] = "Popescu_transfer"
f="~/Downloads/Ts21_reference_scArches_Healthy_Liver_metadata.csv"
t21transfer=fread(f,data.table = F,stringsAsFactors = F)
t21transfer$cellname=unlist(lapply(strsplit(t21transfer$barcode,"-"),function(x) paste(x[1:2],collapse = "-")))
t21transfer = t21transfer[!duplicated(t21transfer[,c("cellname","sample","Predicted")]),] # why??
t21transfer = t21transfer[,c("cellname","sample","Predicted")]
colnames(t21transfer)[3] = "T21_transfer"

df = merge(merge(df.mg,popescu,by=c("sample","cellname")),t21transfer,by=c("sample","cellname"))
table(subset(df,Integrated=="Cycling HSC")$Separate)
table(subset(df,T21_transfer=="Cycling HSCs/MPPs")$Separate)

table(subset(df,Integrated=="Cycling HSC")$T21_transfer)
table(subset(df,T21_transfer=="Cycling HSCs/MPPs")$Integrated)

table(subset(df,T21_transfer=="Cycling HSCs/MPPs" & Integrated=="Cycling HSC")$Separate)


table(subset(df,Integrated=="Cycling HSC")$Separate)

fwrite(df,"~/Downloads/Healthy_liver_annotations.txt",quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
tab = aggregate(df$Separate,by=df[,3:6],length)
tab = tab[order(tab$x,decreasing = T),]
fwrite(tab,"~/Downloads/Healthy_liver_annotations.tab.txt",quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

# subset(df,T21_transfer=="Cycling HSCs/MPPs")


