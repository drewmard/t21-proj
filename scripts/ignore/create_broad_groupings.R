dataf.scvi = fread("~/Downloads/b3_tab.obs.withannot.csv",data.table = F,stringsAsFactors = F)
new_annot_mapping = fread("~/Downloads/new_annot_mapping.csv",data.table = F,stringsAsFactors = F)
dataf.scvi.mg = merge(data.frame(combi_annot_v3=unique(dataf.scvi$combi_annot_v3)),new_annot_mapping,by.x="combi_annot_v3",by.y="combi_annot",all.x = TRUE)
fwrite(dataf.scvi.mg,"~/Downloads/new_annot_mapping.scvi.csv",quote = F,na = NA,sep = ',',row.names = F,col.names = T)

dataf.scvi.mg = merge(data.frame(combi_annot_v3=unique(dataf.scvi$combi_annot_v3)),new_annot_mapping,by.x="combi_annot_v3",by.y="combi_annot",all.x = TRUE)
fwrite(dataf.scvi.mg,"~/Downloads/new_annot_mapping.scvi.csv",quote = F,na = NA,sep = ',',row.names = F,col.names = T)

tmp = merge(data.frame(T21_transfer=unique(refmap$T21_transfer)),df.uniq,by.x="T21_transfer",by.y="leiden_latest",all.x=TRUE)
fwrite(tmp,"~/Downloads/t21_transfer.broad.csv",quote = F,na = NA,sep = ',',row.names = F,col.names = T)
