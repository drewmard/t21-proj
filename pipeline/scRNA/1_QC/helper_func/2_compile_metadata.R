library(data.table)
metadata_names <- fread("~/Downloads/patientMetadata/metadata.csv",data.table = F,stringsAsFactors = F)
metalst <- list()
i=1
for (i in 1:nrow(metadata_names)) {
  metasub = metadata_names[i,]
  f = paste0("~/Downloads/patientMetadata/",metasub[1],"_",metasub[2],"_metadata.csv")
  metalst[[i]] <- fread(f,data.table = F,stringsAsFactors = F)
  metalst[[i]][,"ID1"] <- metasub[1]
  metalst[[i]][,"ID2"] <- metasub[2]
}
metadf = as.data.frame(do.call(rbind,metalst))
fwrite(metadf,"~/Downloads/patientMetadata/scRNA.metadata.txt",row.names = F,col.names = T,quote = F,na = "NA",sep = '\t')

head(metadf)
sum(duplicated(metadf$ID2))

fwrite(as.data.frame(metadf$ID2),"~/Downloads/patientMetadata/samp_ids",row.names = F,col.names = F,quote = F,na = "NA",sep = '\t')
fwrite(as.data.frame(sapply(as.data.frame(metadf$sample),function(x) paste0("/oak/stanford/groups/smontgom/amarder/data/t21/Cellranger/",x))),"~/Downloads/patientMetadata/samp_lst",row.names = F,col.names = F,quote = F,na = "NA",sep = '\t')


