flst = list.files("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names",pattern = "sample.txt")
i=1
y = strsplit(flst,"\\.")

dflst = list()
for (i in 1:length(flst)) {
  f = flst[i]
  organ = y[[i]][1]
  cell = y[[i]][2]
  
  df=fread(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/",f),data.table = F,stringsAsFactors = F)
  df$organ = organ
  df$cell = cell
  
  dflst[[i]] = df
}

dfall = as.data.frame(do.call(rbind,dflst))
dfall$analysis = "Ts21 vs Disomic"
library(dplyr)
dfall = dfall %>% select(analysis,organ,cell,names,everything())
f.out = "~/Downloads/Ts21_vs_Disomic_DE.supp.csv"
fwrite(dfall,f.out,quote = F,na = "NA",sep = ',',row.names = F,col.names = T)


