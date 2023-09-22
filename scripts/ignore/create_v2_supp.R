library(data.table)
df=fread("~/Downloads/cluster_markers.supp.csv",data.table = F,stringsAsFactors = F)
df = subset(df,Rank <= 100)
fwrite(df,"~/Downloads/cluster_markers.supp.v2.csv",quote = F,na = NA,sep = ',',row.names = F,col.names = T)

library(data.table)
df=fread("~/Downloads/Ts21_vs_Disomic_DE.supp.csv",data.table = F,stringsAsFactors = F)
df$Rank = NA
y = paste(df$organ,df$cell)
iter = 0
for (analy in unique(y)) {
  iter = iter+1; print(iter)
  ind = y==analy
  df$Rank[ind] = 1:sum(ind)
}
df.sub = subset(df,Rank <= 100)
fwrite(df.sub,"~/Downloads/Ts21_vs_Disomic_DE.supp.v2.csv",quote = F,na = NA,sep = ',',row.names = F,col.names = T)
