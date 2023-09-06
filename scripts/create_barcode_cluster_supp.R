library(data.table)
# df1 = fread("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/data/meta.10X_Healthy_Femur.umap2d.cells_removed.txt",data.table = F,stringsAsFactors = F)
# df2 = fread("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/data/meta.10X_Healthy_Liver.umap2d.cells_removed.txt",data.table = F,stringsAsFactors = F)
# df3 = fread("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/data/meta.10X_DownSyndrome_Femur.umap2d.cells_removed.txt",data.table = F,stringsAsFactors = F)
# df4 = fread("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/data/meta.10X_DownSyndrome_Liver.umap2d.cells_removed.txt",data.table = F,stringsAsFactors = F)
df1 = fread("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/cellComp/10X_Healthy_Femur.cellComp.csv",data.table = F,stringsAsFactors = F)
df2 = fread("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/cellComp/10X_Healthy_Liver.cellComp.csv",data.table = F,stringsAsFactors = F)
df3 = fread("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/cellComp/10X_DownSyndrome_Femur.cellComp.csv",data.table = F,stringsAsFactors = F)
df4 = fread("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/cellComp/10X_DownSyndrome_Liver.cellComp.csv",data.table = F,stringsAsFactors = F)

colnames(df1)[1]="cell"
colnames(df2)[1]="cell"
colnames(df3)[1]="cell"
colnames(df4)[1]="cell"
colnames(df1)[6]="cell_type"
colnames(df2)[6]="cell_type"
colnames(df3)[6]="cell_type"
colnames(df4)[6]="cell_type"

df1$organ="Femur"
df1$type = "Disomic"
df2$organ="Liver"
df2$type = "Disomic"
df3$organ="Femur"
df3$type = "Ts21"
df4$organ="Liver"
df4$type = "Ts21"

dfall = rbind(df1,df2,df3,df4)
library(dplyr)
dfall = dfall %>% select("cell","patient","sample","sorting","type","organ",everything())
f.out = "~/Downloads/cell_type_labels.supp.csv"
fwrite(dfall,f.out,quote = F,na = NA,sep = ',',row.names = F,col.names = T)


