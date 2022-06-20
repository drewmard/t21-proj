library(data.table)

cell_type="HSCs/MPPs"
cell_type_filename = gsub("/","_",cell_type)

disease_status="DownSyndrome"
f=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_liver_v_femur/",cell_type_filename,'.',disease_status,'.FE.txt')
df <- fread(f,data.table = F,stringsAsFactors = F)
print(aggregate(df$adj.P.Val < 0.1,by=list(df$chromosome_name==21),mean))
df[order(df$adj.P.Val)[1:10],]
subset(df,names=="MYL4")
subset(df,names=="APOC1")
subset(df,names=="GATA1")

disease_status="Healthy"
f=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_liver_v_femur/",cell_type_filename,'.',disease_status,'.FE.txt')
df <- fread(f,data.table = F,stringsAsFactors = F)
print(aggregate(df$adj.P.Val < 0.1,by=list(df$chromosome_name==21),mean))
df[order(df$adj.P.Val)[1:2],]
subset(df,names=="GATA1")

dir="/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_liver_v_femur/"
flist <- list.files(dir)
flist <- flist[grep("FE",flist)]
df.all <- list()
for (i in 1:length(flist)) {
  f=flist[i]
  df=fread(paste0(dir,f),data.table = F,stringsAsFactors = F)
  # df.all[[f]] = df[df$names=="APOC1",]
  df.all[[f]] = df[df$names=="GATA1",]
  # print(f)
  # print(df[df$names=="GATA1",])
  # print("")
}
do.call(rbind,df.all)[,c("AveExpr","logFC","P.Value","adj.P.Val")]
