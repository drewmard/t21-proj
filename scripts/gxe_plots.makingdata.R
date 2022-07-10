# module load R/4.1.2

library(Seurat)
library(data.table)

cell_type_filename="HSCs_MPPs"

iter=0; df <- list()
for (sampletype in c("Liver","Femur")) {
  # disease_status="Healthy"
  for (disease_status in c("Healthy","DownSyndrome")) {
    iter=iter+1
    print(iter)
    f.out <- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/cell_subset/",disease_status,'_',sampletype,'_',cell_type_filename,".rds")
    df[[iter]] <- readRDS(file = f.out)
  }
}

dfcombined = merge(df[[1]],c(df[[2]],df[[3]],df[[4]]))
f.out=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/cell_subset/",cell_type_filename,".all.rds")
saveRDS(dfcombined,f.out)

rm(df)
df <- list()
iter=5
cell_type="Cycling HSCs/MPPs";
cell_type_filename = gsub("/","_",cell_type)
disease_status="DownSyndrome"
sampletype="Liver"
f.out <- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/cell_subset/",disease_status,'_',sampletype,'_',cell_type_filename,".rds")
df[[iter]] <- readRDS(file = f.out)

dfcombined <- merge(dfcombined,df[[5]])
# dfcombined <- NormalizeData(dfcombined,normalization.method="RC",scale.factor = 1e6)

res <- fread("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DEG_list/liver_v_femur.txt",data.table = F,stringsAsFactors = F)
f.out=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/cell_subset/",cell_type_filename,".all_sig_genes.rds")
saveRDS(dfcombined[subset(res,class%in%c("environment-driven","t21-induced","t21-reverted"))$names,],f.out)

res <- fread("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DEG_list/t21_v_healthy.txt",data.table = F,stringsAsFactors = F)
f.out=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/cell_subset/",cell_type_filename,".all_sig_genes.t21.rds")
saveRDS(dfcombined[subset(res,class%in%c("environment-independent","liver-induced","femur-induced"))$names,],f.out)


res <- fread("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DEG_list/liver_v_femur.txt",data.table = F,stringsAsFactors = F)

dfcombined <- AddModuleScore(dfcombined,features=list(subset(res,class=="t21-induced" & logFC.t21 > 0)$names),ctrl=10,name="upreg_t21_induced")
dfcombined <- AddModuleScore(dfcombined,features=list(subset(res,class=="t21-induced" & logFC.t21 < 0)$names),ctrl=10,name="downreg_t21_induced")

dfcombined@meta.data$environment[1:5]
upreg_t21_induced <- apply(apply(dfcombined[["RNA"]][subset(res,class=="t21-induced" & logFC.t21 > 0)$names,],1,scale),1,median)

dfcombined@meta.data[1,]
aggregate(dfcombined$upreg_t21_induced1,by=list(dfcombined$environment,dfcombined$organ),mean)
aggregate(dfcombined$downreg_t21_induced1,by=list(dfcombined$environment,dfcombined$organ),mean)

aggregate(upreg_t21_induced,by=list(dfcombined$environment,dfcombined$organ),median)
aggregate(upreg_t21_induced,by=list(dfcombined$environment,dfcombined$organ),mean)
aggregate(upreg_t21_induced,by=list(dfcombined$environment,dfcombined$organ),median)

env <- factor(dfcombined$environment,levels = c("Healthy","Down Syndrome"))
# summary(lm(upreg_t21_induced~dfcombined$environment*dfcombined$organ))
summary(lm(dfcombined$upreg_t21_induced1~env*dfcombined$organ))
summary(lm(dfcombined$downreg_t21_induced1~env*dfcombined$organ))

summary(lm(unlist(FetchData(dfcombined,vars="GATA1"))~dfcombined$environment*dfcombined$organ))
summary(lm(FetchData(dfcombined,vars="GATA1")[,1]~dfcombined$environment*dfcombined$organ))
summary(lm(FetchData(dfcombined,vars="PVT1")[,1]~dfcombined$environment*dfcombined$organ))

upreg_t21_induced <- apply(apply(dfcombined[["RNA"]][c("PVT1","GATA1"),],1,scale),1,median)



