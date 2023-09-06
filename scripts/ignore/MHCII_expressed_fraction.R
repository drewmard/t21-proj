dflst = list()

gene_name_lst <- c("HLA-DMA","HLA-DPA1","HLA-DRA","HLA-DRB1")
for (i in 1:length(gene_name_lst)) {
  print(i)
  gene_name = gene_name_lst[i]
  x1=mean(dfcombined[["RNA"]][gene_name,dfcombined@meta.data$organ=="Femur" & dfcombined@meta.data$environment=="Healthy"]>0)
  x2=mean(dfcombined[["RNA"]][gene_name,dfcombined@meta.data$organ=="Femur" & dfcombined@meta.data$environment=="Down Syndrome"]>0)
  x3=mean(dfcombined[["RNA"]][gene_name,dfcombined@meta.data$organ=="Liver" & dfcombined@meta.data$environment=="Healthy"]>0)
  x4=mean(dfcombined[["RNA"]][gene_name,dfcombined@meta.data$organ=="Liver" & dfcombined@meta.data$environment=="Down Syndrome"]>0)
  dflst[[i]] = data.frame(gene=gene_name,femur.h=x1,femur.t21=x2,liver.h=x3,liver.t21=x4)
}

do.call(rbind,dflst)
