library(data.table)
library(ggplot2)
fileDir="/oak/stanford/groups/smontgom/amarder/t21_multiome/output/finemap/finemappedtraits_hg19/hg38"
flst = list.files(fileDir)
flst = list.files(fileDir,pattern = "hg38",include.dirs = FALSE)

DeviationRegressionModel <- function(X,Y) {
  X <- as.factor(X)
  Y.i <- tapply(Y, X, median,na.rm=T)
  Z.ij <- abs(Y - Y.i[X])
  res <- summary(lm(Z.ij~X))$coef[2,]
  return(res)
}

traitName = "mono"
f = paste0("~/Documents/Research/t21_multiome/output/scavenge/",traitName,".txt")
trait_mat.save = fread(f,data.table = F,stringsAsFactors = F)
celltypes_all = c("B cells","Early erythroid","Granulocyte progenitors","HSCs","Late erythroid","Mast cells","Megakaryocytes","MEMPs","Neutrophils","NK cells","pDCs","Pro-inflammatory macrophages","Stroma")
res = list()
for (j in (1:length(celltypes_all))) {
  celltype = celltypes_all[j]
  trait_mat.sub = trait_mat.save[trait_mat.save$subclust_v6==celltype,]
  DRM_res = DeviationRegressionModel(X=trait_mat.sub$disease,Y=trait_mat.sub$TRS)
  x1 = trait_mat.sub$TRS[trait_mat.sub$disease=="H"]
  x2 = trait_mat.sub$TRS[trait_mat.sub$disease=="T21"]
  pval_mwu = wilcox.test(x1,x2)$p.value
  pval_t = t.test(x1,x2)$p.value
  res[[j]] = data.frame(cell=celltype,n1=length(x1),n2=length(x2),med1=median(x1),med2=median(x2),mu1=mean(x1),mu2=mean(x2),pval_mwu,pval_t,beta_drm=DRM_res[1],se_drm=DRM_res[2],z_drm=DRM_res[3],pval_drm=DRM_res[4])
}
res.df = as.data.frame(do.call(rbind,res))
res.df$theta = res.df$z_drm/(res.df$med2-res.df$med1)
rownames(res.df)<- NULL
res.df

