library(data.table)

# metadata
meta = fread("/oak/stanford/groups/smontgom/amarder/data/hasaart_et_al_2020_scirep/metadata.csv",data.table = F,stringsAsFactors = F)

# for loop around all files
headdir="/oak/stanford/groups/smontgom/amarder/data/hasaart_et_al_2020_scirep/hg38"
files = list.files(path=headdir,full.names = TRUE)
files_without_string = files[!grepl("unmapp",files)]

iter=0; df=list()
for (inFile in files_without_string) {
  
  iter = iter + 1
  
  # read in mutation data
  tmp=fread(inFile,data.table = F,stringsAsFactors = F)
  
  # add fetus and genotype info
  fetusName=sub("_complete.bed","",sub(paste0(headdir,"/"),"",inFile))
  genotypeInfo=meta[meta$name==fetusName,"genotype"]
  
  tmp$fetus = fetusName
  tmp$genotype = genotypeInfo
  
  # save data
  df[[iter]] = tmp
  
  # end for loop
}

# merge together
data_snv = as.data.frame(do.call(rbind,df))
data_snv$snv = paste(data_snv[,1],data_snv[,2],data_snv[,3],sep = '-')

# subset by genotype
# create a joint preleuk and t21 file
# df.d21 = subset(df.mg,genotype=="D21")
# df.t21 = subset(df.mg,genotype=="T21")
# df.leuk = subset(df.mg,genotype=="T21_MyeloidPreleuk")
# df.t21_leuk = subset(df.mg,genotype%in% c("T21_MyeloidPreleuk","T21"))
# 
# ind = !duplicated(df.d21$snv); sum(!ind); df.d21=df.d21[ind,]
# ind = !duplicated(df.t21$snv); sum(!ind); df.t21=df.t21[ind,]
# ind = !duplicated(df.leuk$snv); sum(!ind); df.leuk=df.leuk[ind,]
# ind = !duplicated(df.t21_leuk$snv); sum(!ind); df.t21_leuk=df.t21_leuk[ind,]

############################################

library(data.table)

atac = fread("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/differential_analysis/pb_de_atac.hsc.txt",data.table = F,stringsAsFactors = F)

# for (k in 1:4) {
# for (k in unique(data_snv$fetus)) {
    
  df.sub=subset(data_snv,genotype=="D21")
  df.sub2=subset(data_snv,genotype=="T21")
  
  snv.df=df.sub
  colnames(snv.df)[1:3] = c("chrom","start","end")
  snv.df2=df.sub2
  colnames(snv.df2)[1:3] = c("chrom","start","end")
  
  # if (k==1) {
  #   id="D21"
  #   df.sub=subset(data_snv,genotype==id)
  # } else if (k==2) {
  #   id="T21"
  #   df.sub=subset(data_snv,genotype==id)
  # } else if (k==3) {
  #   id="T21_MyeloidPreleuk"
  #   df.sub=subset(data_snv,genotype==id)
  # } else if (k==4) {
  #   id="T21_and_T21_MyeloidPreleuk"
  #   df.sub=subset(data_snv,genotype%in%c("T21_MyeloidPreleuk","T21"))
  # }
  # ind = !duplicated(df.sub$snv); sum(!ind); df.sub=df.sub[ind,]
  
  
  # print(paste0(k," - ",id))
  
  y=strsplit(atac$names,"-")
  atac.df = data.frame(chrom=paste0(unlist(lapply(y,function(x)x[[1]]))),
                       start=as.numeric(unlist(lapply(y,function(x)x[[2]])))-1,
                       end=as.numeric(unlist(lapply(y,function(x)x[[3]]))))
  atac.df = cbind(atac.df,atac)
  
  df.mg = as.data.frame(valr::bed_intersect(atac.df,snv.df))
  tab = table(df.mg$names.x)
  tab = as.data.frame(tab)
  colnames(tab)[2]="Freq_d21"
  atac.df$snv_d21 = atac.df$names %in% df.mg$names.x
  atac.df = merge(atac.df,tab,by.x="names",by.y="Var1",all.x=TRUE)
  atac.df$Freq_d21[is.na(atac.df$Freq_d21)] = 0
  
  df.mg = as.data.frame(valr::bed_intersect(atac.df,snv.df2))
  tab = table(df.mg$names.x)
  tab = as.data.frame(tab)
  colnames(tab)[2]="Freq_t21"
  
  atac.df$snv_t21 = atac.df$names %in% df.mg$names.x
  atac.df = merge(atac.df,tab,by.x="names",by.y="Var1",all.x=TRUE)
  atac.df$Freq_t21[is.na(atac.df$Freq_t21)] = 0

  atac.df$span = atac.df$end - atac.df$start
  
  # tmp = subset(atac.df,P.Value > 0.05)
  # sum(tmp$Freq)/sum(tmp$span)
  # tmp = subset(atac.df,P.Value < 0.05 & logFC > 0)
  
  calc_rate = function(dataf,column) {
    return(sum(dataf[,column])/sum(dataf[,"span"]))
  }
  
  bootstrap_snv_rate = function(j) {
    {i=sample(1:nrow(tmp),nrow(tmp),replace = T); tmp2 = tmp[i,]; val = sum(tmp2$Freq)/sum(tmp2$span); 
    return(val)}
  }
  
  bootstrap_ratio = function(j,column) {
    {i=sample(1:nrow(tmp),nrow(tmp),replace = T); tmp2 = tmp[i,]; val = calc_rate(tmp2,column)/background_t21; 
    return(val)}
  }
  
  bootstrap_open_vs_closed = function(j,column) {
    {
      tmp = subset(atac.df,P.Value < 0.05 & logFC > 0)
      i=sample(1:nrow(tmp),nrow(tmp),replace = T); tmp2 = tmp[i,]; val = calc_rate(tmp2,column); 
      tmp = subset(atac.df,P.Value < 0.05 & logFC < 0)
      i=sample(1:nrow(tmp),nrow(tmp),replace = T); tmp2 = tmp[i,]; val2 = calc_rate(tmp2,column); 
      
    return(val/val2)
    }
  }
  
  
  
  # iter = 0; res = list()
  
  # iter=iter+1
  tmp = subset(atac.df,P.Value > 0.05)
  background_d21 = calc_rate(tmp,"Freq_d21")
  background_t21 = calc_rate(tmp,"Freq_t21")
  
  tmp = subset(atac.df,P.Value < 0.05 & logFC > 0)
  val_d21=calc_rate(tmp,"Freq_d21")/background_d21
  val_t21=calc_rate(tmp,"Freq_t21")/background_t21
  val_d21;val_t21
  
  boot.open = unlist(parallel::mclapply(1:100,bootstrap_ratio,column="Freq_t21",mc.cores = 8))
  ci = quantile(boot.open,probs=c(0.025,0.975))
  ci
  
  tmp = subset(atac.df,P.Value < 0.05 & logFC > 0)
  calc_rate(tmp,"Freq_d21")
  calc_rate(tmp,"Freq_t21")
  tmp = subset(atac.df,P.Value < 0.05 & logFC < 0)
  calc_rate(tmp,"Freq_d21")
  calc_rate(tmp,"Freq_t21")
  
  bootstrap_open_vs_closed
  boot.open = unlist(parallel::mclapply(1:100,bootstrap_open_vs_closed,column="Freq_t21",mc.cores = 8))
  ci = quantile(boot.open,probs=c(0.025,0.975))
  ci
  
  
  tmp = subset(atac.df,P.Value < 0.05 & logFC < 0)
  val_d21=calc_rate(tmp,"Freq_d21")/background_d21
  val_t21=calc_rate(tmp,"Freq_t21")/background_t21
  val_d21;val_t21
  
  
  boot.open = unlist(parallel::mclapply(1:100,bootstrap_snv_rate,mc.cores = 8))
  
  
  
  sum(tmp$Freq_d21)/sum(tmp$span)
  val.open = sum(tmp$Freq)/sum(tmp$span)
  boot.open = unlist(parallel::mclapply(1:100,bootstrap_snv_rate,mc.cores = 8))
  ci = quantile(boot.open,probs=c(0.025,0.975))
  res[[iter]] = data.frame(type="Background (P > 0.05)",est=val.open,l=ci[1],h=ci[2])
  
  iter=iter+1
  tmp = subset(atac.df,P.Value < 0.05 & logFC > 0)
  val.open = sum(tmp$Freq)/sum(tmp$span)
  boot.open = unlist(parallel::mclapply(1:100,bootstrap_snv_rate,mc.cores = 8))
  ci = quantile(boot.open,probs=c(0.025,0.975))
  res[[iter]] = data.frame(type="Ts21 open (P < 0.05; LFC > 0)",est=val.open,l=ci[1],h=ci[2])
  
  iter=iter+1
  tmp = subset(atac.df,P.Value < 0.05 & logFC < 0)
  val.open = sum(tmp$Freq)/sum(tmp$span)
  boot.open = unlist(parallel::mclapply(1:100,bootstrap_snv_rate,mc.cores = 8))
  ci = quantile(boot.open,probs=c(0.025,0.975))
  res[[iter]] = data.frame(type="Ts21 closed (P < 0.05; LFC < 0)",est=val.open,l=ci[1],h=ci[2])
  
  iter=iter+1
  tmp = subset(atac.df,adj.P.Val < 0.05 & logFC > 0)
  val.open = sum(tmp$Freq)/sum(tmp$span)
  boot.open = unlist(parallel::mclapply(1:100,bootstrap_snv_rate,mc.cores = 8))
  ci = quantile(boot.open,probs=c(0.025,0.975))
  res[[iter]] = data.frame(type="Ts21 sig open (FDR < 0.05; LFC > 0)",est=val.open,l=ci[1],h=ci[2])
  
  iter=iter+1
  tmp = subset(atac.df,adj.P.Val < 0.05 & logFC < 0)
  val.open = sum(tmp$Freq)/sum(tmp$span)
  boot.open = unlist(parallel::mclapply(1:100,bootstrap_snv_rate,mc.cores = 8))
  ci = quantile(boot.open,probs=c(0.025,0.975))
  res[[iter]] = data.frame(type="Ts21 sig closed (FDR < 0.05; LFC > 0)",est=val.open,l=ci[1],h=ci[2])
  
  out.df = as.data.frame(do.call(rbind,res))
  print(out.df)
  f.out = paste0("/oak/stanford/groups/smontgom/amarder/data/hasaart_et_al_2020_scirep/out/",id,".snv_enrich.txt")
  fwrite(out.df,f.out,quote = F,na = "NA",sep = "\t",row.names = F,col.names = T)
  
}

