library(data.table)
cluster_label="leiden_names"
disease_status = "DownSyndrome"
cell_type_filename = gsub("/","_",cell_type)

for (disease_status in c("Healthy","DownSyndrome")) {
  i=0; for (subset_column in c("patient","sample")) { 
    for (use_pcw in c(FALSE,TRUE)) {
      for (cluster_label in c("leiden_names")) {
        # for (cluster_label in c("combi_annot","leiden_names")) {
        i = i + 1
        # next: what labels should be used?
        if (cluster_label=="combi_annot") {
          cell_type = "HSCs"
          suffix = ".integ1"
        } else if (cluster_label=="leiden_names") {
          cell_type = cell_type
          suffix = ""
        }
        
        age_suffix = ifelse(use_pcw,".sex","")
        f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/pseudobulks/",subset_column,"/DE/",disease_status,'_',cell_type_filename,".de",age_suffix,suffix,".txt")
        res = fread(f,data.table = F,stringsAsFactors = F)
        res = res[,c("names","adj.P.Val")]
        colnames(res)[2] = paste0(colnames(res)[2],".",subset_column,".",use_pcw)
        if (i==1) {
          res.mg = res
        } else {
          res.mg = merge(res.mg,res,by="names")
        }
      }
    }
  }
  f.out=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/pseudobulks/DE.",disease_status,".sex.txt")
  fwrite(res.mg,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
}

table(res.mg$adj.P.Val.sample.FALSE < 0.1,res.mg$adj.P.Val.sample.TRUE < 0.1)
table(res.mg$adj.P.Val.sample.FALSE < 0.1,res.mg$adj.P.Val.sample.TRUE < 0.1)

use_pcw=TRUE
age_suffix = ifelse(use_pcw,".sex","")
f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/pseudobulks/",subset_column,"/DE/",disease_status,'_',cell_type_filename,".de",age_suffix,suffix,".txt")
res1 = fread(f,data.table = F,stringsAsFactors = F)
use_pcw=FALSE
age_suffix = ifelse(use_pcw,".sex","")
f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/pseudobulks/",subset_column,"/DE/",disease_status,'_',cell_type_filename,".de",age_suffix,suffix,".txt")
res2 = fread(f,data.table = F,stringsAsFactors = F)
res.mg = merge(res1,res2,by="names")

table(subset(res.mg,adj.P.Val.y < 0.1)$P.Value.x < 0.05)


