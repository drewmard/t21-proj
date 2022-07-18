# for cpdb sub input

# x=$((1+$(head -n1 10X_DownSyndrome_Liver.norm_count.txt | grep -P '\t' -o | wc -l)))
# echo $x

x2=$((wc -l 10X_DownSyndrome_Liver.meta.txt))
echo $x2

columns=1,$( for (( ii=2; ii<=780300; ii++ )); do echo $ii; done | sort -R | head -108922 | sort -n | tr '\n' ',' | sed 's/,$//' ) 
echo $columns > /oak/stanford/groups/smontgom/amarder/t21-proj/out/full/cpdb_Data/10X_DownSyndrome_Liver.fields.out
cut --fields=${columns} /oak/stanford/groups/smontgom/amarder/t21-proj/out/full/cpdb_Data/10X_DownSyndrome_Liver.norm_count.txt > /oak/stanford/groups/smontgom/amarder/t21-proj/out/full/cpdb_Data/10X_DownSyndrome_Liver.norm_count.sub.txt

##################
# module load R/4.1.2
library(data.table)
meta=fread("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/cpdb_Data/10X_DownSyndrome_Liver.meta.txt",data.table = F,stringsAsFactors = F)
fields=fread("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/cpdb_Data/10X_DownSyndrome_Liver.fields.out",data.table = F,stringsAsFactors = F)
meta.sub <- meta[as.numeric(fields[1,]),]
fwrite(meta.sub,"/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/cpdb_Data/10X_DownSyndrome_Liver.meta.sub.txt",quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
###########