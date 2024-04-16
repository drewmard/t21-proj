library(data.table)
df = fread("~/Downloads/data_snv",data.table = F,stringsAsFactors = F)
df = df[,c(1,3,6)]
colnames(df) = c("#CHROM","POS","ID")
df$ID = gsub("-","_",df$ID)
df$REF = "A"
df$ALT = "G"
fwrite(df,"~/Downloads/data_snv.vcf",quote = F,na = "NA",sep = "\t",row.names = F,col.names = T)


library(data.table)
df = fread("~/Downloads/data_snv",data.table = F,stringsAsFactors = F)
df$id = paste0(substring(df[,1],4),":",df[,3],"-",df[,3])
# cadd = fread("~/Downloads/GRCh38-v1.6_anno_93838d79e66aa5890c99fc20b97e0052.tsv.gz",data.table = F,stringsAsFactors = F)
vep = fread("~/Downloads/mri1f1Kk5dKTZ0HU.txt",data.table = F,stringsAsFactors = F)
vep = vep[!duplicated(vep$Location),]
df2 = merge(df,vep[,c("Location","Consequence","SYMBOL")],by.y="Location",by.x="id")
sum(!(vep$Location %in% df$id))
sum(!(df$id %in% vep$Location))

tab = aggregate(id~V5+Consequence,data=df2,length)#,by=df2[,c("V5","Consequence")],length)
tab$prop = NA
i<-tab$V5=="D21"; tab[i,"prop"] = tab[i,3] / sum(tab[i,3])
i<-tab$V5=="T21"; tab[i,"prop"] = tab[i,3] / sum(tab[i,3])

library(data.table)
disease_status="DownSyndrome"
pseudobulk_covariate="sample"
suffix=""
f = paste0("~/Documents/Research/t21-proj/out/full/pseudobulks/",pseudobulk_covariate,"/DE/",disease_status,'_',"cyc_vs_hsc",".de",suffix,".txt")
atac=fread(f,data.table = F,stringsAsFactors = F)[,1:7]
cre = fread("~/Downloads/GRCh38-cCREs.bed",data.table=F,stringsAsFactors = F)
# cre = fread("/oak/stanford/groups/smontgom/amarder/data/ENCFF394DBM.bed",data.table=F,stringsAsFactors = F)
colnames(cre)[1:3] = c("chrom",'start','end')

df2 = df2[,c(2:4,1,5:ncol(df2))]
colnames(df2)[1:3] = c("chrom",'start','end')
tmp = as.data.frame(valr::bed_intersect(cre,df2))
df2$cre = (df2$id %in% tmp$id.y)
df3 = merge(df2,atac,by.x="SYMBOL",by.y="names",all.x = TRUE)
df3.sub = df3[grepl("intron",df3$Consequence) & !is.na(df3$SYMBOL) & df3$V5=="T21",]
df3.sub = df3[grepl("intron|interg",df3$Consequence) & !is.na(df3$SYMBOL) & !is.na(df3$logFC) & df3$V5!="T21_MyeloidPreleuk",]
# df3.sub = df3[df3$Consequence %in% c("intron_variant","intron_variant,non_coding_transcript_variant") & df3$V5!="T21_MyeloidPreleuk",]
df3.sub = df3[df3$Consequence %in% c("intron_variant","intron_variant,non_coding_transcript_variant") & !is.na(df3$logFC) & df3$V5!="T21_MyeloidPreleuk",]
i <- df3.sub$V5=="T21"
tmp1 = df3.sub[i,]
tmp2 <- df3.sub[!i,]
binom.test(sum(tmp1$cre),nrow(tmp1),sum(tmp2$cre)/nrow(tmp2))

pvec = c()
for (j in 1:10000) {
  # i = sample(1:nrow(tmp1),nrow(tmp1),replace = T)
  # tmp1.boot = tmp1[i,]
  # p = sum(tmp1.boot$cre)/nrow(tmp1.boot)
  # pvec = c(pvec,p)
  i = sample(1:nrow(tmp2),nrow(tmp2),replace = T)
  tmp2.boot = tmp2[i,]
  p = sum(tmp2.boot$cre)/nrow(tmp2.boot)
  pvec = c(pvec,p)
}
library(ggplot2)
g=ggplot(data.frame(pvec),aes(x=pvec)) + geom_histogram(fill='darkblue',col='white',bins=16) + 
  # theme_bw() +
  ggpubr::theme_pubr() + 
  theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5)) + 
  geom_vline(aes(xintercept = sum(tmp1$cre)/nrow(tmp1),color='Down Syndrome'),lty='dashed',lwd=2) +
  geom_vline(aes(xintercept = sum(tmp2$cre)/nrow(tmp2),color='Disomic'),lty='dashed',lwd=2) +
  labs(x="Empirical bootstrap distribution in disomic somatics",
       title="Percentage of intronic SNVs in DS HSC-expressed genes that overlap CREs",
       y='Count') +
  scale_color_manual(name = "Observed external somatic set", values = c(`Down Syndrome` = "yellow3", Disomic = "orange"))
g
pdf("~/Downloads/snv_cre_enrich.pdf",width = 9,height=4)
print(g)
dev.off()
# png("~/Downloads/snv_cre_enrich.png",width = 90*4,height=40*4)
# print(g)
# dev.off()
1 - mean((sum(tmp1$cre)/nrow(tmp1)) > pvec)


binom.test(sum(tmp1$adj.P.Val < 0.05),nrow(tmp1),sum(tmp2$adj.P.Val < 0.05)/nrow(tmp2))

# exact2x2::exact2x2(df3.sub$cre,df3.sub$V5=="T21")
df3.sub = df3[df3$V5=="D21",]
fisher.test(df3.sub$cre,df3.sub$adj.P.Val < 0.01)

fisher.test(df3.sub$cre,df3.sub$V5!="D21")
fisher.test(df3.sub$cre,df3.sub$adj.P.Val < 0.05)
fisher.test(df3.sub$cre,df3.sub$adj.P.Val < 0.01 & df3.sub$logFC > 1)

fisher.test(df2.sub$cre,df2.sub$V5=="T21")


tmp = tmp[!duplicated(tmp$id.y),]
tmp[,c(id.y)]
df3 = merge(df2,tmp[,c("V6.y","V6.x")],by.x="V6",by.y="V6.y",all.x=TRUE)
df3$cre = !is.na(df3$V6.x)
sum(df2$id %in% tmp$id.y)
table(df3$cre)

df.mg = as.data.frame(valr::bed_intersect(cre,snv.df))
print(nrow(df.mg))
print(nrow(snv.df))

fisher.test(df3$V5!="D21",df3$adj.P.Val < 0.05)
fisher.test(df3$V5!="D21",df3$P.Value < 0.05)

table(vep$Consequence)
sum(!(vep$Location %in% paste0(substring(df[,1],4),":",df[,2],"-",df[,2])))
sum(!(vep$Location %in% paste(substring(df[,1],4),df[,3],df[,3],sep = "-")))

table(cadd$Consequence)
