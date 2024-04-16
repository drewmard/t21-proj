library(readxl)
library(ggplot2)
dataf = as.data.frame(read_excel("~/Downloads/41467_2021_21064_MOESM13_ESM.xlsx",sheet = "Table S12.csv",skip=1))
res = fread("~/Downloads/DownSyndrome_HSC_PeakGeneSets/Liver.HSCs_MPPs.sample.txt",data.table = F,stringsAsFactors = F)
df.mg = merge(dataf,res,by.x="Gene",by.y="names")
fisher.test(df.mg$log2FoldChange>0,df.mg$logFC>0)$p.value
cor.test(df.mg$log2FoldChange,df.mg$logFC)
g=ggplot(df.mg,aes(x=log2FoldChange,y=logFC)) + 
  geom_point() + 
  labs(x="Muskens et al: Fetal Liver CD34+ cells",y="Marderstein et al: Fetal Liver HSCs") + 
  theme_bw() +
  ggpubr::theme_pubr() +
  geom_vline(xintercept = 0,col='red',lty='dashed') +
  geom_hline(yintercept = 0,col='red',lty='dashed') +
  geom_abline(slope = 1,intercept = 0,lty='dashed',col='blue') +
  xlim(-10,15) +
  ylim(-3,5) +
  labs(title = "Differential Expression Comparison (Ts21 vs Disomy)") +
  theme(plot.title = element_text(hjust=0.5));g
f="~/Documents/Research/t21-proj/out/full/methylation/de_comp.pdf"
pdf(f,width = 8,height = 4)
print(g)
dev.off()

df.mg = merge(dataf,res,by.x="Gene",by.y="names",all.y=TRUE)
fisher.test(!is.na(df.mg$pvalue),df.mg$P.Value<0.05)
fisher.test(!is.na(df.mg$pvalue),df.mg$adj.P.Val<0.1)

dataf = as.data.frame(read_excel("~/Downloads/41467_2021_21064_MOESM11_ESM.xlsx",sheet = "a_DMRs_overall",skip=1))
dataf$DMRname = paste0(dataf$chrom,"-",dataf$`start (hg19)`,"-",dataf$`end (hg19)`)
head(dataf)
dataf.sub = subset(dataf,`DMR location` %in% c("enhancer",'promoter'))
fisher.test(dataf.sub$DMRCate_meandiff<0,as.numeric(dataf.sub$`Gene expression foldchange (DS vs non-DS FL HSPC)`) > 1)
df.mg = merge(dataf.sub,res,by.x="Associated gene",by.y="names",all.x=TRUE)
fisher.test(df.mg$DMRCate_meandiff<0,as.numeric(df.mg$`Gene expression foldchange (DS vs non-DS FL HSPC)`) > 1)
fisher.test(df.mg$DMRCate_meandiff<0,as.numeric(df.mg$logFC) > 0)

df.mg$DE = ifelse(df.mg$logFC>0,"Upregulated","Downregulated")
df.mg$DE[is.na(df.mg$logFC)] = "Not expressed"
df.mg$Methylation = ifelse(df.mg$DMRCate_meandiff>0,"Hypermethylated","Hypomethylated")
table_data = as.data.frame(table(df.mg$DE,df.mg$Methylation))
library(plyr)
# Create a stacked barplot
table_data = ddply(table_data, "Var2", transform, Freq2 = Freq / sum(Freq))
table_data$Var1 = factor(table_data$Var1,levels = c("Not expressed","Upregulated","Downregulated"))
g=ggplot(table_data, aes(x = Var2, y = Freq2, fill = Var1)) +
  geom_bar(stat = "identity",col='black') +
  labs(title = "Marderstein et al.",
       x = "",
       y = "% of DMRs",
       fill = "") +
  theme_minimal() + 
  scale_fill_viridis_d() +
  theme(plot.title = element_text(hjust=0.5)) +
  theme(panel.grid = element_blank())
pdf("~/Documents/Research/t21-proj/out/full/methylation/meth_exp.marderstein.pdf",width = 5,height=4)
print(g)
dev.off()

table(df.mg$DE,df.mg$DMRCate_meandiff > 0)
df.mg$DE = ifelse(as.numeric(df.mg$`Gene expression foldchange (DS vs non-DS FL HSPC)` > 1),"Upregulated","Downregulated")
df.mg$DE[is.na(as.numeric(df.mg$`Gene expression foldchange (DS vs non-DS FL HSPC)`))] = "Not expressed"

df.mg$Methylation = ifelse(df.mg$DMRCate_meandiff>0,"Hypermethylated","Hypomethylated")
table_data = as.data.frame(table(df.mg$DE,df.mg$Methylation))
library(plyr)
# Create a stacked barplot
table_data = ddply(table_data, "Var2", transform, Freq2 = Freq / sum(Freq))
table_data$Var1 = factor(table_data$Var1,levels = c("Not expressed","Upregulated","Downregulated"))
g=ggplot(table_data, aes(x = Var2, y = Freq2, fill = Var1)) +
  geom_bar(stat = "identity",col='black') +
  labs(title = "Muskens et al.",
       x = "",
       y = "% of DMRs",
       fill = "") +
  theme_minimal() + 
  scale_fill_viridis_d() +
  theme(plot.title = element_text(hjust=0.5)) +
  theme(panel.grid = element_blank())
pdf("~/Documents/Research/t21-proj/out/full/methylation/meth_exp.muskens.pdf",width = 5,height=4)
print(g)
dev.off()

dataf.sub = subset(dataf,!(`DMR location` %in% c("enhancer",'promoter')))
fisher.test(dataf.sub$DMRCate_meandiff<0,as.numeric(dataf.sub$`Gene expression foldchange (DS vs non-DS FL HSPC)`) > 1)
df.mg = merge(dataf.sub,res,by.x="Associated gene",by.y="names",all.x=TRUE)
fisher.test(df.mg$DMRCate_meandiff<0,as.numeric(df.mg$`Gene expression foldchange (DS vs non-DS FL HSPC)`) > 1)
fisher.test(df.mg$DMRCate_meandiff<0,as.numeric(df.mg$logFC) > 0)

df.mg$DE = ifelse(df.mg$logFC>0,"Upregulated","Downregulated")
df.mg$DE[is.na(df.mg$logFC)] = "Not expressed"
df.mg$Methylation = ifelse(df.mg$DMRCate_meandiff>0,"Hypermethylated","Hypomethylated")
table_data = as.data.frame(table(df.mg$DE,df.mg$Methylation))
library(plyr)
# Create a stacked barplot
table_data = ddply(table_data, "Var2", transform, Freq2 = Freq / sum(Freq))
table_data$Var1 = factor(table_data$Var1,levels = c("Not expressed","Upregulated","Downregulated"))
g=ggplot(table_data, aes(x = Var2, y = Freq2, fill = Var1)) +
  geom_bar(stat = "identity",col='black') +
  labs(title = "Marderstein et al.",
       x = "",
       y = "% of DMRs",
       fill = "") +
  theme_minimal() + 
  scale_fill_viridis_d() +
  theme(plot.title = element_text(hjust=0.5)) +
  theme(panel.grid = element_blank())
pdf("~/Documents/Research/t21-proj/out/full/methylation/meth_exp.marderstein.out_prom_enh.pdf",width = 5,height=4)
print(g)
dev.off()

table(df.mg$DE,df.mg$DMRCate_meandiff > 0)
df.mg$DE = ifelse(as.numeric(df.mg$`Gene expression foldchange (DS vs non-DS FL HSPC)` > 1),"Upregulated","Downregulated")
df.mg$DE[is.na(as.numeric(df.mg$`Gene expression foldchange (DS vs non-DS FL HSPC)`))] = "Not expressed"

df.mg$Methylation = ifelse(df.mg$DMRCate_meandiff>0,"Hypermethylated","Hypomethylated")
table_data = as.data.frame(table(df.mg$DE,df.mg$Methylation))
library(plyr)
# Create a stacked barplot
table_data = ddply(table_data, "Var2", transform, Freq2 = Freq / sum(Freq))
table_data$Var1 = factor(table_data$Var1,levels = c("Not expressed","Upregulated","Downregulated"))
g=ggplot(table_data, aes(x = Var2, y = Freq2, fill = Var1)) +
  geom_bar(stat = "identity",col='black') +
  labs(title = "Muskens et al.",
       x = "",
       y = "% of DMRs",
       fill = "") +
  theme_minimal() + 
  scale_fill_viridis_d() +
  theme(plot.title = element_text(hjust=0.5)) +
  theme(panel.grid = element_blank())
pdf("~/Documents/Research/t21-proj/out/full/methylation/meth_exp.muskens.out_prom_enh.pdf",width = 5,height=4)
print(g)
dev.off()


res = fread("~/Downloads/DownSyndrome_HSC_PeakGeneSets/pb_de_atac.hsc.txt",data.table = F,stringsAsFactors = F)
y = strsplit(res$names,"-")
res$chrom = unlist(lapply(y,function(x)x[[1]]))
res$start = as.numeric(unlist(lapply(y,function(x)x[[2]])))
res$end = as.numeric(unlist(lapply(y,function(x)x[[3]])))
# Create a GRanges object
gr <- GRanges(seqnames = res$chrom,
              ranges = IRanges(start = res$start, end = res$end),
              names = res$names)
library(rtracklayer)
chain_file="~/Downloads/hg38ToHg19.over.chain"
chain <- import.chain(chain_file)
res19 <- liftOver(gr, chain)
res19 = unlist(res19)
# Extract relevant columns
res19 <- data.frame(
  chrom = seqnames(res19),
  start = start(res19),
  end = end(res19),
  names = res19$names
)
# dataf.sub = subset(dataf,!(`DMR location` %in% c("enhancer",'promoter')))
library(valr)
dataf.tmp = dataf[,c("chrom","start (hg19)","end (hg19)","DMRname")]
colnames(dataf.tmp)[2:3] = c("start","end")
df.mg = bed_intersect(res19,dataf.tmp,suffix = c("",""))
df.mg = merge(df.mg[,c("names","DMRname")],res,by="names",all.x=TRUE)
df.mg = merge(df.mg,dataf,by="DMRname",all.x=TRUE)
fisher.test(df.mg$logFC > 0,df.mg$DMRCate_meandiff < 0)
fisher.test(df.mg$logFC > 0,df.mg$DMRCate_meandiff < 0)$p.value

df.mg.tmp = subset(df.mg,(`DMR location` %in% c("enhancer",'promoter')))
df.mg.tmp$DE = ifelse(df.mg.tmp$logFC>0,"Increased accessibility","Decreased accessibility")
df.mg.tmp$DE[is.na(df.mg.tmp$logFC)] = "Not accessible"
df.mg.tmp$Methylation = ifelse(df.mg.tmp$DMRCate_meandiff>0,"Hypermethylated","Hypomethylated")
table_data = as.data.frame(table(df.mg.tmp$DE,df.mg.tmp$Methylation))
library(plyr)
# Create a stacked barplot
table_data = ddply(table_data, "Var2", transform, Freq2 = Freq / sum(Freq))
table_data$Var1 = factor(table_data$Var1,levels = c("Not accessible","Increased accessibility","Decreased accessibility"))

g=ggplot(table_data, aes(x = Var2, y = Freq2, fill = Var1)) +
  geom_bar(stat = "identity",col='black') +
  labs(title = "DMRs within prom/enh",
       x = "",
       y = "% of DMRs",
       fill = "") +
  theme_minimal() + 
  scale_fill_viridis_d(option="G") +
  theme(plot.title = element_text(hjust=0.5)) +
  theme(panel.grid = element_blank())
pdf("~/Documents/Research/t21-proj/out/full/methylation/meth_access.pdf",width = 5,height=4)
print(g)
dev.off()

# df.mg.tmp = subset(df.mg,(`DMR location` == c("gene body")))
# df.mg.tmp = subset(df.mg,(`DMR location` == c("intergenic")))
df.mg.tmp = subset(df.mg,!(`DMR location` %in% c("enhancer",'promoter')))
fisher.test(df.mg.tmp$logFC > 0,df.mg.tmp$DMRCate_meandiff < 0)

df.mg.tmp$DE = ifelse(df.mg.tmp$logFC>0,"Increased accessibility","Decreased accessibility")
df.mg.tmp$DE[is.na(df.mg.tmp$logFC)] = "Not accessible"
df.mg.tmp$Methylation = ifelse(df.mg.tmp$DMRCate_meandiff>0,"Hypermethylated","Hypomethylated")
table_data = as.data.frame(table(df.mg.tmp$DE,df.mg.tmp$Methylation))
library(plyr)
# Create a stacked barplot
table_data = ddply(table_data, "Var2", transform, Freq2 = Freq / sum(Freq))
table_data$Var1 = factor(table_data$Var1,levels = c("Not accessible","Increased accessibility","Decreased accessibility"))
g=ggplot(table_data, aes(x = Var2, y = Freq2, fill = Var1)) +
  geom_bar(stat = "identity",col='black') +
  labs(title = "DMRs outside prom/enh",
       x = "",
       y = "% of DMRs",
       fill = "") +
  theme_minimal() + 
  scale_fill_viridis_d(option="G") +
  theme(plot.title = element_text(hjust=0.5)) +
  theme(panel.grid = element_blank())
pdf("~/Documents/Research/t21-proj/out/full/methylation/meth_access.outside_prom_enh.pdf",width = 5,height=4)
print(g)
dev.off()





