library(data.table)
dataf = fread("~/Downloads/a3_tab.obs.withannot.csv",data.table = F,stringsAsFactors = F)
new_annot_mapping = fread("~/Downloads/new_annot_mapping.csv",data.table = F,stringsAsFactors = F)

# dataf = fread("~/Downloads/b3_tab.obs.withannot.csv",data.table = F,stringsAsFactors = F)
# dataf$combi_annot = dataf$combi_annot_v3
# colnames(dataf)[3] = "leiden_latest"
# new_annot_mapping = fread("~/Downloads/new_annot_mapping.scvi.csv",data.table = F,stringsAsFactors = F)
# colnames(new_annot_mapping)[1] = "combi_annot"
# colnames(dataf)[colnames(dataf)=="organ"] <- "Organ"
# colnames(dataf)[colnames(dataf)=="disease_status"] <- "DiseaseStatus"

library(data.table)
i=0; df = list()
for (disease_status in c("DownSyndrome","Healthy")) {
  for (sampletype in c("Femur","Liver")) { 
    i=i+1
    f=paste0("~/Downloads/10X_",disease_status,"_",sampletype,".cellComp.csv")
    # f=paste0("~/Documents/Research/t21-proj/out/full/cellComp/10X_",disease_status,"_",sampletype,".cellComp.csv")
    df[[i]]=fread(f,data.table = F,stringsAsFactors = F)
    colnames(df[[i]])[6] = "leiden_latest"
    df[[i]]$disease=disease_status
    df[[i]]$env = sampletype
  }
}
df.lst=df
df=as.data.frame(do.call(rbind,df.lst))
df.uniq = unique(df[,c("leiden_latest","cell_type_groups")])
# df.uniq[order(df.uniq[,1]),]

dataf2 = merge(dataf,df.uniq,by="leiden_latest",all.x=TRUE)
# there's some drop out of cells that were unknown and therefore got removed
df$patient_sample = paste0(df$patient," ",df$sample)
dataf2$cellname=unlist(lapply(strsplit(dataf2$barcodekey,"-"),function(x) paste(x[1:2],collapse = "-")))
df$cellname=unlist(lapply(strsplit(df$V1,"-"),function(x) paste(x[1:2],collapse = "-")))

dataf3 = merge(dataf2,df[,c("cellname","patient_sample","sorting")],by=c("cellname","patient_sample"),all.x="TRUE")
sum(!is.na(dataf3$sorting))
sum(!is.na(dataf3$cell_type_groups))

library(data.table)
# dont need this step, because reading file at the top of the script now instead:
# new_annot_mapping = fread("~/Downloads/new_annot_mapping.csv",data.table = F,stringsAsFactors = F)
dataf3 = merge(dataf3,new_annot_mapping,by="combi_annot")
dataf3.old = dataf3
rm(dataf2)
rm(df)
rm(dataf)

sum(!(keepers$Patient.ID %in% dataf3$patient_sample))

###########

keepers = read.table(text="Patient ID	Patient	Tissue	Organ	Sorting strategy	Age
15532 T21 15532 A	15532	Down Syndrome	Liver	CD45+	14 pcw
15532 T21 15532 B	15532	Down Syndrome	Liver	CD45+	14 pcw
15559 T21 15559 B	15559	Down Syndrome	Liver	CD45+	12 pcw
15559 T21 15559 C	15559	Down Syndrome	Liver	CD45+	12 pcw
15582 T21 15582 B	15582	Down Syndrome	Liver	CD45+	13 pcw
15646 15646B	15646	Down Syndrome	Liver	CD45+	13 pcw
15669 15669H	15669	Down Syndrome	Liver	CD45+	14 pcw
15712 15712A	15712	Down Syndrome	Liver	CD45+	12 pcw
15724 15724C	15724	Down Syndrome	Liver	CD45+	13 pcw
15734 L15734A	15734	Down Syndrome	Liver	CD45+	14 pcw
15633 L15633B	15633	Healthy	Liver	CD45+	14 pcw
15657 L15657A	15657	Healthy	Liver	CD45+	12 pcw
15781 L15781D	15781	Healthy	Liver	CD45+	13 pcw",sep="\t",header=T)


# dataf3 = subset(dataf3.old,sorting=="CD45+")
dataf3 = subset(dataf3.old,sorting=="CD45+" & patient_sample %in% keepers$Patient.ID)
# tab = aggregate(dataf3$combi_annot,by=list(annot=dataf3$broad_annot,organ=dataf3$Organ,disease=dataf3$DiseaseStatus),length)
tab = aggregate(dataf3$combi_annot,by=list(annot=dataf3$broad_annot,organ=dataf3$Organ,disease=dataf3$DiseaseStatus,patient_sample=dataf3$patient_sample),length)
library(dplyr)
tab <- tab %>%
  group_by(organ, disease,patient_sample) %>%
  mutate(prop = x / sum(x))
tab = as.data.frame(tab)

# tab2 = aggregate(dataf3$combi_annot,by=list(annot=dataf3$cell_type_groups,organ=dataf3$Organ,disease=dataf3$DiseaseStatus),length)
tab2 = aggregate(dataf3$combi_annot,by=list(annot=dataf3$cell_type_groups,organ=dataf3$Organ,disease=dataf3$DiseaseStatus,patient_sample=dataf3$patient_sample),length)
library(dplyr)
tab2 <- tab2 %>%
  group_by(organ, disease,patient_sample) %>%
  mutate(prop = x / sum(x))
tab2 = as.data.frame(tab2)
# tab2

# tab.mg = merge(tab,tab2,by=c("annot","organ","disease"))
# colnames(tab.mg)[4:7] = c("N.int","prop.int","N.sep","prop.sep")
tab.mg = merge(tab,tab2,by=c("annot","organ","disease","patient_sample"))
colnames(tab.mg)[5:8] = c("N.int","prop.int","N.sep","prop.sep")
tab.mg$organ_disease = paste0(tab.mg$organ," ",tab.mg$disease)
library(ggplot2)
# ggplot(tab.mg,aes(x=prop.sep,y=prop.int,col=annot,shape=organ)) + #,fill=disease)) + 
g1=ggplot(tab.mg,aes(x=prop.sep,y=prop.int,col=annot,shape=organ_disease)) + #,fill=disease)) + 
  geom_abline(slope=1,intercept=0,lty='dashed') + 
  geom_point(size=rel(3)) +
  scale_shape_manual(values = c(1,16,2,17),label=c("Ts21","Disomic")) +
  scale_fill_manual(values = c("white", "black")) +
  theme_bw() +
  # ggpubr::theme_pubr() + 
  theme(panel.grid = element_blank()) +
  labs(x="Separate",y="Integrated",shape="Sample Type",col="Cell type") +
  scale_color_brewer(palette="Set2")
cor.test(tab.mg$prop.sep,tab.mg$prop.int)
  
tab = aggregate(dataf3$combi_annot,by=list(annot=dataf3$broad_annot,organ=dataf3$Organ,disease=dataf3$DiseaseStatus,sample=dataf3$patient_sample),length)
library(dplyr)
tab <- tab %>%
  group_by(organ, disease,sample) %>%
  mutate(prop = x / sum(x))
tab = as.data.frame(tab)
cluster_percentages_all = subset(tab,organ=="Liver")
g2<-ggplot(cluster_percentages_all, aes(fill=annot, y=prop*100, x=reorder(sample,disease=="DownSyndrome"),col=disease)) + 
  geom_bar(position="stack", stat="identity") +
  scale_color_manual(values=c("black","orange"),labels=c("Down Syndrome","Healthy"),name="Disease Status") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(hjust=1,angle=60),
        plot.title = element_text(hjust=0.5)) +
  scale_fill_brewer(palette="Set3") +
  labs(x="Patient ID",y="Cell composition (%)",title="CD45+",fill="Cell type");g2
pdf("~/Downloads/abundance_comparison.pdf",width = 14,height=5.5)
print(cowplot::plot_grid(g1,g2,ncol=2))
dev.off()

broad_groups = unique(cluster_percentages_all$annot)
# broad_groups=subset(broad_groups,broad_groups != "Stroma")
res.lst=list()
for (i in 1:length(broad_groups)) {
  print(i)
  tmp = subset(cluster_percentages_all,annot==broad_groups[i])[,c("disease","sample","x","prop")]
  tmp2=subset(unique(cluster_percentages_all[,c("disease","sample")]),!(sample %in% tmp$sample))
  if (nrow(tmp2) > 0) {
    tmp2=data.frame(disease=tmp2$disease,sample=tmp2$sample,x=0,prop=0)
    tmp = rbind(tmp,tmp2)
  }
  tres=t.test(prop~disease,tmp)
  mwu_pval=wilcox.test(prop~disease,tmp)$p.value
  res.lst[[i]] = data.frame(annot=broad_groups[i],mwu_pval,cases=as.numeric(tres$estimate[1]),controls=as.numeric(tres$estimate[2]),t_pval=as.numeric(tres$p.value))
}
res <- as.data.frame(do.call(rbind,res.lst))
res$t_pval.adj <- p.adjust(res$t_pval,method='fdr')
res$mwu_pval.adj <- p.adjust(res$mwu_pval,method='fdr')
res
fwrite(res,"~/Downloads/integrated_abundance_results.cd45.scvi.csv",quote = F,na = "NA",row.names = F,col.names = T,sep = ',')



tab.mg2 = subset(tab.mg,organ=="Liver")
tab.mg2 = merge(subset(tab.mg2,disease=="DownSyndrome"),subset(tab.mg,disease!="DownSyndrome"),by=c("annot"))
library(ggplot2)
ggplot(tab.mg2,aes(x=prop.int.y,y=prop.int.x,col=annot)) + #,fill=disease)) + 
  geom_abline(slope=1,intercept=0,lty='dashed') + 
  geom_point(size=rel(3)) +
  scale_fill_manual(values = c("white", "black")) +
  theme_bw() + theme(panel.grid = element_blank()) +
  labs(x="Control",y="Ts21",title="Integrated") +
  scale_color_brewer(palette="Set2")
ggplot(tab.mg2,aes(x=prop.sep.y,y=prop.sep.x,col=annot)) + #,fill=disease)) + 
  geom_abline(slope=1,intercept=0,lty='dashed') + 
  geom_point(size=rel(3)) +
  scale_fill_manual(values = c("white", "black")) +
  theme_bw() + theme(panel.grid = element_blank()) +
  labs(x="Control",y="Ts21",title="Separate") +
  scale_color_brewer(palette="Set2")


unique(df.uniq[,2])
unique(dataf2$combi_annot)
broad_label="Erythroid"
narrow_labels=c("Erythroid",)
dataf2[dataf2$combi_annot%in%narrow_labels] = broad_label

broad_label="HSC/Progenitors"
narrow_labels=c("Granulocyte progenitors")

broad_label="Myeloid"
narrow_labels=c("")

broad_label="NK/T cells"
narrow_labels=c("")

broad_label="Mast cells"
narrow_labels=c("")

broad_label="Megakaryocytes"
narrow_labels=c("")

broad_label="B cells"
narrow_labels=c("")

broad_label="Stroma"
narrow_labels=c("")




