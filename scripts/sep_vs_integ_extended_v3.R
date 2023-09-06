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


# library(data.table)
# refmap = fread("~/Downloads/scArches_Healthy_Liver_metadata.csv",data.table = F,stringsAsFactors = F)
# refmap$cellname=unlist(lapply(strsplit(refmap$barcode,"-"),function(x) paste(x[1:2],collapse = "-")))
# refmap$patient_sample = paste0(refmap$patient," ",refmap$sample)
# refmap_broad = fread("~/Downloads/ref_broad_map.csv",data.table = F,stringsAsFactors = F)
# refmap = merge(refmap,refmap_broad,by="Predicted")

library(data.table)
refmap = fread("~/Downloads/scArches_Healthy_Liver_metadata.csv",data.table = F,stringsAsFactors = F)
refmap$cellname=unlist(lapply(strsplit(refmap$barcode,"-"),function(x) paste(x[1:2],collapse = "-")))
refmap$patient_sample = paste0(refmap$patient," ",refmap$sample)
refmap_t21 = fread("~/Downloads/Healthy_liver_annotations.txt",data.table = F,stringsAsFactors = F)
refmap_t21 = refmap_t21[,c(1:2,5:6)]
refmap_t21 = merge(refmap[,c("patient","sample","cellname","patient_sample")],refmap_t21,by=c("sample","cellname"))
refmap_broad = fread("~/Downloads/ref_broad_map.csv",data.table = F,stringsAsFactors = F)
refmap = merge(refmap_t21,refmap_broad,by.x="Popescu_transfer",by.y="Predicted")
refmap_broad_t21 = fread("~/Downloads/t21_transfer.broad.csv",data.table = F,stringsAsFactors = F)
refmap = merge(refmap,refmap_broad_t21,by.x="T21_transfer",by.y="T21_transfer")

dataf = fread("~/Downloads/a3_tab.obs.withannot.csv",data.table = F,stringsAsFactors = F)
dataf.scvi = fread("~/Downloads/b3_tab.obs.withannot.csv",data.table = F,stringsAsFactors = F)
dataf = merge(dataf,dataf.scvi[,c("barcodekey","combi_annot_v3")],by="barcodekey")

dataf2 = merge(dataf,df.uniq,by="leiden_latest",all.x=TRUE)
dataf2$cellname=unlist(lapply(strsplit(dataf2$barcodekey,"-"),function(x) paste(x[1:2],collapse = "-")))
df$cellname=unlist(lapply(strsplit(df$V1,"-"),function(x) paste(x[1:2],collapse = "-")))
df$patient_sample = paste0(df$patient," ",df$sample)
dataf3 = merge(dataf2,df[,c("cellname","patient_sample","sorting")],by=c("cellname","patient_sample"),all.x=TRUE)

library(data.table)
new_annot_mapping = fread("~/Downloads/new_annot_mapping.csv",data.table = F,stringsAsFactors = F)
dataf3 = merge(dataf3,new_annot_mapping,by="combi_annot",all.x=TRUE)
new_annot_mapping_scvi = fread("~/Downloads/new_annot_mapping.scvi.csv",data.table = F,stringsAsFactors = F)
colnames(new_annot_mapping_scvi)[2] = "broad_annot_scvi"
dataf3 = merge(dataf3,new_annot_mapping_scvi,by="combi_annot_v3",all.x=TRUE)
dataf3.old = dataf3
rm(dataf2)
rm(df)
rm(dataf)
rm(dataf.scvi)

dataf4=merge(dataf3,refmap,by=c("patient_sample","cellname"),all.x=TRUE)
dataf4.full = dataf4
dataf4 = subset(dataf4.full,Organ=="Liver" & DiseaseStatus=="Healthy" )#& broad_annot != "Stroma" & cell_type_groups !="Stroma" & broad_annot_scvi!="Stroma")
dataf4[dataf4=="Stroma"] = NA

####################

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
dataf3 = subset(dataf4,sorting=="CD45+" & patient_sample %in% keepers$Patient.ID)
# tab = aggregate(dataf4.keep$combi_clust,by=list(annot=dataf4.keep$Ref_Broad,organ=dataf4.keep$Organ,disease=dataf4.keep$DiseaseStatus,patient_sample=dataf4.keep$patient_sample),length)
# library(dplyr)
# tab <- tab %>%
#   group_by(organ, disease,patient_sample) %>%
#   mutate(prop = x / sum(x))
# tab = as.data.frame(tab)
# colnames(tab)[4:5]=paste0(c("N","Prop"),"_refmap")

# tab2 = aggregate(dataf3$combi_annot,by=list(annot=dataf3$cell_type_groups,organ=dataf3$Organ,disease=dataf3$DiseaseStatus),length)
tab = aggregate(dataf3$combi_annot,by=list(annot=dataf3$T21Ref_Broad,organ=dataf3$Organ,disease=dataf3$DiseaseStatus,patient_sample=dataf3$patient_sample),length)
library(dplyr)
tab <- tab %>%
  group_by(organ, disease,patient_sample) %>%
  mutate(prop = x / sum(x))
tab = as.data.frame(tab)

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
g1=ggplot(tab.mg,aes(x=prop.sep,y=prop.int,col=annot)) + #,fill=disease)) + 
  geom_abline(slope=1,intercept=0,lty='dashed') + 
  geom_point(size=rel(3)) +
  # scale_shape_manual(values = c(1,16,2,17),label=c("Ts21","Disomic")) +
  scale_fill_manual(values = c("white", "black")) +
  theme_bw() +
  # ggpubr::theme_pubr() + 
  theme(panel.grid = element_blank()) +
  labs(x="Separate",y="Reference mapping (Ts21)",shape="Sample Type",col="Cell type") +
  scale_color_brewer(palette="Set2"); g1

cor.test(tab.mg$prop.sep,tab.mg$prop.int)
pdf("~/Downloads/abundance_comparison.pdf",width = 14,height=5.5)
print(cowplot::plot_grid(g1,g2,ncol=2))
dev.off()


# dataf3 = subset(dataf3.old,sorting=="CD45+" & patient_sample %in% keepers$Patient.ID)
tab = aggregate(dataf4$combi_clust,by=list(annot=dataf4$broad_annot,organ=dataf4$Organ,disease=dataf4$DiseaseStatus),length)
library(dplyr)
tab <- tab %>%
  group_by(organ, disease) %>%
  mutate(prop = x / sum(x))
tab = as.data.frame(tab)
colnames(tab)[4:5]=paste0(c("N","Prop"),"_harm_integ")

tab2 = aggregate(dataf4$combi_clust,by=list(annot=dataf4$cell_type_groups,organ=dataf4$Organ,disease=dataf4$DiseaseStatus),length)
library(dplyr)
tab2 <- tab2 %>%
  group_by(organ, disease) %>%
  mutate(prop = x / sum(x))
tab2 = as.data.frame(tab2)
colnames(tab2)[4:5]=paste0(c("N","Prop"),"_sep")

tab3 = aggregate(dataf4$combi_clust,by=list(annot=dataf4$Ref_Broad,organ=dataf4$Organ,disease=dataf4$DiseaseStatus),length)
library(dplyr)
tab3 <- tab3 %>%
  group_by(organ, disease) %>%
  mutate(prop = x / sum(x))
tab3 = as.data.frame(tab3)
colnames(tab3)[4:5]=paste0(c("N","Prop"),"_refmap")

tab4 = aggregate(dataf4$combi_clust,by=list(annot=dataf4$broad_annot_scvi,organ=dataf4$Organ,disease=dataf4$DiseaseStatus),length)
library(dplyr)
tab4 <- tab4 %>%
  group_by(organ, disease) %>%
  mutate(prop = x / sum(x))
tab4 = as.data.frame(tab4)
colnames(tab4)[4:5]=paste0(c("N","Prop"),"_scvi_integ")

tab5 = aggregate(dataf4$combi_clust,by=list(annot=dataf4$T21Ref_Broad,organ=dataf4$Organ,disease=dataf4$DiseaseStatus),length)
library(dplyr)
tab5 <- tab5 %>%
  group_by(organ, disease) %>%
  mutate(prop = x / sum(x))
tab5 = as.data.frame(tab5)
colnames(tab5)[4:5]=paste0(c("N","Prop"),"_T21refmap")

tab.mg = merge(merge(merge(merge(tab,tab2,by="annot"),tab3,by='annot'),tab4,by="annot"),tab5,by = "annot")
tab.mg = tab.mg[,c("annot","Prop_sep","Prop_harm_integ","Prop_refmap","Prop_scvi_integ","Prop_T21refmap")]
library(ggplot2)
g1=ggplot(tab.mg,aes(x=Prop_sep,y=Prop_harm_integ,col=annot)) + 
  geom_point(size=rel(3)) + 
  geom_abline(slope=1,intercept=0,col='red',lty='dashed') + 
  ggpubr::theme_pubr() + 
  labs(col='') +
  scale_color_brewer(palette="Set2") +
  labs(x="Separate annotations",
       y="Harmony integration",
       title="Cell type proportions in Healthy Liver") +
  theme(plot.title = element_text(hjust=0.5))

g2=ggplot(tab.mg,aes(x=Prop_sep,y=Prop_refmap,col=annot)) + 
  geom_point(size=rel(3)) + 
  geom_abline(slope=1,intercept=0,col='red',lty='dashed') + 
  ggpubr::theme_pubr() + 
  labs(col='') +
  scale_color_brewer(palette="Set2") +
  labs(x="Separate annotations",
       y="Pulpescu reference - label transfer",
       title="Cell type proportions in Healthy Liver") +
  theme(plot.title = element_text(hjust=0.5))

g3=ggplot(tab.mg,aes(x=Prop_sep,y=Prop_scvi_integ,col=annot)) + 
  geom_point(size=rel(3)) + 
  geom_abline(slope=1,intercept=0,col='red',lty='dashed') + 
  ggpubr::theme_pubr() + 
  labs(col='') +
  scale_color_brewer(palette="Set2") +
  labs(x="Separate annotations",
       y="scVI integration",
       title="Cell type proportions in Healthy Liver") +
  theme(plot.title = element_text(hjust=0.5))

g4=ggplot(tab.mg,aes(x=Prop_sep,y=Prop_T21refmap,col=annot)) + 
  geom_point(size=rel(3)) + 
  geom_abline(slope=1,intercept=0,col='red',lty='dashed') + 
  ggpubr::theme_pubr() + 
  labs(col='') +
  scale_color_brewer(palette="Set2") +
  labs(x="Separate annotations",
       y="Ts21 reference - label transfer",
       title="Cell type proportions in Healthy Liver") +
  theme(plot.title = element_text(hjust=0.5))

pdf("~/Downloads/prop_compare.pdf",width = 24,height=6)
cowplot::plot_grid(g1,g2,g3,g4,ncol=4)
dev.off()



