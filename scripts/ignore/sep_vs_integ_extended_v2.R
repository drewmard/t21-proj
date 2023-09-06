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



