library(data.table)
refmap = fread("~/Downloads/scArches_Healthy_Liver_metadata.csv",data.table = F,stringsAsFactors = F)
refmap$cellname=unlist(lapply(strsplit(refmap$barcode,"-"),function(x) paste(x[1:2],collapse = "-")))
refmap$patient_sample = paste0(refmap$patient," ",refmap$sample)
refmap_broad = fread("~/Downloads/ref_broad_map.csv",data.table = F,stringsAsFactors = F)
refmap = merge(refmap,refmap_broad,by="Predicted")

dataf = fread("~/Downloads/a3_tab.obs.withannot.csv",data.table = F,stringsAsFactors = F)

disease_status="Healthy"
sampletype="Liver"
f=paste0("~/Documents/Research/t21-proj/out/full/cellComp/10X_",disease_status,"_",sampletype,".cellComp.csv")
df=fread(f,data.table = F,stringsAsFactors = F)
disease_status="Healthy"
sampletype="Liver"
colnames(df)[6] = "leiden_latest"
df$disease=disease_status
df$env = sampletype
df.uniq = unique(df[,c("leiden_latest","cell_type_groups")])
df$patient_sample = paste0(df$patient," ",df$sample)

dataf2 = merge(dataf,df.uniq,by="leiden_latest",all.x=TRUE)
dataf2$cellname=unlist(lapply(strsplit(dataf2$barcodekey,"-"),function(x) paste(x[1:2],collapse = "-")))
df$cellname=unlist(lapply(strsplit(df$V1,"-"),function(x) paste(x[1:2],collapse = "-")))
dataf3 = merge(dataf2,df[,c("cellname","patient_sample","sorting")],by=c("cellname","patient_sample"),all.x="TRUE")

library(data.table)
new_annot_mapping = fread("~/Downloads/new_annot_mapping.csv",data.table = F,stringsAsFactors = F)
dataf3 = merge(dataf3,new_annot_mapping,by="combi_annot")
dataf3.old = dataf3
rm(dataf2)
rm(df)
rm(dataf)

dataf4=merge(dataf3,refmap,by=c("patient_sample","cellname"))


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

tab.mg = merge(merge(tab,tab2,by="annot"),tab3,by='annot')

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
       y="Reference query label transfer",
       title="Cell type proportions in Healthy Liver") +
  theme(plot.title = element_text(hjust=0.5))

g3=ggplot(tab.mg,aes(x=Prop_harm_integ,y=Prop_refmap,col=annot)) + 
  geom_point(size=rel(3)) + 
  geom_abline(slope=1,intercept=0,col='red',lty='dashed') + 
  ggpubr::theme_pubr() + 
  labs(col='') +
  scale_color_brewer(palette="Set2") +
  labs(x="Harmony integration",
       y="Reference query label transfer",
       title="Cell type proportions in Healthy Liver") +
  theme(plot.title = element_text(hjust=0.5))

cowplot::plot_grid(g1,g2,g3,ncol=3)

