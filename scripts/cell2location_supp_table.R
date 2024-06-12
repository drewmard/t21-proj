library(data.table)

df1=fread("~/Downloads/abundance_section_celltype.disomy_to_disomy.csv",data.table = F,stringsAsFactors = F)
df2=fread("~/Downloads/abundance_section_celltype.t21_to_disomy.csv",data.table = F,stringsAsFactors = F)
df3=fread("~/Downloads/abundance_section_celltype.public_to_disomy.csv",data.table = F,stringsAsFactors = F)
df4=fread("~/Downloads/abundance_section_celltype.disomy_to_t21.csv",data.table = F,stringsAsFactors = F)
df5=fread("~/Downloads/abundance_section_celltype.t21_to_t21.csv",data.table = F,stringsAsFactors = F)

df.mg = merge(df1,df2,by=c('sample','cell type'),all=TRUE)

df3$`cell type`[df3$`cell type`=="B cell"] = "B cells"
df3$`cell type`[df3$`cell type`=="Early erythroid"] = "Early erythroid cells"
df3$`cell type`[df3$`cell type`=="HSC_MPP"] = "HSCs/MPPs"
df3$`cell type`[df3$`cell type`=="Hepatocyte"] = "Hepatocytes"
df3$`cell type`[df3$`cell type`=="Kupffer Cell"] = "Kupffer Cells"
df3$`cell type`[df3$`cell type`=="Late Erythroid"] = "Late erythroid cells"
df3$`cell type`[df3$`cell type`=="MEMP"] = "MEMPs"
df3$`cell type`[df3$`cell type`=="Megakaryocyte"] = "Megakaryocytes"
df3$`cell type`[df3$`cell type`=="NK"] = "NK cells"
df3$`cell type`[df3$`cell type`=="Pre pro B cell"] = "Pre pro B cells"
df3$`cell type`[df3$`cell type`=="pro-B cell"] = "Pro B cells"

df.mg = merge(df.mg,df3,by=c('sample','cell type'),all=TRUE)
df.mg$status = "Disomy"
df.mg2 = merge(df4,df5,by=c('sample','cell type'),all=TRUE)
df.mg2$`percent_ref-public` = NA
df.mg2$status = "Ts21"

df.mg3 = as.data.frame(rbind(df.mg,df.mg2))

df.mg3$ref_comp_disomy_v_t21 = df.mg3$`percent_ref-disomy` - df.mg3$`percent_ref-t21`
df.mg3$ref_comp_disomy_v_public = df.mg3$`percent_ref-disomy` - df.mg3$`percent_ref-public`
df.mg3$ref_comp_t21_v_public = df.mg3$`percent_ref-t21` - df.mg3$`percent_ref-public`

library(dplyr)
df.mg3 = df.mg3 %>%
  select(status,everything())

f.out = "~/Downloads/cell2location_supp_table.txt"
fwrite(df.mg3,f.out,quote = F,sep = '\t',row.names = F,col.names = T,na = "NA")





