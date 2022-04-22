df.ds <- fread("~/Documents/Research/t21-proj/out/data_small/10X_DownSyndrome_Liver.cellComp.csv",data.table = F,stringsAsFactors = F)
df.healthy <- fread("~/Documents/Research/t21-proj/out/data_small/10X_Healthy_Liver.cellComp.csv",data.table = F,stringsAsFactors = F)
for_cellComp <- fread("~/Documents/Research/t21-proj/out/data_small/for_cellComp.csv",data.table = F,stringsAsFactors = F)

head(for_cellComp)
subset(for_cellComp,Patient=="15532")
subset(df.ds,patient=="15532")[1,]
subset(for_cellComp,Patient=="15582")
unique(subset(df.ds,patient=="15582")[,2:4])

df.ds$Patient_ID <- paste(df.ds$patient,df.ds$sample,sep=" ")
df.healthy$Patient_ID <- paste(df.healthy$patient,df.healthy$sample,sep=" ")

df.ds.sub <- subset(df.ds,Patient_ID %in% for_cellComp[,"Patient ID"])
df.healthy.sub <- subset(df.healthy,Patient_ID %in% for_cellComp[,"Patient ID"])

ds.cd45.counts <- as.data.frame(table(subset(df.ds.sub,sorting=="CD45+")[,c("Patient_ID","cell_type_groups")]))
ds.cd45.total_counts <- as.data.frame(table(subset(df.ds.sub,sorting=="CD45+")[,c("Patient_ID")]))
healthy.cd45.counts <- as.data.frame(table(subset(df.healthy.sub,sorting=="CD45+")[,c("Patient_ID","cell_type_groups")]))
healthy.cd45.total_counts <- as.data.frame(table(subset(df.healthy.sub,sorting=="CD45+")[,c("Patient_ID")]))

ds.cd235a.counts <- as.data.frame(table(subset(df.ds.sub,sorting=="CD235a-")[,c("Patient_ID","cell_type_groups")]))
ds.cd235a.total_counts <- as.data.frame(table(subset(df.ds.sub,sorting=="CD235a-")[,c("Patient_ID")]))
healthy.cd235a.counts <- as.data.frame(table(subset(df.healthy.sub,sorting=="CD235a-")[,c("Patient_ID","cell_type_groups")]))
healthy.cd235a.total_counts <- as.data.frame(table(subset(df.healthy.sub,sorting=="CD235a-")[,c("Patient_ID")]))

convert_to_percentages <- function(counts,total_counts) {
  for (patient.uniq in unique(counts$Patient_ID)) {
    total_num_cells <- subset(total_counts,Var1==patient.uniq)$Freq
    counts$Freq[counts$Patient_ID==patient.uniq] <- counts$Freq[counts$Patient_ID==patient.uniq]/total_num_cells
  }
  return(counts)
}

ds.cd45.counts <- convert_to_percentages(ds.cd45.counts,ds.cd45.total_counts)
healthy.cd45.counts <- convert_to_percentages(healthy.cd45.counts,healthy.cd45.total_counts)
ds.cd235a.counts <- convert_to_percentages(ds.cd235a.counts,ds.cd235a.total_counts)
healthy.cd235a.counts <- convert_to_percentages(healthy.cd235a.counts,healthy.cd235a.total_counts)

library(reshape2)
ds.cd45 <- dcast(ds.cd45.counts,Patient_ID~cell_type_groups,value.var="Freq")
healthy.cd45 <- dcast(healthy.cd45.counts,Patient_ID~cell_type_groups,value.var="Freq")
ds.cd235a <- dcast(ds.cd235a.counts,Patient_ID~cell_type_groups,value.var="Freq")
healthy.cd235a <- dcast(healthy.cd235a.counts,Patient_ID~cell_type_groups,value.var="Freq")


res.lst <- list(); iter=0; for (sorting in c("cd45","cd235a")) {
  if (sorting=="cd45") {
    cases <- ds.cd45
    controls <- healthy.cd45
  } else {
    cases <- ds.cd235a
    controls <- healthy.cd235a
  } 
  for (celltype in colnames(cases)[-1]) { # remove Patient_ID
    iter=iter+1
    tres <- t.test(cases[,celltype],controls[,celltype])
    res.lst[[iter]] <- data.frame(sorting,celltype,cases=as.numeric(tres$estimate[1]),controls=as.numeric(tres$estimate[2]),pval=as.numeric(tres$p.value))
  }
}
as.data.frame(do.call(rbind,res.lst))

ds.cd45$disease_status <- "DownSyndrome"
healthy.cd45$disease_status <- "Healthy"
cd45 <- as.data.frame(rbind(ds.cd45,healthy.cd45))
ggplot(cd45,aes(x=disease_status,y=Erythroid,fill=disease_status)) + 
  theme_bw() + theme(panel.grid=element_blank()) +
  geom_boxplot(outlier.shape=NA) + geom_jitter(width = 0.05) + scale_fill_brewer(palette = "Set3")

cd45.melt <- melt(cd45,id.vars = c("Patient_ID","disease_status"))
g1 <- ggplot(cd45.melt,aes(x=variable,y=100*value,fill=disease_status)) + 
  theme_bw() + 
  theme(panel.grid=element_blank(),
        axis.text.x =element_text(angle=30,hjust=1),
        plot.title = element_text(hjust=0.5)) +
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitterdodge(jitter.width=0.075),alpha=0.5) + 
  scale_fill_brewer(palette = "Set3",name="Disease Status",labels=c("Down Syndrome","Healthy")) +
  labs(x="Cell type",y="Cell composition (%)",title="CD45+ sorting")

ds.cd235a$disease_status <- "DownSyndrome"
healthy.cd235a$disease_status <- "Healthy"
cd235a <- as.data.frame(rbind(ds.cd235a,healthy.cd235a))
cd235a.melt <- melt(cd235a,id.vars = c("Patient_ID","disease_status"))
g2 <- ggplot(cd235a.melt,aes(x=variable,y=100*value,fill=disease_status)) + 
  theme_bw() + 
  theme(panel.grid=element_blank(),
        axis.text.x =element_text(angle=30,hjust=1),
        plot.title = element_text(hjust=0.5)) +
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitterdodge(jitter.width=0.075),alpha=0.5) + 
  scale_fill_brewer(palette = "Set3",name="Disease Status",labels=c("Down Syndrome","Healthy")) +
  labs(x="Cell type",y="Cell composition (%)",title="CD235a- sorting")

library(cowplot)
plot_grid(g1,g2,ncol=2)




