library(data.table)
library(ggplot2)
x=2.4
per_fetus=FALSE

# 1: read in metadata that contains PCW

# 2: identify cell type composition of each sample

for_cellComp <- fread("~/Documents/Research/t21-proj/out/etc/for_cellComp.csv",data.table = F,stringsAsFactors = F)

if (per_fetus) {
  for_cellComp$`Patient ID` <- for_cellComp$Patient
}

convert_to_percentages <- function(counts,total_counts) {
  for (patient.uniq in unique(counts$Patient_ID)) {
    total_num_cells <- subset(total_counts,Var1==patient.uniq)$Freq
    counts$Freq[counts$Patient_ID==patient.uniq] <- counts$Freq[counts$Patient_ID==patient.uniq]/total_num_cells
  }
  return(counts)
}

convert_to_cluster_percentages <- function(counts,total_counts) {
  for (patient.uniq in unique(counts$Patient_ID)) {
    for (broad_cell in unique(counts$cell_type_groups)) {
      total_num_cells <- subset(total_counts,Patient_ID==patient.uniq & cell_type_groups==broad_cell)$Freq
      counts$Freq[counts$Patient_ID==patient.uniq & counts$cell_type_groups==broad_cell] <- counts$Freq[counts$Patient_ID==patient.uniq & counts$cell_type_groups==broad_cell]/total_num_cells
    }
  }
  return(counts)
}

total_counts <- list()
group_counts <- list()
group_percentages <- list()
group_matrix <- list()
cluster_counts <- list()
cluster_percentages <- list()
cluster_matrix <- list()
df <- list()
cluster_to_label_mapping <- list()
all_patients <- list()
for (sampletype in c("Femur","Liver")) { 
  cluster_counts[[sampletype]] <- list()
  cluster_percentages[[sampletype]] <- list()
  cluster_matrix[[sampletype]] <- list()
  
  group_counts[[sampletype]] <- list()
  group_percentages[[sampletype]] <- list()
  group_matrix[[sampletype]] <- list()
  
  total_counts[[sampletype]] <- list()
  df[[sampletype]] <- list()
  all_patients[[sampletype]] <- list()
  cluster_to_label_mapping[[sampletype]] <- list()
  for (disease_status in c("DownSyndrome","Healthy")) {
    cluster_counts[[sampletype]][[disease_status]] <- list()
    cluster_percentages[[sampletype]][[disease_status]] <- list()
    cluster_matrix[[sampletype]][[disease_status]] <- list()
    
    group_counts[[sampletype]][[disease_status]] <- list()
    group_percentages[[sampletype]][[disease_status]] <- list()
    group_matrix[[sampletype]][[disease_status]] <- list()
    
    total_counts[[sampletype]][[disease_status]] <- list()
    cluster_to_label_mapping[[sampletype]][[disease_status]] <- fread(paste0("~/Documents/Research/t21-proj/out/full/cluster_to_label_mapping/10X_",disease_status,"_",sampletype,".cluster_to_label_mapping.csv"),data.table = F,stringsAsFactors = F)[,c(3,4)]
    colnames(cluster_to_label_mapping[[sampletype]][[disease_status]])[1] <- 'cluster'
    
    df[[sampletype]][[disease_status]] <- fread(paste0("~/Documents/Research/t21-proj/out/full/cellComp/10X_",disease_status,"_",sampletype,".cellComp.csv"),data.table = F,stringsAsFactors = F)
    
    df[[sampletype]][[disease_status]]$Patient_ID <- paste(df[[sampletype]][[disease_status]]$patient,df[[sampletype]][[disease_status]]$sample,sep=" ")
    
    if (per_fetus) {
      df[[sampletype]][[disease_status]]$Patient_ID <- df[[sampletype]][[disease_status]]$patient
    }
    
    df[[sampletype]][[disease_status]] <- subset(df[[sampletype]][[disease_status]],Patient_ID %in% for_cellComp[,"Patient ID"])
    df[[sampletype]][[disease_status]] <- subset(df[[sampletype]][[disease_status]],cell_type_groups != "Unknown")
    all_patients[[sampletype]][[disease_status]] <- list()
    
    for (sorting_strategy in c("CD45+","CD235a-")) {
      
      if (sorting_strategy=="CD45+" & sampletype=="Femur") {next}
      
      all_patients[[sampletype]][[disease_status]][[sorting_strategy]] <- unique(subset(df[[sampletype]][[disease_status]],sorting==sorting_strategy)[,"Patient_ID"])
      
      # total_counts[[sampletype]][[disease_status]][[sorting_strategy]] <- as.data.frame(table(subset(df[[sampletype]][[disease_status]],sorting==sorting_strategy)[,c("patient")]))
      # group_counts[[sampletype]][[disease_status]][[sorting_strategy]] <- as.data.frame(table(subset(df[[sampletype]][[disease_status]],sorting==sorting_strategy)[,c("patient","cell_type_groups")]))
      # cluster_counts[[sampletype]][[disease_status]][[sorting_strategy]] <- as.data.frame(table(subset(df[[sampletype]][[disease_status]],sorting==sorting_strategy)[,c("patient",colnames(subset(df[[sampletype]][[disease_status]],sorting==sorting_strategy))[6])]))
      total_counts[[sampletype]][[disease_status]][[sorting_strategy]] <- as.data.frame(table(subset(df[[sampletype]][[disease_status]],sorting==sorting_strategy)[,c("Patient_ID")]))
      group_counts[[sampletype]][[disease_status]][[sorting_strategy]] <- as.data.frame(table(subset(df[[sampletype]][[disease_status]],sorting==sorting_strategy)[,c("Patient_ID","cell_type_groups")]))
      cluster_counts[[sampletype]][[disease_status]][[sorting_strategy]] <- as.data.frame(table(subset(df[[sampletype]][[disease_status]],sorting==sorting_strategy)[,c("Patient_ID",colnames(subset(df[[sampletype]][[disease_status]],sorting==sorting_strategy))[6])]))
      
      colnames(cluster_counts[[sampletype]][[disease_status]][[sorting_strategy]])[2] <- 'cluster'
      cluster_counts[[sampletype]][[disease_status]][[sorting_strategy]] <- merge(cluster_counts[[sampletype]][[disease_status]][[sorting_strategy]],cluster_to_label_mapping[[sampletype]][[disease_status]],by="cluster")
      # cluster_percentages[[sampletype]][[disease_status]][[sorting_strategy]] <- convert_to_cluster_percentages(cluster_counts[[sampletype]][[disease_status]][[sorting_strategy]],group_counts[[sampletype]][[disease_status]][[sorting_strategy]])
      cluster_percentages[[sampletype]][[disease_status]][[sorting_strategy]] <- convert_to_percentages(cluster_counts[[sampletype]][[disease_status]][[sorting_strategy]],total_counts[[sampletype]][[disease_status]][[sorting_strategy]])
      cluster_matrix[[sampletype]][[disease_status]][[sorting_strategy]] <- reshape2::dcast(cluster_percentages[[sampletype]][[disease_status]][[sorting_strategy]],Patient_ID~cluster,value.var="Freq",drop=FALSE)
      
      group_percentages[[sampletype]][[disease_status]][[sorting_strategy]] <- convert_to_percentages(group_counts[[sampletype]][[disease_status]][[sorting_strategy]],total_counts[[sampletype]][[disease_status]][[sorting_strategy]])
      group_matrix[[sampletype]][[disease_status]][[sorting_strategy]] <- reshape2::dcast(group_percentages[[sampletype]][[disease_status]][[sorting_strategy]],Patient_ID~cell_type_groups,value.var="Freq",drop=FALSE)
      
    }
  }
}

sampletype="Liver"
disease_status="Healthy"
sorting_strategy="CD45+"
df1=group_matrix[[sampletype]][[disease_status]][[sorting_strategy]]
# df1$disease_status = disease_status
disease_status="DownSyndrome"
df2=group_matrix[[sampletype]][[disease_status]][[sorting_strategy]]
# df2$disease_status = disease_status
df.mg = merge(for_cellComp,rbind(df1,df2),by.x='Patient ID',by.y="Patient_ID")
df.mg$AgeNum = as.numeric(substring(df.mg$Age,1,2))
cell_vec=colnames(df.mg)[7:14]
g=list();tmp=list()
gall = lapply(1:length(cell_vec),function(i) {
# for (i in 1:length(cell_vec)) {
  cell_type = cell_vec[i]
  tmp[[i]] = df.mg[,c("AgeNum","Age",cell_type)]
  g = ggplot(tmp[[i]],aes(x=as.factor(AgeNum),y=tmp[[i]][,cell_type],fill=Age)) + geom_boxplot() + geom_point() + geom_smooth(method='lm') + 
    theme_bw() + theme(panel.grid = element_blank()) +
    labs(y=paste0(cell_type," (%)"),x="Age (PCW)") +
    scale_fill_brewer(palette="Set2") +guides(fill="none")
  return(g)
}
)

for (i in 1:length(cell_vec)) {
  cell_type = cell_vec[i]
  tmp[[i]] = df.mg[,c("AgeNum","Age",cell_type)]
  # print(kruskal.test(tmp[[i]]$Age,tmp[[i]][,cell_type]))
  # print(wilcox.test(tmp[[i]][,cell_type][tmp[[i]]$AgeNum==12],
  #                   tmp[[i]][,cell_type][tmp[[i]]$AgeNum!=12]))
  print(cor.test(tmp[[i]]$AgeNum,tmp[[i]][,cell_type]))
}

library(cowplot)
pdf("~/Downloads/age_pcw.pdf",width = 11,height=5.5)
print(plot_grid(plotlist=gall,nrow=2,ncol=4))
dev.off()
# dggpubr::ggarrange(g)
# plot_grid(g,nrow = 2,ncol = 4)

