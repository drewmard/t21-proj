library(data.table)
library(ggplot2)

for_cellComp <- fread("~/Documents/Research/t21-proj/out/data_small/for_cellComp.csv",data.table = F,stringsAsFactors = F)

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
for (sampletype in c("Femur","Liver")) { 
  cluster_counts[[sampletype]] <- list()
  cluster_percentages[[sampletype]] <- list()
  cluster_matrix[[sampletype]] <- list()
  
  group_counts[[sampletype]] <- list()
  group_percentages[[sampletype]] <- list()
  group_matrix[[sampletype]] <- list()
  
  total_counts[[sampletype]] <- list()
  df[[sampletype]] <- list()
  cluster_to_label_mapping[[sampletype]] <- list()
  for (disease_status in c("DownSyndrome","Healthy")) {
    cluster_counts[[sampletype]][[disease_status]] <- list()
    cluster_percentages[[sampletype]][[disease_status]] <- list()
    cluster_matrix[[sampletype]][[disease_status]] <- list()
    
    group_counts[[sampletype]][[disease_status]] <- list()
    group_percentages[[sampletype]][[disease_status]] <- list()
    group_matrix[[sampletype]][[disease_status]] <- list()
    
    total_counts[[sampletype]][[disease_status]] <- list()
    cluster_to_label_mapping[[sampletype]][[disease_status]] <- fread(paste0("~/Documents/Research/t21-proj/out/data_small/10X_",disease_status,"_",sampletype,".cluster_to_label_mapping.csv"),data.table = F,stringsAsFactors = F)[,c(3,4)]
    colnames(cluster_to_label_mapping[[sampletype]][[disease_status]])[1] <- 'cluster'
    
    df[[sampletype]][[disease_status]] <- fread(paste0("~/Documents/Research/t21-proj/out/data_small/10X_",disease_status,"_",sampletype,".cellComp.csv"),data.table = F,stringsAsFactors = F)
    df[[sampletype]][[disease_status]]$Patient_ID <- paste(df[[sampletype]][[disease_status]]$patient,df[[sampletype]][[disease_status]]$sample,sep=" ")
    df[[sampletype]][[disease_status]] <- subset(df[[sampletype]][[disease_status]],Patient_ID %in% for_cellComp[,"Patient ID"])

    for (sorting_strategy in c("CD45+","CD235a-")) {
      
      if (sorting_strategy=="CD45+" & sampletype=="Femur") {next}
      total_counts[[sampletype]][[disease_status]][[sorting_strategy]] <- as.data.frame(table(subset(df[[sampletype]][[disease_status]],sorting==sorting_strategy)[,c("Patient_ID")]))
      group_counts[[sampletype]][[disease_status]][[sorting_strategy]] <- as.data.frame(table(subset(df[[sampletype]][[disease_status]],sorting==sorting_strategy)[,c("Patient_ID","cell_type_groups")]))
      cluster_counts[[sampletype]][[disease_status]][[sorting_strategy]] <- as.data.frame(table(subset(df[[sampletype]][[disease_status]],sorting==sorting_strategy)[,c("Patient_ID",colnames(subset(df[[sampletype]][[disease_status]],sorting==sorting_strategy))[6])]))
      
      colnames(cluster_counts[[sampletype]][[disease_status]][[sorting_strategy]])[2] <- 'cluster'
      cluster_counts[[sampletype]][[disease_status]][[sorting_strategy]] <- merge(cluster_counts[[sampletype]][[disease_status]][[sorting_strategy]],cluster_to_label_mapping[[sampletype]][[disease_status]],by="cluster")
      cluster_percentages[[sampletype]][[disease_status]][[sorting_strategy]] <- convert_to_cluster_percentages(cluster_counts[[sampletype]][[disease_status]][[sorting_strategy]],group_counts[[sampletype]][[disease_status]][[sorting_strategy]])
      cluster_matrix[[sampletype]][[disease_status]][[sorting_strategy]] <- dcast(cluster_percentages[[sampletype]][[disease_status]][[sorting_strategy]],Patient_ID~cluster,value.var="Freq")
      
      group_percentages[[sampletype]][[disease_status]][[sorting_strategy]] <- convert_to_percentages(group_counts[[sampletype]][[disease_status]][[sorting_strategy]],total_counts[[sampletype]][[disease_status]][[sorting_strategy]])
      group_matrix[[sampletype]][[disease_status]][[sorting_strategy]] <- dcast(group_percentages[[sampletype]][[disease_status]][[sorting_strategy]],Patient_ID~cell_type_groups,value.var="Freq")
      
    }
  }
}

res <- list()
for (sampletype in c("Femur","Liver")) { 
  res[[sampletype]] <- list()
  for (sorting_strategy in c("CD45+","CD235a-")) {
    if (sorting_strategy=="CD45+" & sampletype=="Femur") {next}
    res.lst <- list(); iter=0; for (celltype in colnames(group_matrix[[sampletype]][[disease_status]][[sorting_strategy]])[-1]) { # remove Patient_ID
    # res.lst <- list(); iter=0; for (celltype in colnames(cluster_matrix[[sampletype]][[disease_status]][[sorting_strategy]])[-1]) { # remove Patient_ID
      iter=iter+1
      # if (celltype %in% colnames(cluster_matrix[[sampletype]][["DownSyndrome"]][[sorting_strategy]]) &
      #     celltype %in% colnames(cluster_matrix[[sampletype]][["Healthy"]][[sorting_strategy]])) {
      tryCatch({
        # tres <- t.test(cluster_matrix[[sampletype]][["DownSyndrome"]][[sorting_strategy]][,celltype],cluster_matrix[[sampletype]][["Healthy"]][[sorting_strategy]][,celltype])
        tres <- t.test(group_matrix[[sampletype]][["DownSyndrome"]][[sorting_strategy]][,celltype],group_matrix[[sampletype]][["Healthy"]][[sorting_strategy]][,celltype])
        ures <- wilcox.test(group_matrix[[sampletype]][["DownSyndrome"]][[sorting_strategy]][,celltype],group_matrix[[sampletype]][["Healthy"]][[sorting_strategy]][,celltype])
        res.lst[[iter]] <- data.frame(sampletype,sorting_strategy,celltype,cases=as.numeric(tres$estimate[1]),controls=as.numeric(tres$estimate[2]),t_pval=as.numeric(tres$p.value),mwu_pval=ures$p.value)
      },error=function(e) {print(paste(sampletype,sorting_strategy,celltype))})
    }
    res[[sampletype]][[sorting_strategy]] <- as.data.frame(do.call(rbind,res.lst))
    res[[sampletype]][[sorting_strategy]]$t_pval.adj <- p.adjust(res[[sampletype]][[sorting_strategy]]$t_pval,method='fdr')
    res[[sampletype]][[sorting_strategy]]$mwu_pval.adj <- p.adjust(res[[sampletype]][[sorting_strategy]]$mwu_pval,method='fdr')
    f.out <- paste0("~/Documents/Research/t21-proj/out/results_etc/cellComp.group.",sampletype,'.',gsub("[[:punct:]]","",sorting_strategy),'.txt')
    fwrite(res[[sampletype]][[sorting_strategy]],file = f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
  }
}

group_matrix.melt <- list()
for (sampletype in c("Femur","Liver")) { 
  group_matrix.melt[[sampletype]] <- list()
  for (sorting_strategy in c("CD45+","CD235a-")) {
    if (sorting_strategy=="CD45+" & sampletype=="Femur") {next}
    
    group_matrix[[sampletype]][["DownSyndrome"]][[sorting_strategy]]$disease_status <- "DownSyndrome"
    group_matrix[[sampletype]][["Healthy"]][[sorting_strategy]]$disease_status <- "Healthy"
    if (sampletype=="Femur" & sorting_strategy=="CD235a-") {colnames(group_matrix[[sampletype]][["Healthy"]][[sorting_strategy]])[colnames(group_matrix[[sampletype]][["Healthy"]][[sorting_strategy]])=="NK cells"] <- "NK/T cells"}
    group_matrix[[sampletype]][[sorting_strategy]] <- as.data.frame(rbind(group_matrix[[sampletype]][["DownSyndrome"]][[sorting_strategy]],group_matrix[[sampletype]][["Healthy"]][[sorting_strategy]]))
    
    group_matrix.melt[[sampletype]][[sorting_strategy]] <- melt(group_matrix[[sampletype]][[sorting_strategy]],id.vars = c("Patient_ID","disease_status"))
    
    system("mkdir -p ~/Documents/Research/t21-proj/out/figures/cellComp/")
    f.out=paste0("~/Documents/Research/t21-proj/out/figures/cellComp/boxplot.group.",sampletype,'.',gsub("[[:punct:]]","",sorting_strategy),'.png')
    png(filename = f.out,width = 4000,height=1800,res=480)
    g1 <- ggplot(group_matrix.melt[[sampletype]][[sorting_strategy]],aes(x=variable,y=100*value,fill=disease_status)) + 
      theme_bw() + 
      theme(panel.grid=element_blank(),
            axis.text.x =element_text(angle=30,hjust=1),
            plot.title = element_text(hjust=0.5)) +
      geom_boxplot(outlier.shape=NA) + 
      geom_point(position=position_jitterdodge(jitter.width=0.075),alpha=0.5) + 
      scale_fill_brewer(palette = "Set3",name="Disease Status",labels=c("Down Syndrome","Healthy")) +
      labs(x="Cell type",y="Cell composition (%)",title=paste0(sorting_strategy))
    print(g1)
    dev.off()
    
    system("mkdir -p ~/Documents/Research/t21-proj/out/figures/cellComp/")
    f.out=paste0("~/Documents/Research/t21-proj/out/figures/cellComp/stacked.group.",sampletype,'.',gsub("[[:punct:]]","",sorting_strategy),'.png')
    png(filename = f.out,width = 3300,height=2600,res=480)
    g2<-ggplot(group_matrix.melt[[sampletype]][[sorting_strategy]], aes(fill=variable, y=value*100, x=Patient_ID,col=disease_status)) + 
      geom_bar(position="stack", stat="identity") +
      scale_color_manual(values=c("black","orange"),labels=c("Down Syndrome","Healthy"),name="Disease Status") +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text.x=element_text(hjust=1,angle=60),
            plot.title = element_text(hjust=0.5)) +
      scale_fill_brewer(palette="Set3",name="Cell type") + 
      labs(x="Patient ID",y="Cell composition (%)",title=sorting_strategy)
    print(g2)
    dev.off()
  }
}

sampletype="Liver"
cluster_percentages_all <- list()
cluster_percentages_all[[sampletype]] <- list()
for (sorting_strategy in c("CD45+","CD235a-")) {
  # sorting_strategy="CD45+"
  cluster_percentages[[sampletype]][["DownSyndrome"]][[sorting_strategy]]$disease_status <- "DownSyndrome"
  cluster_percentages[[sampletype]][["Healthy"]][[sorting_strategy]]$disease_status <- "Healthy"
  
  cluster_percentages_all[[sampletype]][[sorting_strategy]] <- list()
  cell_type_groups <- unique(cluster_percentages[[sampletype]][["Healthy"]][[sorting_strategy]]$cell_type_groups)
  for (cell_type in cell_type_groups) {
    mat1 <- subset(cluster_percentages[[sampletype]][["DownSyndrome"]][[sorting_strategy]],cell_type_groups==cell_type)
    mat2 <- subset(cluster_percentages[[sampletype]][["Healthy"]][[sorting_strategy]],cell_type_groups==cell_type)
    cluster_percentages_all[[sampletype]][[sorting_strategy]][[cell_type]] <- as.data.frame(rbind(mat1,mat2))
    # cluster_matrix_all[[sampletype]][[sorting_strategy]][[cell_type]] <- dcast(group_percentages[[sampletype]][[disease_status]][[sorting_strategy]],Patient_ID~cell_type_groups,value.var="Freq")
    
    system("mkdir -p ~/Documents/Research/t21-proj/out/figures/cellComp/")
    f.out=paste0("~/Documents/Research/t21-proj/out/figures/cellComp/boxplot.cluster.",sampletype,'.',gsub("[[:punct:]]","",sorting_strategy),'.',gsub('\\/','_',cell_type),'.png')
    png(filename = f.out,width = 4000,height=1800,res=480)
    g1 <- ggplot(cluster_percentages_all[[sampletype]][[sorting_strategy]][[cell_type]],aes(x=cluster,y=100*Freq,fill=disease_status)) + 
      theme_bw() + 
      theme(panel.grid=element_blank(),
            axis.text.x =element_text(angle=60,hjust=1),
            plot.title = element_text(hjust=0.5)) +
      geom_boxplot(outlier.shape=NA) + 
      geom_point(position=position_jitterdodge(jitter.width=0.075),alpha=0.5) + 
      scale_fill_brewer(palette = "Set3",name="Disease Status",labels=c("Down Syndrome","Healthy")) +
      labs(x="Cell type",y="Cell composition (%)",title=paste0(sorting_strategy))
    print(g1)
    dev.off()
    
    system("mkdir -p ~/Documents/Research/t21-proj/out/figures/cellComp/")
    f.out=paste0("~/Documents/Research/t21-proj/out/figures/cellComp/stacked.cluster.",sampletype,'.',gsub("[[:punct:]]","",sorting_strategy),'.',gsub('\\/','_',cell_type),'.png')
    png(filename = f.out,width = 3300,height=2600,res=480)
    g2<-ggplot(cluster_percentages_all[[sampletype]][[sorting_strategy]][[cell_type]], aes(fill=cluster, y=Freq*100, x=Patient_ID,col=disease_status)) + 
      geom_bar(position="stack", stat="identity") +
      scale_color_manual(values=c("black","orange"),labels=c("Down Syndrome","Healthy"),name="Disease Status") +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text.x=element_text(hjust=1,angle=60),
            plot.title = element_text(hjust=0.5)) +
      scale_fill_brewer(palette="Set3",name="Cell type") + 
      labs(x="Patient ID",y="Cell composition (%)",title=sorting_strategy)
    print(g2)
    dev.off()
    
  }
}


sampletype="Liver"
sorting_strategy="CD45+"
for (cell_type in c("B cells","HSC/Progenitors","Myeloid")) {
  f.out=paste0("~/Documents/Research/t21-proj/out/figures/cellComp/boxplot.cluster.",sampletype,'.',gsub("[[:punct:]]","",sorting_strategy),'.',gsub('\\/','_',cell_type),'.pdf')
  x=2.4
  pdf(file = f.out,width = 4*x,height=1.8*x)
  g1 <- ggplot(cluster_percentages_all[[sampletype]][[sorting_strategy]][[cell_type]],aes(x=cluster,y=100*Freq,fill=disease_status)) + 
    theme_bw() + 
    theme(panel.grid=element_blank(),
          axis.text.x =element_text(angle=60,hjust=1),
          plot.title = element_text(hjust=0.5)) +
    geom_boxplot(outlier.shape=NA) + 
    geom_point(position=position_jitterdodge(jitter.width=0.075),alpha=0.5) + 
    scale_fill_brewer(palette = "Set3",name="Disease Status",labels=c("Down Syndrome","Healthy")) +
    labs(x="Cell type",y="Cell composition (%)",title=paste0(sorting_strategy))
  print(g1)
  dev.off()
}

sampletype="Liver"
sorting_strategy="CD45+"
f.out=paste0("~/Documents/Research/t21-proj/out/figures/cellComp/boxplot.group.",sampletype,'.',gsub("[[:punct:]]","",sorting_strategy),'.pdf')
x=2.4
pdf(file = f.out,width = 4*x,height=1.8*x)
g1 <- ggplot(group_matrix.melt[[sampletype]][[sorting_strategy]],aes(x=variable,y=100*value,fill=disease_status)) + 
  theme_bw() + 
  theme(panel.grid=element_blank(),
        axis.text.x =element_text(angle=30,hjust=1),
        plot.title = element_text(hjust=0.5)) +
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitterdodge(jitter.width=0.075),alpha=0.5) + 
  scale_fill_brewer(palette = "Set3",name="Disease Status",labels=c("Down Syndrome","Healthy")) +
  labs(x="Cell type",y="Cell composition (%)",title=paste0(sorting_strategy))
print(g1)
dev.off()





# mat1=cluster_matrix[[sampletype]][["DownSyndrome"]][[sorting_strategy]]
# mat2=cluster_matrix[[sampletype]][["Healthy"]][[sorting_strategy]]
# 
# cluster_matrix[[sampletype]][[sorting_strategy]] <- as.data.frame(rbind(mat1,cluster_matrix[[sampletype]][["Healthy"]][[sorting_strategy]]))
# 
# cluster_matrix.melt <- list()
# cluster_matrix.melt[[sampletype]] <- list()
# cluster_matrix.melt[[sampletype]][[sorting_strategy]] <- melt(cluster_matrix[[sampletype]][[sorting_strategy]],id.vars = c("Patient_ID","disease_status"))
# 
# 
# 
# f.out=paste0("~/Documents/Research/t21-proj/out/figures/cellComp/cluster.",sampletype,'.',gsub("[[:punct:]]","",sorting_strategy),'.png')
# png(filename = f.out,width = 4000,height=1800,res=480)
# ggplot(group_matrix.melt[[sampletype]][[sorting_strategy]],aes(x=variable,y=100*value,fill=disease_status)) + 
#   theme_bw() + 
#   theme(panel.grid=element_blank(),
#         axis.text.x =element_text(angle=30,hjust=1),
#         plot.title = element_text(hjust=0.5)) +
#   geom_boxplot(outlier.shape=NA) + 
#   geom_point(position=position_jitterdodge(jitter.width=0.075),alpha=0.5) + 
#   scale_fill_brewer(palette = "Set3",name="Disease Status",labels=c("Down Syndrome","Healthy")) +
#   labs(x="Cell type",y="Cell composition (%)",title=paste0(sorting_strategy))
# dev.off()
# 
