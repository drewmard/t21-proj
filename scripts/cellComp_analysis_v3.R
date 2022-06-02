library(data.table)
library(ggplot2)
x=2.4
per_fetus=FALSE

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
      cluster_percentages[[sampletype]][[disease_status]][[sorting_strategy]] <- convert_to_cluster_percentages(cluster_counts[[sampletype]][[disease_status]][[sorting_strategy]],group_counts[[sampletype]][[disease_status]][[sorting_strategy]])
      cluster_matrix[[sampletype]][[disease_status]][[sorting_strategy]] <- reshape2::dcast(cluster_percentages[[sampletype]][[disease_status]][[sorting_strategy]],Patient_ID~cluster,value.var="Freq",drop=FALSE)
      
      group_percentages[[sampletype]][[disease_status]][[sorting_strategy]] <- convert_to_percentages(group_counts[[sampletype]][[disease_status]][[sorting_strategy]],total_counts[[sampletype]][[disease_status]][[sorting_strategy]])
      group_matrix[[sampletype]][[disease_status]][[sorting_strategy]] <- reshape2::dcast(group_percentages[[sampletype]][[disease_status]][[sorting_strategy]],Patient_ID~cell_type_groups,value.var="Freq",drop=FALSE)
      
    }
  }
}

res <- list()
system(paste0("mkdir -p ~/Documents/Research/t21-proj/out/full/cellComp_res"))
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
    f.out <- paste0("~/Documents/Research/t21-proj/out/full/cellComp_res/cellComp.group.",sampletype,'.',gsub("[[:punct:]]","",sorting_strategy),'.txt')
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
    
    group_matrix.melt[[sampletype]][[sorting_strategy]] <- reshape2::melt(group_matrix[[sampletype]][[sorting_strategy]],id.vars = c("Patient_ID","disease_status"))
    
    system("mkdir -p ~/Documents/Research/t21-proj/out/full/cellComp_fig/")
    f.out=paste0("~/Documents/Research/t21-proj/out/full/cellComp_fig/boxplot.group.",sampletype,'.',gsub("[[:punct:]]","",sorting_strategy),'.pdf')
    # png(filename = f.out,width = 4000,height=1800,res=480)
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
    
    f.out=paste0("~/Documents/Research/t21-proj/out/full/cellComp_fig/stacked.group.",sampletype,'.',gsub("[[:punct:]]","",sorting_strategy),'.pdf')
    # png(filename = f.out,width = 3300,height=2600,res=480)
    pdf(file = f.out,width = 3.3*x,height=2.6*x)
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

cluster_percentages_all <- list()
sampletype="Femur"
sorting_strategy="CD235a-"
for (sampletype in c("Femur","Liver")) { 
  cluster_percentages_all[[sampletype]] <- list()
  for (sorting_strategy in c("CD235a-","CD45+")) {
    
    if (sorting_strategy=="CD45+" & sampletype=="Femur") {next}
    
    # sorting_strategy="CD45+"
    cluster_percentages[[sampletype]][["DownSyndrome"]][[sorting_strategy]]$disease_status <- "DownSyndrome"
    cluster_percentages[[sampletype]][["Healthy"]][[sorting_strategy]]$disease_status <- "Healthy"
    
    cluster_percentages_all[[sampletype]][[sorting_strategy]] <- list()
    cell_type_groups <- unique(cluster_percentages[[sampletype]][["Healthy"]][[sorting_strategy]]$cell_type_groups)
    for (cell_type in cell_type_groups) {
      cell_type_groups <- unique(cluster_percentages[[sampletype]][["Healthy"]][[sorting_strategy]]$cell_type_groups)
      mat1 <- subset(cluster_percentages[[sampletype]][["DownSyndrome"]][[sorting_strategy]],cell_type_groups==cell_type)
      mat2 <- subset(cluster_percentages[[sampletype]][["Healthy"]][[sorting_strategy]],cell_type_groups==cell_type)
      
      cell_missing <- unique(mat2$cluster)[!(unique(mat2$cluster) %in% unique(mat1$cluster))]
      for (cell in cell_missing) {
        tmp <- subset(mat2,cluster==cell)[1,]
        tmp <- tmp[rep(1,length(all_patients[[sampletype]][["DownSyndrome"]][[sorting_strategy]])),]
        tmp$Patient_ID <- all_patients[[sampletype]][["DownSyndrome"]][[sorting_strategy]]
        tmp$disease_status <- "DownSyndrome"
        tmp$Freq <- 0
        mat1 <- as.data.frame(rbind(mat1,tmp))
      }
      cell_missing <- unique(mat1$cluster)[!(unique(mat1$cluster) %in% unique(mat2$cluster))]
      for (cell in cell_missing) {
        tmp <- subset(mat1,cluster==cell)[1,]
        tmp <- tmp[rep(1,length(all_patients[[sampletype]][["Healthy"]][[sorting_strategy]])),]
        tmp$Patient_ID <- all_patients[[sampletype]][["Healthy"]][[sorting_strategy]]
        tmp$disease_status <- "Healthy"
        tmp$Freq <- 0
        mat2 <- as.data.frame(rbind(mat2,tmp))
      }
      cluster_percentages_all[[sampletype]][[sorting_strategy]][[cell_type]] <- as.data.frame(rbind(mat1,mat2))
      
      # f.out=paste0("~/Documents/Research/t21-proj/out/full/cellComp_fig/boxplot.cluster.",sampletype,'.',gsub("[[:punct:]]","",sorting_strategy),'.',gsub('\\/','_',cell_type),'.pdf')
      # # png(filename = f.out,width = 4000,height=1800,res=480)
      # pdf(file = f.out,width = 4*x,height=1.8*x)
      # g1 <- ggplot(cluster_percentages_all[[sampletype]][[sorting_strategy]][[cell_type]],aes(x=cluster,y=100*Freq,fill=disease_status)) + 
      #   theme_bw() + 
      #   theme(panel.grid=element_blank(),
      #         axis.text.x =element_text(angle=60,hjust=1),
      #         plot.title = element_text(hjust=0.5)) +
      #   geom_boxplot(outlier.shape=NA) + 
      #   geom_point(position=position_jitterdodge(jitter.width=0.075),alpha=0.5) + 
      #   scale_fill_brewer(palette = "Set3",name="Disease Status",labels=c("Down Syndrome","Healthy")) +
      #   labs(x="Cell type",y="Cell composition (%)",title=paste0(sorting_strategy))
      # print(g1)
      # dev.off()
      # 
      # f.out=paste0("~/Documents/Research/t21-proj/out/full/cellComp_fig/stacked.cluster.",sampletype,'.',gsub("[[:punct:]]","",sorting_strategy),'.',gsub('\\/','_',cell_type),'.pdf')
      # # png(filename = f.out,width = 3300,height=2600,res=480)
      # pdf(file = f.out,width = 3.3*x,height=2.6*x)
      # g2<-ggplot(cluster_percentages_all[[sampletype]][[sorting_strategy]][[cell_type]], aes(fill=cluster, y=Freq*100, x=Patient_ID,col=disease_status)) + 
      #   geom_bar(position="stack", stat="identity") +
      #   scale_color_manual(values=c("black","orange"),labels=c("Down Syndrome","Healthy"),name="Disease Status") +
      #   theme_bw() +
      #   theme(panel.grid = element_blank(),
      #         axis.text.x=element_text(hjust=1,angle=60),
      #         plot.title = element_text(hjust=0.5)) +
      #   scale_fill_brewer(palette="Set3",name="Cell type") + 
      #   labs(x="Patient ID",y="Cell composition (%)",title=sorting_strategy)
      # print(g2)
      # dev.off()
      
    }
  }
}


sampletype="Liver"
sorting_strategy="CD45+"
for (sampletype in c("Liver","Femur")) { 
  for (sorting_strategy in c("CD235a-","CD45+")) {
    res.lst <- list()
    iter=0; for (cell_type in cell_type_groups) {
      
      # for (cell_type in c("B cells","HSC/Progenitors","Myeloid")) {
      # f.out=paste0("~/Documents/Research/t21-proj/out/full/cellComp_fig/boxplot.cluster.",sampletype,'.',gsub("[[:punct:]]","",sorting_strategy),'.',gsub('\\/','_',cell_type),'.pdf')
      # x=2.4
      # pdf(file = f.out,width = 4*x,height=1.8*x)
      # g1 <- ggplot(cluster_percentages_all[[sampletype]][[sorting_strategy]][[cell_type]],aes(x=cluster,y=100*Freq,fill=disease_status)) + 
      #   theme_bw() + 
      #   theme(panel.grid=element_blank(),
      #         axis.text.x =element_text(angle=60,hjust=1),
      #         plot.title = element_text(hjust=0.5)) +
      #   geom_boxplot(outlier.shape=NA) + 
      #   geom_point(position=position_jitterdodge(jitter.width=0.075),alpha=0.5) + 
      #   scale_fill_brewer(palette = "Set3",name="Disease Status",labels=c("Down Syndrome","Healthy")) +
      #   labs(x="Cell type",y="Cell composition (%)",title=paste0(sorting_strategy))
      # print(g1)
      # dev.off()
      
      cluster.uniq=unique(cluster_percentages_all[[sampletype]][[sorting_strategy]][[cell_type]]$cluster)
      for (clustertype in cluster.uniq) {
        iter=iter+1
        tryCatch({
          y=subset(cluster_percentages_all[[sampletype]][[sorting_strategy]][[cell_type]],cluster==clustertype)
          tres <- t.test(y$Freq~y$disease_status)
          ures <- wilcox.test(y$Freq~y$disease_status)
          res.lst[[iter]] <- data.frame(sampletype,sorting_strategy,cell_type,clustertype,cases=as.numeric(tres$estimate[1]),controls=as.numeric(tres$estimate[2]),t_pval=as.numeric(tres$p.value),mwu_pval=ures$p.value)
        },error=function(e) {print(paste(sampletype,sorting_strategy,cell_type,clustertype))})
      }
    }
    res <- as.data.frame(do.call(rbind,res.lst))
    res$t_fdr.adj <- p.adjust(res$t_pval,method='fdr')
    res$mwu_fdr.adj <- p.adjust(res$mwu_pval,method='fdr')
    f.out=paste0('~/Documents/Research/t21-proj/out/full/cellComp_res/cellComp.cluster.',sampletype,'.',sorting_strategy,'.','cluster','.txt')
    fwrite(res,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
  }
}

sampletype="Liver"
sorting_strategy="CD45+"
f.out=paste0("~/Documents/Research/t21-proj/out/full/cellComp_fig/boxplot.group.",sampletype,'.',gsub("[[:punct:]]","",sorting_strategy),'.pdf')
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

#####################

sampletype="Femur"
cluster_percentages_all[[sampletype]] <- list()
cluster_percentages[[sampletype]][["DownSyndrome"]][[sorting_strategy]]$disease_status <- "DownSyndrome"
cluster_percentages[[sampletype]][["Healthy"]][[sorting_strategy]]$disease_status <- "Healthy"
cluster_percentages_all[[sampletype]][[sorting_strategy]] <- list()
cell_type_groups <- unique(cluster_percentages[[sampletype]][["Healthy"]][[sorting_strategy]]$cell_type_groups)
mat1 <- subset(cluster_percentages[[sampletype]][["DownSyndrome"]][[sorting_strategy]],cell_type_groups==cell_type)
mat2 <- subset(cluster_percentages[[sampletype]][["Healthy"]][[sorting_strategy]],cell_type_groups==cell_type)
cluster_percentages_all[[sampletype]][[sorting_strategy]][[cell_type]] <- as.data.frame(rbind(mat1,mat2))
cell_type="Stroma"
sorting_strategy="CD235a-"
f.out=paste0("~/Documents/Research/t21-proj/out/full/cellComp_fig/stacked.cluster.",sampletype,'.',gsub("[[:punct:]]","",sorting_strategy),'.',gsub('\\/','_',cell_type),'.pdf')
pdf(file = f.out,width = 3.3*x,height=2.9*x)
g2<-ggplot(cluster_percentages_all[[sampletype]][[sorting_strategy]][[cell_type]], aes(fill=cluster, y=Freq*100, x=Patient_ID,col=disease_status)) + 
  geom_bar(position="stack", stat="identity") +
  scale_color_manual(values=c("black","orange"),labels=c("Down Syndrome","Healthy"),name="Disease Status") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(hjust=1,angle=60),
        plot.title = element_text(hjust=0.5)) +
  labs(x="Patient ID",y="Cell composition (%)",title=sorting_strategy,fill="Cell type")
print(g2)
dev.off()

res.lst <- list()
cluster.uniq=unique(cluster_percentages_all[[sampletype]][[sorting_strategy]][[cell_type]]$cluster)
iter=0; for (clustertype in cluster.uniq) {
  iter=iter+1
  tryCatch({
  x=subset(cluster_percentages_all[[sampletype]][[sorting_strategy]][[cell_type]],cluster==clustertype)
  tres <- t.test(x$Freq~x$disease_status)
  ures <- wilcox.test(x$Freq~x$disease_status)
  res.lst[[iter]] <- data.frame(sampletype,sorting_strategy,celltype,clustertype,cases=as.numeric(tres$estimate[1]),controls=as.numeric(tres$estimate[2]),t_pval=as.numeric(tres$p.value),mwu_pval=ures$p.value)
  },error=function(e) {print(paste(sampletype,sorting_strategy,celltype,clustertype))})
}
res <- as.data.frame(do.call(rbind,res.lst))
res$t_pval.adj <- p.adjust(res$t_pval,method='fdr')
res$mwu_pval.adj <- p.adjust(res$mwu_pval,method='fdr')
fwrite(res,'~/Documents/Research/t21-proj/out/full/cellComp_res/cellComp.cluster.Femur.CD235a.txt',quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)



    res.lst <- list(); iter=0; for (celltype in colnames(cluster_matrix[[sampletype]][[disease_status]][[sorting_strategy]])[-1]) { # remove Patient_ID
      # res.lst <- list(); iter=0; for (celltype in colnames(cluster_matrix[[sampletype]][[disease_status]][[sorting_strategy]])[-1]) { # remove Patient_ID
      iter=iter+1
      # if (celltype %in% colnames(cluster_matrix[[sampletype]][["DownSyndrome"]][[sorting_strategy]]) &
      #     celltype %in% colnames(cluster_matrix[[sampletype]][["Healthy"]][[sorting_strategy]])) {
      tryCatch({
        # tres <- t.test(cluster_matrix[[sampletype]][["DownSyndrome"]][[sorting_strategy]][,celltype],cluster_matrix[[sampletype]][["Healthy"]][[sorting_strategy]][,celltype])
        tres <- t.test(group_matrix[[sampletype]][["DownSyndrome"]][[sorting_strategy]][,celltype]~group_matrix[[sampletype]][["Healthy"]][[sorting_strategy]][,celltype])
        ures <- wilcox.test(group_matrix[[sampletype]][["DownSyndrome"]][[sorting_strategy]][,celltype],group_matrix[[sampletype]][["Healthy"]][[sorting_strategy]][,celltype])
        res.lst[[iter]] <- data.frame(sampletype,sorting_strategy,celltype,clustertype,cases=as.numeric(tres$estimate[1]),controls=as.numeric(tres$estimate[2]),t_pval=as.numeric(tres$p.value),mwu_pval=ures$p.value)
      },error=function(e) {print(paste(sampletype,sorting_strategy,celltype))})
    }
    res <- as.data.frame(do.call(rbind,res.lst))
    res$t_pval.adj <- p.adjust(res$t_pval,method='fdr')
    res$mwu_pval.adj <- p.adjust(res$mwu_pval,method='fdr')
    # f.out <- paste0("~/Documents/Research/t21-proj/out/full/cellComp_res/cellComp.group.",sampletype,'.',gsub("[[:punct:]]","",sorting_strategy),'.txt')
    # fwrite(res[[sampletype]][[sorting_strategy]],file = f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
  }
}


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
