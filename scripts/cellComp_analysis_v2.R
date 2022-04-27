library(data.table)
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
        res.lst[[iter]] <- data.frame(sampletype,sorting_strategy,celltype,cases=as.numeric(tres$estimate[1]),controls=as.numeric(tres$estimate[2]),pval=as.numeric(tres$p.value))
      },error=function(e) {print(paste(sampletype,sorting_strategy,celltype))})
    }
    res[[sampletype]][[sorting_strategy]] <- as.data.frame(do.call(rbind,res.lst))
  }
}

