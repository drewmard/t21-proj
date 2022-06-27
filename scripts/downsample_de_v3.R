library(Seurat)
library(data.table)
library(Matrix.utils)
library( "DESeq2" )
library('variancePartition')
library('edgeR')
library('BiocParallel')

cell_type="HSCs/MPPs"
cell_type_filename = gsub("/","_",cell_type)
downsample="Down Syndrome"

disease_status="DownSyndrome"
f.out=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_liver_v_femur/",cell_type_filename,'.',disease_status,'.FE.txt')
full_diff_exp_results <- fread(f.out,data.table = F,stringsAsFactors = F)

geneAnnot=fread("/oak/stanford/groups/smontgom/amarder/t21-proj/out/data_small/gene_annot.txt",data.table = F,stringsAsFactors = F)

iter=0; df <- list()
for (sampletype in c("Liver","Femur")) {
  # disease_status="Healthy"
  for (disease_status in c("Healthy","DownSyndrome")) {
    print(iter)
    iter=iter+1
    f.out <- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/cell_subset/",disease_status,'_',sampletype,'_',cell_type_filename,".rds")
    df[[iter]] <- readRDS(file = f.out)
  }
}

dfcombined = merge(df[[1]],c(df[[2]],df[[3]],df[[4]]))

tmp <- subset(dfcombined@meta.data,leiden_names==cell_type & environment=="Healthy" & organ=="Liver"); 
summary1 <- aggregate(tmp$patient,by=list(patient=tmp$patient,sample=tmp$sample),length)
summary1 = subset(summary1,x >= 10)
tmp <- subset(dfcombined@meta.data,leiden_names==cell_type & environment=="Healthy" & organ=="Femur"); 
summary2 <- aggregate(tmp$patient,by=list(patient=tmp$patient,sample=tmp$sample),length)
summary2 = subset(summary2,x >= 10)

table(dfcombined@meta.data[,c("environment","organ")])

##############################

########################################################################################################

# in this step, we are going to downsample the liver data, for working with that dataset

permNum=1
for (permNum in 1:100) {
  print(permNum)
  set.seed(permNum)
  ind.lst <- c()
  
  tmp <- subset(dfcombined@meta.data,leiden_names==cell_type & environment=="Down Syndrome" & organ=="Liver"); 
  new_summary1 <- aggregate(tmp$patient,by=list(patient=tmp$patient,sample=tmp$sample),length)
  new_summary1.nsamp <- aggregate(tmp$sample,by=list(patient=tmp$patient),function(x) length(unique(x)))
  tmp <- subset(dfcombined@meta.data,leiden_names==cell_type & environment=="Down Syndrome" & organ=="Femur"); 
  new_summary2 <- aggregate(tmp$patient,by=list(patient=tmp$patient,sample=tmp$sample),length)
  new_summary2.nsamp <- aggregate(tmp$sample,by=list(patient=tmp$patient),function(x) length(unique(x)))
  

  # if patient in both femur and liver, and cell count liver > femur, then use X samples from liver
  patient_in_both <- unique(summary1$patient) # [summary1$patient %in% new_summary1$patient])
  iter=0; res <- list()
  for (patient_of_interest1 in patient_in_both) {
    iter=iter+1
    summary1.sub <- subset(summary1,patient==patient_of_interest1)
    potential_patients <- unique(subset(new_summary1.nsamp,x > nrow(summary1.sub))$patient)
    potential_patients <- sample(potential_patients,length(potential_patients),replace = F)
    # patient_of_interest2=potential_patients[1]
    j = 0; identified=FALSE; while (j <= length(potential_patients) & !identified) {
      j = j + 1
      patient_of_interest2 <- potential_patients[j]
      tmp <- subset(new_summary1,patient==patient_of_interest2)
      tmp <- tmp[order(tmp$x,decreasing = T),]
      if (mean(summary1.sub$x < tmp$x[1:nrow(summary1.sub)])==1) {
        res[[iter]] <- data.frame(patient=tmp$patient[1:nrow(summary1.sub)],
                                  sample=tmp$sample[1:nrow(summary1.sub)],
                                  x=summary1.sub$x
        )
        new_summary1 <- new_summary1[!(new_summary1$patient %in% res[[iter]]$patient),]
        new_summary1.nsamp <- new_summary1[!(new_summary1.nsamp$patient %in% res[[iter]]$patient),]
        identified=TRUE
      }
    }
  }
  res.all1 <- as.data.frame(do.call(rbind,res))
  
  patient_in_both <- unique(summary2$patient) # [summary1$patient %in% new_summary1$patient])
  iter=0; res <- list()
  for (patient_of_interest1 in patient_in_both) {
    iter=iter+1
    summary2.sub <- subset(summary2,patient==patient_of_interest1)
    potential_patients <- unique(subset(new_summary2.nsamp,x > nrow(summary2.sub))$patient)
    potential_patients <- sample(potential_patients,length(potential_patients),replace = F)
    # patient_of_interest2=potential_patients[1]
    j = 0; identified=FALSE; while (j <= length(potential_patients) & !identified) {
      j = j + 1
      patient_of_interest2 <- potential_patients[j]
      tmp <- subset(new_summary2,patient==patient_of_interest2)
      tmp <- tmp[order(tmp$x,decreasing = T),]
      unique
      if (mean(summary2.sub$x < tmp$x[1:nrow(summary2.sub)])==1) {
        res[[iter]] <- data.frame(patient=tmp$patient[1:nrow(summary2.sub)],
                                  sample=tmp$sample[1:nrow(summary2.sub)],
                                  x=summary2.sub$x
        )
        new_summary2 <- new_summary2[!(new_summary2$patient %in% res[[iter]]$patient),]
        new_summary2.nsamp <- new_summary2[!(new_summary2.nsamp$patient %in% res[[iter]]$patient),]
        identified=TRUE
      }
    }
  }
  res.all2 <- as.data.frame(do.call(rbind,res))
  
  for (i in 1:nrow(res.all2)) {
    ind <- unname(which(dfcombined$sample==res.all2$sample[i]))
    ind.lst <- c(ind.lst,sample(ind,res.all2$x[i],replace = FALSE))
  }
  ind2.lst <- c()
  for (i in 1:nrow(res.all1)) {
    ind <- unname(which(dfcombined$sample==res.all1$sample[i]))
    ind2.lst <- c(ind2.lst,sample(ind,res.all1$x[i],replace = FALSE))
  }
  
  column_to_use="leiden_names"
  dfsub <- dfcombined[,c(ind.lst,ind2.lst)]
  
  df.aggre <- aggregate.Matrix(
    t(
      GetAssayData(object = dfsub, slot = "counts", assay="RNA")
    ),
    groupings=dfsub$sample,fun="sum")
  
  df.aggre <- t(df.aggre)
  df.aggre <- as.data.frame(df.aggre)
  
  x <- unique(dfsub@meta.data[,c("patient","sample","sorting","environment","organ")])
  rownames(x) <- x[,"sample"]
  
  metadata_to_use <- x[rownames(x) %in% colnames(df.aggre),]
  df.aggre <- as.matrix(df.aggre[,match(rownames(metadata_to_use),colnames(df.aggre))])
  # samples_to_keep <- colnames(df.aggre)
  # df.aggre <- df.aggre[,samples_to_keep,drop=FALSE]
  # metadata_to_use <- metadata_to_use[samples_to_keep,]
  
  # Standard usage of limma/voom
  geneExpr = DGEList( df.aggre[full_diff_exp_results$names,] )
  geneExpr = calcNormFactors( geneExpr )
  
  # Specify parallel processing parameters
  # this is used implicitly by dream() to run in parallel
  param = SnowParam(8, "SOCK", progressbar=TRUE)
  
  # The variable to be tested must be a fixed effect
  form <- ~ organ + sorting #+ patient
  if (length(unique(metadata_to_use$sorting))==1 | length(unique(metadata_to_use$sorting)) > nrow(metadata_to_use)/2) {
    form <- ~ organ
  }
  
  # A positive FC is increased expression in the DS compared to healthy
  # metadata_to_use$environment <- factor(metadata_to_use$environment,levels=c("Healthy","DownSyndrome"))
  metadata_to_use$organ <- factor(metadata_to_use$organ,levels=c("Femur","Liver"))
  metadata_to_use$patient <- factor(metadata_to_use$patient)
  
  # estimate weights using linear mixed model of dream
  vobjDream = voomWithDreamWeights( geneExpr, form, metadata_to_use, BPPARAM=param )
  
  # Fit the dream model on each gene
  # By default, uses the Satterthwaite approximation for the hypothesis test
  fitmm = dream( vobjDream, form, metadata_to_use )
  fitmm = eBayes(fitmm)
  
  # reorganize:
  res.df1 <- as.data.frame(topTable( fitmm, number=Inf,coef = "organLiver" ))
  colnames(res.df1)[colnames(res.df1)=="ID"] <- "names"
  res.df1$disease_status <- disease_status

  # print(subset(res.df1,names%in%c("GATA1","APOC1","MYL4")))
  # mean(res.df1$adj.P.Val < 0.1)
  
  df1 <- merge(res.df1,geneAnnot[,c("hgnc_symbol","chromosome_name","start_position")],by.x="names",by.y="hgnc_symbol")
  df1$chr21 <- factor(ifelse(df1$chromosome_name=='21','Chr 21','Not Chr 21'),levels=c('Not Chr 21','Chr 21'))
  
  # save:
  subset_column='sample'
  system("mkdir -p /oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names_downsamp")
  f.out <- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names_downsamp/","DownSyndrome",".",cell_type_filename,".",subset_column,".perm",permNum,".txt")
  fwrite(df1,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

}




#