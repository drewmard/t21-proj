library(data.table)
library(Seurat)
library(Matrix.utils) # for pseudobulking
library(BiocParallel) # for parallel
library(variancePartition) # for limma voom

sampletype="Femur"
disease_status="Healthy"
for (sampletype in c("Liver","Femur")) {
  for (disease_status in c("Healthy","DownSyndrome")) {
    cell_type = "HSCs/MPPs"
    cell_type_filename = gsub("/","_",cell_type)
    
    print(paste("Reading:",disease_status,sampletype,cell_type))
    f<- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/cell_subset/",disease_status,'_',sampletype,'_',cell_type_filename,".rds")
    dfseurat=readRDS(f)
    
    # Subset matrix down to cell type of interest
    dfsub = subset(dfseurat,leiden_names==cell_type)
    
    for (subset_column in c("patient","sample")) { 
      if (subset_column=="patient") {
        col_of_interest = c("patient","environment","organ","age")
      } else if (subset_column=="sample") {
        col_of_interest = c("patient","sample","sorting","environment","organ","age")
      }
      
      # create pseudobulks:
      df.aggre <- aggregate.Matrix(
        t(
          GetAssayData(object = dfsub, slot = "counts", assay="RNA")
        ),
        groupings=dfsub@meta.data[,subset_column],fun="sum")
      df.aggre <- t(df.aggre)
      df.aggre <- as.data.frame(df.aggre)
      
      # create metadata:
      x <- unique(dfsub@meta.data[,col_of_interest])
      rownames(x) <- as.character(x[,subset_column])
      
      # align metadata and pseudobulks:
      metadata_to_use <- x[rownames(x) %in% colnames(df.aggre),]
      df.aggre <- as.data.frame.matrix(df.aggre[,match(rownames(metadata_to_use),colnames(df.aggre))])
      
      # save:
      f.out = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/pseudobulks/",subset_column,"/pb/",disease_status,'_',sampletype,'_',cell_type_filename,".pb.txt")
      fwrite(df.aggre,f.out,na = "NA",sep = "\t",quote = F,row.names = T,col.names = T)
      
      f.out = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/pseudobulks/",subset_column,"/meta/",disease_status,'_',sampletype,'_',cell_type_filename,".meta.txt")
      fwrite(metadata_to_use,f.out,na = "NA",sep = "\t",quote = F,row.names = T,col.names = T)
    }
  }
}


sampletype="Femur"
disease_status="Healthy"
subset_column="sample"
cell_type = "HSCs/MPPs"
cell_type_filename = gsub("/","_",cell_type)
for (disease_status in c("Healthy","DownSyndrome")) {
  for (subset_column in c("patient","sample")) { 
    
    # disease_status="Healthy"
    print(paste0("Reading data: ",disease_status," ",subset_column,"..."))
    
    sampletype="Liver"
    f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/pseudobulks/",subset_column,"/pb/",disease_status,'_',sampletype,'_',cell_type_filename,".pb.txt")
    df.aggre1<-fread(f,data.table = F,stringsAsFactors = F,header = T)
    f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/pseudobulks/",subset_column,"/meta/",disease_status,'_',sampletype,'_',cell_type_filename,".meta.txt")
    meta1<-fread(f,data.table = F,stringsAsFactors = F,header = T)
    rownames(df.aggre1) = df.aggre1[,1]
    df.aggre1 = df.aggre1[,-1]
    rownames(meta1) = meta1[,1]
    meta1 = meta1[,-1]
    
    sampletype="Femur"
    f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/pseudobulks/",subset_column,"/pb/",disease_status,'_',sampletype,'_',cell_type_filename,".pb.txt")
    df.aggre2<-fread(f,data.table = F,stringsAsFactors = F,header = T)
    f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/pseudobulks/",subset_column,"/meta/",disease_status,'_',sampletype,'_',cell_type_filename,".meta.txt")
    meta2<-fread(f,data.table = F,stringsAsFactors = F,header = T)
    rownames(df.aggre2) = df.aggre2[,1]
    df.aggre2 = df.aggre2[,-1]
    rownames(meta2) = meta2[,1]
    meta2 = meta2[,-1]
    
    # to account for same patients in femur and liver
    rownames(meta1) = paste0(rownames(meta1),".liver")
    rownames(meta2) = paste0(rownames(meta2),".femur")
    colnames(df.aggre1) = paste0(colnames(df.aggre1),".liver")
    colnames(df.aggre2) = paste0(colnames(df.aggre2),".femur")
    
    df.aggre = merge(df.aggre1,df.aggre2,by=0)
    rownames(df.aggre) = df.aggre[,1]
    df.aggre = df.aggre[,-1]
    meta = rbind(meta1,meta2)
    
    # align metadata and pseudobulks:
    metadata_to_use <- meta[rownames(meta) %in% colnames(df.aggre),]
    metadata_to_use$age = as.numeric(substring(metadata_to_use$age,1,2))
    df.aggre <- as.data.frame.matrix(df.aggre[,match(rownames(metadata_to_use),colnames(df.aggre))])
    
    # Specify parallel processing parameters
    # this is used implicitly by dream() to run in parallel
    param = SnowParam(8, "SOCK", progressbar=TRUE)
    
    # The variable to be tested must be a fixed effect
    for (use_pcw in c(TRUE,FALSE)) {
      if (!use_pcw & subset_column=="sample") {
        form <- ~ organ + sorting
      } else if (use_pcw & subset_column=="sample") {
        form <- ~ organ + sorting + age
      } else if (use_pcw & subset_column=="patient") {
        form <- ~ organ  + age
      } else if (!use_pcw & subset_column=="patient") {
        form <- ~ organ
      }
      
      # A positive FC is increased expression in the DS compared to healthy
      metadata_to_use$organ <- factor(metadata_to_use$organ,levels=c("Femur","Liver"))
      if ("sorting" %in% colnames(metadata_to_use)) {
        metadata_to_use$sorting[!(metadata_to_use$sorting %in% c("CD235a-","CD45+"))] <- "Other"
      }
      
      # estimate weights using linear mixed model of dream
      vobjDream = voomWithDreamWeights( df.aggre, form, metadata_to_use, BPPARAM=param )
      
      # Fit the model on each gene
      fitmm = dream( vobjDream, form, metadata_to_use )
      fitmm = eBayes(fitmm)
      
      # reorganize:
      res.df1 <- as.data.frame(topTable( fitmm, number=Inf,coef = "organLiver" ))
      # res.df1 <- as.data.frame(topTable( fitmm, number=Inf,coef = "age" ))
      res.df1$names <- rownames(res.df1)
      res.df1$disease_status <- disease_status
      
      age_suffix = ifelse(use_pcw,".age","")
      dir.create(paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/pseudobulks/",subset_column,"/DE"),showWarnings = FALSE)
      f.out = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/pseudobulks/",subset_column,"/DE/",disease_status,'_',cell_type_filename,".de",age_suffix,".txt")
      fwrite(res.df1,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
    }
  }
}

library(data.table)
disease_status = "Healthy"
for (disease_status in c("Healthy","DownSyndrome")) {
  i=0; for (subset_column in c("patient","sample")) { 
    for (use_pcw in c(TRUE,FALSE)) {
      
      i = i +1
      
      age_suffix = ifelse(use_pcw,".age","")
      f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/pseudobulks/",subset_column,"/DE/",disease_status,'_',cell_type_filename,".de",age_suffix,".txt")
      res = fread(f,data.table = F,stringsAsFactors = F)
      res = res[,c("names","logFC")]
      colnames(res)[2] = paste0(colnames(res)[2],".",subset_column,".",use_pcw)
      if (i==1) {
        res.mg = res
      } else {
        res.mg = merge(res.mg,res,by="names")
      }
    }
  }
  f.out=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/pseudobulks/DE.",disease_status,".txt")
  fwrite(res.mg,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
}


res.save=list(); i =0; for (sampletype in c("DownSyndrome","Healthy")) { 
  i = i+1
  res.mg = fread(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/pseudobulks/DE.",sampletype,".txt"),data.table = F,stringsAsFactors = F)
  # res.mg = fread("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/pseudobulks/DE.DownSyndrome.txt",data.table = F,stringsAsFactors = F)
  # ggplot(res.mg,aes(x=logFC.sample.FALSE,y=logFC.patient.FALSE)) + geom_point() + geom_abline(slope=1,intercept = 0,lty='dashed',col='red') + ggpubr::theme_pubr() + labs(x="LFC (sample-level pseudobulk)",y="LFC (fetus-level pseudobulk)")
  # ggplot(res.mg,aes(x=logFC.sample.TRUE,y=logFC.patient.TRUE)) + geom_point() + geom_abline(slope=1,intercept = 0,lty='dashed',col='red')
  
  colnames(res.mg)[2:5] = c("Fetus (+ age)","Fetus","Sample (+ age)","Sample")
  cor.mat = cor(res.mg[,-1],use='na.or.complete')
  cor.melt = melt(cor.mat)
  
  g <- ggplot(cor.melt,aes(x=Var1,y=Var2,fill=value)) +
    geom_tile(col='black') +
    theme_bw()  +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 60, hjust=1),
          panel.border = element_blank(),
          legend.position = 'none',
          plot.title = element_text(hjust=0.5)
    ) +
    scale_fill_gradient2(low = "blue", mid='white',high = "red",name='Power') +
    geom_text(aes(label=round(value,3))) + 
    labs(x='Differential expression analysis',y='Differential expression analysis',title=sampletype); g
  pdf(paste0("~/Documents/Research/t21-proj/out/full/pseudobulks/lfc_correlation.",sampletype,".pdf"),width=6,height=6)
  print(g)
  dev.off()
  
  colnames(res.mg)[-1] = paste0(colnames(res.mg)[-1]," (",sampletype,")")
  res.save[[i]] = res.mg
}

res.save = merge(res.save[[1]],res.save[[2]],by='names')
cor.mat=cor(res.save[,-1],use='na.or.complete')
cor.melt = melt(cor.mat)
g <- ggplot(cor.melt,aes(x=Var1,y=Var2,fill=scale(value))) +
  geom_tile(col='black') +
  theme_bw()  +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 60, hjust=1),
        panel.border = element_blank(),
        legend.position = 'none',
        plot.title = element_text(hjust=0.5)
  ) +
  scale_fill_gradient2(low = "orange", mid='red',high = "red",name='Power') +
  geom_text(aes(label=round(value,3))) + 
  labs(x='Differential expression analysis',y='Differential expression analysis',"Liver vs Femur"); g
g <- g + scale_x_discrete(limits = rev(levels(cor.melt$Var1)))  # Flip the x-axis and place tick labels at the bottom
pdf(paste0("~/Documents/Research/t21-proj/out/full/pseudobulks/lfc_correlation.","All",".pdf"),width=6,height=6)
print(g)
dev.off()


library(gplots)  # Load the gplots package

heatmap(cor(res.mg[,-1],use='na.or.complete'))

cor(res.mg[,-1],use='na.or.complete')
