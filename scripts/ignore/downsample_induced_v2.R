# load
library(data.table)
library('edgeR')
library(data.table)
library(Seurat)
library(Matrix.utils) # for pseudobulking
library(BiocParallel) # for parallel
library(variancePartition)
# input:
fmeta1 = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/downsample_induced/input/meta1.txt")
fmeta2 = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/downsample_induced/input/meta2.txt")
meta1=fread(fmeta1,data.table = F,stringsAsFactors = F)
meta2=fread(fmeta2,data.table = F,stringsAsFactors = F)
meta1.ordered=meta1[order(meta1$num_cells,decreasing = TRUE),]
meta2.ordered=meta2[order(meta2$num_cells,decreasing = TRUE),]
match_between_conditions = NULL # "sorting" # here, conditions is trisomy vs disomy
match_between_differential_factor = "patient" # and diff factor is cyc vs less cycling

disease_status="DownSyndrome"
sampletype="Liver"
subset_column="sample"

for (k in 1:2) {
  print(paste0("Running for: ",k))
  if (k==1) { 
    pseudobulk_covariate = "sample"
    cluster_label="combi_annot"
    cell_type="HSCs"
    cell_type_filename = gsub("/","_",cell_type)
    f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/downsample_induced/input/",disease_status,".",sampletype,".",cell_type_filename,".",cluster_label,".rds")
    dfseurat1 = readRDS(f)
    cell_type="Cycling HSC"
    cell_type_filename = gsub("/","_",cell_type)
    f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/downsample_induced/input/",disease_status,".",sampletype,".",cell_type_filename,".",cluster_label,".rds")
    dfseurat2 = readRDS(f)
  } else if (k==2) {
    pseudobulk_covariate = "sample"
    cluster_label="leiden_names"
    cell_type="HSCs/MPPs"
    cell_type_filename = gsub("/","_",cell_type)
    f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/downsample_induced/input/",disease_status,".",sampletype,".",cell_type_filename,".",cluster_label,".rds")
    dfseurat1 = readRDS(f)
    cell_type="Cycling HSCs/MPPs"
    cell_type_filename = gsub("/","_",cell_type)
    f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/downsample_induced/input/",disease_status,".",sampletype,".",cell_type_filename,".",cluster_label,".rds")
    dfseurat2 = readRDS(f)
  }
  
  meta1.input = dfseurat1@meta.data 
  meta2.input = dfseurat2@meta.data
  
  get_sampler_df = function(meta1.input,meta1.ordered) { 
    # meta1.input = dfseurat1@meta.data # the new data to be downsampled
    # meta1 is the old data table that specifies number of cells per pseudobulk
    
    # list of potential samples we can use
    tab = aggregate(meta1.input[,1],by=meta1.input[,c(pseudobulk_covariate,match_between_differential_factor,match_between_conditions)],length)
    colnames(tab)[ncol(tab)] = "num_cells"
    tab = subset(tab,num_cells >= 10)
    tab = tab[order(tab$num_cells,decreasing = T),]
    
    keepers = list()
    for (i in 1:nrow(meta1.ordered)) {
      keepers[[i]] = tab[tab$num_cells >= meta1.ordered$num_cells[i],pseudobulk_covariate]
      if (length(keepers[[i]]) == 0) {
        keepers[[i]]
      } 
    }
    # iteration 1
    sample.vec = rep(NA,nrow(meta1.ordered))
    for (i in 1:nrow(meta1.ordered)) {
      if (length(keepers[[i]]) == 0) {
        # sample one of top 10% samples that still remain if no samples are larger than the small dataset sample
        tab.sub = tab[!(tab[,pseudobulk_covariate] %in% sample.vec),]
        sample_to_use = sample(tab.sub[1:ceiling(nrow(tab.sub)*0.1),pseudobulk_covariate],1)
      } else {
        # potential list of samples
        keepers_potential = keepers[[i]][!(keepers[[i]] %in% sample.vec)]
        if (length(keepers_potential)==0) {
          # sample one of top 10% samples that still remain if no samples that remain haven't been used
          tab.sub = tab[!(tab[,pseudobulk_covariate] %in% sample.vec),]
          sample_to_use = sample(tab.sub[1:ceiling(nrow(tab.sub)*0.1),pseudobulk_covariate],1)
        } else {
          # sample one of top 10% samples that still remain if no samples that remain haven't been used
          sample_to_use = sample(keepers_potential,1)
          while (sample_to_use %in% sample.vec) { # this should never happen, I subset down to the keepers that are not already in sample.vec
            print("Blah!")
            sample_to_use = sample(keepers[[i]],1)    # this would be faster if I subset down to the keepers that are already in sample.vec
          }
        }
      }
      sample.vec[i] = as.character(sample_to_use)
    }
    sampler_df = data.frame(sample=as.character(sample.vec),num_cells=meta1.ordered$num_cells)
    return(sampler_df)
  }
  j=1
  
  for (j in 1:100) {
    
    print(paste0("Running: ",j,"/100..."))
    
    set.seed(j)
    meta1.input = dfseurat1@meta.data 
    meta2.input = dfseurat2@meta.data
    sampler_df1 = get_sampler_df(meta1.input,meta1.ordered)
    sampler_df2 = get_sampler_df(meta2.input,meta2.ordered)
    
    get_idx = function(dfseurat1meta.data,sampler_df1) { 
      ind.lst = c()
      for (i in 1:nrow(sampler_df1)) {
        # first, find the indices of all cells from the sample
        ind = unname(which(dfseurat1meta.data[,pseudobulk_covariate]==sampler_df1[i,pseudobulk_covariate]))
        ind.lst <- c(ind.lst,sample(ind,min(sampler_df1[i,"num_cells"],length(ind)),replace = FALSE))
      }
      return(ind.lst)
    }
    ind1 = get_idx(dfseurat1meta.data = dfseurat1@meta.data,sampler_df1)
    ind2 = get_idx(dfseurat1meta.data = dfseurat2@meta.data,sampler_df2)
    df.sub1=dfseurat1[,ind1]
    df.sub2=dfseurat2[,ind2]
    
    get_pseudobulk_and_metadata = function(df.sub1) { 
      seurat_meta = df.sub1@meta.data
      # create pseudobulks:
      df.aggre <- aggregate.Matrix(
        t(
          GetAssayData(object = df.sub1, slot = "counts", assay="RNA")
        ),
        groupings=seurat_meta[,pseudobulk_covariate],fun="sum")
      df.aggre <- t(df.aggre)
      df.aggre <- as.data.frame(df.aggre)
      
      # create metadata:
      col_of_interest = c("patient","sample","sorting","environment","organ","age",cluster_label)
      x <- unique(seurat_meta[,col_of_interest])
      rownames(x) <- as.character(x[,subset_column])
      
      # align metadata and pseudobulks:
      metadata_to_use <- x[rownames(x) %in% colnames(df.aggre),]
      df.aggre <- as.data.frame.matrix(df.aggre[,match(rownames(metadata_to_use),colnames(df.aggre))])
      
      return(list(df.aggre,metadata_to_use))
    }
    
    res = get_pseudobulk_and_metadata(df.sub1); df.aggre1 = res[[1]]; meta1 = res[[2]]
    res = get_pseudobulk_and_metadata(df.sub2); df.aggre2 = res[[1]]; meta2 = res[[2]]
    
    # to account for same patients in femur and liver
    rownames(meta1) = paste0(rownames(meta1),".cyc")
    rownames(meta2) = paste0(rownames(meta2),".hsc")
    colnames(df.aggre1) = paste0(colnames(df.aggre1),".cyc")
    colnames(df.aggre2) = paste0(colnames(df.aggre2),".hsc")
    
    df.aggre = merge(df.aggre1,df.aggre2,by=0)
    rownames(df.aggre) = df.aggre[,1]
    df.aggre = df.aggre[,-1]
    meta1$cond = "cond1"
    meta2$cond = "cond2"
    meta = rbind(meta1,meta2)
    
    # align metadata and pseudobulks:
    metadata_to_use <- meta[rownames(meta) %in% colnames(df.aggre),]
    # metadata_to_use$age = as.numeric(substring(metadata_to_use$age,1,2))
    df.aggre <- as.data.frame.matrix(df.aggre[,match(rownames(metadata_to_use),colnames(df.aggre))])
    
    ########################
    # DE analysis!
    
    # Specify parallel processing parameters
    # this is used implicitly by dream() to run in parallel
    param = SnowParam(8, "SOCK", progressbar=TRUE)
    
    # The variable to be tested must be a fixed effect
    # for (use_pcw in c(TRUE,FALSE)) {
    use_pcw = FALSE
    
    form <- ~ cond + sorting
    
    # A positive FC is increased expression in the DS compared to healthy
    metadata_to_use$cond <- factor(metadata_to_use$cond,levels=c("cond1","cond2"))
    metadata_to_use$sorting = as.character(metadata_to_use$sorting)
    if ("sorting" %in% colnames(metadata_to_use)) {
      metadata_to_use$sorting[!(metadata_to_use$sorting %in% c("CD235a-","CD45+"))] <- "Other"
    }
    metadata_to_use$sorting = as.factor(metadata_to_use$sorting)
    
    # preprocessing:
    geneExpr = DGEList( df.aggre )
    # keep <- filterByExpr(geneExpr, group=metadata_to_use$cond)
    # geneExpr <- geneExpr[keep,]
    geneExpr = calcNormFactors( geneExpr )
    
    # estimate weights using linear mixed model of dream
    vobjDream = voomWithDreamWeights( geneExpr, form, metadata_to_use, BPPARAM=param )
    
    # Fit the model on each gene
    fitmm = dream( vobjDream, form, metadata_to_use )
    fitmm = eBayes(fitmm)
    
    # reorganize:
    res.df1 <- as.data.frame(topTable( fitmm, number=Inf,coef = "condcond2" ))
    # res.df1 <- as.data.frame(topTable( fitmm, number=Inf,coef = "age" ))
    res.df1$names <- rownames(res.df1)
    
    # dir.create(paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/downsample_induced/output"),showWarnings = FALSE)
    f.out = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/downsample_induced/output/",disease_status,".",sampletype,".cyc_vs_hsc.",cluster_label,".perm",j,".txt")
    fwrite(res.df1,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
  }
}


for (k in 1:2) {
  print(paste0("Running for: ",k))
  if (k==1) { 
    cluster_label="combi_annot"
  } else if (k==2) {
    cluster_label="leiden_names"
  }
  
  for (j in 1:100) {
    print(j)
    disease_status="DownSyndrome"
    f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/downsample_induced/output/",disease_status,".",sampletype,".cyc_vs_hsc.",cluster_label,".perm",j,".txt")
    df <- fread(f,data.table=F,stringsAsFactors = F)
    df1.full <- df[!duplicated(df$names),]
    df1.full.lfc.tmp <- df1.full[,c("names","logFC")]
    df1.full.p.tmp <- df1.full[,c("names","P.Value")]
    df1.full.fdr.tmp <- df1.full[,c("names","adj.P.Val")]
    colnames(df1.full.lfc.tmp)[2] <- paste0(colnames(df1.full.lfc.tmp)[2],".",j)
    colnames(df1.full.p.tmp)[2] <- paste0(colnames(df1.full.p.tmp)[2],".",j)
    colnames(df1.full.fdr.tmp)[2] <- paste0(colnames(df1.full.fdr.tmp)[2],".",j)
    if (j==1) {
      df1.full.lfc = df1.full.lfc.tmp
      df1.full.p = df1.full.p.tmp
      df1.full.fdr = df1.full.fdr.tmp
    } else {
      df1.full.lfc = merge(df1.full.lfc,df1.full.lfc.tmp,by='names',all=TRUE)
      df1.full.p = merge(df1.full.p,df1.full.p.tmp,by='names',all=TRUE)
      df1.full.fdr = merge(df1.full.fdr,df1.full.fdr.tmp,by='names',all=TRUE)
    }
  }
  
  
  df1.full.p[is.na(df1.full.p)] = 1
  downsamp_df = data.frame(names=df1.full.p$names,medP=apply(df1.full.p[,!(colnames(df1.full.p) %in% "names")],1,median))
  
  if (cluster_label=="combi_annot") {
    suffix = ".integ1"
  } else {
    suffix = ""
  }
  disease_status="DownSyndrome"
  f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/pseudobulks/",pseudobulk_covariate,"/DE/",disease_status,'_',"cyc_vs_hsc",".de",suffix,".txt")
  df1=fread(f,data.table = F,stringsAsFactors = F)[,1:7]
  colnames(df1)[1:6] = paste0(colnames(df1)[1:6],".t21")
  disease_status="Healthy"
  suffix = ".integ1"
  f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/pseudobulks/",pseudobulk_covariate,"/DE/",disease_status,'_',"cyc_vs_hsc",".de",suffix,".txt")
  df2=fread(f,data.table = F,stringsAsFactors = F)[,1:7]
  colnames(df2)[1:6] = paste0(colnames(df2)[1:6],".h")
  
  df.mg = merge(merge(df1,df2,by="names",all=TRUE),downsamp_df,by="names",all.x=TRUE)
  df.mg$class = "none"
  df.mg$class[df.mg$adj.P.Val.t21 < 0.05 & df.mg$adj.P.Val.h < 0.05] = "cell-type driven"
  df.mg$class[df.mg$adj.P.Val.t21 < 0.05 & (df.mg$adj.P.Val.h > 0.05 | is.na(df.mg$adj.P.Val.h))] = "unknown"
  df.mg$class[df.mg$adj.P.Val.h < 0.05 & (df.mg$adj.P.Val.t21 > 0.05 | is.na(df.mg$adj.P.Val.t21))] = "healthy-specific"
  df.mg$class[df.mg$adj.P.Val.t21 < 0.05 & (df.mg$adj.P.Val.h > 0.05 | is.na(df.mg$adj.P.Val.h)) & df.mg$medP < 0.05] = "Ts21-specific"
  
  f.out = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/downsample_induced/output/",sampletype,".cyc_vs_hsc.",cluster_label,".MedP.txt")
  fwrite(df.mg,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
}


for (k in 1:2) {
  print(paste0("Running for: ",k))
  if (k==1) { 
    cluster_label="combi_annot"
    prefix="integrated_t21"
  } else if (k==2) {
    cluster_label="leiden_names"
    prefix="original_labels"
  }
  
  library("enrichR")
  dbs <- c("GO_Molecular_Function_2021", "GO_Biological_Process_2021","ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")

  f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/downsample_induced/output/",sampletype,".cyc_vs_hsc.",cluster_label,".MedP.txt")
  df.mg = fread(f,data.table = F,stringsAsFactors = F)
  
  gene_lst=subset(df.mg,class=="Ts21-specific" & logFC.t21 > 0)$names
  enriched <- enrichr(gene_lst, dbs)
  y <- enriched[[dbs[1]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); y = (subset(y,geneCt>=3 & Adjusted.P.value<0.1)); y$Set <- dbs[1]; 
  f.out = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/downsample_induced/output/gsea/",prefix,".upreg.",dbs[1],".txt")
  fwrite(y,f.out,quote = F,na = NA,sep = "\t",row.names = F,col.names = T)
  y <- enriched[[dbs[2]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); y = (subset(y,geneCt>=3 & Adjusted.P.value<0.1)); y$Set <- dbs[2]; 
  f.out = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/downsample_induced/output/gsea/",prefix,".upreg.",dbs[2],".txt")
  fwrite(y,f.out,quote = F,na = NA,sep = "\t",row.names = F,col.names = T)
  y <- enriched[[dbs[3]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); y = (subset(y,geneCt>=3 & Adjusted.P.value<0.1)); y$Set <- dbs[3]; 
  f.out = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/downsample_induced/output/gsea/",prefix,".upreg.",dbs[3],".txt")
  fwrite(y,f.out,quote = F,na = NA,sep = "\t",row.names = F,col.names = T)
  
  gene_lst=subset(df.mg,class=="Ts21-specific" & logFC.t21 < 0)$names
  enriched <- enrichr(gene_lst, dbs)
  y <- enriched[[dbs[1]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); y = (subset(y,geneCt>=3 & Adjusted.P.value<0.1)); y$Set <- dbs[1]; 
  f.out = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/downsample_induced/output/gsea/",prefix,".downreg.",dbs[1],".txt")
  fwrite(y,f.out,quote = F,na = NA,sep = "\t",row.names = F,col.names = T)
  y <- enriched[[dbs[2]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); y = (subset(y,geneCt>=3 & Adjusted.P.value<0.1)); y$Set <- dbs[2]; 
  f.out = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/downsample_induced/output/gsea/",prefix,".downreg.",dbs[2],".txt")
  fwrite(y,f.out,quote = F,na = NA,sep = "\t",row.names = F,col.names = T)
  y <- enriched[[dbs[3]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); y = (subset(y,geneCt>=3 & Adjusted.P.value<0.1)); y$Set <- dbs[3]; 
  f.out = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/downsample_induced/output/gsea/",prefix,".downreg.",dbs[3],".txt")
  fwrite(y,f.out,quote = F,na = NA,sep = "\t",row.names = F,col.names = T)
  
  # this fails for when there is no sig GO terms
}



