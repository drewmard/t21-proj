sampletype="Femur"

dbs <- c("GO_Molecular_Function_2018", "GO_Biological_Process_2018","ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")
outdir = "~/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/gsea"
for (db in dbs) {
  system(paste0("mkdir -p ",paste0(outdir,"/",db)))
}

for (sampletype in c("Femur","Liver")) {
  print(sampletype)
  res.df.all.p <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".p.rds"))
  res.df.all.lfc <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".lfc.rds"))
  cell_type_groups <- colnames(res.df.all.p[[1]])[(which(colnames(res.df.all.p[[1]])=="names")+1):(which(colnames(res.df.all.p[[1]])=="chromosome_name")-1)]
  
  for (iter in 2:2) { 
    print(iter)
    df.lfc=res.df.all.lfc[[iter]]
    df.p=res.df.all.p[[iter]]
    
    for (i in 1:length(cell_type_groups)) {
      
      for (subset_to_use in c("not_chr21_up","not_chr21_down","chr21_up")) {
        cell_type = cell_type_groups[i]
        print(paste0("enrichR: ",sampletype," - ",subset_to_use,' (',cell_type,', ',i,'/',length(cell_type_groups),')'))
        if (subset_to_use=="not_chr21_up") {
          x <- df.lfc$names[df.lfc[,cell_type] > 0 & df.p[,cell_type] < 0.1 & df.lfc[,"chr21"]=="Not Chr 21"]
          f.out <- paste0("gsea.not_chr21.up.iter",iter,".",cell_type,".",sampletype,".txt")
        } else if (subset_to_use=="not_chr21_down") {
          x <- df.lfc$names[df.lfc[,cell_type] < 0 & df.p[,cell_type] < 0.1 & df.lfc[,"chr21"]=="Not Chr 21"]
          f.out <- paste0("gsea.not_chr21.down.iter",iter,".",cell_type,".",sampletype,".txt")
        } else if (subset_to_use=="chr21_up") {
          x <- df.lfc$names[df.lfc[,cell_type] > 0 & df.p[,cell_type] < 0.1 & df.lfc[,"chr21"]=="Chr 21"]
          f.out <- paste0("gsea.chr21.up.iter",iter,".",cell_type,".",sampletype,".txt")
        } else if (subset_to_use=="chr21_down") {
          x <- df.lfc$names[df.lfc[,cell_type] < 0 & df.p[,cell_type] < 0.1 & df.lfc[,"chr21"]=="Chr 21"]
          f.out <- paste0("gsea.chr21.down.iter",iter,".",cell_type,".",sampletype,".txt")
        }
        
        x <- x[!is.na(x)]
        df.sub <- df.lfc[df.lfc$names %in% x,c(cell_type,"names")]
        gene_lst <- df.sub[,"names"]

        if (length(gene_lst) > 3) {
          enriched <- enrichr(gene_lst, dbs)
          for (k in 1:length(dbs)) {
            y <- enriched[[dbs[k]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); y <- subset(y,geneCt>3 & Adjusted.P.value<0.1)
            foutpath <- paste0(outdir,"/",dbs[[k]],"/",f.out)
            fwrite(y,foutpath,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
          }
        } else {
          for (k in 1:length(dbs)) {
            y <- data.frame("Term Overlap"=NA,Overlap=NA,P.value=NA,Adjusted.P.value=NA,Old.P.value=NA,Old.Adjusted.P.value=NA,Odds.Ratio=NA,Combined.Score=NA,Genes=NA,geneCt=NA,check.names = F)
            foutpath <- paste0(outdir,"/",dbs[[k]],"/",f.out)
            fwrite(y,foutpath,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
          }
        }
      }
    }
  }
}





###########3


library(data.table)

df.save <- list()
for (sampletype in c("Femur","Liver")) {
  
  for (iter in 2:2) { 
    
    for (db in dbs) {
      
      for (i in 1:length(cell_type_groups)) {
        
        for (subset_to_use in c("not_chr21_up","not_chr21_down","chr21_up","chr21_down")) {
          print(subset_to_use)
          cell_type = cell_type_groups[i]
          print(paste0("GSEA - gene set: ",gene_set,' (',cell_type,', ',i,'/',length(cell_type_groups),')'))
          if (subset_to_use=="not_chr21_up") {
            f.out <- paste0("gsea.not_chr21.up.iter",iter,".",cell_type,".",sampletype,".txt")
          } else if (subset_to_use=="not_chr21_down") {
            f.out <- paste0("gsea.not_chr21.down.iter",iter,".",cell_type,".",sampletype,".txt")
          } else if (subset_to_use=="chr21_up") {
            f.out <- paste0("gsea.chr21.up.iter",iter,".",cell_type,".",sampletype,".txt")
          } else if (subset_to_use=="chr21_down") {
            f.out <- paste0("gsea.chr21.down.iter",iter,".",cell_type,".",sampletype,".txt")
          }
          
          foutpath <- paste0(outdir,"/",db,"/",f.out)
          
          if (!file.exists(foutpath)) { next}
          df <- fread(foutpath,data.table = F,stringsAsFactors = F)
          df.sub <- subset(df,Adjusted.P.value < 0.1)
          
          if (nrow(df.sub) > 0) {
            df.sub$sampletype <- sampletype
            df.sub$iter <- iter
            df.sub$db <- db
            df.sub$cell_type <- cell_type
            df.sub$subset_to_use <- subset_to_use
          }
          
          if (is.null(df.save)) {
            df.save <- df.sub
          }
          df.save <- rbind(df.save,df.sub)
          
        }
      }
    }
  }
}

sort(table(df.save$Term),decreasing = T)[1:30]
table(subset(df.save,iter==2)[,c('subset_to_use','sampletype')])
table(subset(df.save,iter==2)[,c('subset_to_use','sampletype')])

df.save[order(df.save$Adjusted.P.value)[1:15],][,c("Term","Adjusted.P.value","geneCt","db","sampletype","cell_type","subset_to_use")]
df.sub <- subset(df.save,subset_to_use=="not_chr21_down" & db==dbs[1])
(sort(table(df.sub$Term),decreasing=T)[1:10])
df.sub <- subset(df.save,subset_to_use=="not_chr21_down" & db==dbs[2])
(sort(table(df.sub$Term),decreasing=T)[1:10])
df.sub <- subset(df.save,subset_to_use=="not_chr21_down" & db==dbs[3])
(sort(table(df.sub$Term),decreasing=T)[1:10])

df.save[grep("GATA1",df.save$Term),c("Term","Adjusted.P.value","geneCt","db","sampletype","cell_type","subset_to_use")]

df.sub <- subset(df.save,subset_to_use=="not_chr21_down" & db==dbs[3] & cell_type=="Late erythroid cells" & sampletype=="Femur")
dim(df.sub)
df.sub[,c("Term","Adjusted.P.value","geneCt","db","sampletype","cell_type","subset_to_use")]$Term
df.sub[,c("Term","Adjusted.P.value","geneCt","db","sampletype","cell_type","subset_to_use")][1:3,]
subset(df.sub,Term=="GATA1 CHEA")[,c("Term","Adjusted.P.value","geneCt","db","sampletype","cell_type","subset_to_use")]
subset(df.sub,Term=="MYC ENCODE")[,c("Term","Adjusted.P.value","geneCt","db","sampletype","cell_type","subset_to_use")]

df.sub <- subset(df.save,subset_to_use=="not_chr21_up" & db==dbs[3] & cell_type=="Late erythroid cells" & sampletype=="Femur")
dim(df.sub)
df.sub[,c("Term","Adjusted.P.value","geneCt","db","sampletype","cell_type","subset_to_use")]$Term
subset(df.sub,Term=="GATA1 CHEA")[,c("Term","Adjusted.P.value","geneCt","db","sampletype","cell_type","subset_to_use")]
subset(df.sub,Term=="MYC ENCODE")[,c("Term","Adjusted.P.value","geneCt","db","sampletype","cell_type","subset_to_use")]

df.sub1 <- subset(df.save,subset_to_use=="not_chr21_down" & db==dbs[3] & cell_type=="Late erythroid cells" & sampletype=="Femur")
df.sub2 <- subset(df.save,subset_to_use=="not_chr21_up" & db==dbs[3] & cell_type=="Late erythroid cells" & sampletype=="Femur")
df.sub.mg <- merge(df.sub1,df.sub2,by="Term",all=T)
cor.test(-log10(df.sub.mg$Adjusted.P.value.x),y=-log10(df.sub.mg$Adjusted.P.value.y))
df.sub.mg$Adjusted.P.value.x[is.na(df.sub.mg$P.value.x)] <- 1
df.sub.mg$Adjusted.P.value.y[is.na(df.sub.mg$P.value.y)] <- 1

ggplot(df.sub.mg,aes(x=-log10(Adjusted.P.value.x),y=-log10(Adjusted.P.value.y))) + 
  geom_point() + theme_bw() + theme(panel.grid = element_blank()) +
  labs(x='TF enrichment from downreg genes (-log10 FDR)',y='TF enrichment from upreg genes (-log10 FDR)') + 
  geom_smooth(method='loess',span=3,se=F,col='purple')

df.sub <- subset(df.save,subset_to_use=="not_chr21_down" & db==dbs[1] & cell_type=="Late erythroid cells" & sampletype=="Femur")
dim(df.sub)
df.sub <- subset(df.save,subset_to_use=="not_chr21_down" & db==dbs[2] & cell_type=="Late erythroid cells" & sampletype=="Femur")
dim(df.sub)
df.sub[,c("Term","Adjusted.P.value","geneCt","db","sampletype","cell_type","subset_to_use")]#$Term
df.sub[,c("Term","Adjusted.P.value","geneCt","db","sampletype","cell_type","subset_to_use")][1:100,]$Term


