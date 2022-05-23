library(data.table)

df.save <- list()
for (sampletype in c("Femur","Liver")) {

  for (iter in 1:3) { 

    for (gene_set in gene_set_list) {

      for (i in 1:length(cell_type_groups)) {
        
        for (subset_to_use in c("not_chr21_up","not_chr21_down","chr21_up","chr21_down")) {
          print(subset_to_use)
          cell_type = cell_type_groups[i]
          print(paste0("GSEA - gene set: ",gene_set,' (',cell_type,', ',i,'/',length(cell_type_groups),')'))
          if (subset_to_use=="not_chr21_up") {
            f.out <- paste0(outdir,"/",gene_set,".v7.5/gsea.not_chr21.up.iter",iter,".",cell_type,".",sampletype,".txt")
          } else if (subset_to_use=="not_chr21_down") {
            f.out <- paste0(outdir,"/",gene_set,".v7.5/gsea.not_chr21.down.iter",iter,".",cell_type,".",sampletype,".txt")
          } else if (subset_to_use=="chr21_up") {
            f.out <- paste0(outdir,"/",gene_set,".v7.5/gsea.chr21.up.iter",iter,".",cell_type,".",sampletype,".txt")
          } else if (subset_to_use=="chr21_down") {
            f.out <- paste0(outdir,"/",gene_set,".v7.5/gsea.chr21.down.iter",iter,".",cell_type,".",sampletype,".txt")
          }
          
          if (!file.exists(f.out)) { next}
          df <- fread(f.out,data.table = F,stringsAsFactors = F)
          df.sub <- subset(df,padj < 0.1)
          
          if (nrow(df.sub) > 0) {
            df.sub$sampletype <- sampletype
            df.sub$iter <- iter
            df.sub$gene_set <- gene_set
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

x <- paste(df.save$sampletype,
           df.save$iter ,
           df.save$gene_set ,
           df.save$cell_type ,
           df.save$subset_to_use)
x <- paste(df.save$sampletype,
           df.save$gene_set ,
           df.save$cell_type ,
           df.save$subset_to_use)
y <- aggregate(x,by=list(df.save$iter,x),length)
z <- reshape2::dcast(Group.2~Group.1,data=y)
z[is.na(z)] <- 0
head(sort(table(x),decreasing = T))

table(df.save$iter)
table(df.save$subset_to_use)
df.save[df.save$subset_to_use=="chr21_down",]
sort(table(df.save[df.save$subset_to_use=="chr21_up" & df.save$iter==2 & df.save$gene_set=="c2.all",'pathway']),decreasing = T)[1:5]
sort(table(df.save[df.save$subset_to_use=="not_chr21_up" & df.save$iter==2 & df.save$gene_set=="c2.all",'pathway']),decreasing = T)[1:10]
sort(table(df.save[df.save$subset_to_use=="not_chr21_down" & df.save$iter==2 & df.save$gene_set=="c2.all",'pathway']),decreasing = T)[1:10]

df.sub <- df.save[order(df.save$padj),c(1,3,7,9:13)]
subset(df.sub,gene_set!="c8.all")[1:10,]

tab <- sort(table(df.save$pathway[df.save$gene_set=="h.all"]),decreasing = T); tab[1:5]
tab <- sort(table(df.save$pathway[df.save$gene_set=="c8.all"]),decreasing = T); tab[1:5]
tab <- sort(table(df.save$pathway[df.save$gene_set=="c5.go"]),decreasing = T); tab[1:5]
tab <- sort(table(df.save$pathway[df.save$gene_set=="c2.all"]),decreasing = T); tab[1:5]

####

df.sub <- df.save[df.save$cell_type=="Early erythroid cells" & df.save$gene_set=="c2.all",c(1,3,7,9:13)]
df.sub[order(df.sub$padj),][1:5,]
###


iter=1
subset(res.df.all.lfc[[iter]],names=="SOD1")
#SOD1 downregulated in late erythroid cells in femur in analysis 1... but upregulated in analysis 2 and 3

iter=2
subset(res.df.all.lfc[[iter]],names=="NDUFV3")


df.save[df.save$pathway=="DIAZ_CHRONIC_MYELOGENOUS_LEUKEMIA_UP",c(1,3,7,9:13)]
df.save[df.save$pathway=="JAATINEN_HEMATOPOIETIC_STEM_CELL_DN",c(1,3,7,9:13)]
df.save[df.save$pathway=="PILON_KLF1_TARGETS_UP",c(1,3,7,9:13)]
df.save[df.save$pathway=="CAIRO_LIVER_DEVELOPMENT_DN",c(1,3,7,9:13)]
df.save[df.save$pathway=="NABA_SECRETED_FACTORS",] # S100A9 shows up 
df.save[df.save$pathway=="HOUNKPE_HOUSEKEEPING_GENES",c(1,3,7,9:13)]
df.save[df.save$pathway=="FISCHER_DREAM_TARGETS",c(1,3,7,9:13)]


