library(tidyverse)
library(motifbreakR)
library(MotifDb)
library(BiocParallel)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(data.table)

genome <- BSgenome.Hsapiens.UCSC.hg38
bpparam <- MulticoreParam(workers = multicoreWorkers())
print("Loaded genome annotations")

# directories
variants_path = "~/Downloads/rbc_peaks_overlap.bed"
out_name = "~/Downloads/motifbreaker.txt"
temp_name = "~/Downloads/tmp.bed"

# load in variants that show significant allelic imbalance in DHS set and format them 
# create variant ID needed for motifbreakR and drop unused columns 
# for motifbreakR 
variants = fread(variants_path,data.table = F,stringsAsFactors = F)
variants = subset(variants,V5 > 0.2)

print("Variants read")

library(SNPlocs.Hsapiens.dbSNP150.GRCh38)
variants.rs = variants[grep("rs",variants[,4]),]
variants.not_rs = variants[!grep("rs",variants[,4]),]

snps.mb <- snps.from.rsid(rsid = variants.rs[,4],
                          dbSNP = SNPlocs.Hsapiens.dbSNP150.GRCh38 ,
                          search.genome = BSgenome.Hsapiens.UCSC.hg38)
results <- motifbreakR(snpList = snps.mb, filterp = TRUE,
                       pwmList = hocomoco,
                       threshold = 1e-4,
                       method = "ic",
                       BPPARAM = bpparam)
# results2 = calculatePvalue(results)

strong_effects = data.frame(results[results$effect == 'strong'])

head(sort(table(strong_effects$geneSymbol),decreasing = T),20)
table(strong_effects$geneSymbol)["GATA1"]

subset(strong_effects,geneSymbol=="BACH1")
subset(strong_effects,SNP_id=="rs1887428")
subset(variants,variants[,4]=="rs1887428")
subset(pk,peakid=="chr9-4984170-4986204")
subset(res.df.mg,peak=="chr9-4984170-4986204")
subset(de_genes,names=="JAK2")
subset(cyc_HSC_pb,names=="JAK2")

subset(strong_effects,geneSymbol=="GATA2")
subset(variants,variants[,4]=="rs746051334")
subset(pk,peakid=="chr6-37003824-37006559")
subset(res.df.mg,peak=="chr6-37003824-37006559")
subset(de_genes,names=="FGD2")
subset(cyc_HSC_pb,names=="FGD2")

subset(strong_effects,geneSymbol=="GATA1")
subset(strong_effects,SNP_id=="rs746051334")

subset(strong_effects,geneSymbol=="GATA1")
subset(strong_effects,SNP_id=="rs11866877")
subset(variants,variants[,4]=="rs11866877")
subset(pk,peakid=="chr16-119367-120477")
subset(res.df.mg,peak=="chr16-119367-120477")
subset(de_genes,names=="NPRL3")
subset(cyc_HSC_pb,names=="NPRL3")

subset(snps.mb,SNP_id %in% "rs7775698")
subset(results,SNP_id %in% "rs7775698")

subset(strong_effects,SNP_id %in% "rs2075672")
plotMB(results = results, rsid = "rs2075672", effect = "strong")

subset(de_genes,names=="HMGN1")
subset(cyc_HSC_pb,names=="HMGN1")

snps.mb.tmp <- snps.from.rsid(rsid = "rs72928038",
                              dbSNP = SNPlocs.Hsapiens.dbSNP150.GRCh38 ,
                              search.genome = BSgenome.Hsapiens.UCSC.hg38)
results.tmp <- motifbreakR(snpList = snps.mb.tmp, filterp = TRUE,
                           pwmList = hocomoco,
                           threshold = 1e-4,
                           method = "ic",
                           BPPARAM = bpparam)
results.tmp
###########



pk = fread("~/Downloads/rbc.fm_peak.pip_0.2.bed",data.table = F,stringsAsFactors = F)
#

celltype_to_use="HSCs_H"
output_file = paste0("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/scent/out1/",celltype_to_use,".txt")
res.df.h = fread(output_file,data.table = F,stringsAsFactors = F)

celltype_to_use="HSCs_T21"
output_file = paste0("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/scent/out1/",celltype_to_use,".txt")
res.df.t21 = fread(output_file,data.table = F,stringsAsFactors = F)

res.df.mg = merge(res.df.h,res.df.t21,all=TRUE,by=c("gene","peak"))


de_genes = fread("~/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/Liver.HSCs_MPPs.sample.txt",data.table = F,stringsAsFactors = F)
cyc_HSC_pb = fread("~/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/HSC_v_CyclingHSC.txt",data.table = F,stringsAsFactors = F)
subset(de_genes,names=="JAK2")
subset(cyc_HSC_pb,names=="JAK2")

variants.sub = subset(variants.rs,V9 %in%
                        pk[(pk$acc2 - pk$acc3) > 0.05,"peakid"])
results2 = subset(results,SNP_id %in% variants.sub[,4])

results2 = results[variants.sub$V4,]
strong_effects2 = data.frame(results2[results2$effect == 'strong'])

head(sort(table(strong_effects2$geneSymbol),decreasing = T),20)
table(strong_effects2$geneSymbol)["GATA1"]

                      
                      