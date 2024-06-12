library(data.table)
library(dplyr)
df1 = fread("~/Downloads/DownSyndrome_HSC_PeakGeneSets/liver_v_femur.txt",data.table = F,stringsAsFactors = F)
df2 = fread("~/Downloads/DownSyndrome_HSC_PeakGeneSets/Liver.HSCs_MPPs.sample.txt",data.table = F,stringsAsFactors = F)
df3 = fread("~/Downloads/DownSyndrome_HSC_PeakGeneSets/peak_gene.txt",data.table = F,stringsAsFactors = F)
df4 = fread("~/Downloads/DownSyndrome_HSC_PeakGeneSets/DownSyndrome_cyc_vs_hsc.de.txt",data.table = F,stringsAsFactors = F)
df5 = fread("~/Downloads/ts21_disomy_motif_enrichment_results-2.arm.tsv",data.table = F,stringsAsFactors = F)

genes = fread('~/Documents/Research/data/grch38/genes.extended.txt',data.table = F,stringsAsFactors = F)
genes_chr21 = subset(genes,`Chromosome/scaffold name`==21)$`Gene name`

df1.1 = (subset(df1,names %in% genes_chr21 & adj.P.Val.t21 < 0.05)) %>%
  mutate(Analysis = "Liver vs Femur HSCs (Ts21 DE)") %>%
  select(Analysis,
         ID=names,
         effect_size=logFC.t21,
         pval=P.Value.t21,
         fdr=adj.P.Val.t21,
         info=class
  )
df1.2 = (subset(df1,names %in% genes_chr21 & adj.P.Val.h < 0.05)) %>%
  mutate(Analysis = "Liver vs Femur HSCs (Disomy DE)") %>%
  select(Analysis,
         ID=names,
         effect_size=logFC.h,
         pval=P.Value.h,
         fdr=adj.P.Val.h,
         info=class
  )
df2.1 = (subset(df2,names %in% genes_chr21 & adj.P.Val < 0.05)) %>%
  mutate(Analysis = "Ts21 vs Disomy liver HSCs (DE)",
         class=NA) %>%
  select(Analysis,
         ID=names,
         effect_size=logFC,
         pval=P.Value,
         fdr=adj.P.Val,
         info=class
  )
df3.1 = (subset(df3,gene %in% genes_chr21 & fdr_int < 0.2)) %>%
  mutate(Analysis = "SCENT peak-gene (All HSCs - interaction term)",
         class=NA) %>%
  select(Analysis,
         ID=gene_peak,
         effect_size=beta_int,
         pval=boot_basic_p_int,
         fdr=fdr_int,
         info=class
  )
df3.2 = (subset(df3,gene %in% genes_chr21 & fdr_H < 0.2)) %>%
  mutate(Analysis = "SCENT peak-gene (Disomy HSCs)",
         class=NA) %>%
  select(Analysis,
         ID=gene_peak,
         effect_size=beta_H,
         pval=boot_basic_p_H,
         fdr=fdr_H,
         info=class
  )
df3.3 = (subset(df3,gene %in% genes_chr21 & fdr_t21 < 0.2)) %>%
  mutate(Analysis = "SCENT peak-gene (Ts21 HSCs)",
         class=NA) %>%
  select(Analysis,
         ID=gene_peak,
         effect_size=beta_t21,
         pval=boot_basic_p_t21,
         fdr=fdr_t21,
         info=class
  )


df4.1 = (subset(df4,names %in% genes_chr21 & adj.P.Val < 0.05)) %>%
  mutate(Analysis = "Ts21 liver cycling vs less-cycling HSCs (DE)",
         class=NA) %>%
  select(Analysis,
         ID=names,
         effect_size=logFC,
         pval=P.Value,
         fdr=adj.P.Val,
         info=class
  )
df5.1=subset(df5,motif_name %in% c("BACH1","ETS2","Erg","GABPA",genes_chr21) & disomy_negLog10Padj > (-log10(0.05)) & peak_type=="all") %>%
  mutate(Analysis = "Motif enrichment: Ts21 vs Disomy liver HSCs (All Peaks - Disomy enriched)",
         class=NA,
         pval=NA) %>%
  select(Analysis,
         ID=motif_name,
         effect_size=disomy_log2enr,
         pval=pval,
         fdr=disomy_negLog10Padj,
         info=class
  )
df5.2=subset(df5,motif_name %in% c("BACH1","ETS2","Erg","GABPA",genes_chr21) & disomy_negLog10Padj > (-log10(0.05)) & peak_type=="Enhancer") %>%
  mutate(Analysis = "Motif enrichment: Ts21 vs Disomy liver HSCs (Enhancer Peaks - Disomy enriched)",
         class=NA,
         pval=NA) %>%
  select(Analysis,
         ID=motif_name,
         effect_size=disomy_log2enr,
         pval=pval,
         fdr=disomy_negLog10Padj,
         info=class
  )
df5.3=subset(df5,motif_name %in% c("BACH1","ETS2","Erg","GABPA",genes_chr21) & disomy_negLog10Padj > (-log10(0.05)) & peak_type=="Promoter") %>%
  mutate(Analysis = "Motif enrichment: Ts21 vs Disomy liver HSCs (Promoter Peaks - Disomy enriched)",
         class=NA,
         pval=NA) %>%
  select(Analysis,
         ID=motif_name,
         effect_size=disomy_log2enr,
         pval=pval,
         fdr=disomy_negLog10Padj,
         info=class
  )

df5.4=subset(df5,motif_name %in% c("BACH1","ETS2","Erg","GABPA",genes_chr21) & ts21_negLog10Padj > (-log10(0.05)) & peak_type=="all") %>%
  mutate(Analysis = "Motif enrichment: Ts21 vs Disomy liver HSCs (All Peaks - Ts21 enriched)",
         class=NA,
         pval=NA) %>%
  select(Analysis,
         ID=motif_name,
         effect_size=ts21_log2enr,
         pval=pval,
         fdr=ts21_negLog10Padj,
         info=class
  )
df5.5=subset(df5,motif_name %in% c("BACH1","ETS2","Erg","GABPA",genes_chr21) & ts21_negLog10Padj > (-log10(0.05)) & peak_type=="Enhancer") %>%
  mutate(Analysis = "Motif enrichment: Ts21 vs Disomy liver HSCs (Enhancer Peaks - Ts21 enriched)",
         class=NA,
         pval=NA) %>%
  select(Analysis,
         ID=motif_name,
         effect_size=ts21_log2enr,
         pval=pval,
         fdr=ts21_negLog10Padj,
         info=class
  )
df5.6=subset(df5,motif_name %in% c("BACH1","ETS2","Erg","GABPA",genes_chr21) & ts21_negLog10Padj > (-log10(0.05)) & peak_type=="Promoter") %>%
  mutate(Analysis = "Motif enrichment: Ts21 vs Disomy liver HSCs (Promoter Peaks - Ts21 enriched)",
         class=NA,
         pval=NA) %>%
  select(Analysis,
         ID=motif_name,
         effect_size=ts21_log2enr,
         pval=pval,
         fdr=ts21_negLog10Padj,
         info=class
  )

chr21_result_summary = rbind(df1.1,
      df1.2,
      df2.1,
      df3.1,
      df3.2,
      df3.3,
      df4.1,
      df5.1,
      df5.2,
      df5.3,
      df5.4,
      df5.5,
      df5.6) %>% as.data.frame()

f.out = "~/Downloads/chr21_result_summary.txt"
fwrite(chr21_result_summary,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)




