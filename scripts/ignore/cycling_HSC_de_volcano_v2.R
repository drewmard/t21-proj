library(data.table)
library(ggplot2)
library(dplyr)

df = fread("~/Documents/Research/t21-proj/out/full/pseudobulks/sample/DE/DownSyndrome_cyc_vs_hsc.de.txt",data.table = F,stringsAsFactors = F)

g=ggplot(df,aes(x=logFC,y=-log10(P.Value),col=adj.P.Val < 0.1)) + geom_point(size=rel(0.5)) +
  ggpubr::theme_pubr() + labs(x="LFC",y=expression("-log"[10]~italic(P))) +
  scale_color_manual(values=c("#E4E4EB","#6AD7FC")) +
  guides(col=F)
pdf("~/Downloads/DownSyndrome_cyc_vs_hsc.de.pdf",width=7,height=3)
print(g)
dev.off()

# gene_lst = c("DNA2","LONP1","MGME1","RNF26","POLA1","MRE11","PIN1")#,"IFNAR1","MAVS","FUNDC1","NBR1","VDAG1","")
# s="MRE11;DHX9;PRKDC;CGAS;DDX41;GPATCH3;HERC5;DHX33;IRAK1;G3BP1;POLR2E;UFD1;POLR2F;POLR2H;POLR2K;XRCC6;XRCC5;RNF26;POLA1;ITCH;MAVS;NPLOC4;POLR3A;POLR3B;PIN1;POLR3G;POLR3K"
# gene_lst = as.character(unlist(strsplit(s,";")))
df$color = as.numeric(df$adj.P.Val < 0.1)
mito_lst = c("DNA2","LONP1","MGME1")
# s = "BECN1;VCP;TOMM40;VPS4A;PIK3R4;ATG101;FUNDC1;UBQLN1;EPG5;MFN2;ATG5;VTI1B;TOMM70;ATG3;CSNK2A1;VPS33A;VTA1;VPS37C;YOD1;SNF8;STAM;VPS37B;WIPI2;DYNLL1;ATG13;PRKAB1;DYNLL2;HGS;ATG16L1;PLAA;ATG4B;CHMP4A;VDAC1;PGAM5;VPS25;TOMM5;ARL8B"
# macroautophagy_lst = as.character(unlist(strsplit(s,";")))
# subset(df,names %in% macroautophagy_lst)[,]
# subset(df,names %in% macroautophagy_lst & names %in% c("VDAC1","ULK1","FUNDC1","NBR1"))[,]
macroautophagy_lst = c("ATG13","ATG3","YOD1","VDAC1")
ifn_lst = c("RNF26","POLA1","MRE11","PIN1")#,"IFNAR1","MAVS","FUNDC1","NBR1","VDAG1","")
# df$color[df$names %in% gene_lst] = 2
df$color[df$names %in% mito_lst] = 2
df$color[df$names %in% macroautophagy_lst] = 3
df$color[df$names %in% ifn_lst] = 4

s="SFPQ;UBQLN1;HTRA2;AKT1;P4HB;SOD2;SIRT1;PPIA"
oxidative_lst = as.character(unlist(strsplit(s,";")))
subset(df,names %in% oxidative_lst)[1:10,]
oxidative_lst = c("SFPQ","UBQLN1","SOD2")
s="BUB1B;CCP110;CDC20;PSMD8;PSMD9;CDC23;PSMD7;PSMD2;CDC27;PSMD3;CEP250;PSMD1;KNTC1;NEK2;FBXO5;PRKACA;TMEM14B;CDK5RAP2;CEP135;ANAPC7;CEP131;PSME3;ZW10;PSME4;ANAPC5;ANAPC1;ANAPC2;KLHL18;PSMD10;PSMD12;ANAPC15;PSMD11;PSMD14;PSMD13;CUL3;CUL1;HMMR;ANAPC10;ANAPC11;ZFP36L1;CCNB1;UBB;PRKAR2B;CEP70;UBC;CEP72;CEP192;CEP76;CEP78;PLK4;ODF2;UBE2C;PLK1;MAD2L1BP;CDC6;ANLN;TPX2;UBE2S;CDK2;CDK1;CEP57;YWHAE;CEP164;KIF14;ACTR1A;PCM1;CNTRL;HSP90AA1;TUBB;HAUS4;HAUS3;HAUS6;HAUS5;DYNLL1;TUBG1;CKAP5;HAUS2;HAUS1;DDB1;PSMA3;PSMA4;CEP63;MAPRE1;DCTN1;SSNA1;DCTN3;PSMA7;AURKA;PSMB6;PSMB7;FZR1;PSMB5;HAUS8;PSMB2;PSMB3;E2F1;BUB3;PCNT;DYNC1H1;CEP152;NDE1;TUBB4B;CENPE;NEDD1;CENPF;PSMC5;PSMC6;PSMC3;PSMC4;CENPJ;PSMC1;ALMS1;PSMC2;CEP41;MAD2L1"
cycling_lst = as.character(unlist(strsplit(s,";")))
s="RB1;BTG3;INTS13;CDKN1B;BUB1B;HNRNPU;SRA1;MKI67;CKS1B;CDC20;PSMD8;PSMD9;CDC23;MAEA;PSMD7;PSMD2;CDC27;PSMD3;PSMD1;RCC1;NEK2;FBXO5;CDK5RAP3;DLGAP5;PRMT5;BORA;ANAPC7;SCRIB;CDC25C;SIRT1;CDC25B;THAP1;DDB1;PSMA3;PSMA4;PSME3;PSME4;CKS2;ANAPC5;TTLL12;KIF20B;ANAPC1;ANAPC2;PSMD10;PSMD12;ANAPC15;PSMD11;PSMD14;PSMD13;CUL1;IQGAP1;PKMYT1;TTL;ANAPC10;PSMA7;ANAPC11;PSMB6;SDCBP;PSMB7;CCNB1;FZR1;PSMB5;PSMB2;PSMB3;BUB3;EIF4E;UBE2C;PLK1;GBF1;CDK11B;PSMC5;PSMC6;PSMC3;RPS6KB1;UBE2S;PSMC4;PSMC1;PSMC2;CDK1;MDM2;PIN1;MAD2L1"
cycling_lst2 = as.character(unlist(strsplit(s,";")))
# subset(df,names %in% cycling_lst & names %in% cycling_lst2)[1:10,]
cycling_lst3 = c("MKI67","CCNB1","PLK1","CDK1","NEK2")
df$color[df$names %in% cycling_lst3] = 5
df$color[df$names %in% oxidative_lst] = 6

df$color = as.factor(df$color)

df <- df %>% mutate(color_label = case_when(
  color == 0 ~ "Not sig",    # Change color == 1 to "Label1"
  color == 1 ~ "FDR < 0.1",    # Change color == 1 to "Label1"
  color == 2 ~ "Mitochondrial dysfunction",    # Change color == 2 to "Label2"
  color == 3 ~ "Macroautophagy",    
  color == 4 ~ "Type I interferon signalling pathway",
  color == 5 ~ "Cell cycling",    
  color == 6 ~ "Oxidative stress-induced\nintrinsic apoptotic signaling")   )

df$color_label <- factor(df$color_label, unique(df[order(df$color),c("color","color_label")])[,2])

  
# subset(df,color==2)
# fwrite(subset(df,color==2),"~/Downloads/regulation_of_type_I_interferon_production.csv",quote = F,na = "NA",sep = ",",row.names = F,col.names = T)
g=ggplot(df %>% arrange(color),aes(x=logFC,y=-log10(P.Value),col=color_label)) + geom_point() +
  ggpubr::theme_pubr() + labs(x="LFC",y=expression("-log"[10]~italic(P)),col="") +
  scale_color_manual(values=c("#E4E4EB","#6AD7FC","red","orange","purple","#EBEB78","#F7BAFF")) +
  theme(legend.text = element_text(size = 8));g #+
  # guides(col=F);g

# Subset the data where color > 1
subset_df <- df %>% filter(color %in% c(2,3,4,5,6))

# Position adjustment for the text
position_adj <- position_nudge(x = -0.4, y = 0.4)

# Add text to the plot
g_with_text <- g + geom_text(data = subset_df, aes(label = names), position = position_adj,size=2.5) + guides(col=guide_legend(override.aes = aes(label = "")))

# Display the plot
print(g_with_text)

pdf("~/Downloads/DownSyndrome_cyc_vs_hsc.de.pdf",width=7,height=3.5)
print(g)
dev.off()

pdf("~/Downloads/DownSyndrome_cyc_vs_hsc.de.labels.pdf",width=7,height=3.5)
print(g_with_text)
dev.off()


