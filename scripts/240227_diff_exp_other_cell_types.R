library(ggplot2)
library(reshape2)
cell_types = c("")

df = readRDS("~/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.Liver.lfc.rds")[[1]]
df[is.na(df)] = 0

cor.mat=cor(subset(df,chr21=="Not Chr 21")[,2:22],use='pairwise.complete.obs')
# Perform hierarchical clustering on the correlation matrix
d <- as.dist(1 - cor.mat)  # convert correlation matrix to distance matrix
hc <- hclust(d, method = "complete")  # perform hierarchical clustering

# Reorder the correlation matrix based on clustering
ordered_cor.mat <- cor.mat[hc$order, hc$order]

# Convert to long format
melted_cor <- melt(ordered_cor.mat)

# Plot heatmap with ggplot2
g=ggplot(melted_cor, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  coord_fixed() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Correlation Matrix between non-chr21 gene LFC",x="Cell type",y="Cell type"); g
f.out = "~/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/corr.Z4a.pdf"
pdf(f.out,width = 8,height=8)
print(g)
dev.off()


# Plot the reordered correlation matrix
library(corrplot)
corrplot(ordered_cor.mat, method = "color", type = "upper", tl.col = "black", tl.srt = 45)

distances = dist(apply(df[,2:22],2,scale))
dfall = readRDS("~/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.Liver.p.rds")[[1]]
df=dfall
thres=0.05
df[,2:22][dfall[,2:22]>=thres] = 0
df[,2:22][dfall[,2:22]<thres] = 1
df[,2:22][is.na(dfall[,2:22])] = 0
distances = dist(df)
# Perform hierarchical clustering
hc <- hclust(distances)
# Clusters
k <- 10
clusters <- cutree(hc, k)

df$cluster=clusters
table(df$cluster)
aggregate((df[,2:22]),by=list(df$cluster),mean,na.rm=T)

aggregate(abs(df[,2:22]),by=list(df$cluster),median)

subset(df,cluster==4)[1:5,]



oxidative_stress_genes <- c(
  "SOD1", "SOD2", "CAT", "GPX1", "GPX4", 
  "GSTP1", "NFE2L2", "HMOX1", "NQO1", "GSR",
  "TXN", "TXNRD1", "TXNRD2", "PRDX1", "PRDX2",
  "PRDX3", "PRDX4", "PRDX5", "PRDX6", "GCLC",
  "GCLM", "GLRX", "GLRX2", "GSS", "MGST1",
  "MGST3", "MGST2", "ALDH2", "ALDH3A1", "ALDH9A1",
  "ALOX5", "ALOX15", "ALOX12", "ALOXE3", "NRF1",
  "NRF2", "NRF3", "KEAP1", "ARE", "HO-1",
  "G6PD", "PRMT1", "NOS1", "NOS2", "NOS3",
  "PRDX1", "PRDX2", "PRDX3", "PRDX4", "PRDX5"
)

subset(df,names %in% oxidative_stress_genes)
aggregate(subset(df,names %in% oxidative_stress_genes),by=list(subset(df,names %in% oxidative_stress_genes)$cluster),mean,na.rm=T)

subset(df,names %in% oxidative_stress_genes)[,"HSCs_MPPs"]

# Plot dendrogram
plot(hc, main = "Dendrogram of Hierarchical Clustering")


library(data.table)
library(enrichR)
cell_type = "HSCs_MPPs"
res.lst = list()
for (cell_type in c("HSCs_MPPs","Pro B cells","Pre pro B cells","pDCs","NK progenitors","NK cells","Monocyte progenitors","MEMPs","Megakaryocytes","Mast cells","Late erythroid cells","Kupffer cells","Inflammatory macrophages","Granulocyte progenitors","Early erythroid cells","cDC2")) {
# for (cell_type in c("Monocyte progenitors","MEMPs","Megakaryocytes","Mast cells","Late erythroid cells","Kupffer cells","Inflammatory macrophages","Granulocyte progenitors","Early erythroid cells","cDC2")) {
  de_res = fread(paste0("~/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/Liver.",cell_type,".sample.txt"),data.table = F,stringsAsFactors = F)
  gene_lst <- subset(de_res,adj.P.Val < 0.05 & logFC > 0)$names
  # gene_lst <- subset(MEMP,adj.P.Val < 0.05)$names
  
  dbs <- c("GO_Molecular_Function_2018", "GO_Biological_Process_2018","ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")
  enriched <- enrichr(gene_lst, dbs)
  # y <- enriched[[dbs[1]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); y <- subset(y,geneCt>3 & Adjusted.P.value<0.1);y
  y <- enriched[[dbs[2]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); y <- subset(y,geneCt>3 & Adjusted.P.value<0.1);y
  if (nrow(y) > 0) {
    tmp = data.frame(cell_type,term=y$Term)
  } else {
    tmp = data.frame(cell_type=character(),term=character())
  }
  res.lst[[cell_type]] = tmp
}

res.df = do.call(rbind,res.lst)

# Assuming your data frame is named 'res.df'
library(tidyr)
library(dplyr)

# Create a binary indicator column (set to 1) for presence of each term-cell_type combination
presence_df <- res.df %>%
  mutate(presence = 1) %>%
  spread(cell_type, presence, fill = 0)
  # select(-term)  # Remove the original 'term' column
rownames(presence_df) = presence_df$term
presence_df$term <- NULL  # Remove the 'term' column after setting row names

# Compute Jaccard distance between rows (terms) of the binary matrix
jaccard_dist <- dist(as.matrix(presence_df), method = "binary")

# Perform hierarchical clustering based on Jaccard distance
hc <- hclust(jaccard_dist, method = "complete")

# Reorder the correlation matrix based on clustering
presence_df.ordered <- presence_df[hc$order,]

# Compute Jaccard distance between rows (terms) of the binary matrix
jaccard_dist <- dist(t(as.matrix(presence_df)), method = "binary")

# Perform hierarchical clustering based on Jaccard distance
hc <- hclust(jaccard_dist, method = "complete")

# Reorder the correlation matrix based on clustering
presence_df.ordered <- presence_df.ordered[,hc$order]

# Convert to long format 
presence_df.ordered$term = rownames(presence_df)
melted_cor <- melt(presence_df.ordered,by='term')

# Plot heatmap with ggplot2
g=ggplot(melted_cor, aes(term, variable, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  # theme_bw() +
  theme_minimal() +
  # theme(panel.grid = element_blank()) +
  # coord_fixed() +
  # theme(axis.text.x = element_blank()) +
  coord_flip() +
  theme(axis.text.y = element_blank()) +
  labs(title = "",x="GO Term",y="Cell type") +
  guides(fill="none"); g
f.out = "~/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/go.Z4b.pdf"
pdf(f.out,width = 15.5,height=5)
print(g)
dev.off()

presence_df$total = apply(presence_df,1,sum)
presence_df[order(presence_df$total,decreasing = T),][1:50,]
subset(presence_df,HSCs_MPPs==1 & total<=3)
presence_df["pattern recognition receptor signaling pathway (GO:0002221)",]
data.frame(colnames(presence_df),apply(presence_df,2,sum))

mean(presence_df$total[colnames(presence_df)[1]],na.rm = T)

mean(presence_df$total[which(presence_df[,colnames(presence_df)[1]]==1)])
mean(presence_df$total[which(presence_df[,colnames(presence_df)[2]]==1)])
median(presence_df$total[which(presence_df[,colnames(presence_df)[3]]==1)])

presence_df[which(presence_df[,colnames(presence_df)[3]]==1),]
presence_df["regulation of type I interferon-mediated signaling pathway (GO:0060338)",]



library(data.table)
library(enrichR)
cell_type = "HSCs_MPPs"
res.lst = list()
for (cell_type in c("HSCs_MPPs","Pro B cells","Pre pro B cells","pDCs","NK progenitors","NK cells","Monocyte progenitors","MEMPs","Megakaryocytes","Mast cells","Late erythroid cells","Kupffer cells","Inflammatory macrophages","Granulocyte progenitors","Early erythroid cells","cDC2")) {
  # for (cell_type in c("Monocyte progenitors","MEMPs","Megakaryocytes","Mast cells","Late erythroid cells","Kupffer cells","Inflammatory macrophages","Granulocyte progenitors","Early erythroid cells","cDC2")) {
  de_res = fread(paste0("~/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/Liver.",cell_type,".sample.txt"),data.table = F,stringsAsFactors = F)
  gene_lst <- subset(de_res,adj.P.Val < 0.05 & logFC > 0)$names
  # gene_lst <- subset(MEMP,adj.P.Val < 0.05)$names
  
  dbs <- c("GO_Molecular_Function_2018", "GO_Biological_Process_2018","ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")
  enriched <- enrichr(gene_lst, dbs)
  # y <- enriched[[dbs[1]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); y <- subset(y,geneCt>3 & Adjusted.P.value<0.1);y
  y <- enriched[[dbs[2]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); y <- subset(y,geneCt>3 & Adjusted.P.value<0.1);y
  if (nrow(y) > 0) {
    tmp = data.frame(cell_type,term=y$Term,P.value=y$P.value,Adjusted.P.value=y$Adjusted.P.value,Genes=y$Genes)
  } else {
    tmp = data.frame(cell_type=character(),term=character(),P.value=character(),Adjusted.P.value=character(),Genes=character())
  }
  res.lst[[cell_type]] = tmp
}
res.df = as.data.frame(do.call(rbind,res.lst))

f.out = paste0("~/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/Liver.sample.enrichr_GO_BP.txt")
fwrite(res.df,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

res.lst = list()
for (cell_type in c("HSCs_MPPs","Pro B cells","Pre pro B cells","pDCs","NK progenitors","NK cells","Monocyte progenitors","MEMPs","Megakaryocytes","Mast cells","Late erythroid cells","Kupffer cells","Inflammatory macrophages","Granulocyte progenitors","Early erythroid cells","cDC2")) {
  # for (cell_type in c("Monocyte progenitors","MEMPs","Megakaryocytes","Mast cells","Late erythroid cells","Kupffer cells","Inflammatory macrophages","Granulocyte progenitors","Early erythroid cells","cDC2")) {
  de_res = fread(paste0("~/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/Liver.",cell_type,".sample.txt"),data.table = F,stringsAsFactors = F)
  gene_lst <- subset(de_res,adj.P.Val < 0.05 & logFC < 0)$names
  # gene_lst <- subset(MEMP,adj.P.Val < 0.05)$names
  
  dbs <- c("GO_Molecular_Function_2018", "GO_Biological_Process_2018","ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")
  enriched <- enrichr(gene_lst, dbs)
  # y <- enriched[[dbs[1]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); y <- subset(y,geneCt>3 & Adjusted.P.value<0.1);y
  y <- enriched[[dbs[2]]]; y$geneCt <- as.numeric(gsub("\\/[0-9]*","",y[,'Overlap'],fixed=F)); y <- subset(y,geneCt>3 & Adjusted.P.value<0.1);y
  if (nrow(y) > 0) {
    tmp = data.frame(cell_type,term=y$Term,P.value=y$P.value,Adjusted.P.value=y$Adjusted.P.value,Genes=y$Genes)
  } else {
    tmp = data.frame(cell_type=character(),term=character(),P.value=character(),Adjusted.P.value=character(),Genes=character())
  }
  res.lst[[cell_type]] = tmp
}
res.df = as.data.frame(do.call(rbind,res.lst))

f.out = paste0("~/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/Liver.sample.enrichr_GO_BP.down.txt")
fwrite(res.df,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)



library(data.table)
HSC = fread("~/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/Liver.HSCs_MPPs.sample.txt",data.table = F,stringsAsFactors = F)
MEMP = fread("~/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/Liver.MEMPs.sample.txt",data.table = F,stringsAsFactors = F)
EarlyEry = fread("~/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/Liver.Early erythroid cells.sample.txt",data.table = F,stringsAsFactors = F)
LateEry = fread("~/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/Liver.Late erythroid cells.sample.txt",data.table = F,stringsAsFactors = F)
Megakary = fread("~/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/Liver.Megakaryocytes.sample.txt",data.table = F,stringsAsFactors = F)
Mast = fread("~/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/Liver.Mast cells.sample.txt",data.table = F,stringsAsFactors = F)
# subset(Mast,names=="GATA1")
# HSC = fread("~/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/Femur.HSCs_MPPs.sample.txt",data.table = F,stringsAsFactors = F)
# MEMP = fread("~/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/Femur.MEMPs.sample.txt",data.table = F,stringsAsFactors = F)
# EarlyEry = fread("~/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/Femur.Early erythroid cells.sample.txt",data.table = F,stringsAsFactors = F)
# LateEry = fread("~/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/Femur.Late erythroid cells.sample.txt",data.table = F,stringsAsFactors = F)
# Megakary = fread("~/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/Femur.Megakaryocytes.sample.txt",data.table = F,stringsAsFactors = F)

# col_to_use="logFC"
col_to_use="P.Value"
df.mg = merge(HSC[,c("names",col_to_use)],
              MEMP[,c("names",col_to_use)],
              by="names",all = TRUE)
colnames(df.mg)[2:3] = c("HSC","MEMP")
df.mg = merge(df.mg,EarlyEry[,c("names",col_to_use)],by="names",all = TRUE); colnames(df.mg)[ncol(df.mg)] = "EarlyEry"
df.mg = merge(df.mg,LateEry[,c("names",col_to_use)],by="names",all = TRUE); colnames(df.mg)[ncol(df.mg)] = "LateEry"
df.mg = merge(df.mg,Megakary[,c("names",col_to_use)],by="names",all = TRUE); colnames(df.mg)[ncol(df.mg)] = "Megakary"
subset(df.mg,names=="GATA1")
subset(df.mg,names=="GATA2")
subset(df.mg,names=="RUNX1")
subset(df.mg,names=="SOD1")

######################

library(data.table)
df1=fread("~/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/Liver.sample.enrichr_GO_BP.down.txt",data.table = F,stringsAsFactors = F)
df2=fread("~/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/Liver.sample.enrichr_GO_BP.txt",data.table = F,stringsAsFactors = F)


df1$set = "down"
df2$set = "up"

fwrite(rbind(df1,df2),"~/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/Liver.sample.enrichr_GO_BP.all.txt",quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

