library(data.table)
dir="/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/cpdb_Results"
# disease_status="Healthy"
disease_status="DownSyndrome"
envir="Liver"
dfmu = fread(paste0(dir,"/",disease_status,"_",envir,"/means.txt"),data.table = F,stringsAsFactors = F)
dfp = fread(paste0(dir,"/",disease_status,"_",envir,"/pvalues.txt"),data.table = F,stringsAsFactors = F)

# Get the number of cell types
# N_cellTypes = sqrt(ncol(dfmu)-11) # may be wrong

keep <- !(grepl("HLA",dfmu$interacting_pair) | grepl("complex",dfmu$interacting_pair) | dfmu$gene_a=="" | dfmu$gene_b=="" | apply(dfmu[,c("receptor_a","receptor_b")],1,sum) != 1)
# removing interacting_pairs like: 12oxoLeukotrieneB4_byPTGR1_LTB4R
# sometimes both genes not annotated as a receptor: SFRP5 or WNT6
dfmu <- dfmu[keep,]

# keep <- !(grepl("HLA",dfp$interacting_pair) | grepl("complex",dfp$interacting_pair)) # same as dfmu
dfp <- dfp[keep,]

rng <- 12:ncol(dfp)

dfmu.melt <- reshape2::melt(dfmu[,c(5,6,8,9,rng)],id.vars=c("gene_a","gene_b","receptor_a","receptor_b"))
x = strsplit(as.character(dfmu.melt$variable),"\\|")
dfmu.melt$cell1 = as.character(unlist(lapply(x,function(y){y[1]})))
dfmu.melt$cell2 = as.character(unlist(lapply(x,function(y){y[2]})))
dfmu.melt = subset(dfmu.melt,!(cell1==cell2))
dfmu.melt$Ligand = ifelse(!dfmu.melt$receptor_b,dfmu.melt$gene_b,dfmu.melt$gene_a)
dfmu.melt$Receptor = ifelse(dfmu.melt$receptor_b,dfmu.melt$gene_b,dfmu.melt$gene_a)
dfmu.melt$Sender = ifelse(!dfmu.melt$receptor_b,dfmu.melt$cell2,dfmu.melt$cell1)
dfmu.melt$Receiver = ifelse(dfmu.melt$receptor_b,dfmu.melt$cell2,dfmu.melt$cell1)
dfmu.melt = dfmu.melt[,c("Ligand","Receptor","Sender","Receiver","value")]
colnames(dfmu.melt)[5] <- "MeanExpr"

# # check to make sure it worked:
# subset(dfmu.melt,Ligand %in% c("CCL4","SLC7A1") & Receptor %in% c("CCL4","SLC7A1") &
#          Sender%in%c("Activated stellate cells","B cells") & Receiver%in%c("Activated stellate cells","B cells") )

dfp.melt <- reshape2::melt(dfp[,c(5,6,8,9,rng)],id.vars=c("gene_a","gene_b","receptor_a","receptor_b"))
x = strsplit(as.character(dfp.melt$variable),"\\|")
dfp.melt$cell1 = as.character(unlist(lapply(x,function(y){y[1]})))
dfp.melt$cell2 = as.character(unlist(lapply(x,function(y){y[2]})))
dfp.melt = subset(dfp.melt,!(cell1==cell2))
dfp.melt$Ligand = ifelse(!dfp.melt$receptor_b,dfp.melt$gene_b,dfp.melt$gene_a)
dfp.melt$Receptor = ifelse(dfp.melt$receptor_b,dfp.melt$gene_b,dfp.melt$gene_a)
dfp.melt$Sender = ifelse(!dfp.melt$receptor_b,dfp.melt$cell2,dfp.melt$cell1)
dfp.melt$Receiver = ifelse(dfp.melt$receptor_b,dfp.melt$cell2,dfp.melt$cell1)
dfp.melt = dfp.melt[,c("Ligand","Receptor","Sender","Receiver","value")]
colnames(dfp.melt)[5] <- "P"

df.mg <- merge(dfmu.melt,dfp.melt,by=c("Ligand","Receptor","Sender","Receiver"))
df.mg <- df.mg[order(df.mg$P,-1*df.mg$MeanExpr),]
df.mg <- subset(df.mg,P < 0.05)

f.out = paste0(dir,"/",disease_status,"_",envir,"/",disease_status,"_",envir,".processed_output.txt")
fwrite(df.mg,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)


table(df.mg$P==1)

dfmu.melt[100:105,]


x = strsplit(colnames(dfp)[rng],"\\|")
z <- lapply(x,function(y) {y[1]==y[2]})
remove <- which(unlist(z))
rng <- rng[-remove];



dfmu <- dfmu[,c(2,rng)]
dfp <- dfp[,c(2,rng)]


unique(unlist(strsplit(colnames(dfp)[12:ncol(dfp)],"\\|")))
columns = ["interacting_pair"] + mean.columns[-N_cellTypes**2:].tolist()

df_means = mean[columns]
df_pvals = pval[columns]

df_means.sort_values(by="interacting_pair", ascending=True, inplace=True)
df_means.reset_index(inplace=True, drop=True)

df_pvals.sort_values(by="interacting_pair", ascending=True, inplace=True)
df_pvals.reset_index(inplace=True, drop=True)

keep_idx = list()
for i in df_pvals["interacting_pair"].tolist():
  tmp = df_pvals[df_pvals["interacting_pair"]==i]
pvalues = tmp.T[tmp.index.tolist()[0]].tolist()[1:]
if all(pvalues) > 0.01:
  pass
else:
  keep_idx.append(tmp.index.tolist()[0])

df_means = df_means[df_means.index.isin(keep_idx)]
df_pvals = df_pvals[df_pvals.index.isin(keep_idx)]

dfs = list()
for col in df_means.columns[1:]:
  tmp = df_means[["interacting_pair"]]
tmp["interacting_group"] = col
tmp["color"] = df_means[col].tolist()
tmp["size"] = df_pvals[col].tolist()
dfs.append(tmp)

df = pd.concat(dfs, axis=0)
df = df.rename(columns={"interacting_group": "Cell group",
  "interacting_pair": "L-R pair",
  "size":"Corrected p-value",
  "color":"log(1+expression)"})
corrected = 0.01/len(df) # We will be conservative and apply the Bonferroni correction

df["log(1+expression)"] = round(np.log1p(df["log(1+expression)"]), 3)
df = df[df["Corrected p-value"]<=corrected]
df = df[df["log(1+expression)"]>=1.0]
df["tissue"] = g
df["Combination"] = df["L-R pair"]+": "+df["Cell group"]

dfs_final.append(df)
df_sorted = df.sort_values(by="log(1+expression)",ascending=False)
df_sorted.reset_index(drop=True, inplace=True)
del df_sorted["Combination"]
df_sorted.to_csv(out_dir+"Results/csvs/10pc/{}_filtered_pairs.csv".format(g))
'''
    df = pd.concat(dfs_final)
    df_ds = df[df["tissue"].str.contains("DownSyndrome")]
    df_h = df[df["tissue"].str.contains("Healthy")]
    combs_ds, combs_h = df_ds["Combination"].tolist(), df_h["Combination"].tolist()

    combs_diff = [x for x in combs_ds if x not in combs_h]
    df_ds = df_ds[df_ds["Combination"].isin(combs_diff)]
    del df_ds["Combination"]
    df_ds.sort_values(by="log(1+expression)",ascending=False, inplace=True)
    df_ds.reset_index(drop=True, inplace=True)
    df_ds.to_csv("Results/csvs/DS_only_Healthy_all_used_{}_filtered_pairs.csv".format("_".join(list(g)[0].split("_")[1:])))

