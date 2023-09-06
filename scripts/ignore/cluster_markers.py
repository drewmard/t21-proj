import scanpy as sc
import pandas as pd

for disease_status in ["DownSyndrome","Healthy"]:
    for sampletype in ["Femur","Liver"]:


# Load your single-cell RNA-seq object
direc = "/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/"
f = direc + '10X_'+disease_status+'_'+sampletype+'.umap2d.cells_removed.h5ad'
adata = sc.read(f)

a=[columnName[8:] for columnName in adata.obs.columns if 'leiden_v' in columnName]
colName="leiden_v" + str(max([int(x) for x in a if x.isdigit()]))
print("\n * Cell type column name to use: " + colName)

# adata.obs[colName]

sc.tl.rank_genes_groups(adata, groupby=colName, method='wilcoxon')

# dedf = sc.get.rank_genes_groups_df(adata, key="rank_genes_groups_leiden",group="0")
# dedf

# dedf = sc.get.rank_genes_groups_df(adata, key="rank_genes_groups_annotations",group="Neutrophils")
# dedf

# Get a list of unique group names
group_names = adata.obs[colName].unique()
# group_names = adata.obs["numerical_labels"].unique()

mapping = adata.obs.drop_duplicates(subset=['numerical_labels', colName])

# adata.uns["rank_genes_groups"]

# Loop through each group and retrieve the rank genes groups DataFrame
for group_name in group_names:
	rank_genes_groups_df = sc.get.rank_genes_groups_df(adata, group=group_name, key="rank_genes_groups_leiden_filtered")
	subset_df = mapping[mapping['numerical_labels'] == group_name]
	value = subset_df.iloc[0][colName]
	value = value.replace("/","_")
	print(f"Rank genes groups for group: {group_name} {value}")
	rank_genes_groups_df['numerical_labels'] = group_name
	rank_genes_groups_df['cell_type'] = value
	rank_genes_groups_df['Rank'] = range(1,len(rank_genes_groups_df)+1)
	rank_genes_groups_df['type'] = disease_status
	rank_genes_groups_df['organ'] = sampletype
	outfile = "/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/marker_genes/" + disease_status + "." + sampletype + "." + group_name + "." + value + ".csv"
	rank_genes_groups_df.to_csv(outfile, index=False)




group_name="B cells"
rank_genes_groups_df = sc.get.rank_genes_groups_df(adata, group=group_name, key="rank_genes_groups")
subset_df = mapping[mapping['numerical_labels'] == group_name]
value = subset_df.iloc[0][colName]
value = value.replace("/","_")
print(f"Rank genes groups for group: {group_name} {value}")
rank_genes_groups_df['numerical_labels'] = group_name
rank_genes_groups_df['cell_type'] = value
rank_genes_groups_df['Rank'] = range(1,len(rank_genes_groups_df)+1)
rank_genes_groups_df['type'] = disease_status
rank_genes_groups_df['organ'] = sampletype
outfile = "/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/marker_genes/" + disease_status + "." + sampletype + "." + group_name + "." + value + ".csv"
rank_genes_groups_df.to_csv(outfile, index=False)




