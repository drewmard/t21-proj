import pandas as pd
import scanpy as sc
import os
import numpy as np

headDir="/oak/stanford/groups/smontgom/amarder/t21-proj"
# headDir = "/oak/stanford/groups/smontgom/amarder/t21-proj/out/subset/data"
# suffix = ".umap.subset.cells_removed"
subDirec="data"
suffixDirec="full"

# preprocess 
outdir = headDir + "/out/full/cpdb_Data"
os.system("mkdir -p " + outdir)

disease_status = "Healthy"
env = "Femur"

for disease_status in ["DownSyndrome","Healthy"]:
	for env in ["Liver","Femur"]:

		savepath_counts = outdir + "/10X_" + disease_status + "_" + env + ".norm_count.txt"
		savepath_meta = outdir + "/10X_" + disease_status + "_" + env + ".meta.txt"

		# data after filtering and normalising
		# adata = sc.read(headDir + "/10X_" + disease_status + "_" + env + suffix + ".h5ad")

		fout="10X_" + disease_status + "_" + env + ".umap3d.cells_removed.h5ad"
		foutpath=headDir + "/out/" + suffixDirec + "/" + subDirec + "/" + fout
		print("\n * Reading Seurat data from file..." + foutpath)
		adata=sc.read_h5ad(foutpath)
		print("Done.")

		# we recommend using the normalised non-log transformed data - you can save adata.X before log-transforming in adata.norm for example
		adata.norm = sc.pp.normalize_total(adata, target_sum=1e6, inplace=False)["X"]
		df_expr_matrix = adata.norm.T
		# np.savetxt(savepath_counts, df_expr_matrix, delimiter='\t')

		df_expr_matrix = pd.DataFrame(df_expr_matrix.toarray())
		# Set cell ids as columns
		df_expr_matrix.columns = adata.obs.index
		# Genes should be either Ensembl IDs or gene names
		df_expr_matrix.set_index(adata.raw.var.index, inplace=True) 
		df_expr_matrix.to_csv(savepath_counts,sep='\t')

		# generating meta file
		a=[columnName[8:] for columnName in adata.obs.columns if 'leiden_v' in columnName]
		colName="leiden_v" + str(max([int(x) for x in a if x.isdigit()]))

		df_meta = pd.DataFrame(data={'Cell': list(adata.obs.index), 
		                         'cell_type': list(adata.obs[colName])})
		df_meta.set_index('Cell',inplace=True)
		df_meta.to_csv(savepath_meta, sep='\t')



