#!/usr/bin/env python

import scanpy as sc

for disease_status in ["DownSyndrome","Healthy"]:
	for sampletype in ["Femur","Liver"]:
		headdir="/oak/stanford/groups/smontgom/amarder/t21-proj"
		disease_status="Healthy"
		sampletype="Liver"

		print("1")
		fout="10X_" + disease_status + "_" + sampletype + ".umap.h5ad"
		foutpath=headdir + "/out/data/" + fout
		print("\n * Reading in data..." + foutpath)
		adata=sc.read_h5ad(foutpath)

		print("2")

		adata_subset=adata[:10000,]

		print("3")

		fout="10X_" + disease_status + "_" + sampletype + ".umap.subset.h5ad"
		foutpath=headdir + "/out/data/" + fout
		print("\n * Saving data to file..." + foutpath)
		adata_subset.write(foutpath)

		print("4")

