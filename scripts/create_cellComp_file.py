#!/usr/bin/env python

print("Go")

import scanpy as sc
headdir="/oak/stanford/groups/smontgom/amarder/t21-proj"
disease_status="Healthy"
sampletype="Liver"
suffix=""
for disease_status in ["DownSyndrome","Healthy"]:
	for sampletype in ["Femur","Liver"]:

		print(disease_status + " " + sampletype)

		if disease_status=="Healthy" and sampletype=="Femur":
		    cell_types_to_remove=["0 (to remove)",
		                          "Osteoblasts,4",
		                          "Megakaryocytes,3"]
		    colName="leiden_v7"
		elif disease_status=="DownSyndrome" and sampletype=="Femur":
		    cell_types_to_remove=["Odd PTPRC+ cells",
		                          "34,0",
		                          "Odd NK cells",
		                          "Pre pro B cells,4"]
		    colName="leiden_v12"
		elif disease_status=="Healthy" and sampletype=="Liver":
		    cell_types_to_remove=["38,0",
		                          "31,0",
		                          "Erythroid cells,2,1",
		                          "34",
		                          "Megakaryocytes,2,0",
		                          "36,1",
		                          "36,0,0",
		                          "36,0,1",
		                          "38,1",
		                          "Megakaryocytes,2,1"]
		    colName="leiden_v7"
		elif disease_status=="DownSyndrome" and sampletype=="Liver":
		    cell_types_to_remove=["do not remove any cells"]
		    colName="leiden_v10"
		else:
		    print("Error!")
		    break



		fout="10X_" + disease_status + "_" + sampletype + ".umap"+suffix+".cells_removed.h5ad"
		foutpath=headdir + "/out/data/" + fout

		print("\n * Reading in data..." + foutpath)
		adata=sc.read_h5ad(foutpath)
		# print(adata.obs[["patient","sample","sorting"]].drop_duplicates().sort_values("patient"))

		foutpath=headdir + "/out/data_small/" + "10X_" + disease_status + "_" + sampletype + ".cellComp.csv"
		adata.obs.loc[:,["patient","sample","sorting","numerical_labels",colName,"cell_type_groups"]].to_csv(foutpath)





print("Done!")

