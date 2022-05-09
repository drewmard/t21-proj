#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import scanpy as sc
import pandas as pd
import numpy as np

print("Done.")


# In[1]:


headdir="/oak/stanford/groups/smontgom/amarder/t21-proj"
# disease_status="Healthy"
# sampletype="Liver"
suffix="subset"

for sampletype in ["Liver","Femur"]:
    # suffix="subset"
    suffixDirec=""
    if suffix=="subset":
        suffixDirec="subset"
    else:
        suffixDirec="full"
    os.system("mkdir -p "+headdir + "/out/" + suffixDirec)

    def return_fileName(headdir,disease_status,sampletype,suffix):
        foutpath=""
        if suffix=="subset":
            fout="10X_" + disease_status + "_" + sampletype + ".umap.subset.cells_removed.h5ad"
            foutpath=headdir + "/out/subset/data/" + fout
        else:
            fout="10X_" + disease_status + "_" + sampletype + ".umap2d.cells_removed.h5ad"
            foutpath="/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data" + "/" + fout
        return foutpath

    foutpath=return_fileName(headdir,"Healthy",sampletype,suffix)
    print("\n * Reading in data..." + foutpath)
    adata=sc.read_h5ad(foutpath)

    foutpath=return_fileName(headdir,"DownSyndrome",sampletype,suffix)
    print("\n * Reading in data..." + foutpath)
    adata2=sc.read_h5ad(foutpath)

    # disease_status="DownSyndrome"

    print("Data read completed.")


    # In[2]:





    # In[13]:


    # pd.set_option('display.max_columns', None)
    # print(adata3.obs)
    # pd.set_option('display.max_columns', 20)


    # In[102]:


    a=[columnName[8:] for columnName in adata.obs.columns if 'leiden_v' in columnName]
    colName1="leiden_v" + str(max([int(x) for x in a if x.isdigit()]))
    print("\n * Cell type column name to use: " + colName1)
    adata.obs["leiden_names"]=adata.obs[colName1]

    a=[columnName[8:] for columnName in adata2.obs.columns if 'leiden_v' in columnName]
    colName2="leiden_v" + str(max([int(x) for x in a if x.isdigit()]))
    print("\n * Cell type column name to use: " + colName2)
    adata2.obs["leiden_names"]=adata2.obs[colName2]

    # adata=adata[:1000,]
    # adata2=adata2[:1000,]
    healthy_cells = np.unique(adata.obs["cell_type_groups"])
    ds_cells = np.unique(adata2.obs["cell_type_groups"])

    clusters_for_DE = [cell_type for cell_type in healthy_cells if cell_type in ds_cells]
    subDirec="DE_cell_type_groups"
    os.system("mkdir -p "+headdir + "/out/" + suffixDirec +"/" + subDirec)

    # adata3=adata.concatenate(adata2)


    # In[115]:


    # healthy_cells = np.unique(adata.obs[colName1])
    # ds_cells = np.unique(adata2.obs[colName2])
    # clusters_for_DE = [cell_type for cell_type in healthy_cells if cell_type in ds_cells]

    # cell_type=clusters_for_DE[0]
    for cell_type in clusters_for_DE:

        ind=np.where(adata.obs[["cell_type_groups"]]==cell_type)[0]
        ind2=np.where(adata2.obs[["cell_type_groups"]]==cell_type)[0]
        adata4=adata[ind,].concatenate(adata2[ind2,]) #in cell_types_of_interest

        # ind=np.where(adata3.obs[["cell_type_groups"]]==cell_type)[0]
        # adata4=adata3[ind,] #in cell_types_of_interest
        sc.tl.rank_genes_groups(adata4, 'environment', method='wilcoxon')

        cell_type_filename = cell_type.replace("/","_")
        filename_out=headdir + "/out/" + suffixDirec +"/" + subDirec + "/" + cell_type_filename + ".txt"

        sc.get.rank_genes_groups_df(adata4,group='Down Syndrome').to_csv(filename_out,na_rep='NA',index=False)

    ##########################################
    # # random sampling:
    # adata4.obs["randNumCol"] = np.random.permutation(adata4.obs["environment"].values)

    # # random 50/50 labels
    # import random
    # adata4.obs["randNumCol"] = random.choices(["Down Syndrome","Healthy"],k=adata4.shape[0])

    # sc.tl.rank_genes_groups(adata4, 'randNumCol', method='wilcoxon')
    # sc.get.rank_genes_groups_df(adata4,group="Down Syndrome")


    # In[169]:


    healthy_cells = np.unique(adata.obs[colName1])
    ds_cells = np.unique(adata2.obs[colName2])
    clusters_for_DE = [cell_type for cell_type in healthy_cells if cell_type in ds_cells]
    subDirec="DE_leiden_names"
    os.system("mkdir -p "+headdir + "/out/" + suffixDirec +"/" + subDirec)

    # cell_type=clusters_for_DE[0]
    for cell_type in clusters_for_DE:
        
        ind=np.where(adata.obs[["leiden_names"]]==cell_type)[0]
        ind2=np.where(adata2.obs[["leiden_names"]]==cell_type)[0]
        adata4=adata[ind,].concatenate(adata2[ind2,]) #in cell_types_of_interest

        # ind=np.where(adata3.obs[["leiden_names"]]==cell_type)[0]
        # adata4=adata3[ind,] #in cell_types_of_interest
        sc.tl.rank_genes_groups(adata4, 'environment', method='wilcoxon')

        cell_type_filename = cell_type.replace("/","_")
        filename_out=headdir + "/out/" + suffixDirec +"/" + subDirec + "/" + cell_type_filename + ".txt"

        sc.get.rank_genes_groups_df(adata4,group='Down Syndrome').to_csv(filename_out,na_rep='NA',index=False)

    ##########################################
    # # random sampling:
    # adata4.obs["randNumCol"] = np.random.permutation(adata4.obs["environment"].values)

    # # random 50/50 labels
    # import random
    # adata4.obs["randNumCol"] = random.choices(["Down Syndrome","Healthy"],k=adata4.shape[0])

    # sc.tl.rank_genes_groups(adata4, 'randNumCol', method='wilcoxon')
    # sc.get.rank_genes_groups_df(adata4,group="Down Syndrome")


    # In[ ]:


    # enr_res = gseapy.enrichr(gene_list=glist,
    #                      organism='Human',
    #                      gene_sets='GO_Biological_Process_2018',
    #                      description='pathway',
    #                      cutoff = 0.5)

