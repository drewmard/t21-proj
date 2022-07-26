        if suffix=="subset":
            fout="10X_" + disease_status + "_" + sampletype + ".umap.subset.cells_removed.h5ad"
            foutpath=headdir + "/out/data/" + fout
        else:
            fout="10X_" + disease_status + "_" + sampletype + ".h5ad"
            foutpath="/oak/stanford/groups/smontgom/amarder/data/t21/ScanpyObjects" + "/" + fout

        print("\n * Reading in data..." + foutpath)
        adata=sc.read_h5ad(foutpath)

        print("Data read completed.")

        # In[44]:

        a=[columnName[8:] for columnName in adata.obs.columns if 'leiden_v' in columnName]
        colName="leiden_v" + str(max([int(x) for x in a if x.isdigit()]))
        print("\n * Cell type column name to use: " + colName)

        if disease_status=="Healthy" and sampletype=="Femur":    
            cell_types_to_remove=["0 (to remove)", 
                                  "Osteoblasts,4 (to remove)", 
                                  "Megakaryocytes,3 (to remove)"]
        elif disease_status=="DownSyndrome" and sampletype=="Femur":
            cell_types_to_remove=["Odd PTPRC+ cells (to remove)", 
                                  "34,0 (to remove)", 
                                  "Odd NK cells (to remove)", 
                                  "Pre pro B cells,4 (to remove)"]
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
        elif disease_status=="DownSyndrome" and sampletype=="Liver":
            cell_types_to_remove=["do not remove any cells"]
        else:
            print("Error!")
            break

        print("\n * Indexing the scanpy object...")
        adata=adata[~adata.obs[colName].isin(cell_types_to_remove)]

        # In[48]:

        print("\n * Grouping clusters within broad cell types...")
        print("\n * Creating dictionary...")
        list_of_cell_types=np.unique(adata.obs[[colName]])
        celltypeDict={}
        broad_cell_type="Erythroid"; old_cell_types="Erythroid"
        celltypeDict[broad_cell_type] = [s for s in list_of_cell_types if old_cell_types.upper() in s.upper()]
        broad_cell_type="Myeloid"; old_cell_types="Neutrophils"
        celltypeDict[broad_cell_type] = [s for s in list_of_cell_types if old_cell_types.upper() in s.upper()]
        broad_cell_type="Myeloid"; old_cell_types="Macrophages"
        celltypeDict[broad_cell_type] += [s for s in list_of_cell_types if old_cell_types.upper() in s.upper()]
        broad_cell_type="Myeloid"; old_cell_types="pDCs"
        celltypeDict[broad_cell_type] += [s for s in list_of_cell_types if old_cell_types.upper() in s.upper()]
        broad_cell_type="Myeloid"; old_cell_types="cDC2"
        celltypeDict[broad_cell_type] += [s for s in list_of_cell_types if old_cell_types.upper() in s.upper()]
        broad_cell_type="Myeloid"; old_cell_types="Kupffer"
        celltypeDict[broad_cell_type] += [s for s in list_of_cell_types if old_cell_types.upper() in s.upper()]
        broad_cell_type="Myeloid"; old_cell_types="Monocyte"
        celltypeDict[broad_cell_type] += [s for s in list_of_cell_types if old_cell_types.upper() in s.upper()]
        broad_cell_type="Myeloid"; old_cell_types="Osteoclasts"
        celltypeDict[broad_cell_type] += [s for s in list_of_cell_types if old_cell_types.upper() in s.upper()]
        broad_cell_type="HSC/Progenitors"; old_cell_types="HSC"
        celltypeDict[broad_cell_type] = [s for s in list_of_cell_types if old_cell_types.upper() in s.upper()]
        broad_cell_type="HSC/Progenitors"; old_cell_types="MPP"
        celltypeDict[broad_cell_type] += [s for s in list_of_cell_types if old_cell_types.upper() in s.upper()]
        broad_cell_type="HSC/Progenitors"; old_cell_types="MEMP"
        celltypeDict[broad_cell_type] += [s for s in list_of_cell_types if old_cell_types.upper() in s.upper()]
        broad_cell_type="HSC/Progenitors"; old_cell_types="Granulocyte"
        celltypeDict[broad_cell_type] += [s for s in list_of_cell_types if old_cell_types.upper() in s.upper()]
        broad_cell_type="B cells"; old_cell_types="B cells"
        celltypeDict[broad_cell_type] = [s for s in list_of_cell_types if old_cell_types.upper() in s.upper()]
        broad_cell_type="Mast cells"; old_cell_types="Mast cells"
        celltypeDict[broad_cell_type] = [s for s in list_of_cell_types if old_cell_types.upper() in s.upper()]
        broad_cell_type="Megakaryocytes"; old_cell_types="Megakaryocytes"
        celltypeDict[broad_cell_type] = [s for s in list_of_cell_types if old_cell_types.upper() in s.upper()]
        if disease_status=="DownSyndrome" and sampletype=="Femur":
            broad_cell_type="NK/T cells"; old_cell_types="NK "
            celltypeDict[broad_cell_type] = [s for s in list_of_cell_types if old_cell_types.upper() in s.upper()]
            broad_cell_type="NK/T cells"; old_cell_types="lymphoid"
            celltypeDict[broad_cell_type] += [s for s in list_of_cell_types if old_cell_types.upper() in s.upper()]
        else:
            broad_cell_type="NK/T cells"; old_cell_types="NK "
            celltypeDict[broad_cell_type] = [s for s in list_of_cell_types if old_cell_types.upper() in s.upper()]
        if disease_status=="DownSyndrome" and sampletype=="Femur":
            celltypeDict["Unknown"] = ["Unknown"]
        elif disease_status=="DownSyndrome" and sampletype=="Liver":
            celltypeDict["Unknown"] = ["X"]
        elif disease_status=="Healthy" and sampletype=="Femur":
            celltypeDict["Unknown"] = ["Unknown 1","Unknown 2"]

        list_of_used_cell_types=[]
        for sublist in list(celltypeDict.values()):
            for item in sublist:
                list_of_used_cell_types.append(item)

        list_of_not_used_cell_types=list(set(list_of_cell_types) - set(list_of_used_cell_types))

        celltypeDict["Stroma"] = list_of_not_used_cell_types

        print(celltypeDict)

        print("\n * Saving new cell type labels to 'cell_type_groups' column in the metadata...")
        adata.obs["cell_type_groups"] = np.nan
        def create_cell_type_groups(row):
            word=row[colName] 
            topiclist = [topic for topic in celltypeDict if word in celltypeDict[topic]]
            return(topiclist[0])

        print("\n * About to run create_cell_type_groups...")
        adata.obs["cell_type_groups"] = adata.obs.apply(create_cell_type_groups, axis=1)
        print("\n * create_cell_type_groups running complete.")
        print("\n * New labels created.")

        print("\n * Mapping to numeric labels...")
        le = preprocessing.LabelEncoder()
        le.fit(adata.obs.loc[:,colName])
        adata.obs['numerical_labels'] = le.transform(adata.obs.loc[:,colName]).astype(str)
        print("\n * Numeric labels created.")

