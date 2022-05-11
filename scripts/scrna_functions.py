import os
import gzip
import itertools
import umap
import scanpy as sc
import numpy as np
import pandas as pd
#import scrublet as scr
import seaborn as sns
import plotly as py
import plotly.graph_objs as go
import json

from scipy.io import mmwrite
from scipy.sparse import issparse, csr_matrix
from pandas.api.types import is_numeric_dtype
from matplotlib import pylab as plt
from IPython.display import display

from plotly.subplots import make_subplots

class scRNA_functions:
    def __init__(self):

        print(" * Initialising ...")
        py.offline.init_notebook_mode(connected=True)
        sc.settings.verbosity = 0

        self.__annotationMatrix = ''



# ****************************************** Private functions ******************************************
    # Finding the ribosomial genes in the provided annotation matrix
    def __findRibosomialGenes(self, verbose=True):

        rps_genes = self.__annotationMatrix.id[self.__annotationMatrix.symbol.str.startswith('RPS')]
        rpl_genes = self.__annotationMatrix.id[self.__annotationMatrix.symbol.str.startswith('RPL')]

        genes = rps_genes.append(rpl_genes)

        if verbose:
            for rib in genes:
                print(rib, end=" ")
            print("\n")

        return genes


    # Findind the mitochondrial genes in the provided annotation matrix
    def __findMitochondrialGenes(self, verbose=True):

        genes = self.__annotationMatrix.id[self.__annotationMatrix.symbol.str.startswith('MT-')]

        if verbose:
            for mt in genes:
                print(mt, end=" ")
            print("\n")

            for mt in self.__annotationMatrix.symbol[self.__annotationMatrix.symbol.str.match('MT-')]:
                print(mt, end=" ")
            print("\n")

            print("\n")

        return genes


    # Coverting the palette into hexadecimal format
    def __rgb2hex(self, palette):

        colors = []

        if isinstance(palette[0][0], float):
            for i in range(len(palette)):
                r = int(255*palette[i][0])
                g = int(255*palette[i][1])
                b = int(255*palette[i][2])
                colors.append('#%02x%02x%02x' % (r,g,b))
        else:
            for i in range(len(palette)):
                r = int(palette[i][0])
                g = int(palette[i][1])
                b = int(palette[i][2])
                colors.append('#%02x%02x%02x' % (r,g,b))

        return colors


    # Defining and use a simple function to label the plot in axes coordinates
    def __label(self, x, color, label):
        ax = plt.gca()
        ax.text(.6, .2, label, fontweight="bold",
                ha="left", va="center", transform=ax.transAxes)


        # Plotting a single ridge plot
    def __plotSingleRidgePlot(self, df, palette, height=0.75, pdf=None, title=""):

        sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

        g = sns.FacetGrid(df, row=df.columns[1], hue=df.columns[1], aspect=15, height=height, palette=palette, sharey=False)

        # Draw the densities in a few steps
        g.map(sns.kdeplot, df.columns[0], clip_on=False, shade=True, alpha=1, lw=1.5, bw=.2)
        g.map(sns.kdeplot, df.columns[0], clip_on=False, color="w", lw=2, bw=.2)
        g.map(plt.axhline, y=0, lw=2, clip_on=False)

        g.map(self.__label, df.columns[0])

        # Set the subplots to overlap
        g.fig.subplots_adjust(hspace=-.25)

        # Remove axes details that don't play well with overlap
        g.set_titles(title)
        g.set(yticks=[])
        g.set_xlabels(df.columns[0], fontweight="bold")
        #g.despine(bottom=True, left=True)
        plt.tight_layout()
        if pdf is not None:
            pdf.savefig()
        else:
            plt.show(block=False)
        plt.close()


    def __normaliseGroupsDataFrame(self, df):

        df_norm = df.T.copy()

        column_list = []
        sum_list    = []

        columns = df_norm.columns
        print(columns)
        for col in df_norm.columns:
            column_list.append(df_norm[col].values)
            sum_list.append(np.sum(df_norm[col]))

        print(np.sum(df.sum())/2)

        avg = np.sum(df.sum())/2 # Divide by 2 for the DS vs. healthy analysis, but I should generalise this
        new_column = []
        for idx, col in enumerate(df_norm.columns):

            df_norm[col] = df_norm[col].astype("category")
            df_norm[col] = np.round(avg*(column_list[idx]/sum_list[idx]), decimals=0).astype(int)
        display(df_norm)

        return df_norm.T

    # Coverting the palette into hexadecimal format
    def __rgb2hex(self, palette):

        colors = []

        if isinstance(palette[0][0], float):
            for i in range(len(palette)):
                r = int(255*palette[i][0])
                g = int(255*palette[i][1])
                b = int(255*palette[i][2])
                colors.append('#%02x%02x%02x' % (r,g,b))
        else:
            for i in range(len(palette)):
                r = int(palette[i][0])
                g = int(palette[i][1])
                b = int(palette[i][2])
                colors.append('#%02x%02x%02x' % (r,g,b))

        return colors

    def  __barplotGroupsDataFrame(self, df, palette, title="", show_df=False):

        groups = []
        for c in df.columns:
            groups.append(df[c].values)

        if show_df:
            display(df)

        fig = plt.figure(figsize=(12, 5))
        sns.set(font_scale=1.5)
        sns.set_style("white")

        categories = df.index.tolist()
        p = []
        for i, g in enumerate(groups):
            if i > 0:
                p.append(plt.bar(categories, g, 0.8, bottom=df[df.columns[:i]].sum(axis=1), color=palette[i]))
            else:
                p.append(plt.bar(categories, g, 0.8, bottom=0, color=palette[i]))
        sns.despine(offset=10, trim=False)

        plt.xticks(rotation=90)

        plt.legend(tuple(p), tuple(df.columns.tolist()), bbox_to_anchor=(1, 1))
        plt.ylabel('Number of cells')

        plt.title(title)

        plt.show(block=False)

        # Normalisation of barplot in range [0,100]
        new_col     = []
        new_col_val = []
        for idx, c in enumerate(df.columns):

            new_col_val.append(df[c] / df.sum(axis=1) * 100)
            new_col.append(c+" [0,100]")

        df_new = df.copy()
        for idx, col in enumerate(new_col):
            df_new[col] = new_col_val[idx]
        df_new = df_new[new_col]

        if show_df:
            display(df_new)

        for c in df_new.columns:
            groups.append(df_new[c].values)

        groups = []
        for c in df_new.columns:
            groups.append(df_new[c].values)

        fig = plt.figure(figsize=(12, 5))
        sns.set(font_scale=1.5)
        sns.set_style("white")

        categories = df_new.index.tolist()
        p = []
        for i, g in enumerate(groups):
            if i > 0:
                p.append(plt.bar(categories, g, 0.8, bottom=df_new[df_new.columns[:i]].sum(axis=1), color=palette[i]))
            else:
                p.append(plt.bar(categories, g, 0.8, bottom=0, color=palette[i]))
        sns.despine(offset=10, trim=False)

        plt.xticks(rotation=90)

        plt.legend(tuple(p), tuple(df_new.columns.tolist()), bbox_to_anchor=(1, 1))
        plt.ylabel('Number of cells [0, 100]')

        plt.title(title)

        plt.show(block=False)



# ****************************************** Loading functions ******************************************
    # Reading the provided annotation matrix
    def readAnnotationMatrix(self, path, delimiter='\t', verbose=True):

        if verbose:
            print(" * Reading Annotation Matrix, path=%s\n"%path)

        self.__annotationMatrix = pd.read_csv(path, delimiter=delimiter)


    # Reading the excel file with the marker genes
    def readMarkerListFoetal(self, markerListPath, marker_genes, key, variable="id", sheet_name=None):

        print(" * Reading the given marker list")

        df    = pd.read_excel(markerListPath, sheet_name=sheet_name)
        genes = df[variable].tolist()

        print(" * Adding %s marker genes to %s"%(sheet_name, key))
        marker_genes[key] = genes


    # Loading the gene expression matrix provided by Sanger pipeline and building the scanpy object
    def loadSmartSeqSanger(self, path, delimiter='\t', mit=True, ercc=True, rib=True, verbose=True,
            min_genes=0, min_cells=0, title=""):

        if verbose:
            print(" * Loading sample, path=%s"%path)

        rawData = pd.read_csv(path, delimiter=delimiter)
        rawData = rawData.iloc[0:-4:,:]
        indeces = rawData.iloc[:,0]

        names   = {}
        for i,idx in enumerate(indeces):
            names[i] = idx

        cols = {}
        for i,col in enumerate(rawData.columns):
            cols[col] = 'X' + col

        rawData.rename(index=names, inplace=True)
        rawData       = rawData.drop(columns=['Unnamed: 0'])
        rawData.rename(columns=cols, inplace=True)
        countsRawData = rawData[rawData.index.str.match('ENSG0')]

        smartSeq = sc.AnnData(countsRawData.T)

        if min_genes > 0:
            sc.pp.filter_cells(smartSeq, min_genes=min_genes)

        if min_cells > 0:
            sc.pp.filter_genes(smartSeq, min_cells=min_cells)

        smartSeq.obs['sample'] = title
        smartSeq.obs['sample'] = smartSeq.obs['sample'].astype('category')

        smartSeq.obs['n_counts']   = smartSeq.X.sum(1)
        smartSeq.obs['log_counts'] = np.log(smartSeq.obs['n_counts'])
        smartSeq.obs['n_genes']    = (smartSeq.X > 0).sum(1)

        if mit:
            try:
                if verbose:
                    print(" * Calculating the mitochondrial percentage")

                mitGenes = self.__findMitochondrialGenes(verbose=False)
                smartSeq.obs['percent_mito'] = np.sum(rawData.loc[mitGenes].T, axis=1) / smartSeq.obs['n_counts']

            except:
                print(" * Warning: the Annotation Matrix must be provided to calculate the mitochondrial percentage!")

        if ercc:
            erccGenes = rawData.index[rawData.index.str.startswith('ERCC-')]
            smartSeq.obs['ercc_content'] = np.sum(rawData.loc[erccGenes].T, axis=1) / np.sum(rawData.T, axis=1)

            if verbose:
                print(" * Calculating the ercc content")

        if rib:
            try:
                if verbose:
                    print(" * Calculating the ribosomial percentage")

                ribGenes = self.__findRibosomialGenes(verbose=False)
                smartSeq.obs['percent_rib'] = np.sum(rawData.loc[ribGenes].T, axis=1) / smartSeq.obs['n_counts']

            except:
                print(" * Warning: the Annotation Matrix must be provided to calculate the ribosomial percentage!")

        if verbose:
            print(" * Initial SmartSeq2 Object: %d genes across %d single cells"%(smartSeq.n_vars, smartSeq.n_obs))
            print("\n")

        return smartSeq


    # Loading the gene expression matrix provided by CRUCK pipeline and building the scanpy object
    def loadSmartSeqCRUK(self, path, delimiter='\t', mit=True, ercc=True, rib=True, verbose=True,
            min_genes=0, min_cells=0, title=""):

        if verbose:
            print(" * Loading sample, path=%s"%path)

        rawData = pd.read_csv(path, delimiter=delimiter)
        rawData = rawData.iloc[0:-4:,:]
        rawData = rawData[rawData.columns.drop(list(rawData.filter(regex='lostreads')))]
        indeces = rawData.iloc[:,0]

        names   = {}
        for i,idx in enumerate(indeces):
            names[i] = idx

        cols = {}
        for i,col in enumerate(rawData.columns):
            split    = col.split(".")[0:2]
            if len(split) == 1:
                cols[col] = col
            else:
                split[1]  = split[1].replace("_", "-")
                cols[col] = split[0]+split[1]

        rawData.rename(index=names, inplace=True)
        rawData       = rawData.drop(columns=['Unnamed: 0'])

        rawData.rename(columns=cols, inplace=True)
        countsRawData = rawData[rawData.index.str.match('ENSG0')]

        smartSeq = sc.AnnData(countsRawData.T)

        if min_genes > 0:
            sc.pp.filter_cells(smartSeq, min_genes=min_genes)

        if min_cells > 0:
            sc.pp.filter_genes(smartSeq, min_cells=min_cells)

        smartSeq.obs['sample'] = title
        smartSeq.obs['sample'] = smartSeq.obs['sample'].astype('category')

        smartSeq.obs['n_counts']   = smartSeq.X.sum(1)
        smartSeq.obs['log_counts'] = np.log(smartSeq.obs['n_counts'])
        smartSeq.obs['n_genes']    = (smartSeq.X > 0).sum(1)

        if mit:
            try:
                if verbose:
                    print(" * Calculating the mitochondrial percentage")

                mitGenes = self.__findMitochondrialGenes(verbose=False)
                smartSeq.obs['percent_mito'] = np.sum(rawData.loc[mitGenes].T, axis=1) / smartSeq.obs['n_counts']

            except:
                print(" * Warning: the Annotation Matrix must be provided to calculate the mitochondrial percentage!")

        if ercc:
            erccGenes = rawData.index[rawData.index.str.startswith('ERCC-')]
            smartSeq.obs['ercc_content'] = np.sum(rawData.loc[erccGenes].T, axis=1) / np.sum(rawData.T, axis=1)

            if verbose:
                print(" * Calculating the ercc content")

        if rib:
            try:
                if verbose:
                    print(" * Calculating the ribosomial percentage")

                ribGenes = self.__findRibosomialGenes(verbose=False)
                smartSeq.obs['percent_rib'] = np.sum(rawData.loc[ribGenes].T, axis=1) / smartSeq.obs['n_counts']

            except:
                print(" * Warning: the Annotation Matrix must be provided to calculate the ribosomial percentage!")

        if verbose:
            print(" * Initial SmartSeq2 Object: %d genes across %d single cells"%(smartSeq.n_vars, smartSeq.n_obs))
            print("\n")

        return smartSeq

    # Loading the gene expression matrix provided by Sanger pipeline and building the scanpy object
    def load10XSanger(self, path, min_genes=10, min_cells=10, title="", patient="", verbose=True, n_reads=False, dublets=True, folder_mtx="filtered", chromosome_dir=None):

        if folder_mtx != "raw" and folder_mtx != "Foetal_all" and folder_mtx != "filtered" and folder_mtx != "EmptyDrops" and folder_mtx != "Healthy":
            print("Error, select row of filtered matrix folder")
            return

        if verbose:
            print(" * Loading sample, path=%s"%path)

        if folder_mtx=="EmptyDrops" or folder_mtx=="Healthy":
            adata = sc.read_10x_mtx(path+"/outputEmptyDrops", make_unique=dublets)
        elif folder_mtx=="Foetal_all":
            adata = sc.read_10x_mtx(path+"/GRCh38", make_unique=dublets)
        else:
            adata = sc.read_10x_mtx(path+"/%s_feature_bc_matrix"%folder_mtx, make_unique=dublets)

        if min_genes > 0:
            sc.pp.filter_cells(adata, min_genes=min_genes)

        if min_cells > 0:
            sc.pp.filter_genes(adata, min_cells=min_cells)

        adata.X = adata.X.toarray()

        # Keep cells that express more than 1 UMI counts in more than 10 single cells
        adata = adata[(adata.X > 1).sum(1) > 10]

        adata.obs['patient'] = patient
        adata.obs['sample'] = title
        #adata.obs['sample'] = adata.obs['sample'].astype('category')

        #adata.obs['environment'] = 'DownSyndrome'
        #if folder_mtx=="Healthy":
        #    adata.obs['environment'] = 'Healthy'
        #if folder_mtx=="Foetal_all":
        #    del adata.obs["environment"]
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        adata.var["rb"] = (adata.var_names.str.startswith("RPS")) | (adata.var_names.str.startswith("RPL"))

        # Add chromosome numbers if a directory to the JSON is specified
        if chromosome_dir is not None:
            environment = "healthy" if folder_mtx=="Healthy" else "downsyndrome"
            with open(chromosome_dir+'chromosome_numbers_{0}.json'.format(environment)) as json_file:
                chromosomes = json.load(json_file)
            sample = path.split("/")[-1]
            chromo_by_sample = chromosomes[sample] # Dictionary of dictionaries across all samples
            adata.var["chromosome"] = np.array([chromo_by_sample[ID] if ID in chromo_by_sample.keys() else 'Not considered' for ID in adata.var["gene_ids"]])

        sc.pp.calculate_qc_metrics(adata,
                qc_vars=["mt", "rb"],
                percent_top=None,
                log1p=False,
                inplace=True)



        if n_reads:
            try:
                reads            = pd.read_csv(path+"/reads_per_barcode", header=None, sep="\s")
                reads            = reads.rename(columns={0: "n_reads", 1: "barcode"})
                reads["barcode"] = reads["barcode"]+"-1"
                reads            = reads[reads["barcode"].isin(adata.obs_names)]
                reads            = reads.set_index("barcode")
                reads            = reads.to_dict()

                n_reads = []
                for i in adata.obs_names:
                    n_reads.append(reads["n_reads"][i])
                adata.obs["total_reads"] = n_reads
            except:
                print(" * Warning, use samtools to generate 'reads_per_barcode'")

        if verbose:
            print(" * Initial 10X Object: %d genes across %d single cells"%(adata.n_vars, adata.n_obs))
            print("\n")

        return adata


    def __checkCoExpression(self, adata, test, use_raw, verbose):

        if verbose:
            if len(test) == 1:
                print(" * Analysing gene", test, "...")
            else:
                print(" * Analysing genes", test, "...")

        if use_raw == True:
            X = adata.raw.X
        else:
            X = adata.X

        if issparse(X):
            X = X.todense()

        X = np.array(X)

        geneExp = pd.DataFrame(data    = X[:, adata.var.index.isin(test)],
                columns = test,
                index   = adata.obs.index)

        geneExp = geneExp > 0

        selected = [True]*geneExp.shape[0]
        for col in geneExp.columns:
            selected = selected & geneExp[col]

        geneExp = geneExp[selected]

        return 100*geneExp.shape[0]/adata.n_obs, selected



# ****************************************** Computation functions ******************************************
    # Getting the annotation matrix as dataframe
    def getAnnotationMatrix(self):
        return self.__annotationMatrix


    # Mapping Ensembl IDs to the corresponding gene names
    def mappingEnsemblToAnnotated(self, adata, verbose=True, raw_data=True):

        if verbose:
            print(" * Applying the mapping based on the provided Annotation Matrix ...")

        IDs = list(set(self.__annotationMatrix.id).intersection(set(adata.var.index.tolist())))
        adata.var["Ensembl"] = adata.var.index.tolist()

        filtered = self.__annotationMatrix[self.__annotationMatrix.id.isin(IDs)]
        names    = pd.Series(filtered.symbol.values, index=filtered.id).to_dict()

        adata.var.rename(index=names, inplace=True)

        if verbose:
            print(" * Gene names renamed from ensembl to annotated gene name")

        if raw_data:
            try:
                adata.raw.var.rename(index=names, inplace=True)

                if verbose:
                    print(" * Gene names renamed (raw data) from ensembl to annotated gene name")

            except:
                pass


    # Mapping gene names to the corresponding Ensembl IDs
    def mappingAnnotatedToEnsembl(self, adata, verbose=True, raw_data=True):

        if verbose:
            print(" * Applying the mapping based on the provided Annotation Matrix ...")

        IDs = list(set(self.__annotationMatrix.symbol).intersection(set(adata.var.index.tolist())))

        filtered = self.__annotationMatrix[self.__annotationMatrix.symbol.isin(IDs)]
        names    = pd.Series(filtered.id.values, index=filtered.symbol).to_dict()

        if verbose:
            print(" * Gene names renamed from annotated to ensembl gene name")

        if raw_data:
            try:

                adata.raw.var.rename(index=names, inplace=True)

                if verbose:
                    print(" * Gene names renamed (raw data) from ensembl to annotated gene name")

            except:
                pass


    # Calculating the cell cycle scores
    def CellCycleScoring(self, path, adata, verbose=True):

        cell_cycle_genes = [x.strip() for x in open(path)]

        s_genes   = cell_cycle_genes[:43]
        g2m_genes = cell_cycle_genes[43:]


        if verbose:
            print(" * Searching s genes")
            print(self.__annotationMatrix)

        s_genesIDs = []
        for gene in s_genes:
            A = self.__annotationMatrix.id[self.__annotationMatrix["symbol"] == gene]
            if len(A.values) > 0:
                for i in range(0, len(A.values)):
                    s_genesIDs.append(A.values[i])

        if verbose:
            print(" * Searching g2m genes")

        g2m_genesIDs = []
        for gene in g2m_genes:
            A = self.__annotationMatrix.id[self.__annotationMatrix["symbol"] == gene]
            if len(A.values) > 0:
                for i in range(0, len(A.values)):
                    g2m_genesIDs.append(A.values[i])

        sc.tl.score_genes_cell_cycle(adata, s_genes=s_genesIDs, g2m_genes=g2m_genesIDs)


        if verbose:
            print(" * Cell cycle phase calculated")


    # Calculating the cell cycle scores
    def CellCycleScoring10x(self, path, adata, verbose=True):

        cell_cycle_genes = [x.strip() for x in open(path)]

        s_genes   = cell_cycle_genes[:43]
        g2m_genes = cell_cycle_genes[43:]

        sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

        if verbose:
            print(" * Cell cycle phase calculated")


    # Computing UMAP using the umap-learn package
    def computeUMAP(self, dataset1=None, dataset2=None, n_neighbors=5, fitAll=False):

        ump  = umap.UMAP(n_neighbors=n_neighbors, random_state=42)

        if (dataset1 is not None) and (dataset2 is not None):

            if fitAll:
                print(" * Fitting on all datasets")
                allData = np.append(dataset1, dataset2, axis=0)
                ump.fit(allData)
            else:
                print(" * Fitting only on dataset1")
                ump.fit(dataset1)

            return ump.transform(dataset1), ump.transform(dataset2)

        elif dataset1 is not None:
            print(" * Only dataset1")
            ump.fit(dataset1)
            return ump.transform(dataset1)

        elif dataset2 is not None:
            print(" * Only dataset2")
            ump.fit(dataset2)
            return ump.transform(dataset2)

    def mergeGroups(self, adata, newName=None, labels=None, group_by=None, verbose=False):

        toExit = False

        if labels is None:
            print(" * Error! None for 'labels' is invalid. Please provide a dictionary with old and new categories")
            toExit = True

        # Checking if the cluster is not None
        if group_by is None:
            print(" * Error! None for 'group_by' is an invalide column of .obs, please provide one of the column in .obs")
            toExit = True

        # Checking if the desidered group is in the .obs structure
        if group_by not in adata.obs.columns.tolist():
            print(" * Error! %s for 'group_by' is an invalide column of .obs, please provide one of the column in .obs"%group_by)
            toExit = True

        if toExit:
            return

        if newName is None:
            newName = group_by + "_merged"

        if verbose:
            print(" * Merging groups and renaming (column %s) ... "%newName)

        adata.obs[newName] = adata.obs[group_by]

        for label in labels:

            old  = adata.obs[newName][adata.obs[group_by] == label]
            new  = str(labels[label])

            if verbose:
                print(" \t* Old label %s - new label %s"%(label, new))

            adata.obs[newName].replace(old.values.tolist(), new, inplace=True)
            #adata.obs[newName].replace(old, new, inplace=True)

        adata.obs[newName] = adata.obs[newName].astype('category')

        if verbose:
            print(" *  Groups merged")


    # Merging clusters using a user defined dictionary
    def mergeClusters(self, adata, newName=None, labels=None, cluster=None, verbose=False):

        toExit = False
        # Checking if the cluster is not None
        if cluster is None:
            print(" * Error! None for 'cluster' is an invalide column of .obs, please provide one of the column in .obs")
            toExit = True

        # Checking if the desidered group is in the .obs structure
        if cluster not in adata.obs.columns.tolist():
            print(" * Error! %s for 'cluster' is an invalide column of .obs, please provide one of the column in .obs"%cluster)
            toExit = True

        if toExit:
            return

        self.mergeGroups(adata, newName = newName, labels=labels, group_by=cluster, verbose=verbose)


    # Looking at the expression of genes provided by the user cluster by cluster
    def checkMarkerGenesPerCluster(self, adata, marker_cluster, clusterName=None, group_by=None, sortMarkerGenes=True):

        toExit = False

        # Checking if the group_by is not None
        if group_by is None:
            print(" * Warning! None for group_by is an invalide column of .obs, please provide one of the column in .obs")
            toExit = True

        # Checking if the desidered group is in the .obs structure
        if group_by not in adata.obs.columns.tolist():
            print(" * Error! %s for group_by is an invalide column of .obs, please provide one of the column in .obs"%group_by)
            toExit = True

        if toExit:
            return

        clusters = adata.obs[group_by].cat.categories

        matrix = np.zeros((len(marker_cluster), len(adata.obs[group_by].cat.categories)))

        for row, gene in enumerate(marker_cluster):
            for col, name in enumerate(clusters):
                cl            = adata[adata.obs[group_by] == name]
                geneExp_cl    = cl.raw.X[:, cl.var.index == gene]
                geneExp_cl    = np.array(geneExp_cl[geneExp_cl > 0]).flatten()
                percentage_cl = 100*(len(geneExp_cl)/cl.n_obs)

                matrix[row][col] = percentage_cl

        genes = pd.DataFrame(data=matrix, index=marker_cluster,  columns=clusters)

        if sortMarkerGenes:
            if clusterName is not None:
                genes = genes.sort_values(by=[clusterName], ascending=False)

        for gene in genes.index:
            print("* Gene %s is expressed in:"%('\033[1m'+gene+'\033[0m'))

            g = genes[genes.index == gene]
            g = g.sort_values(by=gene, ascending=False, axis=1)

            for c in g.columns:

                cl         = adata[adata.obs[group_by] == c]
                geneExp_cl = cl.raw.X[:, cl.var.index == gene]
                geneExp_cl = np.array(geneExp_cl[geneExp_cl > 0]).flatten()
                n_cells_cl = len(geneExp_cl)

                if c == clusterName:
                    print('\033[1m'+"\t -> %3.2f%% (%d out of %d) of cells of cluster %s"%(g[c], n_cells_cl, cl.n_obs, c)+'\033[0m')
                else:
                    print("\t -> %3.2f%% (%d out of %d) of cells of cluster %s"%(g[c], n_cells_cl, cl.n_obs, c))

            display(g)


    def run_subclustering(self, adata, algorithm="kmeans", restrict_to=None, n_clusters=[2], space="X_umap", key_added="subcluster"):

        groupby = restrict_to[0]

        if isinstance(restrict_to[1], list):
            pass
        else:
            restrict_to[1] = [restrict_to[1]]

        if isinstance(n_clusters, list):
            pass
        else:
            n_clusters = [n_clusters]

        cluster_names = {}
        for idx,cluster in enumerate(restrict_to[1]):

            if algorithm =="kmeans":
                from sklearn.cluster import KMeans

                print("* Applying KMeans to %s cluster using %d n_clusters"%(cluster, n_clusters[idx]))
                X      = np.array(adata[adata.obs[groupby]==cluster].obsm[space])
                kmeans = KMeans(n_clusters=n_clusters[idx], random_state=42).fit(X)
                labels = kmeans.labels_

            indeces = adata[adata.obs[groupby]==cluster].obs.index.tolist()
            values  = [cluster+",%s"%i for i in labels]
            cluster_names.update(dict(zip(indeces, values)))

        indeces = adata[np.logical_not(adata.obs[groupby].isin(restrict_to[1]))].obs.index.tolist()
        values  = adata[np.logical_not(adata.obs[groupby].isin(restrict_to[1]))].obs[groupby].tolist()
        cluster_names.update(dict(zip(indeces, values)))

        # if key_added is already in the dataframe, remove it and add it again!
        if key_added in adata.obs.columns.tolist():
            adata.obs.drop(columns=[key_added], inplace=True)

        df = pd.DataFrame.from_dict(cluster_names, orient='index', columns=[key_added])
        adata.obs = pd.concat([adata.obs, df], axis=1, sort=False)

        return adata


    # Removing from the expression matrix possible doublets
    '''
    def scrublet_doublet_removal_10X(self, adata):
        scrub = scr.Scrublet(adata.X, expected_doublet_rate=0.06)
        adata.obs['doublet_scores'], adata.obs['predicted_doublets'] = scrub.scrub_doublets(min_counts=2,
                min_cells=3,
                min_gene_variability_pctl=85,
                n_prin_comps=30)
        adata.obs['predicted_doublets'] = adata.obs.predicted_doublets.astype('category')

        try:
            print('Number of doublets: %i' %(adata.obs['predicted_doublets'].value_counts()[1]))
        except:
            print('Number of doublets: 0')

        return adata
    '''
    def checkCoexpression(self,
            adata,
            genes1     = None,
            genes2     = None,
            use_raw    = True,
            space      = "X_umap",
            components = ['1,2'],
            height     = 4,
            width      = 4,
            best_shape = False,
            plot       = True,
            heatmap    = False,
            verbose    = True,
            expression = False,
            pdf        = None
            ):

        toExit = False
        if genes1 is None and genes2 is None:
            print("* Warning! Provide two list of genes (genes1 and genes2) or a list of genes (either genes1 or genes2)")
            toExit = True

        if space is None and plot:
            print(" * Warning! None for space is an invalide space, please provide one of the space in .obsm")
            toExit = True

        if space not in list(adata.obsm) and plot:
            print(" * Error! '%s' space is not available, please calculate the '%s' space and save it in .obsm['%s']"%(space,
                space,
                space))
            print("\t * The available spaces are:", list(adata.obsm))
            toExit = True

        if use_raw:
            try:
                adata.raw.X
            except:
                print(" * Error! '.raw' is not available. Either set 'use_raw = False' or create '.raw'")
                toExit = True

        genes = genes1 + genes2
        genes = list(set(genes).difference(set(adata.var.index)))

        if len(genes) > 0:
            print(" * Error! Genes", genes, "are not expressed in the provided dataset")
            toExit = True

        if plot:
            try:
                ax_label = space.split("_")[1]
            except:
                ax_label = space

            components = components[0].split(",")

        if toExit:
            return None

        if plot:
            try:
                adata.obsm[space][:, [int(x) - 1 for x in components]]
            except:
                print(" * Error! The components", components, "are not available for space '%s'"%space)
                return None

        if genes1 is None or genes2 is None:
            if genes1 is None:
                test = genes2
            else:
                test = genes1

            if len(test) == 1 and expression == False:
                return None

            if plot:
                f, axs = plt.subplots(1, 1, figsize=(width, height))
                sns.despine(offset=10, trim=False)
                sns.set(font_scale=1)
                sns.set_style("white")

            result, selected = self.__checkCoExpression(adata, test, use_raw, verbose)

            if verbose:
                if len(test) == 1:
                    print("\t The gene is expressed in the %.2f%% of cells"%result)
                    print()
                else:
                    print("\t The genes are coexpressed in the %.2f%% of cells"%result)
                    print()

            if plot:
                sns.set_style("white")
                df = pd.DataFrame(data    = adata.obsm[space][:, [int(x) - 1 for x in components]],
                        index   = adata.obs.index,
                        columns = [ax_label.upper() + "_" + x for x in components])

                zeros = df[np.logical_not(selected)]
                ones  = df[selected]

                if zeros.shape[0] > 0:
                    sns.scatterplot(zeros[zeros.columns.tolist()[0]],
                            zeros[zeros.columns.tolist()[1]],
                            c=["#808080"], ax=axs)
                    if ones.shape[0] > 0:
                        sns.scatterplot(ones[ones.columns.tolist()[0]],
                                ones[ones.columns.tolist()[1]],
                                c=["#4C0013"], ax=axs)

                        axs.set_title(test)

                plt.tight_layout()
                if pdf is not None:
                    pdf.savefig()
                else:
                    plt.show(block=False)
                plt.close("all")

                return None

        else:
            if len(genes1) == 1 and len(genes2) == 1 and genes1 == genes2 and expression == False:
                return None

            tests = [genes1, genes2]
            tests = list(itertools.product(*tests))

            if best_shape:
                if len(genes1) < len(genes2):
                    n_cols = len(genes1)
                    n_rows = len(genes2)
                else:
                    n_cols = len(genes2)
                    n_rows = len(genes1)
            else:
                n_cols = len(genes1)
                n_rows = len(genes2)

            if plot:
                f, axs = plt.subplots(n_rows, n_cols, figsize=(width*n_cols, height*n_rows))
                sns.despine(offset=10, trim=False)
                sns.set(font_scale=1)
                sns.set_style("white")

            results = pd.DataFrame(data    = 0,
                    index   = genes1,
                    columns = genes2).astype(float)

            id_cells = pd.DataFrame(data    = 0,
                    index   = genes1,
                    columns = genes2)

            for r in range(n_rows):
                for c in range(n_cols):

                    idx  = r*n_cols+c
                    test = tests[idx]

                    if test[0] == test[1]:
                        test = [test[0]]
                    else:
                        test = list(test)

                    if len(test) == 1 and expression == False:
                        continue

                    result, selected = self.__checkCoExpression(adata, test, use_raw, verbose)

                    if verbose:
                        if len(test) == 1:
                            print("\t The gene is expressed in %.2f%% of cells"%result)
                            print()
                        else:
                            print("\t The genes are coexpressed in %.2f%% of cells"%result)
                            print()

                    if len(test) == 1:
                        results[test[0]][test[0]]  = result
                        id_cells[test[0]][test[0]] = adata.obs.index[selected].tolist()
                    else:
                        results[test[1]][test[0]]  = result
                        id_cells[test[1]][test[0]] = adata.obs.index[selected].tolist()

                    if plot:
                        sns.set_style("white")

                        df = pd.DataFrame(data    = adata.obsm[space][:, [int(x) - 1 for x in components]],
                                index   = adata.obs.index,
                                columns = [ax_label.upper() + "_" + x for x in components])

                        zeros = df[np.logical_not(selected)]
                        ones  = df[selected]

                        if len(tests) == 1:
                            if zeros.shape[0] > 0:
                                sns.scatterplot(zeros[zeros.columns.tolist()[0]],
                                        zeros[zeros.columns.tolist()[1]],
                                        c=["#808080"], ax=axs)
                                if ones.shape[0] > 0:
                                    sns.scatterplot(ones[ones.columns.tolist()[0]],
                                            ones[ones.columns.tolist()[1]],
                                            c=["#4C0013"], ax=axs)

                                    axs.set_title(test)
                        else:
                            if zeros.shape[0] > 0:
                                sns.scatterplot(zeros[zeros.columns.tolist()[0]],
                                        zeros[zeros.columns.tolist()[1]],
                                        c=["#808080"], ax=axs[r,c])
                                if ones.shape[0] > 0:
                                    sns.scatterplot(ones[ones.columns.tolist()[0]],
                                            ones[ones.columns.tolist()[1]],
                                            c=["#4C0013"], ax=axs[r,c])

                                    axs[r,c].set_title(test)

            if plot:
                plt.tight_layout()
                if pdf is not None:
                    pdf.savefig()
                else:
                    plt.show(block=False)
                plt.close("all")

            if heatmap:
                sns.heatmap(results, annot=True, fmt=".2f", linewidths=.5, vmin=0, vmax=100, cmap="viridis")
                plt.tight_layout()
                if pdf is not None:
                    pdf.savefig()
                else:
                    plt.show(block=False)
                plt.close("all")

            return results, id_cells


    def write10x(self, adata, path="custom_feature_bc_matrix", layer=None, overwrite=True, compression=True):

        if overwrite == False & os.path.exists(path):
            print("The folder '%s' already exists! If you want to overwrite it, please set 'overwrite=True'"%path)
            return
        else:
            if not os.path.exists(path):
                os.makedirs(path)

        if layer is None:
            counts = adata.X.T
        else:
            try:
                counts = adata.layers[layer].T
            except:
                print("The layer '%s' does not exist!"%layer)
                return

        # MATRIX
        print("Writing 'matrix.mtx' to folder '%s'..."%path)
        if issparse(counts) == False:
            counts = csr_matrix(counts)

        mmwrite(path+os.sep+"matrix.mtx", counts, field="integer")

        if compression:
            with open(path+os.sep+"matrix.mtx", 'rb') as f_in:
                f_out = gzip.open(path+os.sep+"matrix.mtx.gz", 'wb')
                f_out.write(f_in.read())
                f_out.close()
            os.remove(path+os.sep+"matrix.mtx")


        # BARCODES
        print("Writing 'barcodes.tsv.gz' to folder '%s'..."%path)
        barcodes = adata.obs.index.tolist()
        np.savetxt(path+os.sep+"barcodes.tsv", barcodes, fmt="%s")

        if compression:
            with open(path+os.sep+"barcodes.tsv", 'rb') as f_in:
                f_out = gzip.open(path+os.sep+"barcodes.tsv.gz", 'wb')
                f_out.write(f_in.read())
                f_out.close()
            os.remove(path+os.sep+"barcodes.tsv")


        # FEATURES
        print("Writing 'features.tsv' to folder '%s'..."%path)
        gene_ids      = adata.var["gene_ids"].tolist()
        gene_names    = adata.var.index.tolist()
        feature_types = adata.var["feature_types"].tolist()
        features      = list(zip(gene_ids, gene_names, feature_types))
        np.savetxt(path+os.sep+"features.tsv", features, fmt="%s", delimiter='\t')

        if compression:
            with open(path+os.sep+"features.tsv", 'rb') as f_in:
                f_out = gzip.open(path+os.sep+"features.tsv.gz", 'wb')
                f_out.write(f_in.read())
                f_out.close()
            os.remove(path+os.sep+"features.tsv")


        print("Succefully written files!")



# ****************************************** Metadata functions ******************************************
    # Loading the provided metadata
    def loadMetadata(self, path, adata, delimiter=",", verbose=True, rows=10):

        adataCopy = adata.copy()

        if verbose:
            print(" * Loading metadata, path=%s"%path)

        metadata = pd.read_csv(path, delimiter=delimiter)
        metadata.rename(index=metadata.CELL_NAME, inplace=True)
        metadata = metadata.drop(columns=['Unnamed: 0'])

        adataCopy.obs = pd.concat([adata.obs, metadata.loc[adata.obs.index.tolist()]], axis=1)
        adataCopy.obs.sort_index(axis=1, inplace=True)

        if verbose:
            display(adataCopy.obs.head(rows))

        return adataCopy


    # Assigning the origin to cells
    def assignOrigin(self, adata, namesIn=[], namesOut=[], verbose=True):

        adataCopy = adata.copy()
        adataCopy.obs['origin'] = adataCopy.obs.COMMON_NAME

        for idx,name in enumerate(namesIn):
            names = adataCopy.obs.origin[adataCopy.obs.origin.str.contains(name)]
            adataCopy.obs['origin'].replace(names, namesOut[idx], inplace=True)
        adataCopy.obs['origin'] = adataCopy.obs.origin.astype('category')

        if verbose:
            print(adataCopy.obs.origin.cat.categories.tolist())

        return adataCopy


    # Assigning the gate to cells
    def assignGate(self, adata, namesIn=[], namesOut=[], verbose=True):

        adataCopy = adata.copy()
        adataCopy.obs['gate'] = adataCopy.obs.COMMON_NAME

        for idx,name in enumerate(namesIn):
            names = adataCopy.obs.gate[adataCopy.obs.gate.str.contains(name)]
            adataCopy.obs['gate'].replace(names, namesOut[idx], inplace=True)
        adataCopy.obs['gate'] = adataCopy.obs.gate.astype('category')

        if verbose:
            print(adataCopy.obs.gate.cat.categories.tolist())

        return adataCopy


    # Assigning the cell type to cells
    def assignCellType(self, adata, namesIn=[], namesOut=[], verbose=True):

        adataCopy = adata.copy()
        adataCopy.obs['cell_type'] = adataCopy.obs.COMMON_NAME

        for idx,name in enumerate(namesIn):
            names = adataCopy.obs.cell_type[adataCopy.obs.cell_type.str.contains(name)]
            adataCopy.obs['cell_type'].replace(names, namesOut[idx], inplace=True)
        adataCopy.obs['cell_type'] = adataCopy.obs.cell_type.astype('category')

        if verbose:
            print(adataCopy.obs.cell_type.cat.categories.tolist())

        return adataCopy



# ****************************************** Plotting functions using the same style ******************************************
    # Plotting a group of violin plots based on a variable
    def plotViolinVariableGroup(self, adata, variable, group_by=None, log=False, cut=0, pointSize=4,
            width=8, height=8, ax=None, use_raw=False, rotation=None, pdf=None, title=""):

        toExit = False

        # Checking if the group_by is not None
        if variable is None:
            print(" * Warning! None for variable is an invalide column of .obs, please provide one of the column in .obs")
            toExit = True

        # Checking if the desidered group is in the .obs structure
        if group_by is not None and group_by not in adata.obs.columns.tolist():
            print(" * Error! %s for group_by is an invalide column of .obs, please provide one of the column in .obs"%group_by)
            toExit = True

        if toExit:
            return

        if ax is not None:
            sc.pl.violin(adata, variable, groupby=group_by, log=log, jitter=0.4, size=pointSize,
                    multi_panel=False, ax=ax, show=False, cut=0, use_raw=use_raw, rotation=rotation)

        else:

            f, axs = plt.subplots(1,1,figsize=(width,height))
            sns.set(font_scale=1.5)
            sns.set_style("white")
            print("Making violin plots")
            sc.pl.violin(adata, variable, groupby=group_by, log=log, jitter=0.4, size=0,
                    multi_panel=False, ax=axs, show=False, cut=0, use_raw=use_raw, rotation=rotation)

            #sns.despine(offset=10, trim=False)
            plt.tight_layout()
            plt.title(title)
            if pdf is not None:
                pdf.savefig()
            else:
                plt.show(block=False)
            plt.close()

    # Plotting a scatter plot based on two variables (x and y)
    def plotScatter(self, adata, variable, x="n_genes", y="ercc_content", palette=["gray"],
            pointSize=150, width=8, height=8, ax=None, colormap=None, pdf=None):

        if ax is not None:
            if colormap is None:
                sc.pl.scatter(adata, x=x, y=y, color=variable, palette=palette,
                        size=pointSize, ax=ax, show=False)
            else:
                sc.pl.scatter(adata, x=x, y=y, color=variable, color_map=colormap,
                        size=pointSize, ax=ax, show=False)


        else:

            plt.rcParams['figure.figsize']=(10,10)
            sns.set(font_scale=1.5)
            sns.set_style("white")

            if colormap is None:
                sc.pl.scatter(adata, x=x, y=y, color=variable, palette=palette,
                        size=pointSize, show=False)
            else:
                sc.pl.scatter(adata, x=x, y=y, color=variable, color_map=colormap,
                        size=pointSize, show=False)

                #sns.despine(offset=10, trim=False)
            plt.tight_layout()
            if pdf is not None:
                pdf.savefig()
            else:
                plt.show(block=False)
            plt.close()


    # Plotting a histogram plot based on a variable
    def plotHistrogram(self, variable, kde=False, bins=50, color="black", width=8, height=8, ax=None):

        if ax is not None:

            sns.distplot(variable, bins=bins, kde=kde, color=color)

        else:

            f, axs = plt.subplots(1,1,figsize=(width,height))
            sns.set(font_scale=1.5)
            sns.set_style("white")

            sns.distplot(variable, bins=bins, kde=kde, ax=axs, color=color)

            sns.despine(offset=10, trim=False)
            plt.show(block=False)


    # Plotting a scatter plot for each cluster with the name of the top n marker genes
    def plotGenesRankGroups(self, adata, n_genes=20, key="rank_genes_groups", ncols=2, fontsize=12, pdf=None):

        sns.set(font_scale=1.5)
        sns.set_style("white")

        sc.pl.rank_genes_groups(adata,
                n_genes=n_genes,
                ncols=ncols,
                key=key,
                fontsize=fontsize,
                show=False)

        sns.despine(trim=False)
        if pdf is not None:
            pdf.savefig()
            plt.close()
        #else:
        #    plt.show(block=False)


    # Plotting two UMAP components colouring the cells based on a variable
    def plotUMAP(self, adata, variable=None, pointSize=150, width=8, height=8, palette=None,
            components='1,2', ax=None, colormap=None, use_raw=False, use_title=None, pdf=None):

        toExit = False

        # Checking if the variable is not None
        if variable is None:
            print(" * Warning! None for variable is an invalide column of .obs, please provide one of the column in .obs")
            toExit = True

        if toExit:
            return

        # Checking and palette to be used
        if palette is None:

            try:
                palette = adata.uns[variable+"_colors"]
            except:
                palette = sns.color_palette("tab20c")
        else:
            if isinstance(palette[0], tuple):
                palette = self.__rgb2hex(palette)

        if use_title is None:
            title = variable
        else:
            title = use_title

        # Using the provided axis
        if ax is not None:

            if variable is None:
                sc.pl.umap(adata, size=pointSize, palette=palette, components=components, ax=ax, show=False, use_raw=use_raw, title=title)

            else:
                if colormap is None:
                    sc.pl.umap(adata, color=variable, size=pointSize, palette=palette, components=components, ax=ax, show=False, use_raw=use_raw, title=title)
                else:
                    sc.pl.umap(adata, color=variable, size=pointSize, cmap=colormap, components=components, ax=ax, show=False, use_raw=use_raw, title=title)

        # Create a figure to plot the UMAP
        else:
            f, axs = plt.subplots(1,1,figsize=(width,height))
            sns.set(font_scale=1.5)
            sns.set_style("white")

            if variable is None:
                sc.pl.umap(adata, size=pointSize, palette=palette, components=components, ax=axs, show=False, use_raw=use_raw, title=title)

            else:
                if colormap is None:
                    sc.pl.umap(adata, color=variable, size=pointSize, palette=palette, components=components, ax=axs, show=False, use_raw=use_raw, title=title)
                else:
                    sc.pl.umap(adata, color=variable, size=pointSize, cmap=colormap, components=components, ax=axs, show=False, use_raw=use_raw, title=title)

            #sns.despine(offset=10, trim=False)
            plt.tight_layout()
            if pdf is not None:
                pdf.savefig()
            else:
                plt.show(block=False)
            plt.close()

    def plot3DContinuous(self, adata, space=None, title="", scale=None, poinSize=6, showgrid=True, write_to=None):

        toExit = False

        # Checking if the desidered space is not None
        if space is None:
            print(" * Warning! None for space is an invalide space, please provide one of the space in .obsm")
            toExit = True

        if scale is None:
            print(" * Warning! Please provide values for a relavant colour scale")
            toExit = True

        # Checking if the desidered space is in the .obsm structure
        if space not in list(adata.obsm):
            print(" * Error! %s space is not available, please calculate the %s space and save it in .obsm['%s']"%(space, space, space))
            print(" * The available are:", list(adata.obsm))
            toExit = True

        if toExit:
            return

        try:
            ax_label = space.split("_")[1]
        except:
            ax_label = space

        listGroups = []

        # Retrieve the first 3 components
        data  = adata.obsm[space]
        continuous = adata.obs[scale].values

        x = data[:,0]
        y = data[:,1]
        z = data[:,2]

        trace = go.Scatter3d(x=x,
                y=y,
                z=z,
                mode='markers',
                marker=dict(size=poinSize, color=continuous, colorbar=dict(title=scale), colorscale="Viridis"),
                name=""
                )
        listGroups.append(trace)

        # Setting plotly layout
        layout = go.Layout(scene = dict(xaxis = dict(title=ax_label.upper() + " 1",
            showbackground=False,
            showgrid=showgrid,
            gridcolor="rgb(218,218,218)",
            zeroline=False,
            showline=True,
            ticks='',
            showticklabels=False
            ),
            yaxis = dict(title=ax_label.upper() + " 2",
                showbackground=False,
                gridcolor="rgb(218,218,218)",
                showgrid=showgrid,
                zeroline=False,
                showline=True,
                ticks='',
                showticklabels=False
                ),
            zaxis = dict(title=ax_label.upper() + " 3",
                showbackground=False,
                gridcolor="rgb(218,218,218)",
                showgrid=showgrid,
                zeroline=False,
                showline=True,
                ticks='',
                showticklabels=False
                )
            ),
            height=600,
            legend=dict(orientation="h", y=.05),
            margin=dict(r=0, b=0, l=0, t=0)
            )
        fig = go.Figure(data=listGroups, layout=layout)
        if write_to is not None:
            fig.write_image(write_to+".png") # Save a regular file
            fig.write_html(write_to+".html") # Save an HTML file
        else:
            py.offline.iplot(fig, filename=title)

    # Plotting the data in a 3D space using plotly library. The components must be in the .obsm structure.
    def plot3D(self, adata, group_by=None, space=None, title="",
            poinSize=6, showgrid=True, palette=None, write_to=None):

        toExit = False

        # Checking if the desidered space is not None
        if space is None:
            print(" * Warning! None for space is an invalide space, please provide one of the space in .obsm")
            toExit = True

        # Checking if the group_by is not None
        if group_by is None:
            print(" * Warning! None for group_by is an invalide column of .obs, please provide one of the column in .obs")
            toExit = True

        # Checking if the desidered group is in the .obs structure
        if group_by not in adata.obs.columns.tolist():
            print(" * Error! %s for group_by is an invalide column of .obs, please provide one of the column in .obs"%group_by)
            toExit = True

        # Checking if the desidered space is in the .obsm structure
        if space not in list(adata.obsm):
            print(" * Error! %s space is not available, please calculate the %s space and save it in .obsm['%s']"%(space, space, space))
            print(" * The available are:", list(adata.obsm))
            toExit = True

        if toExit:
            return

        try:
            ax_label = space.split("_")[1]
        except:
            ax_label = space

        #keys  = adata.obs[group_by].cat.categories.tolist()
        keys  = adata.obs[group_by].unique().tolist()
        print(keys)

        # Checking and palette to be used
        if palette is None:

            try:
                palette = adata.uns[group_by+"_colors"]
            except:
                palette = sns.color_palette("tab20c", n_colors=len(keys))
        else:
            if isinstance(palette[0], tuple):
                palette = self.__rgb2hex(palette)



        listGroups = []

        # Retrieve the first 3 components
        for idx,key in enumerate(keys):

            data  = adata[adata.obs[group_by] == str(key)]
            data  = data.obsm[space]

            x = data[:,0]
            y = data[:,1]
            z = data[:,2]

            trace = go.Scatter3d(x=x,
                    y=y,
                    z=z,
                    mode='markers',
                    marker=dict(size=poinSize, color=palette[idx%len(palette)]),
                    name=str(key)
                    )
            listGroups.append(trace)

        # Setting plotly layout
        layout = go.Layout(scene = dict(xaxis = dict(title=ax_label.upper() + " 1",
            showbackground=False,
            showgrid=showgrid,
            gridcolor="rgb(218,218,218)",
            zeroline=False,
            showline=True,
            ticks='',
            showticklabels=False
            ),
            yaxis = dict(title=ax_label.upper() + " 2",
                showbackground=False,
                gridcolor="rgb(218,218,218)",
                showgrid=showgrid,
                zeroline=False,
                showline=True,
                ticks='',
                showticklabels=False
                ),
            zaxis = dict(title=ax_label.upper() + " 3",
                showbackground=False,
                gridcolor="rgb(218,218,218)",
                showgrid=showgrid,
                zeroline=False,
                showline=True,
                ticks='',
                showticklabels=False
                )
            ),
            height=600,
            legend=dict(orientation="h", y=.05),
            margin=dict(r=0, b=0, l=0, t=0)
            )
        fig = go.Figure(data=listGroups, layout=layout)
        if write_to is not None:
            fig.write_image(write_to+".png") # Save a regular file
            fig.write_html(write_to+".html") # Save an HTML file
        else:
            py.offline.iplot(fig, filename=title)


    # Plotting pie-charts using plotly
    def plotPieCharts(self, adata, variable=None, group_by=None, cols=2, palette=None, normalise=True, show_df=True, return_df=False, write_to=False):

        toExit = False

        # Checking if the group_by is not None
        if group_by is None:
            print(" * Warning! None for group_by is an invalide column of .obs, please provide one of the column in .obs")
            toExit = True

        # Checking if the variable is not None
        if variable is None:
            print(" * Warning! None for variable is an invalide column of .obs, please provide one of the column in .obs")
            toExit = True

        # Checking if the desidered group is in the .obs structure
        if group_by not in adata.obs.columns.tolist():
            print(" * Error! %s for group_by is an invalide column of .obs, please provide one of the column in .obs"%group_by)
            toExit = True

        # Checking if the desidered variable is in the .obs structure
        if variable not in adata.obs.columns.tolist():
            print(" * Error! %s for variable is an invalide column of .obs, please provide one of the column in .obs"%variable)
            toExit = True

        if toExit:
            return

        #keys      = adata.obs[group_by].unique()
        #variables = adata.obs[variable].unique()
        keys      = adata.obs[group_by].cat.categories.tolist()
        variables = adata.obs[variable].cat.categories.tolist()

        # Checking and palette to be used
        if palette is None:

            try:
                palette = adata.uns[variable+"_colors"]
            except:
                palette = sns.color_palette("tab20c", n_colors=len(variables))
        else:
            if isinstance(palette[0], tuple):
                palette = self.__rgb2hex(palette)

        if len(keys) > 1:
            rows = int(len(keys)/cols)

            if rows*cols < len(keys):
                rows += 1

        else:
            rows = 1
            cols = 1

        # Create subplots: use 'domain' type for Pie subplot
        fig = make_subplots(rows=rows,
                            cols=cols,
                            specs=[[{'type':'domain'}]*cols]*rows,
                            subplot_titles=keys)


        colorsDict = dict(zip(variables, palette))

        dfs = []
        for key in keys:
            reduced = adata[adata.obs[group_by]==key]
            dfs.append(pd.DataFrame(reduced.obs.groupby([variable]).size(), columns=[key]))

        df = pd.concat(dfs, axis=1)
        df.fillna(0, inplace=True)
        #del df.index.name

        display(df)
        if normalise:
            print("* Applying normalisation among groups in %s" %group_by)
            df = self.__normaliseGroupsDataFrame(df)

        if show_df:
            display(df)

        columns_df = df.columns.tolist()
        for r in range(0, rows):
            for c in range(0, cols):
                idx = r*cols + c
                if idx >= len(keys):
                    break

                if palette is None:
                    fig.add_trace(go.Pie(labels=df.index.tolist(),
                                        values=df[columns_df[idx]]),
                                r+1, c+1)
                else:
                    colors = []

                    for v in df.index.tolist():
                        colors.append(colorsDict[v])

                    fig.add_trace(go.Pie(labels=df.index.tolist(),
                                        values=df[columns_df[idx]],
                                        marker_colors=colors),
                                r+1, c+1)


        # Use `hole` to create a donut-like pie chart
        fig.update_traces(hole=.4, hoverinfo="label+value")

        if 500*cols > 1000:
            width = 1000
        else:
            width = 500*cols


        fig.update_layout(
            height=rows*600, width=width,
            title_text=variable,
            legend=dict(orientation="h", x=0., y=0.))

        if write_to is not None:
            fig.write_image(write_to+".png") # Save a regular file
            fig.write_html(write_to+".html") # Save an HTML file
        else:
            py.offline.iplot(fig, filename=title)

        if return_df:
            return df

    def plotBarPlotGroup(self, adata, variable=None, group_by=None, normalise=False, multi_barplot=False, palette=None, show_df=True, return_df=False):

        toExit = False

        # Checking if the group_by is not None
        if group_by is None:
            print(" * Warning! None for group_by is an invalide column of .obs, please provide one of the column in .obs")
            toExit = True

        # Checking if the variable is not None
        if variable is None:
            print(" * Warning! None for variable is an invalide column of .obs, please provide one of the column in .obs")
            toExit = True

        # Checking if the desidered group is in the .obs structure
        if group_by not in adata.obs.columns.tolist():
            print(" * Error! %s for group_by is an invalide column of .obs, please provide one of the column in .obs"%group_by)
            toExit = True

        # Checking if the desidered variable is in the .obs structure
        if variable not in adata.obs.columns.tolist():
            print(" * Error! %s for variable is an invalide column of .obs, please provide one of the column in .obs"%variable)
            toExit = True

        if toExit:
            return

        keys      = adata.obs[group_by].cat.categories.tolist()
        variables = adata.obs[variable].cat.categories.tolist()

        # Checking and palette to be used
        if palette is None:
            try:
                palette = adata.uns[variable+"_colors"]
            except:
                palette = sns.color_palette("tab20c", n_colors=len(keys))
        else:
            if isinstance(palette[0], tuple):
                palette = self.__rgb2hex(palette)

        dfs = []
        for key in keys:
            reduced = adata[adata.obs[group_by]==key]
            dfs.append(pd.DataFrame(reduced.obs.groupby([variable]).size(), columns=[key]))

        df = pd.concat(dfs, axis=1)
        df.fillna(0, inplace=True)
        del df.index.name

        if normalise:
            print("* Applying normalisation among groups in %s" %group_by)
            df = self.__normaliseGroupsDataFrame(df)

        if multi_barplot:
            for col in df.columns:
                df_norm_tmp = pd.DataFrame(df_norm[col])
                self.__barplotGroupsDataFrame(df_norm_tmp, palette=palette, title="Normalised barplot",show_df=show_df)
        else:
            self.__barplotGroupsDataFrame(df, palette=palette, title="Normalised barplot",show_df=show_df)

        if return_df:
            return df


    def plotBarPlotFeatures(self, adata, group_by1=None, group_by2=None, features=[], show_df=False, height=20, width=60, palette=None, pdf=None):

        toExit = False

        # Checking if the group_by1 is not None
        if group_by1 is None:
            print(" * Error! None for 'group_by1' is an invalide column of .obs, please provide one of the column in .obs")
            toExit = True
        # Checking if the desidered group_by1 is in the .obs structure
        if group_by1 not in adata.obs.columns.tolist():
            print(" * Error! %s for 'group_by1' is an invalide column of .obs, please provide one of the column in .obs"%group_by1)
            toExit = True
        # Checking if the desidered group_by2 is in the .obs structure
        if group_by2 not in adata.obs.columns.tolist() and group_by2 is not None:
            print(" * Error! %s for 'group_by2' is an invalide column of .obs, please provide one of the column in .obs"%group_by2)
            toExit = True
        # Checking if the desidered features are numeric and in the .obs structure
        for f in features:
            if f not in adata.obs.columns.tolist():
                print(" * Error! %s for 'feature' is an invalide column of .obs, please provide one of the column in .obs"%f)
                toExit = True
            elif is_numeric_dtype(adata.obs[f]) == False:
                print(" * Error! %s for 'feature' is not a numeric feature"%f)
                toExit = True
        if toExit:
            return

        # Checking and palette to be used
        if palette is None:

            try:
                palette = adata.uns[group_by1+"_colors"]
            except:
                palette = sns.color_palette("tab20c")
        else:
            if isinstance(palette[0], tuple):
                palette = self.__rgb2hex(palette)


        if group_by2 != None and group_by1 != None:

            dfs_median_sample = []

            for cluster in adata.obs[group_by2].cat.categories:
                reduced = adata[adata.obs[group_by2]==cluster]

                dict_samples = {}
                for sample in reduced.obs[group_by1].cat.categories:
                    reduced2 = reduced[reduced.obs[group_by1]==sample]

                    dict_samples[sample] = [reduced2.n_obs,
                            sample,
                            cluster] + [np.median(reduced2.obs[f]) for f in features]

                    df_sample = pd.DataFrame.from_dict(dict_samples, orient="index")
                df_sample.rename({0:"Number of cell",
                    1:group_by1,
                    2:group_by2}, axis="columns", inplace=True)

                for idx, f in enumerate(features):
                    df_sample.rename({idx+3:f}, axis="columns", inplace=True)

                dfs_median_sample.append(df_sample)

            dfs = pd.concat(dfs_median_sample)
            dfs[group_by1] = dfs.index.astype("category")
            dfs.reset_index(drop=True, inplace=True)

            dfs[group_by1].cat.reorder_categories(adata.obs[group_by1].cat.categories, inplace=True)

            if show_df:
                display(dfs)

            f, axs = plt.subplots(1,len(features)+1,figsize=(width,height))
            sns.set(font_scale=2.6)
            sns.set_style("white")

            sns.barplot(y=group_by2, x="Number of cell", hue=group_by1, data=dfs, ax=axs[0], orient="h", palette=palette)
            for idx, f in enumerate(features):
                sns.barplot(y=group_by2, x=f, hue=group_by1, data=dfs, ax=axs[idx+1], orient="h", palette=palette)
                axs[idx+1].axvline(x=np.median(dfs[f]), linewidth=4, color='b')
                axs[idx+1].axvline(x=np.median(dfs[f])-np.std(dfs[f]),
                        linewidth=4, color='r', linestyle="--")
                axs[idx+1].axvline(x=np.median(dfs[f])+np.std(dfs[f]),
                        linewidth=4, color='r', linestyle="--")

                for idx, f in enumerate(features):
                    axs[idx].get_legend().remove()

            if len(features) > 0:
                axs[-1].legend(loc='best')
                axs[-1].legend(frameon=False)
            else:
                axs.legend(loc='best')
                axs.legend(frameon=False)

            #sns.despine(offset=10, trim=False)
            plt.tight_layout()
            if pdf is not None:
                pdf.savefig()
            else:
                plt.show(block=False)
            plt.close("all")

        else:

            dfs_median_sample = []

            dict_samples = {}
            for sample in adata.obs[group_by1].cat.categories:
                reduced = adata[adata.obs[group_by1]==sample]

                dict_samples[sample] = [reduced.n_obs,
                        sample] + [np.median(reduced.obs[f]) for f in features]

                df_sample = pd.DataFrame.from_dict(dict_samples, orient="index")
                df_sample.rename({0:"Number of cell", 1:group_by1}, axis="columns", inplace=True)

            for idx, f in enumerate(features):
                df_sample.rename({idx+2:f}, axis="columns", inplace=True)

            df_sample[group_by1] = df_sample.index.astype("category")
            df_sample.reset_index(drop=True, inplace=True)

            df_sample[group_by1].cat.reorder_categories(adata.obs[group_by1].cat.categories, inplace=True)

            if show_df:
                display(df_sample)

            f, axs = plt.subplots(1,len(features)+1,figsize=(width,height))
            sns.set(font_scale=2.3)
            sns.set_style("white")

            sns.barplot(y=group_by1, x="Number of cell", data=df_sample, ax=axs[0], orient="h", palette=palette)

            for idx, f in enumerate(features):
                sns.barplot(y=group_by1, x=f, data=df_sample, ax=axs[idx+1], orient="h", palette=palette)
                axs[idx+1].axvline(x=np.median(df_sample[f]), linewidth=4, color='b')
                axs[idx+1].axvline(x=np.median(df_sample[f])-np.std(df_sample[f]),
                        linewidth=4, color='r', linestyle="--")
                axs[idx+1].axvline(x=np.median(df_sample[f])+np.std(df_sample[f]),
                        linewidth=4, color='r', linestyle="--")

                #sns.despine(offset=10, trim=False)
            plt.tight_layout()
            if pdf is not None:
                pdf.savefig()
            else:
                plt.show(block=False)
            plt.close("all")


    def plotBoxPlotGroup(self, adata, group_by1=None, group_by2=None, height=10, width=30, palette=None):

        toExit = False

        # Checking if the group_by1 is not None
        if group_by1 is None:
            print(" * Error! None for 'group_by1' is an invalide column of .obs, please provide one of the column in .obs")
            toExit = True
        # Checking if the desidered group_by1 is in the .obs structure
        if group_by1 not in adata.obs.columns.tolist():
            print(" * Error! %s for 'group_by1' is an invalide column of .obs, please provide one of the column in .obs"%group_by1)
            toExit = True

        # Checking if the group_by1 is not None
        if group_by2 is None:
            print(" * Error! None for 'group_by2' is an invalide column of .obs, please provide one of the column in .obs")
            toExit = True
        # Checking if the desidered group_by1 is in the .obs structure
        if group_by2 not in adata.obs.columns.tolist():
            print(" * Error! %s for 'group_by2' is an invalide column of .obs, please provide one of the column in .obs"%group_by2)
            toExit = True

        # Checking if the group_by2 is numeric
        if is_numeric_dtype(adata.obs[group_by2]) == False:
            print(" * Error! %s for 'group_by2' is not a numeric column in .obs"%group_by2)
            toExit = True
        if toExit:
            return

        # Checking and palette to be used
        if palette is None:

            try:
                palette = adata.uns[group_by1+"_colors"]
            except:
                palette = sns.color_palette("tab20c")
        else:
            if isinstance(palette[0], tuple):
                palette = self.__rgb2hex(palette)

        f, axs = plt.subplots(1,1,figsize=(width,height))
        sns.despine(offset=10, trim=False)
        sns.set(font_scale=1.5)

        sns.set_style("white")
        sns.boxplot(x       = group_by1,
                y       = group_by2,
                data    = adata.obs,
                palette = palette,
                ax      = axs)

        plt.setp(axs.get_xticklabels(), rotation=90)
        plt.tight_layout()



    def PCA_ElbowPlot(self, adata, n_pcs=50, log=True, width=16, height=8, pdf=None):

        sns.set(rc={'figure.figsize':(width,height)})
        sns.set(font_scale=1.5)
        sns.set_style("white")

        sc.pl.pca_variance_ratio(adata, n_pcs=n_pcs, log=log, show=False)

        #sns.despine(offset=10, trim=False)
        if pdf is not None:
            pdf.savefig()
        else:
            plt.show(block=False)
        plt.close()

    # Plotting the HVGs
    def plotHVGs(self, adata, width=8, height=8, pdf=None):

        sns.set(rc={'figure.figsize':(width,height)})
        sns.set(font_scale=1.5)
        sns.set_style("white")

        sc.pl.highly_variable_genes(adata, show=False)

        #sns.despine(offset=10, trim=False)
        if pdf is not None:
            pdf.savefig()
        else:
            plt.show(block=False)
        plt.close()

    # Plotting a scatter plot with two series
    def plotScatterDatasets(self, dataset1=None, dataset2=None, axis=None, label1="Real", label2="Synthetic", title="", width=8, height=8):

        if axis is None:
            f, axis = plt.subplots(1,1,figsize=(width,height))
            sns.set(font_scale=1.5)
            sns.set_style("white")

        if dataset1 is not None:
            axis.scatter(dataset1[:,0], dataset1[:,1], s=36, label=label1)

        if dataset2 is not None:
            axis.scatter(dataset2[:,0], dataset2[:,1], s=12, label=label2)

        axis.legend()
        axis.set_title(title)

        if axis is None:
            sns.despine(offset=10, trim=False)
            plt.tight_layout()
            plt.show()


    # Plotting ridge plots grouping the data
    def ridge_plot(self, adata, genes=[], group_by=None, use_raw=True, height=0.75, palette=None, pdf=None, title=""):

        toExit = False

        # Checking if the group_by is not None
        if group_by is None:
            print(" * Warning! None for group_by is an invalide column of .obs, please provide one of the column in .obs")
            toExit = True

        # Checking if the desidered group is in the .obs structure
        if group_by not in adata.obs.columns.tolist():
            print(" * Error! %s for group_by is an invalide column of .obs, please provide one of the column in .obs"%group_by)
            toExit = True

        if toExit:
            return

        if palette is None:

            try:
                palette = adata.uns[group_by+"_colors"]
            except:
                palette = sns.color_palette("tab20c")
        else:
            if isinstance(palette[0], tuple):
                palette = self.__rgb2hex(palette)

        if not isinstance(genes, list):
            genes = [genes]

        for gene in genes:

            if use_raw:
                genes_exp = adata.raw.X[:, adata.var.index==gene]
            else:
                genes_exp = adata.X[:, adata.var.index==gene]
            try:
                genes_exp = np.array(genes_exp.todense())
            except:
                pass

            dataframe = pd.DataFrame(data=genes_exp, index=adata.obs.index, columns=[gene])
            dataframe[group_by] = adata.obs[group_by]

            dataframe2 = dataframe.groupby([group_by]).sum()
            l = dataframe2[dataframe2[dataframe2.columns[0]]!=0].index.tolist()
            dataframe = dataframe[dataframe[group_by].isin(l)]
            dataframe[group_by] = dataframe[group_by].cat.remove_unused_categories()
            del dataframe2

            self.__plotSingleRidgePlot(dataframe, palette, height=height, pdf=pdf, title=title)
