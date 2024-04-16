# /home/amarder/anaconda3/envs/minimal_env/bin/python

import warnings
warnings.simplefilter("ignore")
import sys
import os

import scanpy as sc
import numpy as np
import pandas as pd
import seaborn as sns

import matplotlib
matplotlib.use("Agg") # Want a non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as mpdf

from math import log
from scipy import sparse
from matplotlib import cm
from matplotlib import colors
from scipy.sparse import issparse
from scipy.spatial import distance

# Global colour settings
myColors = ['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe',
            '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080',
            '#307D7E', '#000000', '#DDEFFF', '#000035', '#7B4F4B', '#A1C299', '#300018', '#C2FF99', '#0AA6D8', '#013349',
            '#00846F', '#8CD0FF', '#3B9700', '#04F757', '#C8A1A1', '#1E6E00', '#DFFB71', '#868E7E', '#513A01', '#CCAA35',
            '#800080', '#DAA520', '#1E90FF', '#3CB371', '#9370DB', '#8FBC8F', '#00FF7F', '#0000CD', '#556B2F', '#FF00FF',
            '#CD853F', '#6B8E23', '#008000', '#6495ED', '#00FF00', '#DC143C', '#FFFF00', '#00FFFF', '#FF4500', '#4169E1',
            '#48D1CC', '#191970', '#9ACD32', '#FFA500', '#00FA9A', '#2E8B57', '#40E0D0', '#D2691E', '#66CDAA', '#FFEFD5',
            '#20B2AA', '#FF0000', '#EEE8AA', '#BDB76B', '#E9967A', '#AFEEEE', '#000080', '#FF8C00', '#B22222', '#5F9EA0',
            '#ADFF2F', '#FFE4B5', '#7B68EE', '#7FFFD4', '#0000FF', '#BA55D3', '#90EE90', '#FFDAB9', '#6A5ACD', '#8B0000',
            '#8A2BE2', '#CD5C5C', '#F08080', '#228B22', '#FFD700', '#006400', '#98FB98', '#00CED1', '#00008B', '#9400D3',
            '#9932CC', '#4B0082', '#F0E68C', '#483D8B', '#008B8B', '#8B008B', '#4682B4']
myColorTissues = ['#ffe119', '#3cb44b', '#e6194b']

# Define a nice colour map for gene expression
colors2 = plt.cm.Reds(np.linspace(0, 1, 128))
colors3 = plt.cm.Greys_r(np.linspace(0.7,0.8,20))
colorsComb = np.vstack([colors3, colors2])
mymap = colors.LinearSegmentedColormap.from_list('my_colormap', colorsComb)

sns.set(font_scale=1.5)
sns.set_style("white")

idtype="h"
idtype="ds"
outDir = "/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h_v2"
f = "/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h_v2/Multiome.ATAC_only2.h.h5ad"
annotations="subclust_v6"
df = sc.read(f)
root_node="HSCs" # for DS
threshold=0.35

sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(3, 3), facecolor='white')  # low dpi (dots per inch) yields small inline figures

# # Trajectory analysis with PAGA graphs and a force-directed graph embedding
adata_paga  = df.copy()
adata_paga = adata_paga[~adata_paga.obs[annotations].isin(['Stroma', 'No markers'])]

sc.pp.neighbors(adata_paga,
            n_pcs=adata_paga.obsm["X_harmony"].shape[1],
            use_rep="X_harmony",
            knn=True,
            random_state=42,
            method='umap',
            metric='euclidean')

sc.tl.paga(adata_paga, groups=annotations)

os.chdir(outDir)
f_out=".ATAC.PAGA_raw_graph_harmony.png"
sc.pl.paga(adata_paga, threshold=threshold, node_size_scale=0.8, edge_width_scale=0.3, fontsize=5, save=f_out)

# Recompute the embedding using PAGA initialisation
kwargs = {"maxiter": 1000}
sc.tl.draw_graph(adata_paga,
             init_pos  = 'paga',
             root      = np.flatnonzero(adata_paga.obs[annotations]==root_node)[0],
             adjacency = adata_paga.uns['neighbors']['connectivities'],
             **kwargs
             )

###########


# Now plot the combined PAGA graph + embedding
f, axs = plt.subplots(1,1,figsize=(25,30))
f_out2=".ATAC.PAGA_trajectory_embedding_harmony.png"
sc.pl.draw_graph(adata_paga,
             color=annotations,
             palette=myColors,
             legend_loc='right margin',
             ax=axs,
             show=False,
             size=300,
             legend_fontsize='large',
             save=f_out2
        )

plt.tight_layout()
plt.close("all")

# Compute the diffusion map
adata_paga.uns['iroot'] = np.flatnonzero(adata_paga.obs[annotations]==root_node)[0]
sc.tl.diffmap(adata_paga, n_comps=50)

# Investigate the effect of applying HARMONY here
# sc.pp.neighbors(adata_paga,
#             n_pcs=adata_paga.obsm["X_diffmap"].shape[1],
#             use_rep="X_diffmap",
#             knn=True,
#             random_state=42,
#             method='umap',
#             metric='euclidean')
sc.tl.dpt(adata_paga)

f, axs = plt.subplots(1,1,figsize=(25,30))
f_out3=".ATAC.PAGA_diffusion_pseudotime_harmony.png".format(idtype)
# sc.pl.draw_graph(adata_paga,
#      color="dpt_pseudotime",
#      legend_loc="right margin",
#      ax=axs,
#      show=False,
#      size=300,
#      legend_fontsize='large',
#      save=f_out3
# )
sc.pl.draw_graph(adata_paga,
     color="dpt_pseudotime",
     legend_loc="right margin",
     ax=axs,
     show=False,
     size=300,
     legend_fontsize='large',
     save=f_out3
)
plt.tight_layout()
plt.close("all")

#######

# Plot the average pseudotimes
times = adata_paga.obs.groupby(by=annotations).mean()["dpt_pseudotime"].sort_values()
#uncs = adata_paga.obs.groupby(by=annotations).std()["dpt_pseudotime"].reindex(times.index) # Standard deviation
uncs = 2*1.96*adata_paga.obs.groupby(by=annotations).sem()["dpt_pseudotime"].reindex(times.index) # 95 % confidence from the standard error in the mean
f, axs = plt.subplots(1,1,figsize=(15,10))
f_out4="{0}/figures/ATAC.PAGA_Mean_Pseudotimes_harmony.pdf".format(outDir)

pdf = mpdf.PdfPages(f_out4)
times.plot.barh(yerr=uncs)
plt.ylabel("Cell type")
plt.xlabel("Mean pseudotime")
plt.title("PAGA")
plt.tight_layout()
pdf.savefig()
plt.close("all")
pdf.close()

# PAGA connection plot following comparison
f_out5="{0}/figures/ATAC.PAGA_connections_harmony.pdf".format(outDir)
# pdf = mpdf.PdfPages(".dfcombined.{0}.post_cluster_peaks.v2.res_0.6.chromvar_PAGA_connections_harmony.pdf")
pdf = mpdf.PdfPages(f_out5)

f, axs = plt.subplots(1,1,figsize=(16,16))
sc.pl.draw_graph(adata_paga, size=100, legend_fontsize="large", legend_loc="right margin", ax=axs, show=False)
sc.pl.paga(adata_paga,
   pos              = adata_paga.uns['paga']['pos'],
   show             = False,
   node_size_power  = 0.5,
   ax               = axs,
   threshold        = threshold,
   node_size_scale  = 0.8,
   edge_width_scale = 0.3,
   fontsize         = 5,
   text_kwds        = {'alpha':1})

plt.tight_layout()
pdf.savefig()
plt.close("all")

pdf.close()

f_out=".ATAC.PAGA_raw_graph_harmony2.png"
sc.pl.paga(adata_paga, threshold=threshold, node_size_scale=0.8, edge_width_scale=0.3, fontsize=5, save=f_out)


f_df_out = "/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h_v2/Multiome.ATAC_only2.h.paga.h5ad"

adata_paga.write(f_df_out, compression="gzip")





