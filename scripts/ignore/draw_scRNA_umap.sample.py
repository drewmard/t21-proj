print("Beginning script...")
import scanpy as sc
import sys
import seaborn as sns
from matplotlib import pylab as plt
import os
import numpy as np
import pandas as pd
from sklearn import preprocessing
from palettable.colorbrewer.qualitative import Set3_8,Set3_9

print("Module loading completed.")

headdir="/oak/stanford/groups/smontgom/amarder/t21-proj"
disease_status="Healthy"
sampletype="Liver"
suffix="" # ".subset"
make_new_3D_UMAP=False

# initialize:
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

suffixDirec=""
if suffix=="subset":
    suffixDirec="subset"
else:
    suffixDirec="full"
os.system("mkdir -p "+headdir + "/out/" + suffixDirec)


for disease_status in ["DownSyndrome","Healthy"]:
      for sampletype in ["Femur","Liver"]:

            # if disease_status=="Healthy" and sampletype=="Liver":
            #       colors_to_use = Set3_8.mpl_colors
            # else:
            #       colors_to_use = Set3_9.mpl_colors

            if suffix=="subset":
                  fout="10X_" + disease_status + "_" + sampletype + ".umap.subset.cells_removed.h5ad"
                  foutpath=headdir + "/out/data/" + fout
            else:
                  fout="10X_" + disease_status + "_" + sampletype + ".h5ad"
                  foutpath="/oak/stanford/groups/smontgom/amarder/data/t21/ScanpyObjects" + "/" + fout

            print("\n * Reading in data..." + foutpath)
            adata=sc.read_h5ad(foutpath)

            subDirec="data"
            os.system("mkdir -p "+headdir + "/out/" + suffixDirec +"/" + subDirec)
            fout="10X_" + disease_status + "_" + sampletype + ".umap3d.cells_removed.h5ad"
            foutpath=headdir + "/out/" + suffixDirec + "/" + subDirec + "/" + fout
            print("\n * Saving 3d data to file..." + foutpath)
            adata_3d=sc.read_h5ad(foutpath)
            print("Done.")
            fout="10X_" + disease_status + "_" + sampletype + ".umap2d.cells_removed.h5ad"
            foutpath=headdir + "/out/" + suffixDirec + "/" + subDirec + "/" + fout
            print("\n * Saving 2d data to file..." + foutpath)
            adata=sc.read_h5ad(foutpath)
            print("Done.")
            print("\n * Script completed.")

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

            # f, axs = plt.subplots(1,1,figsize=(20,16))
            # sns.set(font_scale=1.5)
            # sns.set_style("white")

            subDirec="umap_cluster"
            os.system("mkdir -p "+headdir + "/out/" + suffixDirec +"/" + subDirec)

            fplotout=headdir + "/out/" + suffixDirec + "/" + subDirec + "/" + "10X_"+disease_status+"_"+sampletype+".umap"+".sample.pdf"
            # os.remove(fplotout)
            print("\n * Plotting & saving UMAP (numerical labels)... " + fplotout)
            f, axs = plt.subplots(1,1,figsize=(26,26))
            sns.set(font_scale=2)
            sns.set_style("white")
            new_plot=sc.pl.umap(adata, color="sample", size=150, palette=myColors, components='1,2', ax=axs, show=False, use_raw=False, title=disease_status + ' ' + sampletype,legend_loc="on data")
            plt.tight_layout()
            plt.savefig(fplotout)
            plt.show()
            plt.close()
            print("\n * Plot saved.")

            ##############################

            for components_to_use in ["1,2","1,3","2,3"]:

                  Name_components_to_use = components_to_use.replace(",",".")
                  fplotout=headdir + "/out/" + suffixDirec + "/" + subDirec + "/" + "10X_"+disease_status+"_"+sampletype+".umap."+Name_components_to_use+".sample.pdf"
                  print("\n * Plotting & saving UMAP..." + fplotout)
                  f, axs = plt.subplots(1,1,figsize=(26,26))
                  sns.set(font_scale=2)
                  sns.set_style("white")
                  sc.pl.umap(adata_3d, color="sample", size=150, palette=myColors, components=[components_to_use], ax=axs, show=False, use_raw=False, title=disease_status + ' ' + sampletype,legend_loc="on data")
                  plt.tight_layout()
                  plt.savefig(fplotout)
                  plt.show()
                  plt.close()
                  print("\n * Plot saved.")

                  fplotout=headdir + "/out/" + suffixDirec + "/" + subDirec + "/" + "10X_"+disease_status+"_"+sampletype+".umap."+Name_components_to_use+".sample.no_legend.pdf"
                  print("\n * Plotting & saving UMAP (no legend)..." + fplotout)
                  f, axs = plt.subplots(1,1,figsize=(26,26))
                  sns.set(font_scale=2)
                  sns.set_style("white")
                  new_plot=sc.pl.umap(adata_3d, color="sample", size=150, palette=myColors, components=[components_to_use], ax=axs, show=False, use_raw=False, title=disease_status + ' ' + sampletype,legend_loc="None")
                  plt.tight_layout()
                  plt.savefig(fplotout)
                  plt.show()
                  plt.close()
                  print("\n * Plot saved.")

print("\n * Script complete...")

# plt.tight_layout()
# pdf.savefig()
# plt.close()

