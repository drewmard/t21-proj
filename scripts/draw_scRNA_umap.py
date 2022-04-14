import scanpy as sc
import sys
import seaborn as sns
from matplotlib import pylab as plt
import os
import matplotlib.backends.backend_pdf as mpdf


sys.path.append("/oak/stanford/groups/smontgom/amarder/t21_download/Functions/")
direc="/oak/stanford/groups/smontgom/amarder/data/t21/ScanpyObjects"
headdir="/oak/stanford/groups/smontgom/amarder/t21"
# Get the global settings
#from global_settings import global_settings
#global_settings()

# Get the bespoke analysis functions
from scRNA_functions import scRNA_functions
fc = scRNA_functions()


print("\n * Reading in data...")
for 10X_DownSyndrome_Femur.h5ad 10X_DownSyndrome_Liver.h5ad  10X_Healthy_Femur.h5ad  10X_Healthy_Liver.h5ad
for disease_status in ["DownSyndrome","Healthy"]:
      for sampletype in ["Femur","Liver"]:
            f="10X_" + disease_status "_" + sampletype + ".h5ad"
            fpath=direc + "/" + f

            print("\n * Reading in data...")
            adata=sc.read_h5ad(fpath)
            print("\n * Computing UMAP (2 components)...")
            sc.tl.umap(adata, random_state=10, n_components=2, init_pos='random')

            fout="10X_" + disease_status "_" + sampletype + ".umap.h5ad"
            foutpath=direc + "/" + f
            print("\n * Saving data to file..." + foutpath)
            adata.write(foutpath)

            # print("\n * Reading in data...")
            # fpath="/oak/stanford/groups/smontgom/amarder/t21/out/data/10X_Healthy_Liver.umap.h5ad"
            # adata=sc.read_h5ad(fpath)

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

            os.chdir(headdir + "/out/figures")
            pdf = mpdf.PdfPages("10X_"+disease_status+"_"+sampletype+".umap.pdf")
            print("\n * Writing UMAP...")
            # sc.pl.umap(adata,color="leiden",palette=myColors,save="10X_Healthy_Liver.umap.png")
            fc.plotUMAP(adata, variable="leiden", palette=myColors, pdf=pdf,width=30,height=16)
            print("\n * Script complete...")
            pdf.close()

# plt.tight_layout()
# pdf.savefig()
# plt.close()

