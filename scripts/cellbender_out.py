# conda activate ddqc

# import scanpy as sc
import pegasus as pg
import os
import anndata
import pandas as pd

# Specify the path to your file
file_path = "~/tmp/ts21_cellbender_download.txt"

# Read the file into a Pandas DataFrame
metadata = pd.read_csv(file_path, delimiter='\t')

# Specify the base directory
base_directory = "/oak/stanford/groups/smontgom/amarder/data/t21/PostAmbient_Nov2023"

# disease="Down Syndrome"
# organ = "Liver"
disease="Healthy"
sampletype = "Femur"

for disease in ['Down Syndrome','Healthy']:
	for sampletype in ["Liver","Femur"]:

		# Filter rows based on conditions
		filtered_df = metadata[(metadata['environment'] == disease) & (metadata['organ'] == sampletype)]

		# Extract CellrangerFile values from the filtered DataFrame
		cellranger_files = filtered_df['CellrangerFile'].tolist()
		cellbender_h5 = [base_directory + "/" + path + "/corrected_filtered.h5" for path in cellranger_files]

		# Create a dictionary
		sample_dict = {'Sample': cellranger_files, 'Location': cellbender_h5}

		# Read file & aggregate into a single pegasus object
		dataset = pg.aggregate_matrices(sample_dict)

		# Convert to anndata (scanpy) and save!
		output_file = "/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/cellbender."+disease+"."+sampletype+".h5ad"
		dataset.to_anndata().write_h5ad(output_file)
