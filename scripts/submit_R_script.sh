#!/usr/bin/bash

# example run:
# sbatch --account=smontgom --partition=batch --time=24:00:00 --mem=356G --nodes=1 --ntasks=1 --cpus-per-task=1 submit_R_script.sh anndata_to_seurat.R
# sbatch --account=smontgom --partition=batch --time=24:00:00 --mem=356G --nodes=1 --ntasks=1 --cpus-per-task=1 submit_R_script.sh pseudobulk_de.R

module load R/4.0

echo "Start"

echo "Beginning script $1"
# run script present as the first argument
Rscript $1
