#!/usr/bin/bash

# # # sbatch --account=smontgom --partition=batch --time=7-1:00:00 --mem=256G --nodes=1 --ntasks=1 --cpus-per-task=4 cpdb_t21_liver.sh

module load R/4.1.2

Rscript /oak/stanford/groups/smontgom/amarder/t21-proj/scripts/create_cpdb_input.R

/oak/stanford/groups/smontgom/amarder/t21-proj/scripts/run_cellphonedb.sh DownSyndrome Liver false


