#!/bin/bash

# # # sbatch --account=smontgom --partition=batch --time=3-1:00:00 --mem=128G --nodes=1 --ntasks=1 /oak/stanford/groups/smontgom/amarder/t21-proj/scripts/RNA_cluster.peaks.sh DS_Multiome_h

module load R/4.1.2
Rscript /oak/stanford/groups/smontgom/amarder/t21-proj/scripts/RNA_cluster.peaks.R $1