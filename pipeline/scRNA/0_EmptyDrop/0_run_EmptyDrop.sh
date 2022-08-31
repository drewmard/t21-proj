#!/bin/bash

# # # sbatch --account=smontgom --partition=batch --time=3-1:00:00 --mem=128G --nodes=1 --ntasks=1 /oak/stanford/groups/smontgom/amarder/t21-proj/scripts/RNA_cluster.peaks.sh 2 2 DS_Multiome_h

i=$1

module load R/4.1.2
echo "Rscript /Users/andrewmarderstein/Documents/Research/t21-proj/pipeline/scRNA/0_EmptyDrop/0_run_EmptyDrop.R $i"
Rscript /Users/andrewmarderstein/Documents/Research/t21-proj/pipeline/scRNA/0_EmptyDrop/0_run_EmptyDrop.R $i