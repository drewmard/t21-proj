#!/bin/bash

# # # sbatch --account=smontgom --partition=batch --time=3-1:00:00 --mem=128G --nodes=1 --ntasks=1 /oak/stanford/groups/smontgom/amarder/t21-proj/scripts/RNA_cluster.peaks.sh 0 3 DS_Multiome_h

module load R/4.1.2
dir=/oak/stanford/groups/smontgom/amarder/neuro-variants
start=$1
end=$2
dataset=$3
echo "Rscript /oak/stanford/groups/smontgom/amarder/t21-proj/scripts/RNA_cluster.peaks.R $dir $start $end $dataset"
Rscript /oak/stanford/groups/smontgom/amarder/t21-proj/scripts/RNA_cluster.peaks.R $dir $start $end $dataset