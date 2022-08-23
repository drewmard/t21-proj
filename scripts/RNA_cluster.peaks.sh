#!/bin/bash

# # # sbatch --account=smontgom --partition=batch --time=3-1:00:00 --mem=128G --nodes=1 --ntasks=1 /oak/stanford/groups/smontgom/amarder/t21-proj/scripts/RNA_cluster.peaks.sh DS_Multiome_h 0 3

module load R/4.1.2
dir="/oak/stanford/groups/smontgom/amarder/neuro-variants"
dataset=$1
start=$2
end=$3
Rscript /oak/stanford/groups/smontgom/amarder/t21-proj/scripts/RNA_cluster.peaks.R $dir $dataset $start $end