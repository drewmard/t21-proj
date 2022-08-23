#!/bin/bash

# # # sbatch --account=smontgom --partition=batch --time=7-1:00:00 --mem=128G --nodes=1 --ntasks=1 /oak/stanford/groups/smontgom/amarder/t21-proj/RNA_cluster.peaks.sh $dataset

dataset=DS_Multiome_h

Rscript RNA_cluster.peaks.R