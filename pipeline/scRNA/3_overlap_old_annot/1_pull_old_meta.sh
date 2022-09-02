#!/bin/bash

# # # sbatch --account=smontgom --partition=batch --time=3-1:00:00 --mem=256G --nodes=1 --ntasks=1 /oak/stanford/groups/smontgom/amarder/t21-proj/pipeline/scRNA/3_overlap_old_annot/1_pull_old_meta.sh

module load R/4.1.2
Rscript /oak/stanford/groups/smontgom/amarder/t21-proj/pipeline/scRNA/3_overlap_old_annot/1_pull_old_meta.R