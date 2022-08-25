#!/bin/bash

# # # sbatch --account=smontgom --partition=batch --time=7-1:00:00 --mem=256G --nodes=1 --ntasks=1 /oak/stanford/groups/smontgom/amarder/t21-proj/scripts/FigR/figr.sh DS_Multiome_h

module load R/4.1.2
dir=/oak/stanford/groups/smontgom/amarder/neuro-variants
dataset=$1
echo $dataset
echo "Rscript /oak/stanford/groups/smontgom/amarder/t21-proj/scripts/FigR/figr.R"
Rscript /oak/stanford/groups/smontgom/amarder/t21-proj/FigR/figr.R 