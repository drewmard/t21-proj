#!/bin/bash

# # sbatch --account=smontgom --partition=batch --time=2-1:00:00 --mem=256G --nodes=1 --ntasks=1 /oak/stanford/groups/smontgom/amarder/t21-proj/scripts/multiome/run_chromvar.sh

module load R/4.1.2
path_to_script=/oak/stanford/groups/smontgom/amarder/t21-proj/scripts/multiome/chromvar.R
Rscript $path_to_script
