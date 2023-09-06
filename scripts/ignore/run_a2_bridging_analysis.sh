#!/bin/sh

# # # sbatch --account=smontgom --partition=batch --time=48:00:00 --mem=256G --nodes=1 --ntasks=1 --cpus-per-task=1 run_a2_bridging_analysis.sh
Rscript /oak/stanford/groups/smontgom/amarder/t21-proj/scripts/a2_bridging_analysis.R
