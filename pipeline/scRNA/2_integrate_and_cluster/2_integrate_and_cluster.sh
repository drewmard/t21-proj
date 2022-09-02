#!/bin/bash

# # # sbatch --account=smontgom --partition=batch --time=3-1:00:00 --mem=224G --nodes=1 --ntasks=1 /oak/stanford/groups/smontgom/amarder/t21-proj/pipeline/scRNA/2_integrate_and_cluster/2_integrate_and_cluster.sh

module load R/4.1.2
Rscript /oak/stanford/groups/smontgom/amarder/t21-proj/pipeline/scRNA/2_integrate_and_cluster/2_integrate_and_cluster.R