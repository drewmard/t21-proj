#!/bin/bash

for i in {3..3}; do

# sbatch --account=default --partition=interactive --time=1-1:00:00 --mem=4G --nodes=1 --ntasks=1 /oak/stanford/groups/smontgom/amarder/t21-proj/pipeline/scRNA/0_EmptyDrop/0_run_EmptyDrop.sh $i
sbatch --account=smontgom --partition=batch --time=1-0:00:00 --mem=4G --nodes=1 --ntasks=1 /oak/stanford/groups/smontgom/amarder/t21-proj/pipeline/scRNA/0_EmptyDrop/0_run_EmptyDrop.sh $i

done
