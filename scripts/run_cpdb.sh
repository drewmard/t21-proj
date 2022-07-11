#!/usr/bin/bash

# example run:
# srun --account=smontgom --partition=batch --time=1:00:00 --mem=356G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash
# srun --account=default --partition=interactive --time=24:00:00 --mem=128G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash
# sbatch --account=default --partition=interactive --time=24:00:00 --mem=128G --nodes=1 --ntasks=1 --cpus-per-task=1 run_script.sh diff_exp.py
# sbatch --account=smontgom --partition=batch --time=24:00:00 --mem=256G --nodes=1 --ntasks=1 --cpus-per-task=1 run_cellphonedb.sh DownSyndrome Liver true


echo "Start"
source /oak/stanford/groups/smontgom/amarder/bin/conda_init.sh 
conda activate cpdb

echo "Beginning script $1"
# run script present as the first argument
$1 $2 $3
