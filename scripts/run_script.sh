#!/usr/bin/bash

# example run:
# sbatch --account=smontgom --partition=batch --time=24:00:00 --mem=256G --nodes=1 --ntasks=1 --cpus-per-task=1 run_script.sh subset_data.py
# sbatch --account=smontgom --partition=batch --time=24:00:00 --mem=356G --nodes=1 --ntasks=1 --cpus-per-task=1 run_script.sh plots_v4.py

echo "Start"
source /oak/stanford/groups/smontgom/amarder/bin/conda_init.sh 
conda activate minimal_env

echo "Beginning script $1"
# run script present as the first argument
python -u $1
