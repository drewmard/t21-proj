#!/usr/bin/bash

# sbatch --account=smontgom --partition=batch --time=7-1:00:00 --mem=256G --nodes=1 --ntasks=1 --cpus-per-task=4 run_cellphonedb.sh Healthy Liver false
# sbatch --account=smontgom --partition=batch --time=7-1:00:00 --mem=256G --nodes=1 --ntasks=1 --cpus-per-task=4 run_cellphonedb.sh DownSyndrome Liver false
# sbatch --account=smontgom --partition=batch --time=7-1:00:00 --mem=500G --nodes=1 --ntasks=1 --cpus-per-task=1 /oak/stanford/groups/smontgom/amarder/t21-proj/scripts/run_cellphonedb.sh DownSyndrome Liver true
# --nodes=1 --ntasks=1 --cpus-per-task=1

echo "Amount of memory used: ${SLURM_MEM_PER_NODE}"
echo "Amount of CPUs/cores requested: ${SLURM_CPUS_PER_TASK}"

echo "Start:"

source /oak/stanford/groups/smontgom/amarder/bin/conda_init.sh
conda activate cpdb

disease_status=Healthy
env=Femur

disease_status=$1
envir=$2
toSubsample=$3

numCells=108922 # Healthy Liver cell count

headDir=/oak/stanford/groups/smontgom/amarder/t21-proj

# middleDir=out/subset/data/cpdb_Data
# outmiddleDir=out/subset/results/cpdb

middleDir=out/full/cpdb_Data
outmiddleDir=out/full/cpdb_Results

dir=$headDir/$middleDir
outDir=$headDir/$outmiddleDir

meta=${dir}/10X_${disease_status}_${envir}.meta.txt
count=${dir}/10X_${disease_status}_${envir}.norm_count.txt

envir=Liver
disease_status=DownSyndrome
if [[ "$envir" == "Liver" && "$disease_status" == "DownSyndrome" ]]; then
meta=${dir}/10X_${disease_status}_${envir}.meta.sub.txt
count=${dir}/10X_${disease_status}_${envir}.norm_count.sub.csv
fi

mkdir -p $outDir
if [ "$toSubsample" = true ]; then

echo "Running with subsampling..."
cellphonedb method statistical_analysis $meta $count --counts-data hgnc_symbol --output-path $outDir --project-name ${disease_status}_${envir} --threshold 0.1 --threads ${SLURM_CPUS_PER_TASK} --subsampling --subsampling-log false --subsampling-num-cells $numCells

else

echo "Count file: $count"

echo "Running with absolutely no subsampling whatsoever..."
cellphonedb method statistical_analysis $meta $count --counts-data hgnc_symbol --output-path $outDir --project-name ${disease_status}_${envir} --threshold 0.1 --threads ${SLURM_CPUS_PER_TASK}

fi


