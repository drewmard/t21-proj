toSubsample=false
numCells=10000

headDir=/oak/stanford/groups/smontgom/amarder/t21-proj

# middleDir=out/subset/data/cpdb_Data
# outmiddleDir=out/subset/results/cpdb

middleDir=/out/full/cpdb_Data
outmiddleDir=out/full/cpdb_Results

dir=$headDir/$middleDir
outDir=$headDir/$outmiddleDir

disease_status=Healthy
env=Femur
meta=${dir}/10X_${disease_status}_${env}.meta.txt
count=${dir}/10X_${disease_status}_${env}.norm_count.txt

mkdir -p $outDir
if [ "$toSubsample" = true ]; then

cellphonedb method statistical_analysis $meta $count --counts-data hgnc_symbol --output-path $outDir --project-name ${disease_status}_${env} --threshold 0.1 --threads 64 --subsampling --subsampling-log false --subsampling-num-cells $numCells

else

cellphonedb method statistical_analysis $meta $count --counts-data hgnc_symbol --output-path $outDir --project-name ${disease_status}_${env} --threshold 0.1 --threads 64

fi


