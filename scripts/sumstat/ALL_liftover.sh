echo "Running liftOver..."
HEADDIR=/oak/stanford/groups/smontgom/amarder/t21-proj/scripts
bedFile=$HEADDIR/ALL_finemap.hg19.bed
hg19ToHg38chain=/oak/stanford/groups/smontgom/amarder/bin/chain/hg19ToHg38.over.chain.gz
source /oak/stanford/groups/smontgom/amarder/bin/conda_init.sh
conda activate liftover
liftOver $bedFile $hg19ToHg38chain $HEADDIR/ALL_finemap.hg38.bed $HEADDIR/ALL_finemap.unmapp.bed
conda deactivate 
