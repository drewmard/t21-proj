library(data.table)
df = fread("/oak/stanford/groups/smontgom/amarder/t21-proj/scripts/ALL_finemap.csv",data.table = F,stringsAsFactors = F)
bed = data.frame(df$chromosome,df$pos_hg19-1,df$pos_hg19,df$liftID,df$PP)
f.out = "/oak/stanford/groups/smontgom/amarder/t21-proj/scripts/ALL_finemap.hg19.bed"
fwrite(bed,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = F)

system('
echo "Running liftOver..."
HEADDIR=/oak/stanford/groups/smontgom/amarder/t21-proj/scripts
bedFile=$HEADDIR/ALL_finemap.hg19.bed
hg19ToHg38chain=/oak/stanford/groups/smontgom/amarder/bin/chain/hg19ToHg38.over.chain.gz
source /oak/stanford/groups/smontgom/amarder/bin/conda_init.sh
conda activate liftover
liftOver $bedFile $hg19ToHg38chain $HEADDIR/ALL_finemap.hg38.bed $HEADDIR/ALL_finemap.unmapp.bed
conda deactivate 
')

f.in = "/oak/stanford/groups/smontgom/amarder/t21-proj/scripts/ALL_finemap.hg38.bed"
df = fread(f.in,data.table = F,stringsAsFactors = F)
head(df)

#########


library(data.table)
df = fread("/home/amarder/R/x86_64-pc-linux-gnu-library/4.1/SCAVENGE/extdata/mono.PP001.bed",data.table = F,stringsAsFactors = F)
df = df[,c(1:5)]
df[,2] = df[,2]-1
f.out = "/home/amarder/R/x86_64-pc-linux-gnu-library/4.1/SCAVENGE/extdata/mono.PP001.hg19.bed"
fwrite(df,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = F)

system('
echo "Running liftOver..."
HEADDIR=/home/amarder/R/x86_64-pc-linux-gnu-library/4.1/SCAVENGE/extdata
bedFile=$HEADDIR/mono.PP001.hg19.bed
hg19ToHg38chain=/oak/stanford/groups/smontgom/amarder/bin/chain/hg19ToHg38.over.chain.gz
source /oak/stanford/groups/smontgom/amarder/bin/conda_init.sh
conda activate liftover
liftOver $bedFile $hg19ToHg38chain $HEADDIR/mono.PP001.hg38.bed $HEADDIR/mono.PP001.unmapp.bed
conda deactivate 
')


