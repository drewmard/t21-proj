library(data.table)

library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# f=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data_pb_leiden/",sampletype,".pb.",cell_type_filename,".",subset_column,".txt")
# df.aggre <- fread(f,data.table = F,stringsAsFactors = F)

annot <- getBM(attributes = c('hgnc_symbol', 'chromosome_name',
                              'start_position', 'end_position'),
               filters = 'hgnc_symbol', 
               values = rownames(df.aggre), 
               mart = ensembl)
annot1 <- annot[annot$chromosome_name %in% c(as.character(seq(1,22)),"X","Y"),]
fileOut <- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/data_small/gene_annot.txt")
fwrite(annot1,fileOut,na = "NA",sep = '\t',quote = F,row.names = F,col.names = T)
