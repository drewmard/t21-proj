# load
library(data.table)

# input:
fmeta1 = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/downsample_induced/input/meta1.txt")
fmeta2 = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/downsample_induced/input/meta2.txt")
meta1=fread(fmeta1,data.table = F,stringsAsFactors = F)
meta2=fread(fmeta2,data.table = F,stringsAsFactors = F)
match_between_conditions = "sorting" # here, conditions is trisomy vs disomy
match_between_differential_factor = "patient" # and diff factor is cyc vs less cycling
pseudobulk_covariate = "sample"

# preprocess
# meta1$diff_factor = "diff1"
# meta2$diff_factor = "diff2"

# need to find data that has enough cells to match on

# create a list by the "match_between_differential_factor"
# meta1.input = dfseurat1@meta.data
tab = aggregate(meta1.input[,1],by=meta1.input[,c(pseudobulk_covariate,match_between_differential_factor,match_between_conditions)],length)
colnames(tab)[ncol(tab)] = "num_cells"
tab = subset(tab,num_cells >= 10)
meta1.input.lst <- split(tab, tab$patient)

# and then lapply through list to find all that contain all "match_between_conditions" and more cells in each one
keepers = list()
for (match1 in unique(meta1[,match_between_differential_factor])) {
  
  relevant_idx = meta1[,match_between_differential_factor]==match1
  meta1.sub = meta1[relevant_idx,]
  
  x = meta1.input.lst[[1]]
  
  if (!is.null(match_between_conditions)) {
    keepers[[as.character(match1)]]=lapply(meta1.input.lst,function(x) {
      meta1.sub.input = merge(meta1.sub,x,by=match_between_conditions)
      meta1.sub.input = subset(meta1.sub.input,num_cells.y > num_cells.x)
      
      # if one of the "match_between_conditions" is missing, then can't use this one!
      if(any(!(meta1.sub[,match_between_conditions] %in% meta1.sub.input[,match_between_conditions]))) {
        return(data.frame())
      } else {
        return(meta1.sub.input)
      }
    })
  } else {
    keepers[[as.character(match1)]] = lapply(meta1.input.lst,function(x) {
      meta1.sub$rand = "rand"
      x$rand = "rand"
      meta1.sub.input = merge(meta1.sub,x,by="rand")
      meta1.sub.input = subset(meta1.sub.input,num_cells.y > num_cells.x)
      # if one of the "match_between_conditions" is missing, then can't use this fetus!
      if(nrow(meta1.sub.input) < nrow(meta1.sub)) {
        return(data.frame())
      } else {
        return(meta1.sub.input)
      }
    })
  }
  y = keepers[[as.character(match1)]]
  keepers[[as.character(match1)]] = c()
  keepers[[as.character(match1)]] = unique(y[,paste0(match_between_differential_factor,".y")])
}
keepers.df = do.call(rbind,keepers)

# this is a simple merge to a tmp
# check nrow
# and then check that num_cells is greater
rbind(meta1,meta2)





# patient, sample, sorting, original_label, integrated_label
meta1 = subset(meta2.full,integrated_label=="HSCs/MPPs")[,c("patient","sample","sorting","original_label","combi_annot")]

# give a data.frame of:
# covariates of interest to match on, e.g. sorting
table

# if sample has both cycling HSCs and HSCs w/ cell count >= 10, then consider those samples

# can set additional restrictions that are dataset-specific:
# 1
# 2
# 3

########################################
# now go through each sample in the Healthy:
# goal: if there are samples that can only be used for the largest samples, then sample from those first 
# sort by min(celltype1,celltype2) size

# find potential matching samples (save to a list)

# next: sample for each sample, starting from the largest
# and remove those used samples from the list for all subsequent samples

# end output: healthy sample, Ts21 sample, celltype1 count X1, celltype2 count X2
# if celltype2 not included, then skip (leave NA)

# now sample cell type 1 from dfseurat1 and cell type 2 from m=dfseurat2:
for (i in 1:nrow(end_output)) {
  # first, find the indices of all cells from Ts21 sample
  ind = unname(which(dfseurat1$sample==dfseurat2$sample))
  ind.lst <- c(ind.lst,sample(ind,res.all2$x[i],replace = FALSE))
}

# run DE analysis:
# note, don't remove based on gene expression since not interested in multiple testing corrections
# just interested in how gene expression changes

######################################################

# condition1 induced
# condition2 induced
# condition-independent
# unknown
# none

######################################################

# how does this compare to interaction?
# problem if low power to detect interactions
# e.g. low expression or variability in one condition
# then interaction tests don't work
# can test this by splitting large samples in 2
# and then setting a 10-100% effect of condition on 1-30% of genes
# and seeing how this method compares to interactions




