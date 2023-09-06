library(data.table)
library(glmmTMB)

df = fread("~/Downloads/mitosox_exp_230719_2.csv",data.table = F,stringsAsFactors = F)
df$Sample = as.factor(df$Sample)
cell.uniq = unique(df$Cell_type)
res.lst=list()

# 
# df = fread("~/Downloads/mtg_exp_230719_2.csv",data.table = F,stringsAsFactors = F)
# df$Sample = as.factor(df$Sample)
# cell.uniq = unique(df$Cell_type)
# res.lst=list()
# 
aggre = aggregate(df$MitoSOX,df[,c("Sample","Status","Age","Cell_type")],mean)
library(ggplot2)
ggplot(aggre,aes(x=Status,y=x)) + geom_violin() + geom_point() + facet_grid(cols=vars(Cell_type)) + ggpubr::theme_pubr() + labs(y="MitoSOX")
ggplot(df,aes(x=Status,y=MitoSOX,col=Sample)) + geom_boxplot() + facet_grid(cols=vars(Cell_type)) + theme_bw()

for (i in 1:length(cell.uniq)) {
  print(i)
  cell = cell.uniq[i]
  df.sub = subset(df,Cell_type==cell)
  
  library(lme4)
  df.sub$Status = factor(df.sub$Status,levels = c("Healthy","DS"))
  min_val = min(df.sub$MitoSOX)
  if (min_val < 0) {
    df.sub$MitoSOX = df.sub$MitoSOX + abs(min_val) + 1
  }
  df.sub$MitoSOX = round(df.sub$MitoSOX)
  # print(summary(lm(MitoSOX~Status + Age + (Sample),df.sub)))
# }
  # hist(df.sub$MitoSOX)
  # model = lmer(MitoSOX ~ Status + (1 | Sample), data = df.sub)
  # summary(model)
  # 
  # model <- glmer(MitoSOX ~ Status + (1 | Sample), data = df.sub)#, family = glmmTMB::nbinom1)
  # summary(model)
  # 
  # model <- glmer(MitoSOX ~ Status + (1 | Sample), data = df.sub, family = glmmTMB::nbinom1)
  # summary(model)
  
  # model <- glmmTMB(MitoSOX~Status + scale(Age) + (1 | Sample), df.sub, family = gaussian, REML = FALSE)
  # summary(model)
  
  df.sub$MitoSOX = rntransform(df.sub$MitoSOX)
  
  model <- glmmTMB(MitoSOX~Status + Age + (1 | Sample), df.sub, family = gaussian, REML = FALSE)
  
  # model <- glmmTMB(MitoSOX~Status + Age + (1 | Sample), df.sub, family = nbinom1, REML = FALSE)
  tab=coef(summary(model))[[1]]
  res = tab["StatusDS",]
  # res
  res.lst[[i]] = res
}

res.df = data.frame(do.call(rbind,res.lst))
colnames(res.df) = c('est','se','stat','pval')
res.df$cell = cell.uniq
fwrite(res.df,"~/Downloads/mitosox_results.txt",quote = F,na = "NA",sep = "\t",row.names = F,col.names = T)

