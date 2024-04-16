library(data.table)
library(reshape2)
meta = fread("~/Downloads/celltype_align.csv",data.table = F,stringsAsFactors = F)
cormat = fread("~/Downloads/CorrMat.public_vs_ours.csv",data.table = F,stringsAsFactors = F)
# rownames(cormat)=cormat[,1];cormat = cormat[,-1]
cormelt = melt(cormat,id.vars="V1",variable.name = "column")
cormelt$id = paste(cormelt[,1],cormelt[,2])
colnames(cormelt)[3] = "Public"
cormelt1 = subset(cormelt,id %in% paste(meta$Public,meta$Disomy))

cormat = fread("~/Downloads/CorrMat.disomySections.2condi.csv",data.table = F,stringsAsFactors = F)
# rownames(cormat)=cormat[,1];cormat = cormat[,-1]
cormelt = melt(cormat,id.vars="V1",variable.name = "column")
cormelt$id = paste(cormelt[,1],cormelt[,2])
colnames(cormelt)[3] = "Disomy"
cormelt2 = subset(cormelt,id %in% paste(meta$Disomy,meta$Disomy))

cormat = fread("~/Downloads/CorrMat.t21Sections.2condi.csv",data.table = F,stringsAsFactors = F)
# rownames(cormat)=cormat[,1];cormat = cormat[,-1]
cormelt = melt(cormat,id.vars="V1",variable.name = "column")
cormelt$id = paste(cormelt[,1],cormelt[,2])
colnames(cormelt)[3] = "Trisomy"
cormelt3 = subset(cormelt,id %in% paste(meta$Trisomy,meta$Disomy))

cormat = fread("~/Downloads/CorrMat.disjoint_disomy_sample_sets.csv",data.table = F,stringsAsFactors = F)
# rownames(cormat)=cormat[,1];cormat = cormat[,-1]
cormelt = melt(cormat,id.vars="V1",variable.name = "column")
cormelt$id = paste(cormelt[,1],cormelt[,2])
colnames(cormelt)[3] = "DisomySep"
cormelt4 = subset(cormelt,id %in% paste(meta$Disomy,meta$Disomy))

cormat = fread("~/Downloads/CorrMat.disjoint_t21_sample_sets.csv",data.table = F,stringsAsFactors = F)
# rownames(cormat)=cormat[,1];cormat = cormat[,-1]
cormelt = melt(cormat,id.vars="V1",variable.name = "column")
cormelt$id = paste(cormelt[,1],cormelt[,2])
colnames(cormelt)[3] = "TrisomySep"
cormelt5 = subset(cormelt,id %in% paste(meta$Trisomy,meta$Trisomy))

cormat = fread("~/Downloads/CorrMat.allSections.2condi.csv",data.table = F,stringsAsFactors = F)
# rownames(cormat)=cormat[,1];cormat = cormat[,-1]
cormelt = melt(cormat,id.vars="V1",variable.name = "column")
cormelt$id = paste(cormelt[,1],cormelt[,2])
colnames(cormelt)[3] = "Trisomy_v_Disomy.all"
cormelt6 = subset(cormelt,id %in% paste(meta$Disomy,meta$Trisomy))

res = merge(merge(cormelt2,cormelt3,by='column',all = T),cormelt1,by="column",all = T)
t.test(res$Public,res$Disomy)
t.test(res$Disomy,res$Trisomy)
ggplot(res,aes())
boxplot_data <- data.frame(
  Group = rep(c("Disomy", "Trisomy"), each = length(res$Disomy)),
  Value = c(res$Disomy, res$Trisomy)
)

# Create boxplot
library(ggplot2)
boxplot <- ggplot(boxplot_data, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot() +
  labs(title = "",
       x = "Sections",
       y = "Correlations") +
  ggpubr::theme_pubr() +
  guides(fill="none")

pdf(paste0("~/Documents/Research/t21-proj/out/full/visium_plots/","disomy_v_trisomy",".pdf"),width = 6,height=5)
print(boxplot)
dev.off()

t.test(res$Trisomy,res$Public)

res2 = merge(merge(res,cormelt4,by='column',all = T),cormelt5,by="column",all = T)

t.test(res2$Disomy,res2$DisomySep)
t.test(res2$DisomySep,res2$TrisomySep)

boxplot_data <- data.frame(
  Group = rep(c("Disomy (All) v Trisomy", "Disomy 1 v Disomy 2"), times = c(length(res2$Disomy),length(res2$DisomySep))),
  Value = c(res2$Disomy, res2$DisomySep)
)

# Create boxplot
library(ggplot2)
boxplot <- ggplot(boxplot_data, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot() +
  labs(title = "",
       x = "References",
       y = "Correlations") +
  ggpubr::theme_pubr() +
  guides(fill="none")

t_test_result = t.test(Value~Group,boxplot_data)
  if(t_test_result$p.value < 0.05) {
    boxplot <- boxplot + 
      geom_text(aes(x = 1.5, y = max(boxplot_data$Value,na.rm = T), label = "p < 0.05"), color = "black")
  }

pdf(paste0("~/Documents/Research/t21-proj/out/full/visium_plots/","disomy_disjoint_compare",".pdf"),width = 6,height=5)
print(boxplot)
dev.off()




boxplot_data <- data.frame(
  Group = rep(c("Disomy 1 v Disomy 2", "Trisomy 1 v Trisomy 2"), times = c(length(res2$Disomy),length(res2$DisomySep))),
  Value = c(res2$DisomySep, res2$TrisomySep)
)

# Create boxplot
library(ggplot2)
boxplot <- ggplot(boxplot_data, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot() +
  labs(title = "",
       x = "References",
       y = "Correlations") +
  ggpubr::theme_pubr() +
  guides(fill="none")

t_test_result = t.test(Value~Group,boxplot_data)
if(t_test_result$p.value < 0.05) {
  boxplot <- boxplot + 
    geom_text(aes(x = 1.5, y = max(boxplot_data$Value,na.rm = T), label = "p < 0.05"), color = "black")
}
boxplot

pdf(paste0("~/Documents/Research/t21-proj/out/full/visium_plots/","all_disjoint_compare",".pdf"),width = 6,height=5)
print(boxplot)
dev.off()


res3 = merge(res2,cormelt6,by='column',all = T)
t.test(res3$Public,res3$Disomy)
t.test(res3$Public,res3$Trisomy_v_Disomy.all)
t.test(res3$Public,res3$DisomySep)


boxplot_data <- data.frame(
  Group = rep(c("Public v Disomy","Trisomy v Disomy", "Disomy 1 v Disomy 2"), times = c(length(res2$Disomy),length(res2$DisomySep),length(res3$DisomySep))),
  Value = c(res3$Public, res3$Disomy,res3$DisomySep)
)

# Create boxplot
library(ggplot2)
boxplot <- ggplot(boxplot_data, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot() +
  labs(title = "",
       x = "References",
       y = "Correlations") +
  ggpubr::theme_pubr() +
  guides(fill="none")
pdf(paste0("~/Documents/Research/t21-proj/out/full/visium_plots/","all_datasets_compared",".pdf"),width = 6,height=5)
print(boxplot)
dev.off()

boxplot

