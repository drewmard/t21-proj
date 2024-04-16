df1=fread("~/Downloads/cell_type_frequencies.pivot.separated.csv",data.table = F,stringsAsFactors = F)
df2=fread("~/Downloads/cell_type_frequencies.pivot.combined.csv",data.table = F,stringsAsFactors = F)

df1$Erythroid=df1$`Late erythroid cells`+df1$`Early erythroid cells`+df1$`Cycling erythroid cells`
df2$Erythroid=df2$`Late erythroid cells`+df2$`Early erythroid cells`+df2$`Cycling erythroid cells`
t.test(df2$Erythroid~df2$environment)
t.test(df1$Erythroid~df1$environment)
t.test(df1$`Hepatic stellate cells`~df1$environment)
t.test(df1$Hepatocytes~df1$environment)
t.test(df2$Hepatocytes~df2$environment)
t.test(df2$`Hepatic stellate cells`~df2$environment)

cor.test(df1[,"Late erythroid cells"],df2[,"Late erythroid cells"])
cor.test(df1[,"B cells"],df2[,"B cells"])
cor.test(df1[,"Megakaryocytes"],df2[,"Megakaryocytes"])
plot(df1[,"Megakaryocytes"],df2[,"Megakaryocytes"])
cor.test(df1[,"Megakaryocytes"],df2[,"Megakaryocytes"])

val = cor(t(df1[,3:ncol(df2)]),t(df2[1,3:ncol(df2)]))

t.test(val~df1$environment)

df1[12,]
df2[12,]

melted_cor = cor(df1[,3:ncol(df1)],df2[,3:ncol(df2)]) %>% melt()
melted_cor = subset(melted_cor,Var1==Var2)
fwrite(subset(melted_cor,Var1==Var2)[,-2],"~/Downloads/combined_vs_separate_corr.txt",quote = F,na = "NA",sep = '\t',row.names = F,col.names = F)

cor(df1[,4:5],df2[,4:5]) %>% melt()

cor.test(df1[,"B cells"],df2[,"B cells"])

plot(df1[,"Early erythroid cells"],df2[,"Early erythroid cells"])
