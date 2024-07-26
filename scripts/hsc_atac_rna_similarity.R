library(data.table)
f = "~/Downloads/Cell_mutual_information_ATAC_RNA_similarity (1).csv"
df=fread(f,data.table = F,stringsAsFactors = F)
wilcox.test(pointwise_mutual_information~Condition,df)$p.value
# t.test(pointwise_mutual_information~Condition,df)

f = "~/Downloads/Megakaryocyte_Differentiation_Probability_raw_data.csv"
df=fread(f,data.table = F,stringsAsFactors = F)
binom.test(124,1791+124,40/(40+2348))$p.value
# chisq.test(matrix(c(124,1791,40,2348),2,2))$p.value

f = "~/Downloads/Erythroid_Differentiation_Probability_raw_data.csv"
df=fread(f,data.table = F,stringsAsFactors = F)
binom.test(396,1519+396,138/(2250+138))$p.value
# chisq.test(matrix(c(138,2250,396,1519),2,2))$p.value
