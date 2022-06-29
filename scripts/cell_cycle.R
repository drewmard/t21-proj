# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

dfcombined <- CellCycleScoring(dfcombined, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

df1 = subset(dfcombined@meta.data[,c("Phase","S.Score","G2M.Score","environment","organ")],organ=="Liver" & environment=="Healthy")
table(df1$Phase)/nrow(df1)
df2 = subset(dfcombined@meta.data[,c("Phase","S.Score","G2M.Score","environment","organ")],organ=="Liver" & environment=="Down Syndrome")
table(df2$Phase)/nrow(df2)
df3 = subset(dfcombined@meta.data[,c("Phase","S.Score","G2M.Score","environment","organ")],organ=="Femur" & environment=="Healthy")
table(df3$Phase)/nrow(df3)
df4 = subset(dfcombined@meta.data[,c("Phase","S.Score","G2M.Score","environment","organ")],organ=="Femur" & environment=="Down Syndrome")
table(df4$Phase)/nrow(df4)

phase_to_test="S"
prop.test(c(sum(df1$Phase==phase_to_test),sum(df2$Phase==phase_to_test)),
          c(nrow(df1),nrow(df2))
          )
prop.test(c(sum(df3$Phase==phase_to_test),sum(df4$Phase==phase_to_test)),
          c(nrow(df3),nrow(df4))
)
phase_to_test="G1"
prop.test(c(sum(df1$Phase==phase_to_test),sum(df3$Phase==phase_to_test)),
          c(nrow(df1),nrow(df3))
)

df1 = subset(dfcombined@meta.data[,c("PHASE","S.Score","G2M.Score","environment","organ")],organ=="Liver" & environment=="Healthy")
table(df1$PHASE)/nrow(df1)
df2 = subset(dfcombined@meta.data[,c("PHASE","S.Score","G2M.Score","environment","organ")],organ=="Liver" & environment=="Down Syndrome")
table(df2$PHASE)/nrow(df2)
df3 = subset(dfcombined@meta.data[,c("PHASE","S.Score","G2M.Score","environment","organ")],organ=="Femur" & environment=="Healthy")
table(df3$PHASE)/nrow(df3)
df4 = subset(dfcombined@meta.data[,c("PHASE","S.Score","G2M.Score","environment","organ")],organ=="Femur" & environment=="Down Syndrome")
table(df4$PHASE)/nrow(df4)
phase_to_test="G2M+S"
prop.test(c(sum(df1$PHASE==phase_to_test),sum(df2$PHASE==phase_to_test)),
          c(nrow(df1),nrow(df2))
)
phase_to_test="G2M+S"
prop.test(c(sum(df1$PHASE==phase_to_test),sum(df3$PHASE==phase_to_test)),
          c(nrow(df1),nrow(df2))
)

summary(lm(as.numeric(PHASE=="G2M+S")~environment+organ,data=dfcombined@meta.data))
summary(lm(as.numeric(Phase!="G1")~environment+organ,data=dfcombined@meta.data))

cor(dfcombined$G2M.Score,dfcombined$G2M_score)
cor(dfcombined$S.Score,dfcombined$S_score)

wilcox.test(df1$S.Score,df2$S.Score)
wilcox.test(df1$S.Score,df3$S.Score)
wilcox.test(df2$S.Score,df4$S.Score)
wilcox.test(df3$S.Score,df4$S.Score)

wilcox.test(df1$G2M.Score,df2$G2M.Score)
wilcox.test(df1$G2M.Score,df3$G2M.Score)
wilcox.test(df2$G2M.Score,df4$G2M.Score)
wilcox.test(df3$G2M.Score,df4$G2M.Score)
