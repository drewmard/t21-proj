library(tls)

sampletype="Liver"
res.df.all.p <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".p.rds"))
res.df.all.lfc <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".lfc.rds"))
cell_type_groups <- colnames(res.df.all.p[[1]])[(which(colnames(res.df.all.p[[1]])=="names")+1):(which(colnames(res.df.all.p[[1]])=="chromosome_name")-1)]

for (cell_type in cell_type_groups) {
  keep_genes <- res.df.all.p[[2]][!is.na(res.df.all.p[[2]][,cell_type]),'names']
  res.df.all.p[[1]][!(res.df.all.p[[1]][,'names']%in%keep_genes),cell_type] <- NA
  res.df.all.lfc[[1]][!(res.df.all.lfc[[1]][,'names']%in%keep_genes),cell_type] <- NA
}

res.df.all.lfc.liver <- res.df.all.lfc[[2]]
df1 <- res.df.all.lfc[[2]][,c("names","HSCs_MPPs")]
colnames(df1)[2] <- sampletype

sampletype="Femur"
res.df.all.p <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".p.rds"))
res.df.all.lfc <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".lfc.rds"))
cell_type_groups <- colnames(res.df.all.p[[1]])[(which(colnames(res.df.all.p[[1]])=="names")+1):(which(colnames(res.df.all.p[[1]])=="chromosome_name")-1)]

for (cell_type in cell_type_groups) {
  keep_genes <- res.df.all.p[[2]][!is.na(res.df.all.p[[2]][,cell_type]),'names']
  res.df.all.p[[1]][!(res.df.all.p[[1]][,'names']%in%keep_genes),cell_type] <- NA
  res.df.all.lfc[[1]][!(res.df.all.lfc[[1]][,'names']%in%keep_genes),cell_type] <- NA
}

df2 <- res.df.all.lfc[[2]][,c("names","chromosome_name","HSCs_MPPs")]
colnames(df2)[3] <- sampletype

df.mg <- merge(df1,df2,by='names')
# ggplot(subset(df.mg,chromosome_name==21),aes(x=Liver,y=Femur)) + 
ggplot(subset(df.mg,chromosome_name==21),aes(x=Femur,y=Liver)) +
  geom_point() + 
  geom_abline(slope=1,intercept = 0) + 
  geom_smooth(method='lm') + 
  theme_bw()

tls(y1~x1 - 1,data=data.frame(x1,y1))

summary(lm(Femur~Liver,data=df.mg))
summary(lm(Femur~Liver,data=subset(df.mg,chromosome_name==21)))
summary(lm(Femur~Liver,data=subset(df.mg,chromosome_name!=21)))
summary(lm(Liver~Femur,data=subset(df.mg,chromosome_name==21)))

dim(subset(df.mg,chromosome_name==21))
subset(res.df.all.p[[2]],names=="COL6A1")

t.test(subset(df.mg,chromosome_name==21)$Femur,
       subset(df.mg,chromosome_name==21)$Liver,paired = T)
t.test(abs(df.mg$Femur),
       abs(df.mg$Liver),paired = T)
t.test(abs(subset(df.mg,chromosome_name==21)$Femur),
       abs(subset(df.mg,chromosome_name==21)$Liver),paired = T)
t.test(abs(subset(df.mg,chromosome_name!=21)$Femur),
       abs(subset(df.mg,chromosome_name!=21)$Liver),paired = T)

# wilcox.test(abs(subset(df.mg,chromosome_name==21)$Femur),
#                 abs(subset(df.mg,chromosome_name==21)$Liver),paired = T)

tab <- lapply(res.df.all.p,function(x) {aggregate(x[,cell_type_groups]<0.1,by=list(x$chr21),mean,na.rm=T)})
tab <- lapply(tab,function(x) {rownames(x) <- x[,1]; x <- x[,-1]; x <- as.data.frame(t(x))})
tab <- do.call(cbind,tab)

colnames(tab)[1:2] <- paste0(colnames(tab)[1:2]," 1")
colnames(tab)[3:4] <- paste0(colnames(tab)[3:4]," 2")
colnames(tab)[5:6] <- paste0(colnames(tab)[5:6]," 3")

