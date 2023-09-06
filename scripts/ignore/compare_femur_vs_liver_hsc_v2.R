library(tls)
library(ggplot2)

tls_bootstrap <- function() {
  Niter=1000; tls.boot <- rep(NA,Niter)
  for (iter in 1:Niter) {
    M=nrow(df.tmp)
    df.tmp2 <- df.tmp[sample(1:M,M,replace = T),]
    tls.boot[iter] <- tls(Liver~Femur - 1,data=df.tmp2)$coefficient
  }
  tls.boot <- sort(tls.boot)
  ci_low=tls.boot[25]; ci_high=tls.boot[975]
  return(c(ci_low,ci_high))
}

################

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
res.df.all.p.liver <- res.df.all.p[[2]]
cell_type_groups.liver <- cell_type_groups

sampletype="Femur"
res.df.all.p <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".p.rds"))
res.df.all.lfc <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".lfc.rds"))
cell_type_groups <- colnames(res.df.all.p[[1]])[(which(colnames(res.df.all.p[[1]])=="names")+1):(which(colnames(res.df.all.p[[1]])=="chromosome_name")-1)]

for (cell_type in cell_type_groups) {
  keep_genes <- res.df.all.p[[2]][!is.na(res.df.all.p[[2]][,cell_type]),'names']
  res.df.all.p[[1]][!(res.df.all.p[[1]][,'names']%in%keep_genes),cell_type] <- NA
  res.df.all.lfc[[1]][!(res.df.all.lfc[[1]][,'names']%in%keep_genes),cell_type] <- NA
}

res.df.all.lfc.femur <- res.df.all.lfc[[2]]
res.df.all.p.femur <- res.df.all.p[[2]]
cell_type_groups.femur <- cell_type_groups

cell_type_groups <- cell_type_groups.femur[cell_type_groups.femur %in% cell_type_groups.liver]

df.res <- list()
# for (i in 1:4) {
for (i in 1:length(cell_type_groups)) {
  print(i)
  cell_type=cell_type_groups[i]
  print(cell_type)
  df1 <- res.df.all.lfc.liver[,c("names",cell_type)]; df1$P.liver <- res.df.all.p.liver[,c(cell_type)]
  df2 <- res.df.all.lfc.femur[,c("names","chromosome_name",cell_type)]; df2$P.femur <- res.df.all.p.femur[,c(cell_type)]
  colnames(df1)[2] <- "Liver"
  colnames(df2)[3] <- "Femur"
  
  df.mg <- merge(df1,df2,by='names')
  # df.tmp <- subset(df.mg,!is.na(Femur) & !is.na(Liver) & chromosome_name==21)
  # tls(Liver~Femur - 1,data=df.tmp)
  # df.tmp <- subset(df.mg,!is.na(Femur) & !is.na(Liver) & chromosome_name==21)
  # tls(Liver~Femur - 1,data=df.tmp)
  
  set.seed(031995)
  df.tmp <- subset(df.mg,!is.na(Femur) & !is.na(Liver) & chromosome_name==21 & P.liver < 0.25)
  ci <- tls_bootstrap()
  res1 <- data.frame(chr="chr21",cell_type,coef=tls(Liver~Femur - 1,data=df.tmp)$coefficient,ci_low=ci[1],ci_hi=ci[2])
  
  df.tmp <- subset(df.mg,!is.na(Femur) & !is.na(Liver) & chromosome_name!=21 & P.liver < 0.25)
  coef <- tls(Liver~Femur - 1,data=df.tmp)$coefficient
  ci <- tls_bootstrap()
  res2 <- data.frame(chr="not chr21",cell_type,coef=coef,ci_low=ci[1],ci_hi=ci[2])
  
  df.res[[cell_type]] <- rbind(res1,res2)
}
df.res.all <- do.call(rbind,df.res); rownames(df.res.all) <- NULL; df.res.all

tmp <- subset(df.res.all,chr=="chr21")
mean(tmp$coef)
mean((tmp$ci_low > 1 & tmp$ci_hi > 1) | (tmp$ci_low < 1 & tmp$ci_hi < 1))
mean((tmp$ci_low > 1 & tmp$ci_hi > 1))
mean((tmp$ci_low < 1 & tmp$ci_hi < 1))

tmp <- subset(df.res.all,chr=="not chr21")
mean(tmp$coef)

mean((tmp$ci_low > 1 & tmp$ci_hi > 1) | (tmp$ci_low < 1 & tmp$ci_hi < 1))
mean((tmp$ci_low > 1 & tmp$ci_hi > 1))
mean((tmp$ci_low < 1 & tmp$ci_hi < 1))
subset(tmp,ci_low > 1 & ci_hi > 1)

ci <- tls_bootstrap()

library(ggplot2)
ggplot(df.tmp,aes(x=Femur,y=Liver)) + geom_smooth() + geom_point()
tls(Liver~Femur - 1,data=df.tmp)
tls(Liver~Femur - 1,data=subset(df.tmp,P.femur < 0.25))

summary(lm(Liver~Femur - 1,data=df.tmp))
summary(lm(Liver~Femur,data=df.tmp))


print(i)
cell_type="Late erythroid cells"
cell_type="HSCs_MPPs"
print(cell_type)
df1 <- res.df.all.lfc.liver[,c("names",cell_type)]; df1$P.liver <- res.df.all.p.liver[,c(cell_type)]
df2 <- res.df.all.lfc.femur[,c("names","chromosome_name",cell_type)]; df2$P.femur <- res.df.all.p.femur[,c(cell_type)]
colnames(df1)[2] <- "Liver"
colnames(df2)[3] <- "Femur"

df.mg <- merge(df1,df2,by='names')
# df.tmp <- subset(df.mg,!is.na(Femur) & !is.na(Liver) & chromosome_name==21)
# tls(Liver~Femur - 1,data=df.tmp)
# df.tmp <- subset(df.mg,!is.na(Femur) & !is.na(Liver) & chromosome_name==21)
# tls(Liver~Femur - 1,data=df.tmp)

set.seed(031995)
df.tmp <- subset(df.mg,!is.na(Femur) & !is.na(Liver) & chromosome_name==21 & P.liver < 0.25)
coef <- tls(Liver~Femur - 1,data=df.tmp)$coefficient; coef
df.tmp <- subset(df.mg,!is.na(Femur) & !is.na(Liver) & chromosome_name!=21 & P.liver < 0.25)
coef <- tls(Liver~Femur - 1,data=df.tmp)$coefficient; coef
library("deming")
?deming
# res <- deming(Liver~Femur,data=df.tmp)
# unname(res$coefficients[2])
tls(Liver~Femur - 1,data=df.tmp)$coefficient
lm(Liver~Femur - 1,data=df.tmp)$coefficient
plot(Liver,Femur,data=df.tmp)

# ggplot(subset(df.mg,chromosome_name==21),aes(x=Liver,y=Femur)) + 
df.tmp <- subset(df.mg,!is.na(Femur) & !is.na(Liver) & chromosome_name!=21 & P.liver < 0.25)
df.tmp <- subset(df.mg,!is.na(Femur) & !is.na(Liver) & chromosome_name!=21)

ggplot(df.tmp,aes(x=Femur,y=Liver)) +
  geom_point() + 
  geom_abline(slope=1,intercept = 0,lty='dashed',col='red') + 
  # geom_abline(slope=tls(Liver~Femur - 1,data=df.tmp)$coefficient,intercept=0,col='blue',size=1) + 
  # geom_abline(slope=lm(Liver~Femur - 1,data=df.tmp)$coefficient,intercept=0,col='purple',size=1) + 
  theme_bw() + labs(x="logFC (femur)",y='logFC (liver)') + geom_smooth(method='lm',se=F)
df.mg[order(df.mg$Liver,decreasing = T)[1:3],]
subset(df.mg,names=="MEG3")
subset(df.mg,names=="GAS6")
subset(df.mg,Liver > Femur & Femur > 0 & P.liver < 0.05 & P.femur < 0.05)[1,]
subset(df.mg,Femur > Liver & Femur > 0 & P.liver < 0.05 & P.femur < 0.05)[1,]
subset(df.mg,Femur < Liver & Femur < 0 & Liver < 0  & P.liver < 0.05 & P.femur < 0.05)[1,]
subset(df.mg,Femur < Liver & Femur < 0 & Liver < 0  & P.liver < 0.05)[1,]
tmp <- subset(df.mg,sign(Femur) == sign(Liver) & P.liver < 0.5 & P.femur < 0.5)

tmp <- subset(df.mg,(P.liver < 0.05 | P.femur < 0.05) & chromosome_name==21)
table(sign(tmp$Femur),sign(tmp$Liver))
table(sign(tmp$Femur),sign(tmp$Liver))

mean(abs(tmp$Liver) > abs(tmp$Femur),na.rm = T)
which.max(abs(tmp$Liver) - abs(tmp$Femur))
which.max(abs(tmp$Femur) - abs(tmp$Liver))
tmp[c(230,353),]

tmp <- subset(df.mg,Femur > 0 & Liver > 0)
mean(df.mg$Femur < df.mg$Liver,na.rm = T)

tmp <- subset(df.mg,Femur < Liver & Femur < 0 & Liver < 0 & P.liver < 0.5)
nrow(tmp)
tmp <- subset(df.mg,Femur > Liver & Femur < 0 & Liver < 0 & P.liver < 0.5)
nrow(tmp)

subset(P.liver < 0.05)

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

