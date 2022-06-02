sampletype="Liver"
res.df.all.p <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".p.rds"))
res.df.all.lfc <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".lfc.rds"))
cell_type_groups <- colnames(res.df.all.p[[1]])[(which(colnames(res.df.all.p[[1]])=="names")+1):(which(colnames(res.df.all.p[[1]])=="chromosome_name")-1)]
iter=2
df1.lfc=res.df.all.lfc[[iter]]
df1.p=res.df.all.p[[iter]]

x <- x[!is.na(x)]
df.sub <- df.lfc[df.lfc$names %in% x,c(cell_type,"names")]
gene_lst <- df.sub[,"names"]

df.obs <- list()
chrNum=22
for (chrNum in seq(1:22)) {
  print(chrNum)
  df1.lfc.sub <- subset(df1.lfc,chromosome_name==chrNum)
  df1.p.sub <- subset(df1.p,chromosome_name==chrNum)
  df1.lfc.sort <- df1.lfc.sub[order(df1.lfc.sub$chromosome_name,df1.lfc.sub$start_position),]
  df1.p.sort <- df1.p.sub[order(df1.lfc.sub$chromosome_name,df1.lfc.sub$start_position),]
  # df1.lfc.sort$Megakaryocytes <- rep(c(-1,1),length.out=nrow(df1.lfc.sub))
  flip_ct.lst <- list()
  j=1
  for (j in 1:length(cell_type_groups)) {
    cell_type = cell_type_groups[j]
    flip_ct = 0; prev <- NULL
    i=2
    for (i in 1:nrow(df1.lfc.sort)) {
      curr <- df1.lfc.sort[i,cell_type]
      curr.p <- df1.p.sort[i,cell_type]
      if (is.null(curr.p) | is.na(curr.p)) {next}
      if (i==1 | is.null(prev)) {
        prev <- curr
      } else {
        if (sign(curr)!=sign(prev)) {
          flip_ct = flip_ct + 1
          prev <- curr
        }
      }
    }
    flip_ct.lst[[cell_type]] <- flip_ct
  }
  df.obs[[chrNum]] <- as.data.frame(flip_ct.lst)
  df.obs[[chrNum]]$chr <- chrNum
}
df.obs.all <- as.data.frame(do.call(rbind,df.obs))
rownames(df.obs.all) <- NULL

###########

df.perm.chr <- list()
for (chrNum in 6:22) {
  print(paste0("chr: ",chrNum))
  df1.lfc.sub <- subset(df1.lfc,chromosome_name==chrNum)#[1:20,]
  df1.p.sub <- subset(df1.p,chromosome_name==chrNum)#[1:20,]
  df.perm <- list()
  for (k in 1:10) {
    print(paste0("iter: ",k))
    idx <- sample(1:nrow(df1.lfc.sub),size = nrow(df1.lfc.sub),replace = F)
    df1.lfc.sort <- df1.lfc.sub[idx,]
    df1.p.sort <- df1.p.sub[idx,]
    flip_ct.lst <- list()
    # for (j in 8:8) {
    for (j in 1:length(cell_type_groups)) {
      cell_type = cell_type_groups[j]
      flip_ct = 0; prev <- NULL
      for (i in 1:nrow(df1.lfc.sort)) {
        # for (i in 1:10) {
        curr <- df1.lfc.sort[i,cell_type]
        curr.p <- df1.p.sort[i,cell_type]
        if (is.null(curr.p) | is.na(curr.p)) {next}
        if (i==1 | is.null(prev)) {
          prev <- curr
        } else {
          if (sign(curr)!=sign(prev)) {
            flip_ct = flip_ct + 1
            prev <- curr
          }
        }
      }
      flip_ct.lst[[cell_type]] <- flip_ct
    }
    df.perm[[k]] <- as.data.frame(flip_ct.lst)
  }
  df.perm.chr[[chrNum]] <- do.call(rbind,df.perm)
}


p.perm.df <- list()
for (chrNum in 1:22) {
  m = ncol(df.perm.chr[[chrNum]])
  p.perm <- rep(NA,m)
  for (j in 1:m) {
    cell_type = colnames(df.perm.chr[[chrNum]])[j]
    p.perm[j] <- mean(df.obs.all[df.obs.all$chr==chrNum,cell_type] < df.perm.chr[[chrNum]][,cell_type])
  }
  p.perm.df[[chrNum]] <- p.perm
}
p.perm.df.all <- as.data.frame(do.call(rbind,p.perm.df))
colnames(p.perm.df.all) <- colnames(df.perm.chr[[1]])
p.perm.df.all$chr <- 1:nrow(p.perm.df.all)
subset(p.perm.df.all,chr==21)

apply(p.perm.df.all,1,function(x) mean(x <= 0.1))


#######

res.df.all.p <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".p.rds"))
res.df.all.lfc <- readRDS(paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/res.",sampletype,".lfc.rds"))
cell_type_groups <- colnames(res.df.all.p[[1]])[(which(colnames(res.df.all.p[[1]])=="names")+1):(which(colnames(res.df.all.p[[1]])=="chromosome_name")-1)]

for (cell_type in cell_type_groups) {
  keep_genes <- res.df.all.p[[2]][!is.na(res.df.all.p[[2]][,cell_type]),'names']
  res.df.all.p[[1]][!(res.df.all.p[[1]][,'names']%in%keep_genes),cell_type] <- NA
  res.df.all.lfc[[1]][!(res.df.all.lfc[[1]][,'names']%in%keep_genes),cell_type] <- NA
}

tab <- lapply(res.df.all.p,function(x) {aggregate(x[,cell_type_groups]<0.1,by=list(x$chr21),mean,na.rm=T)})
tab <- lapply(tab,function(x) {rownames(x) <- x[,1]; x <- x[,-1]; x <- as.data.frame(t(x))})
tab <- do.call(cbind,tab)

colnames(tab)[1:2] <- paste0(colnames(tab)[1:2]," 1")
colnames(tab)[3:4] <- paste0(colnames(tab)[3:4]," 2")
colnames(tab)[5:6] <- paste0(colnames(tab)[5:6]," 3")

tab$celltype <- rownames(tab)






